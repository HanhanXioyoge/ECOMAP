# sabio_pipeline.py
# -*- coding: utf-8 -*-
"""
SABIO pipeline (with batching, checkpoint-resume, and 2-phase execution):
Phase p1: SABIO lookups (name -> SabioCompoundID; details -> ChebiID/PubChemID/InChI/KEGG)
          - Writes an incremental CSV after each batch (checkpoint-friendly).
          - ChebiID and KeggCompoundID keep multiple values (deduplicated) joined by ';'.
          - InChIKey, MetaNetXID left empty for p1 (to be filled in p2).
Phase p2: Enrichment (MetaNetX and PubChem fallsbacks)
          - Rule 1: If InChIKey missing (or IK & InChI both missing) but SabioCompoundID exists,
                    query MetaNetX by SABIO xref (identifiers.org/sabiork.compound:<ID>)
                    to obtain MetaNetXID, InChI, InChIKey and backfill empties.
          - Rule 2: For rows with InChI present, map InChI -> (MetaNetXID, InChIKey) via MNX.
          - Rule 3: If IK still missing but PubChem CID exists, fetch IK/InChI via PubChem.
          - Writes incremental CSV per batch (checkpoint-friendly), can be resumed.

CLI examples:
  Phase 1:
    python sabio_pipeline.py sabio_substrate_test.csv --phase p1 --sabio_substrate_p1.csv --batch-size 5 --resume
  Phase 2:
    python sabio_pipeline.py sabio_substrate_p1.csv --phase p2 --out sabio_substrate_final.csv --batch-size 5 --resume

Python 3.9+ recommended.
"""

import os
import sys
import time
import json
import re
import argparse
from typing import List, Dict, Iterable, Optional, Set, Any, Tuple

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

# -------------------------- Endpoints --------------------------
SABIO_SYN_URL  = "https://sabiork.h-its.org/sabioRestWebServices/searchCompoundSynonyms"
SABIO_DET_URL  = "https://sabiork.h-its.org/sabioRestWebServices/searchCompoundDetails"

PUG_REST       = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MNX_SPARQL     = "https://rdf.metanetx.org/sparql/"

HEADERS_SPARQL_JSON = {"Accept": "application/sparql-results+json"}
POST_FORM = {"Content-Type": "application/x-www-form-urlencoded"}

# -------------------------- Session w/ retries --------------------------
def make_session() -> requests.Session:
    """Create a requests Session with retries and a helpful UA."""
    s = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    s.headers.update({
        "User-Agent": "sabio-pipeline/2.0 (+https://example.org)",
        "Accept": "*/*",
    })
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.mount("http://",  HTTPAdapter(max_retries=retries))
    return s

SESSION = make_session()

# -------------------------- Utils --------------------------
def chunk_list(items: List, n: int) -> Iterable[List]:
    for i in range(0, len(items), n):
        yield items[i:i+n]

def chunk_by_chars(str_items: List[str], max_chars: int) -> Iterable[List[str]]:
    """Group strings so total characters per group stay under max_chars (for SPARQL VALUES)."""
    bucket, length = [], 0
    for s in str_items:
        add = len(s) + 3
        if length + add > max_chars and bucket:
            yield bucket
            bucket, length = [s], add
        else:
            bucket.append(s)
            length += add
    if bucket:
        yield bucket

def uniq_keep_order(seq: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s); out.append(s)
    return out

def safe_first(seq: Optional[Iterable[str]]) -> str:
    if not seq:
        return ""
    for x in seq:
        if x:
            return x
    return ""

def is_nullish(s: str) -> bool:
    """Treat literal 'null'/'none'/'na' etc. as empty."""
    if s is None:
        return True
    t = str(s).strip().lower()
    return t in {"null", "none", "nan", "na", "n/a", "-"}

def null_to_empty(s: str) -> str:
    return "" if is_nullish(s) else str(s)

def normalize_inchi(s: str) -> str:
    """Strip optional 'InChI=' prefix; keep '1S/...', '1/...', etc."""
    s = null_to_empty(s).replace("\ufeff", "")
    if s.lower().startswith("inchi="):
        s = s[6:]
    return s.strip()

def normalize_inchikey(s: str) -> str:
    """Strip 'InChIKey=' prefix, remove spaces, uppercase."""
    s = null_to_empty(s).replace("\ufeff", "")
    if s.lower().startswith("inchikey="):
        s = s[9:]
    return s.replace(" ", "").upper()

def split_ids(s: str) -> List[str]:
    s = null_to_empty(s)
    if not s:
        return []
    toks = re.split(r"[;,\|/()\[\]{}]+|\s+", s.strip())
    return [t for t in (tok.strip() for tok in toks) if t]

def normalize_chebi_multi(s: str) -> str:
    """Normalize to 'CHEBI:<digits>'; keep multiple, deduplicated, joined by ';'."""
    out: List[str] = []
    for t in split_ids(s):
        last = t.split(":")[-1]
        if re.fullmatch(r"\d+", last):
            out.append(f"CHEBI:{last}")
        else:
            digits = re.sub(r"[^\d]", "", last)
            if digits:
                out.append(f"CHEBI:{digits}")
    return ";".join(uniq_keep_order(out))

def normalize_kegg_cpd_multi(s: str) -> str:
    """Keep only KEGG COMPOUND IDs 'Cxxxxx'; multiple allowed, joined by ';'."""
    out: List[str] = []
    for t in split_ids(s):
        kid = t.split(":")[-1].upper()
        if re.fullmatch(r"C\d{5}", kid):
            out.append(kid)
    return ";".join(uniq_keep_order(out))

def normalize_pubchem(s: str) -> str:
    """Digits only, or empty."""
    s = null_to_empty(s)
    return re.sub(r"\D+", "", s)

def parse_sabio_tsv(text: str) -> List[Dict[str, str]]:
    """
    Parse SABIO responses (text/tsv). Assumes header row is first line with tab separators.
    Convert literal 'null' etc. to empty strings.
    """
    rows: List[Dict[str, str]] = []
    if not text:
        return rows
    lines = [ln for ln in text.splitlines() if ln.strip() != ""]
    if not lines:
        return rows
    header = [h.strip() for h in lines[0].split("\t")]
    for ln in lines[1:]:
        cols = [null_to_empty(c.strip()) for c in ln.split("\t")]
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))
        elif len(cols) > len(header):
            cols = cols[:len(header)]
        rows.append({header[i]: cols[i] for i in range(len(header))})
    return rows

# -------------------------- I/O helpers (checkpoint-resume) --------------------------
OUT_COLS = ["substrate","InChI","InChIKey","SabioCompoundID","ChebiID","PubChemID","KeggCompoundID","MetaNetXID"]

def ensure_outfile_header(path: str):
    """Create file with header if not exists."""
    if not os.path.exists(path):
        pd.DataFrame(columns=OUT_COLS).to_csv(path, index=False)

def load_done_keys(path: str, key_col: str = "substrate") -> Set[str]:
    """Load already processed keys from an existing output file to support resume."""
    if not os.path.exists(path):
        return set()
    try:
        # Use low-memory friendly reading
        done = pd.read_csv(path, usecols=[key_col], dtype=str)
        return set((null_to_empty(x).strip() for x in done[key_col].tolist()))
    except Exception:
        return set()

def append_rows(path: str, df: pd.DataFrame):
    """Append rows to CSV (without rewriting header)."""
    # Write header iff file does not exist or is empty
    write_header = (not os.path.exists(path)) or os.path.getsize(path) == 0
    df.to_csv(path, mode="a", index=False, header=write_header)

# -------------------------- SABIO: name -> SabioCompoundID --------------------------
def sabio_name_to_id(name: str, timeout: int = 45) -> str:
    """
    Query SABIO 'searchCompoundSynonyms' with CompoundName; pick best SabioCompoundID.
    Heuristics:
      - If any row has Name (case-insensitive) exactly equal to input, prefer its ID.
      - Else pick the most frequent SabioCompoundID in the result.
      - Else first row's ID.
    """
    name = (name or "").strip()
    if not name:
        return ""
    params = {"CompoundName": name, "fields[]": ["SabioCompoundID", "Name", "NameType"]}
    try:
        r = SESSION.post(SABIO_SYN_URL, params=params, timeout=timeout)
        r.raise_for_status()
        rows = parse_sabio_tsv(r.text)
        if not rows:
            return ""
        id_counts: Dict[str, int] = {}
        exact_ids: List[str] = []
        for row in rows:
            sid = null_to_empty(row.get("SabioCompoundID"))
            nm  = null_to_empty(row.get("Name"))
            if sid:
                id_counts[sid] = id_counts.get(sid, 0) + 1
                if nm.lower() == name.lower():
                    exact_ids.append(sid)
        if exact_ids:
            return exact_ids[0]
        if id_counts:
            return sorted(id_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        return null_to_empty(rows[0].get("SabioCompoundID"))
    except Exception:
        return ""

def sabio_names_to_ids(names: List[str], sleep: float = 0.08) -> Dict[str, str]:
    """Map each (lowercased) name -> SabioCompoundID (single best)."""
    out: Dict[str, str] = {}
    for nm in names:
        key = (nm or "").strip().lower()
        if not key:
            out[key] = ""
            continue
        sid = sabio_name_to_id(nm)
        out[key] = sid
        time.sleep(sleep)  # gentle throttling
    return out

# -------------------------- SABIO: details by SabioCompoundID --------------------------
def sabio_details_by_ids(ids: List[str], timeout: int = 60, sleep: float = 0.08) -> Dict[str, Dict[str, str]]:
    """
    For each SabioCompoundID, call searchCompoundDetails to get:
    ChebiID, PubChemID, InChI, KeggCompoundID, SabioCompoundID.
    Returns: id -> dict(fields) (normalized; 'null' coerced to empty)
    """
    out: Dict[str, Dict[str, str]] = {}
    uniq_ids = uniq_keep_order([s for s in ids if (s or "").strip()])
    for sid in uniq_ids:
        params = {"SabioCompoundID": sid,
                  "fields[]": ["ChebiID","PubChemID","InChI","SabioCompoundID","KeggCompoundID"]}
        try:
            r = SESSION.post(SABIO_DET_URL, params=params, timeout=timeout)
            r.raise_for_status()
            rows = parse_sabio_tsv(r.text)

            acc = {"ChebiID": [], "PubChemID": [], "InChI": [], "KeggCompoundID": []}
            for row in rows:
                for k in acc:
                    v = null_to_empty(row.get(k, ""))
                    if v:
                        acc[k].append(v)

            raw_chebi = " ; ".join(acc["ChebiID"])
            raw_pubch = " ; ".join(acc["PubChemID"])
            raw_inchi = " ; ".join(acc["InChI"])
            raw_kegg  = " ; ".join(acc["KeggCompoundID"])

            out[sid] = {
                "ChebiID":        normalize_chebi_multi(raw_chebi),
                "PubChemID":      normalize_pubchem(raw_pubch),
                "InChI":          normalize_inchi(raw_inchi),
                "SabioCompoundID": sid,
                "KeggCompoundID": normalize_kegg_cpd_multi(raw_kegg),
            }
        except Exception:
            out[sid] = {"ChebiID":"", "PubChemID":"", "InChI":"", "SabioCompoundID":sid, "KeggCompoundID":""}
        time.sleep(sleep)
    return out

# -------------------------- MetaNetX: by SABIO xref (sid) --------------------------
def mnx_from_sabio_ids(sabio_ids: List[str],
                       batch_max: int = 200,
                       max_chars: int = 22000,
                       sleep: float = 0.1,
                       timeout: int = 60) -> Dict[str, Tuple[str, str, str]]:
    """
    Query MetaNetX by SABIO chemXref:
      IRI: https://identifiers.org/sabiork.compound:<ID>
    Fetch: MNXM id, mnx:inchi, mnx:inchikey.
    Returns: { SabioCompoundID -> (MNXMxxxxx, InChI(without prefix), InChIKey) }
    """
    iris: List[str] = []
    iri2sid: Dict[str, str] = {}
    for sid in uniq_keep_order([s for s in sabio_ids if (s or "").strip()]):
        iri = f"https://identifiers.org/sabiork.compound:{sid}"  # equals to 'sabiork.compound:<sid>'
        iris.append(iri)
        iri2sid[iri] = sid

    out: Dict[str, Tuple[str, str, str]] = {}

    for group in chunk_by_chars(iris, max_chars=max_chars):
        for sub in chunk_list(group, batch_max):
            values = " ".join(f"<{s}>" for s in sub)
            query = f"""
PREFIX mnx: <https://rdf.metanetx.org/schema/>
SELECT ?xref ?met ?inchi ?ik WHERE {{
  VALUES ?xref {{ {values} }}
  ?met a mnx:CHEM ; mnx:chemXref ?xref .
  OPTIONAL {{ ?met mnx:inchi    ?inchi }}
  OPTIONAL {{ ?met mnx:inchikey ?ik    }}
}}
"""
            try:
                r = SESSION.post(
                    MNX_SPARQL,
                    data={"query": query},
                    headers={**HEADERS_SPARQL_JSON, **POST_FORM},
                    timeout=timeout
                )
                r.raise_for_status()
                js = r.json()
                for b in js.get("results", {}).get("bindings", []):
                    xref = b.get("xref", {}).get("value", "")
                    met  = b.get("met",  {}).get("value", "")
                    inch = normalize_inchi(b.get("inchi", {}).get("value", ""))
                    ikey = normalize_inchikey(b.get("ik", {}).get("value", ""))
                    if not xref or not met:
                        continue
                    tail = met.replace("#","/").split("/")[-1]
                    if not tail.startswith("MNXM"):
                        continue
                    sid = iri2sid.get(xref, "")
                    if sid and sid not in out:
                        out[sid] = (tail, inch, ikey)
                time.sleep(sleep)
            except Exception as e:
                print(f"[WARN] MNX (SABIO xref) batch failed: {e}", file=sys.stderr)
    return out

# -------------------------- MetaNetX: by InChI --------------------------
def mnx_from_inchis(inchis: List[str],
                    batch_max: int = 150,
                    max_chars: int = 12000,
                    sleep: float = 0.1,
                    timeout: int = 60) -> Dict[str, Tuple[str, str]]:
    """
    Map InChI strings to (MetaNetXID, InChIKey) via mnx:inchi and mnx:inchikey.
    Submit both '1S/...' and 'InChI=1S/...' forms in VALUES for robust matching.
    Returns: { normalized_inchi_without_prefix -> (MNXMxxxxx, InChIKey) }
    """
    def forms(i: str) -> List[str]:
        i = normalize_inchi(i)
        return [i, f"InChI={i}"] if i else []

    cleaned = uniq_keep_order([normalize_inchi(i) for i in inchis if normalize_inchi(i)])
    out: Dict[str, Tuple[str, str]] = {}

    expanded = []
    for i in cleaned:
        expanded.extend(forms(i))

    for group in chunk_by_chars(expanded, max_chars=max_chars):
        for sub in chunk_list(group, batch_max):
            values = " ".join(json.dumps(s) for s in sub)
            query = f"""
PREFIX mnx: <https://rdf.metanetx.org/schema/>
SELECT ?i ?met ?ik WHERE {{
  VALUES ?i {{ {values} }}
  ?met a mnx:CHEM ; mnx:inchi ?i .
  OPTIONAL {{ ?met mnx:inchikey ?ik }}
}}
"""
            try:
                r = SESSION.post(
                    MNX_SPARQL,
                    data={"query": query},
                    headers={**HEADERS_SPARQL_JSON, **POST_FORM},
                    timeout=timeout
                )
                r.raise_for_status()
                js = r.json()
                for b in js.get("results", {}).get("bindings", []):
                    i_raw = b.get("i", {}).get("value", "")
                    met   = b.get("met", {}).get("value", "")
                    ik    = normalize_inchikey(b.get("ik", {}).get("value", ""))
                    if not i_raw or not met:
                        continue
                    tail = met.replace("#","/").split("/")[-1]
                    if not tail.startswith("MNXM"):
                        continue
                    key = normalize_inchi(i_raw)
                    if key and key not in out:
                        out[key] = (tail, ik)
                time.sleep(sleep)
            except Exception as e:
                print(f"[WARN] MNX (InChI) batch failed: {e}", file=sys.stderr)
    return out

# -------------------------- PubChem properties (fallback for IK/InChI) --------------------------
def pubchem_props_for_cids(cids: List[str], props: List[str] = ["InChI","InChIKey"],
                           batch: int = 50, timeout: int = 45, sleep: float = 0.06) -> Dict[str, Dict[str, str]]:
    """
    Batch-fetch PubChem properties for many CIDs.
    Returns: cid -> {prop: value}
    """
    res: Dict[str, Dict[str, str]] = {}
    cids = uniq_keep_order([normalize_pubchem(c) for c in cids if normalize_pubchem(c)])
    if not cids:
        return res
    prop_str = ",".join(props)
    for group in chunk_list(cids, batch):
        ids = ",".join(group)
        url = f"{PUG_REST}/compound/cid/{ids}/property/{prop_str}/JSON"
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            js = r.json()
            for item in js.get("PropertyTable", {}).get("Properties", []):
                cid = str(item.get("CID",""))
                if not cid:
                    continue
                rec = res.setdefault(cid, {})
                for p in props:
                    val = null_to_empty(item.get(p, ""))
                    if p == "InChI":
                        val = normalize_inchi(val)
                    if p == "InChIKey":
                        val = normalize_inchikey(val)
                    if val and not rec.get(p):
                        rec[p] = val
        except Exception:
            pass
        time.sleep(sleep)
    return res

# -------------------------- Batch workers for phases --------------------------
def process_batch_p1(df_batch: pd.DataFrame) -> pd.DataFrame:
    """
    Phase 1 batch worker:
    - name -> SabioCompoundID
    - details by SabioCompoundID -> ChebiID/PubChemID/InChI/KeggCompoundID (multi kept)
    - InChIKey, MetaNetXID left empty for phase 1
    """
    dfb = df_batch.copy()
    dfb["substrate"] = dfb["substrate"].astype(str).fillna("").map(lambda s: s.strip())
    names = dfb["substrate"].tolist()
    lc = [(n or "").strip().lower() for n in names]

    name2sid = sabio_names_to_ids(names)
    sids = [name2sid.get(k, "") for k in lc]
    dfb["SabioCompoundID"] = sids

    uniq_ids = uniq_keep_order([s for s in sids if s])
    id2det = sabio_details_by_ids(uniq_ids)

    det_inchi, det_pubchem, det_chebi, det_kegg = [], [], [], []
    for sid in sids:
        rec = id2det.get(sid, {}) if sid else {}
        det_inchi.append(rec.get("InChI",""))
        det_pubchem.append(rec.get("PubChemID",""))
        det_chebi.append(rec.get("ChebiID",""))
        det_kegg.append(rec.get("KeggCompoundID",""))

    dfb["InChI"]          = [normalize_inchi(x) for x in det_inchi]
    dfb["PubChemID"]      = [normalize_pubchem(x) for x in det_pubchem]
    dfb["ChebiID"]        = [normalize_chebi_multi(x) for x in det_chebi]
    dfb["KeggCompoundID"] = [normalize_kegg_cpd_multi(x) for x in det_kegg]
    dfb["InChIKey"]       = ""
    dfb["MetaNetXID"]     = ""

    return dfb[OUT_COLS].fillna("")

def process_batch_p2(df_batch: pd.DataFrame) -> pd.DataFrame:
    """
    Phase 2 batch worker (enrichment + fallbacks):
    - If IK missing (or IK & InChI both missing) and have SabioCompoundID: MNX via SABIO xref -> fill MNX/IK/InChI
    - If InChI present: MNX via InChI -> fill MNX/IK
    - If IK still missing but CID exists: PubChem properties -> fill IK/InChI
    """
    dfb = df_batch.copy()
    for c in OUT_COLS:
        if c not in dfb.columns:
            dfb[c] = ""
    dfb = dfb[OUT_COLS].fillna("")
    dfb.reset_index(drop=True, inplace=True)

    # --- MNX via SABIO xref (identifiers.org/sabiork.compound:<ID>) ---
    need_idx_sid = [
        i for i, (ik, inch, sid) in enumerate(zip(dfb["InChIKey"], dfb["InChI"], dfb["SabioCompoundID"]))
        if str(sid).strip() and (not str(ik).strip() or (not str(ik).strip() and not str(inch).strip()))
    ]
    if need_idx_sid:
        need_sids = uniq_keep_order([str(dfb.at[i, "SabioCompoundID"]).strip() for i in need_idx_sid])
        sid2mnx = mnx_from_sabio_ids(need_sids)
        for i in need_idx_sid:
            sid = str(dfb.at[i, "SabioCompoundID"]).strip()
            tup = sid2mnx.get(sid)
            if not tup:
                continue
            mnx_id, inch, ik = tup
            if inch and not str(dfb.at[i, "InChI"]).strip():
                dfb.at[i, "InChI"] = normalize_inchi(inch)
            if ik and not str(dfb.at[i, "InChIKey"]).strip():
                dfb.at[i, "InChIKey"] = normalize_inchikey(ik)
            if mnx_id and not str(dfb.at[i, "MetaNetXID"]).strip():
                dfb.at[i, "MetaNetXID"] = mnx_id

    # --- MNX via InChI (fills MNX & IK) ---
    inchis = [normalize_inchi(i) for i in dfb["InChI"].tolist()]
    mnx_map = mnx_from_inchis(inchis)
    for i, inc in enumerate(inchis):
        if not inc:
            continue
        tup = mnx_map.get(inc)
        if not tup:
            continue
        mnx_id, ik = tup
        if mnx_id and not str(dfb.at[i, "MetaNetXID"]).strip():
            dfb.at[i, "MetaNetXID"] = mnx_id
        if ik and not str(dfb.at[i, "InChIKey"]).strip():
            dfb.at[i, "InChIKey"] = normalize_inchikey(ik)

    # --- PubChem properties fallback (IK/InChI) ---
    need_ik_idx = [i for i, (ik, cid) in enumerate(zip(dfb["InChIKey"], dfb["PubChemID"])) if (not str(ik).strip()) and str(cid).strip()]
    if need_ik_idx:
        need_cids = uniq_keep_order([str(dfb.at[i, "PubChemID"]).strip() for i in need_ik_idx])
        cid_props = pubchem_props_for_cids(need_cids, props=["InChIKey","InChI"])
        for i in need_ik_idx:
            cid = str(dfb.at[i, "PubChemID"]).strip()
            if not cid:
                continue
            ik  = cid_props.get(cid, {}).get("InChIKey", "")
            inc = cid_props.get(cid, {}).get("InChI", "")
            if ik and not str(dfb.at[i, "InChIKey"]).strip():
                dfb.at[i, "InChIKey"] = normalize_inchikey(ik)
            if inc and not str(dfb.at[i, "InChI"]).strip():
                dfb.at[i, "InChI"] = normalize_inchi(inc)

    return dfb[OUT_COLS].fillna("")

# -------------------------- Phase runners (with resume & incremental write) --------------------------
def run_phase_p1(input_csv: str, out_csv: str, batch_size: int, resume: bool):
    df = pd.read_csv(input_csv, dtype=str)
    if "substrate" not in df.columns:
        raise ValueError("Phase p1 requires input CSV with column 'substrate'.")

    # Prepare output (header if needed)
    ensure_outfile_header(out_csv)

    # Resume: skip already processed substrates
    done = load_done_keys(out_csv, key_col="substrate") if resume else set()
    print(f"[P1] Total rows: {len(df)}; Already done (resume): {len(done)}")

    # Build todo list preserving order
    todo_idx = [i for i, s in enumerate(df["substrate"].astype(str).tolist()) if (null_to_empty(s).strip() not in done)]
    total_batches = (len(todo_idx) + batch_size - 1) // batch_size if todo_idx else 0

    for bi, start in enumerate(range(0, len(todo_idx), batch_size), start=1):
        sel = todo_idx[start:start+batch_size]
        df_batch = df.iloc[sel][["substrate"]].copy()
        try:
            out_batch = process_batch_p1(df_batch)
        except Exception as e:
            print(f"[P1][ERROR] Batch {bi}/{total_batches} failed: {e}", file=sys.stderr)
            continue
        append_rows(out_csv, out_batch)
        print(f"[P1][BATCH {bi}/{total_batches}] Wrote {len(out_batch)} rows. Progress: {bi}/{total_batches}")

def run_phase_p2(input_csv: str, out_csv: str, batch_size: int, resume: bool):
    df = pd.read_csv(input_csv, dtype=str)
    # Validate required columns
    missing = [c for c in OUT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Phase p2 requires input CSV with columns: {OUT_COLS}. Missing: {missing}")

    ensure_outfile_header(out_csv)

    done = load_done_keys(out_csv, key_col="substrate") if resume else set()
    print(f"[P2] Total rows: {len(df)}; Already done (resume): {len(done)}")

    todo_idx = [i for i, s in enumerate(df["substrate"].astype(str).tolist()) if (null_to_empty(s).strip() not in done)]
    total_batches = (len(todo_idx) + batch_size - 1) // batch_size if todo_idx else 0

    for bi, start in enumerate(range(0, len(todo_idx), batch_size), start=1):
        sel = todo_idx[start:start+batch_size]
        df_batch = df.iloc[sel].copy()
        try:
            out_batch = process_batch_p2(df_batch)
        except Exception as e:
            print(f"[P2][ERROR] Batch {bi}/{total_batches} failed: {e}", file=sys.stderr)
            continue
        append_rows(out_csv, out_batch)
        print(f"[P2][BATCH {bi}/{total_batches}] Wrote {len(out_batch)} rows. Progress: {bi}/{total_batches}")

# -------------------------- CLI --------------------------
def build_argparser():
    p = argparse.ArgumentParser(description="SABIO substrates → IDs pipeline (batched, resumable, 2 phases)")
    p.add_argument("input_csv", help="Input CSV. p1 requires 'substrate'. p2 requires p1 columns.")
    p.add_argument("--phase", choices=["p1", "p2"], required=True,
                   help="p1: SABIO lookups. p2: MNX+PubChem enrichment.")
    p.add_argument("--out", required=True, help="Output CSV for this phase. Will append in batches.")
    p.add_argument("--batch-size", type=int, default=5, help="Batch size (default: 5).")
    p.add_argument("--resume", action="store_true", help="Skip substrates already present in --out.")
    return p

def main():
    ap = build_argparser()
    args = ap.parse_args()

    if args.phase == "p1":
        run_phase_p1(args.input_csv, args.out, args.batch_size, args.resume)
    elif args.phase == "p2":
        run_phase_p2(args.input_csv, args.out, args.batch_size, args.resume)
    else:
        raise ValueError("Unknown phase")

if __name__ == "__main__":
    main()
