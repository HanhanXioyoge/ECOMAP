# brenda_pipeline.py
# -*- coding: utf-8 -*-
"""
Batched, checkpointed, multi-stage pipeline for mapping BRENDA substrate names to IDs.

Flow (column ownership contract):
  Stage A (BRENDA)  : writes InChI, InChIKey
  Stage B (PubChem) : writes PubChemID
  Stage C (Xrefs)   : writes ChebiID (multi ';'), KeggCompoundID (multi ';')
  Stage D (MetaNetX): writes MetaNetXID  (via InChIKey, then fallback via ChEBI/KEGG xref)

Artifacts layout under --workdir (auto-created):
  workdir/
    manifest.json
    snapshot.csv                    (substrate,row_uid)
    batches/batch_0001.uids         (row_uid list per batch)
    cache/
      name2struct.json              (lower_name -> {InChI, InChIKey})
      ikey2cid.json                 (InChIKey -> CID)
      cid2xrefs.json                (CID -> {"ChEBI":[...], "KEGG":[...]}
      ikey2mnx.json                 (InChIKey -> MNXMxxxxx)
      xref2mnx.json                 (xref IRI -> MNXMxxxxx)
    stageA_name2struct/
      batch_0001.ok.csv             (row_uid,InChI,InChIKey)
      ...
    stageB_ikey2cid/
      batch_0001.ok.csv             (row_uid,PubChemID)
      ...
    stageC_cid2xrefs/
      batch_0001.ok.csv             (row_uid,ChebiID,KeggCompoundID)
      ...
    stageD_struct2mnx/
      batch_0001.ok.csv             (row_uid,MetaNetXID)
      ...
    final/
      final.csv                     (merged snapshot + A+B+C+D)

Example:
  python brenda_pipeline.py run all brenda_substrate.csv --workdir D:\project\ECOMAP\ECOMAP\scripts\database\brenda --batch-size 5 --resume
  python brenda_pipeline.py merge brenda_substrate.csv --workdir D:\project\ECOMAP\ECOMAP\scripts\database\brenda
"""

import os
import io
import re
import sys
import json
import time
import hashlib
import argparse
import shutil
from typing import List, Dict, Iterable, Optional, Set, Any

import requests
import pandas as pd
from requests.adapters import HTTPAdapter, Retry

# ----------------------------- Endpoints -----------------------------
BRENDA_SPARQL = "https://sparql.dsmz.de/api/brenda"
MNX_SPARQL    = "https://rdf.metanetx.org/sparql/"
KEGG_API      = "https://rest.kegg.jp"
PUG_REST      = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW      = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"

HEADERS_SPARQL_JSON = {"Accept": "application/sparql-results+json"}
POST_FORM = {"Content-Type": "application/x-www-form-urlencoded"}

# ----------------------------- Session & Regex -----------------------------
def make_session() -> requests.Session:
    """HTTP session with retries/backoff and a clear UA."""
    s = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    s.headers.update({
        "User-Agent": "brenda-pipeline/2.1 (+local)",
        "Accept": "*/*",
    })
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.mount("http://",  HTTPAdapter(max_retries=retries))
    return s

SESSION = make_session()
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

# ----------------------------- Small utils -----------------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def atomic_write_text(path: str, text: str):
    """Write text atomically by first writing to .tmp then os.replace."""
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    os.replace(tmp, path)

def atomic_write_df_csv(path: str, df: pd.DataFrame):
    """Write CSV atomically to avoid partial files."""
    tmp = path + ".tmp"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)

def read_json(path: str) -> dict:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def write_json_atomic(path: str, obj: dict):
    atomic_write_text(path, json.dumps(obj, ensure_ascii=False, indent=2))

def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def sha1_hex(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()

def chunk_list(items: List, n: int) -> Iterable[List]:
    for i in range(0, len(items), n):
        yield items[i:i+n]

def chunk_by_chars(str_items: List[str], max_chars: int) -> Iterable[List[str]]:
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

def normalize_inchi(s: str) -> str:
    s = (s or "").strip().replace("\ufeff", "")
    if s.lower().startswith("inchi="):
        s = s[6:]
    return s

def normalize_inchikey(s: str) -> str:
    s = (s or "").strip().replace("\ufeff", "")
    if s.lower().startswith("inchikey="):
        s = s[9:]
    return s.replace(" ", "").upper()

def is_valid_inchikey(ikey: str) -> bool:
    ik = normalize_inchikey(ikey)
    return bool(_INCHIKEY_RE.match(ik))

def _walk_strings(obj: Any) -> Iterable[str]:
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _walk_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk_strings(v)

# ----------------------------- Stage A: BRENDA name → InChI/InChIKey -----------------------------
def brenda_names_to_inchi_keys(names: List[str],
                               batch_max_names: int = 40,
                               max_chars: int = 8000,
                               sleep: float = 0.15,
                               timeout: int = 60) -> Dict[str, Dict[str, str]]:
    """
    Query DSMZ BRENDA SPARQL to fetch InChI/InChIKey by exact label/synonym (lowercased).
    Returns: lowercased name -> {'InChI': str, 'InChIKey': str}
    """
    lcnames = [ (s or "").strip().lower() for s in names ]
    uniq_names = uniq_keep_order([s for s in lcnames if s])
    out: Dict[str, Dict[str, str]] = { s: {'InChI': '', 'InChIKey': ''} for s in uniq_names }

    def build_query(sub_names: List[str], include_synonyms: bool) -> str:
        values = " ".join(json.dumps(s) for s in sub_names)
        syn_block = """
  OPTIONAL { ?chem d3o:hasSynonym ?syn . }
  BIND(IF(BOUND(?syn), LCASE(STR(?syn)), "") AS ?syn_lc)
  FILTER(?label_lc = ?q || (?syn_lc != "" && ?syn_lc = ?q))
""" if include_synonyms else "  FILTER(?label_lc = ?q)\n"
        return f"""
PREFIX d3o:  <https://purl.dsmz.de/schema/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT DISTINCT ?q ?inchi ?ikey WHERE {{
  VALUES ?q {{ {values} }}
  ?chem a d3o:ChemicalMaterial ;
        rdfs:label ?label .
  OPTIONAL {{
    ?chem d3o:hasStructure ?s .
    OPTIONAL {{ ?s d3o:hasInChI ?inchi . }}
    OPTIONAL {{ ?s d3o:hasInChIKey ?ikey . }}
  }}
  BIND(LCASE(STR(?label)) AS ?label_lc)
{syn_block}
}}
"""

    # prepare work items
    work: List[List[str]] = []
    for group in chunk_by_chars(uniq_names, max_chars=max_chars):
        work.extend(list(chunk_list(group, batch_max_names)))

    def post_query(sub_names: List[str], include_synonyms: bool,
                   retries: int = 3, base_sleep: float = 0.5) -> Optional[dict]:
        q = build_query(sub_names, include_synonyms)
        for attempt in range(retries + 1):
            try:
                r = SESSION.post(BRENDA_SPARQL, data={"query": q},
                                 headers={**HEADERS_SPARQL_JSON, **POST_FORM},
                                 timeout=timeout)
                r.raise_for_status()
                return r.json()
            except requests.HTTPError as e:
                code = getattr(e.response, "status_code", None)
                # backoff on server/ratelimit; otherwise bubble
                if code and (code == 429 or 500 <= code < 600):
                    time.sleep(base_sleep * (2 ** attempt))
                    continue
                raise
            except Exception:
                time.sleep(base_sleep * (2 ** attempt))
        return None

    idx = 0
    while idx < len(work):
        sub = work[idx]; idx += 1
        data = post_query(sub, include_synonyms=True)
        if data is None:
            data = post_query(sub, include_synonyms=False)
        if data is None:
            if len(sub) > 1:
                mid = len(sub)//2
                work.insert(idx, sub[mid:])
                work.insert(idx, sub[:mid])
                print(f"[INFO] Stage A: split heavy batch ({len(sub)})", file=sys.stderr)
                continue
            else:
                print(f"[WARN] Stage A: failed for name {sub[0]!r}", file=sys.stderr)
                continue
        try:
            for b in data.get("results", {}).get("bindings", []):
                q  = b.get("q", {}).get("value", "").strip().lower()
                i  = b.get("inchi", {}).get("value", "")
                ik = b.get("ikey", {}).get("value", "")
                if q in out:
                    if i and not out[q]['InChI']:     out[q]['InChI']    = i
                    if ik and not out[q]['InChIKey']: out[q]['InChIKey'] = ik
        except Exception as ex:
            print(f"[WARN] Stage A: parse error: {ex}", file=sys.stderr)
        time.sleep(sleep)
    return out

# ----------------------------- Stage B: InChIKey → PubChem CID -----------------------------
def inchikey_to_pubchem_cid(ikey: str, timeout: int = 30) -> str:
    """Return first PubChem CID for an InChIKey; '' if none."""
    ikey = normalize_inchikey(ikey)
    if not ikey:
        return ""
    url = f"{PUG_REST}/compound/inchikey/{requests.utils.quote(ikey)}/cids/JSON"
    try:
        r = SESSION.get(url, timeout=timeout)
        r.raise_for_status()
        js = r.json()
        cids = js.get("IdentifierList", {}).get("CID", [])
        return str(cids[0]) if cids else ""
    except Exception:
        return ""

# ----------------------------- Stage C: CID → (ChEBI[], KEGG[]) -----------------------------
def _extract_sourceids_case(js: dict, accept_names: List[str]) -> list:
    out = []
    info_list = js.get("InformationList", {}).get("Information", [])
    accepts = [a.lower() for a in accept_names]
    for info in info_list:
        src = str(info.get("SourceName", "")).strip().lower()
        if src in accepts:
            vals = info.get("SourceID", []) or info.get("SourceIDList", [])
            if isinstance(vals, list):
                out.extend([str(v) for v in vals if v])
            elif isinstance(vals, str):
                out.append(vals)
    # dedup keep order
    seen = set(); dedup = []
    for v in out:
        if v not in seen:
            seen.add(v); dedup.append(v)
    return dedup

def pubchem_xrefs_sourceid_multi(cid: str, sources: List[str], timeout: int = 30) -> List[str]:
    acc: List[str] = []
    for src in sources:
        url = f"{PUG_REST}/compound/cid/{cid}/xrefs/SourceID/JSON?source={requests.utils.quote(src)}"
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            acc.extend(_extract_sourceids_case(r.json(), [src]))
        except Exception:
            pass
        time.sleep(0.03)
    # de-dup keep order
    return uniq_keep_order(acc)

def pubchem_pugview_xrefs(cid: str, timeout: int = 30) -> Dict[str, List[str]]:
    out = {'KEGG': [], 'ChEBI': []}
    if not cid:
        return out
    url = f"{PUG_VIEW}/{cid}/JSON"
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        js = r.json()
        texts = " | ".join(_walk_strings(js))
        chebi_hits = re.findall(r"CHEBI:\s*\d+", texts, flags=re.IGNORECASE)
        kegg_hits  = re.findall(r"\b[CD]\d{5}\b", texts)
        chebi_norm = ["CHEBI:" + re.sub(r"[^\d]", "", h.split(":")[-1]) for h in chebi_hits]
        kegg_norm  = [h.upper() for h in kegg_hits]
        out['ChEBI'] = uniq_keep_order(chebi_norm)
        out['KEGG']  = uniq_keep_order(kegg_norm)
    except Exception:
        pass
    return out

def collect_kegg_chebi_for_cids(cids: List[str]) -> Dict[str, Dict[str, List[str]]]:
    result: Dict[str, Dict[str, List[str]]] = {cid: {'KEGG': [], 'ChEBI': []} for cid in cids}
    for cid in cids:
        if not cid:
            continue
        chebi = pubchem_xrefs_sourceid_multi(cid, ["ChEBI"])
        kegg  = pubchem_xrefs_sourceid_multi(cid, ["KEGG", "KEGG COMPOUND", "KEGG DRUG"])
        if not chebi or not kegg:
            pv = pubchem_pugview_xrefs(cid)
            if not chebi: chebi = pv.get("ChEBI", [])
            if not kegg:  kegg  = pv.get("KEGG", [])
        result[cid]['ChEBI'] = uniq_keep_order(chebi)
        result[cid]['KEGG']  = uniq_keep_order(kegg)
        time.sleep(0.03)
    return result

def kegg_conv_from_pubchem(cids: List[str], batch: int = 20, timeout: int = 30) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = { cid: [] for cid in cids }
    def _one_batch(group: List[str]):
        if not group: return
        dbentries = "+".join([f"pubchem:{cid}" for cid in group])
        url = f"{KEGG_API}/conv/compound/{dbentries}"
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            for line in r.text.strip().splitlines():
                if not line.strip(): continue
                left, right = line.split("\t")
                if right.startswith("pubchem:"):
                    cid = right.split(":",1)[1]
                    if left.startswith("cpd:"):
                        out.setdefault(cid, []).append(left.split(":",1)[1])
        except requests.HTTPError as e:
            code = getattr(e.response, "status_code", None)
            if code in (413, 414) or (code and 500 <= code < 600):
                if len(group) > 1:
                    mid = len(group)//2
                    _one_batch(group[:mid]); _one_batch(group[mid:])
        except Exception:
            pass
        time.sleep(0.05)
    for group in chunk_list([c for c in cids if c], batch):
        _one_batch(group)
    for cid in list(out.keys()):
        out[cid] = uniq_keep_order(out[cid])
    return out

# ----------------------------- Stage D helpers -----------------------------
def mnx_from_inchikeys(ikeys: List[str],
                       batch_max: int = 120,
                       max_chars: int = 10000,
                       sleep: float = 0.1,
                       timeout: int = 60) -> Dict[str, str]:
    clean = uniq_keep_order([normalize_inchikey(x) for x in ikeys if normalize_inchikey(x)])
    out: Dict[str, str] = {}
    def build_query(sub: List[str]) -> str:
        values = " ".join(json.dumps(s) for s in sub)
        return f"""
PREFIX mnx: <https://rdf.metanetx.org/schema/>
SELECT ?ik ?met WHERE {{
  VALUES ?ik {{ {values} }}
  ?met a mnx:CHEM ; mnx:inchikey ?ik .
}}
"""
    def post_query(sub: List[str], retries: int = 3, base_sleep: float = 0.5) -> Optional[dict]:
        q = build_query(sub)
        for attempt in range(retries + 1):
            try:
                resp = SESSION.post(MNX_SPARQL, data={"query": q},
                                    headers={**HEADERS_SPARQL_JSON, **POST_FORM},
                                    timeout=timeout)
                resp.raise_for_status()
                return resp.json()
            except requests.HTTPError as e:
                code = getattr(e.response, "status_code", None)
                if code in (413, 414) or (code and 500 <= code < 600):
                    time.sleep(base_sleep * (2 ** attempt)); continue
                else:
                    raise
            except Exception:
                time.sleep(base_sleep * (2 ** attempt))
        return None
    work: List[List[str]] = []
    for group in chunk_by_chars(clean, max_chars=max_chars):
        work.extend(list(chunk_list(group, batch_max)))
    idx = 0
    while idx < len(work):
        sub = work[idx]; idx += 1
        data = post_query(sub)
        if data is None:
            if len(sub) > 1:
                mid = len(sub)//2
                work.insert(idx, sub[mid:]); work.insert(idx, sub[:mid])
                print("[INFO] Stage D: split heavy batch", file=sys.stderr)
                continue
            else:
                print(f"[WARN] Stage D: failed for IK {sub[0]}", file=sys.stderr)
                continue
        try:
            for b in data.get("results", {}).get("bindings", []):
                ik  = b.get("ik", {}).get("value", "")
                met = b.get("met", {}).get("value", "")
                tail = met.replace("#","/").split("/")[-1] if met else ""
                if ik and tail.startswith("MNXM"):
                    out[ik.upper()] = tail
        except Exception as ex:
            print(f"[WARN] Stage D parse error: {ex}", file=sys.stderr)
        time.sleep(sleep)
    return out

def _split_semicolon(s: str) -> List[str]:
    """Split a semicolon-delimited string into tokens (trim & drop empties)."""
    toks = [t.strip() for t in str(s or "").split(";")]
    return [t for t in toks if t]

def build_xref_iris_from_ids(chebi_multi: str, kegg_multi: str) -> List[str]:
    """
    Build MetaNetX-accepted chemXref IRIs from multi-valued IDs.
      - ChEBI: identifiers.org + OBO PURL (http/https)
      - KEGG : identifiers.org/kegg.compound (http/https)
    Input format examples:
      chebi_multi = "CHEBI:15377;CHEBI:15378"
      kegg_multi  = "C00001;C00006"
    """
    iris: List[str] = []
    # ChEBI: CHEBI:<digits>
    for raw in _split_semicolon(chebi_multi):
        last = raw.split(":")[-1]
        digits = re.sub(r"[^\d]", "", last)
        if not digits:
            continue
        iris.append(f"https://identifiers.org/chebi:{digits}")
        iris.append(f"http://purl.obolibrary.org/obo/CHEBI_{digits}")
        iris.append(f"https://purl.obolibrary.org/obo/CHEBI_{digits}")
        iris.append(f"http://identifiers.org/chebi:{digits}")  # add http identifiers as well
    # KEGG COMPOUND: Cxxxxx
    for raw in _split_semicolon(kegg_multi):
        kid = raw.split(":")[-1].upper()
        if re.fullmatch(r"C\d{5}", kid):
            iris.append(f"https://identifiers.org/kegg.compound:{kid}")
            iris.append(f"http://identifiers.org/kegg.compound:{kid}")
    return uniq_keep_order(iris)

def mnx_from_xref_iris(iris: List[str],
                       batch_max: int = 200,
                       max_chars: int = 22000,
                       sleep: float = 0.08,
                       timeout: int = 60) -> Dict[str, str]:
    """
    Map many chemXref IRIs -> MNXM via MetaNetX SPARQL.
    Returns: { iri -> MNXMxxxxx }
    """
    clean = uniq_keep_order([u for u in iris if u])
    out: Dict[str, str] = {}
    for group in chunk_by_chars(clean, max_chars=max_chars):
        for sub in chunk_list(group, batch_max):
            values = " ".join(f"<{s}>" for s in sub)
            query = f"""
PREFIX mnx: <https://rdf.metanetx.org/schema/>
SELECT ?xref ?met WHERE {{
  VALUES ?xref {{ {values} }}
  ?met a mnx:CHEM ; mnx:chemXref ?xref .
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
                    x = b.get("xref", {}).get("value", "")
                    m = b.get("met",  {}).get("value", "")
                    tail = m.replace("#","/").split("/")[-1] if m else ""
                    if x and tail.startswith("MNXM"):
                        out[x] = tail
                time.sleep(sleep)
            except Exception as e:
                print(f"[WARN] MNX xref batch failed: {e}", file=sys.stderr)
    return out

# ----------------------------- Snapshot / Batches / Manifest -----------------------------
def build_snapshot(input_csv: str, workdir: str) -> pd.DataFrame:
    """Create snapshot.csv with row_uid if not exists; return snapshot df."""
    ensure_dir(workdir)
    snap_path = os.path.join(workdir, "snapshot.csv")
    if os.path.exists(snap_path):
        return pd.read_csv(snap_path, dtype=str).fillna("")
    # load input
    src = pd.read_csv(input_csv, dtype=str).fillna("")
    if "substrate" not in src.columns:
        raise ValueError("Input CSV must contain a 'substrate' column.")
    src = src[["substrate"]].copy()
    # stable row_uid: sha1(index|substrate)
    uids = []
    for i, row in src.iterrows():
        uids.append(sha1_hex(f"{i}|{row['substrate']}"))
    snap = pd.DataFrame({"row_uid": uids, "substrate": src["substrate"].astype(str).map(lambda s: s.strip())})
    atomic_write_df_csv(snap_path, snap)
    return snap

def build_manifest_and_batches(input_csv: str, workdir: str, batch_size: int):
    """Create manifest.json and batches/*.uids if not present."""
    ensure_dir(workdir)
    manifest_path = os.path.join(workdir, "manifest.json")
    batches_dir   = os.path.join(workdir, "batches")
    ensure_dir(batches_dir)

    snap = build_snapshot(input_csv, workdir)
    if os.path.exists(manifest_path) and len(os.listdir(batches_dir)) > 0:
        return  # already present

    # batch split by order
    n = len(snap)
    total_batches = (n + batch_size - 1) // batch_size
    for bi, start in enumerate(range(0, n, batch_size), start=1):
        end = min(start + batch_size, n)
        uids = snap["row_uid"].iloc[start:end].tolist()
        with open(os.path.join(batches_dir, f"batch_{bi:04d}.uids"), "w", encoding="utf-8") as f:
            for u in uids:
                f.write(u + "\n")

    manifest = {
        "input_csv": os.path.abspath(input_csv),
        "snapshot_sha256": file_sha256(os.path.join(workdir, "snapshot.csv")),
        "batch_size": batch_size,
        "total_rows": n,
        "total_batches": total_batches,
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    write_json_atomic(manifest_path, manifest)

def list_batches(workdir: str) -> List[str]:
    batches_dir = os.path.join(workdir, "batches")
    ensure_dir(batches_dir)
    files = [fn for fn in os.listdir(batches_dir) if fn.endswith(".uids")]
    files.sort()
    return [os.path.join(batches_dir, fn) for fn in files]

# ----------------------------- Cache helpers -----------------------------
def cache_path(workdir: str, name: str) -> str:
    p = os.path.join(workdir, "cache")
    ensure_dir(p)
    return os.path.join(p, name)

def load_cache(workdir: str) -> Dict[str, Dict]:
    return {
        "name2struct": read_json(cache_path(workdir, "name2struct.json")),
        "ikey2cid"   : read_json(cache_path(workdir, "ikey2cid.json")),
        "cid2xrefs"  : read_json(cache_path(workdir, "cid2xrefs.json")),
        "ikey2mnx"   : read_json(cache_path(workdir, "ikey2mnx.json")),
        "xref2mnx"   : read_json(cache_path(workdir, "xref2mnx.json")),   # NEW
    }

def persist_cache(workdir: str, caches: Dict[str, Dict]):
    write_json_atomic(cache_path(workdir, "name2struct.json"), caches.get("name2struct", {}))
    write_json_atomic(cache_path(workdir, "ikey2cid.json"),    caches.get("ikey2cid", {}))
    write_json_atomic(cache_path(workdir, "cid2xrefs.json"),   caches.get("cid2xrefs", {}))
    write_json_atomic(cache_path(workdir, "ikey2mnx.json"),    caches.get("ikey2mnx", {}))
    write_json_atomic(cache_path(workdir, "xref2mnx.json"),    caches.get("xref2mnx", {}))  # NEW

# ----------------------------- Stage runners (per-batch) -----------------------------
def run_stage_A_batch(workdir: str, batch_uid_file: str, resume: bool = False):
    out_dir = os.path.join(workdir, "stageA_name2struct")
    ensure_dir(out_dir)
    batch_id = os.path.splitext(os.path.basename(batch_uid_file))[0].replace(".uids","")
    ok_path  = os.path.join(out_dir, f"{batch_id}.ok.csv")
    if resume and os.path.exists(ok_path):
        return  # already done

    snap = pd.read_csv(os.path.join(workdir, "snapshot.csv"), dtype=str).fillna("")
    uids = [ln.strip() for ln in open(batch_uid_file,"r",encoding="utf-8").read().splitlines() if ln.strip()]
    df = snap[snap["row_uid"].isin(uids)].copy()
    names = df["substrate"].astype(str).map(lambda s: s.strip()).tolist()

    caches = load_cache(workdir)
    name2struct = caches["name2struct"]

    # separate names: cached vs need-fetch
    need = []
    for nm in names:
        key = (nm or "").strip().lower()
        if not key or key in name2struct:
            continue
        need.append(nm)

    # query BRENDA for missing
    if need:
        res = brenda_names_to_inchi_keys(need, batch_max_names=40, max_chars=8000, sleep=0.15, timeout=60)
        # merge into cache
        for nm in need:
            k = (nm or "").strip().lower()
            rec = res.get(k, {'InChI':'','InChIKey':''})
            name2struct[k] = {"InChI": rec.get("InChI",""), "InChIKey": rec.get("InChIKey","")}

    # build output
    InChI = []
    InChIKey = []
    for nm in names:
        k = (nm or "").strip().lower()
        rec = name2struct.get(k, {})
        InChI.append(normalize_inchi(rec.get("InChI","")))
        InChIKey.append(normalize_inchikey(rec.get("InChIKey","")))

    out_df = pd.DataFrame({"row_uid": uids, "InChI": InChI, "InChIKey": InChIKey})
    atomic_write_df_csv(ok_path, out_df)
    caches["name2struct"] = name2struct
    persist_cache(workdir, caches)

def run_stage_B_batch(workdir: str, batch_uid_file: str, resume: bool = False):
    out_dir = os.path.join(workdir, "stageB_ikey2cid")
    ensure_dir(out_dir)
    batch_id = os.path.splitext(os.path.basename(batch_uid_file))[0].replace(".uids","")
    ok_path  = os.path.join(out_dir, f"{batch_id}.ok.csv")
    if resume and os.path.exists(ok_path):
        return

    # dependency: stage A batch output
    stageA = pd.read_csv(os.path.join(workdir, "stageA_name2struct", f"{batch_id}.ok.csv"), dtype=str).fillna("")
    uids = stageA["row_uid"].tolist()
    iks  = stageA["InChIKey"].tolist()

    caches = load_cache(workdir)
    ikey2cid = caches["ikey2cid"]

    cids = []
    for ik in iks:
        key = normalize_inchikey(ik)
        if not key:
            cids.append("")
            continue
        if key in ikey2cid:
            cids.append(ikey2cid[key])
            continue
        cid = inchikey_to_pubchem_cid(key)
        ikey2cid[key] = cid or ""
        cids.append(ikey2cid[key])
        time.sleep(0.02)

    out_df = pd.DataFrame({"row_uid": uids, "PubChemID": cids})
    atomic_write_df_csv(ok_path, out_df)
    caches["ikey2cid"] = ikey2cid
    persist_cache(workdir, caches)

def run_stage_C_batch(workdir: str, batch_uid_file: str, resume: bool = False):
    out_dir = os.path.join(workdir, "stageC_cid2xrefs")
    ensure_dir(out_dir)
    batch_id = os.path.splitext(os.path.basename(batch_uid_file))[0].replace(".uids","")
    ok_path  = os.path.join(out_dir, f"{batch_id}.ok.csv")
    if resume and os.path.exists(ok_path):
        return

    stageB = pd.read_csv(os.path.join(workdir, "stageB_ikey2cid", f"{batch_id}.ok.csv"), dtype=str).fillna("")
    uids = stageB["row_uid"].tolist()
    cids = [str(x).strip() for x in stageB["PubChemID"].tolist()]

    caches = load_cache(workdir)
    cid2xrefs = caches["cid2xrefs"]

    uniq = uniq_keep_order([c for c in cids if c])
    need = [c for c in uniq if c not in cid2xrefs]

    if need:
        # xrefs via PUG REST & fallback; then KEGG conv fallback if KEGG list empty
        xrefs = collect_kegg_chebi_for_cids(need)
        miss = [cid for cid in need if not xrefs.get(cid,{}).get("KEGG")]
        if miss:
            conv = kegg_conv_from_pubchem(miss, batch=20)
            for cid, lst in conv.items():
                had = xrefs.setdefault(cid, {}).get("KEGG", [])
                xrefs[cid]["KEGG"] = uniq_keep_order((had or []) + (lst or []))
        for cid in need:
            cid2xrefs[cid] = xrefs.get(cid, {"KEGG":[], "ChEBI":[]})

    # normalize multi-values and join with ';'
    kegg_out, chebi_out = [], []
    for cid in cids:
        info = cid2xrefs.get(cid, {}) if cid else {}
        # ChEBI
        cheb_raw = (info.get("ChEBI", []) or [])
        seen_c = set(); cheb_norm = []
        for v in cheb_raw:
            tail = str(v).split(":")[-1]
            digits = re.sub(r"[^\d]", "", tail)
            if digits:
                val = f"CHEBI:{digits}"
                if val not in seen_c:
                    seen_c.add(val); cheb_norm.append(val)
        chebi_out.append(";".join(cheb_norm) if cheb_norm else "")
        # KEGG (compound only)
        kegg_raw = (info.get("KEGG", []) or [])
        seen_k = set(); kegg_norm = []
        for v in kegg_raw:
            kid = str(v).split(":")[-1].upper()
            if re.fullmatch(r"C\d{5}", kid) and kid not in seen_k:
                seen_k.add(kid); kegg_norm.append(kid)
        kegg_out.append(";".join(kegg_norm) if kegg_norm else "")

    out_df = pd.DataFrame({"row_uid": uids,
                           "ChebiID": chebi_out,
                           "KeggCompoundID": kegg_out})
    atomic_write_df_csv(ok_path, out_df)
    caches["cid2xrefs"] = cid2xrefs
    persist_cache(workdir, caches)

def run_stage_D_batch(workdir: str, batch_uid_file: str, resume: bool = False):
    """
    Stage D:
      1) InChIKey → MetaNetX (primary)
      2) Fallback: for rows still missing MNX, use ChebiID / KeggCompoundID (from Stage C)
         to build xref IRIs and query MetaNetX via mnx:chemXref.
    """
    out_dir = os.path.join(workdir, "stageD_struct2mnx")
    ensure_dir(out_dir)
    batch_id = os.path.splitext(os.path.basename(batch_uid_file))[0].replace(".uids","")
    ok_path  = os.path.join(out_dir, f"{batch_id}.ok.csv")
    if resume and os.path.exists(ok_path):
        return

    # Stage A (structures)
    stageA = pd.read_csv(os.path.join(workdir, "stageA_name2struct", f"{batch_id}.ok.csv"), dtype=str).fillna("")
    uids = stageA["row_uid"].tolist()
    iks  = stageA["InChIKey"].tolist()

    # Stage C (xref ids) - optional (if not ready, fallback won't run but IK step will)
    stageC_path = os.path.join(workdir, "stageC_cid2xrefs", f"{batch_id}.ok.csv")
    stageC = None
    if os.path.exists(stageC_path):
        stageC = pd.read_csv(stageC_path, dtype=str).fillna("")
        # align by row_uid
        stageC = stageC.set_index("row_uid").reindex(uids).reset_index()
    else:
        print(f"[INFO] Stage D: Stage C output not found for {batch_id}, xref fallback will be skipped.", file=sys.stderr)

    caches = load_cache(workdir)
    ikey2mnx = caches["ikey2mnx"]
    xref2mnx = caches["xref2mnx"]

    # -------- 1) IK → MNX (primary) --------
    need_iks = uniq_keep_order([normalize_inchikey(x) for x in iks if normalize_inchikey(x) and normalize_inchikey(x) not in ikey2mnx])
    if need_iks:
        map_new = mnx_from_inchikeys(need_iks)
        for ik, mnx in map_new.items():
            ikey2mnx[ik] = mnx or ""

    # current MNX assignment from IK
    mnx_out = [ikey2mnx.get(normalize_inchikey(x), "") for x in iks]
    filled_by_ik = sum(1 for v in mnx_out if v)
    print(f"[D] IK→MNX filled: {filled_by_ik}/{len(uids)}", file=sys.stderr)

    # -------- 2) Fallback via ChEBI/KEGG xref --------
    if stageC is not None:
        # find rows still missing MNX and having xrefs
        miss_idx = [i for i, v in enumerate(mnx_out) if not v]
        chebi_list = stageC["ChebiID"].tolist() if "ChebiID" in stageC.columns else [""] * len(uids)
        kegg_list  = stageC["KeggCompoundID"].tolist() if "KeggCompoundID" in stageC.columns else [""] * len(uids)

        # collect IRIs (preserve per-row order: ChEBI first, then KEGG)
        per_row_iris: List[List[str]] = []
        all_iris: List[str] = []
        for i in range(len(uids)):
            iris = build_xref_iris_from_ids(str(chebi_list[i]), str(kegg_list[i]))
            per_row_iris.append(iris)
            all_iris.extend(iris)
        all_iris = uniq_keep_order(all_iris)

        # use cache first
        need_iris = [iri for iri in all_iris if iri not in xref2mnx]
        if need_iris:
            iri2mnx_new = mnx_from_xref_iris(need_iris)
            for iri in need_iris:
                xref2mnx[iri] = iri2mnx_new.get(iri, "")

        # assign for missing rows: prefer first mapped IRI (order: ChEBI -> KEGG)
        filled_by_xref = 0
        for i in miss_idx:
            for iri in per_row_iris[i]:
                mnx = xref2mnx.get(iri, "")
                if mnx:
                    mnx_out[i] = mnx
                    filled_by_xref += 1
                    break
        print(f"[D] XREF→MNX filled: {filled_by_xref}/{len(miss_idx)} (from Chebi/KEGG)", file=sys.stderr)
    else:
        print(f"[D] XREF→MNX skipped (Stage C missing).", file=sys.stderr)

    out_df = pd.DataFrame({"row_uid": uids, "MetaNetXID": mnx_out})
    atomic_write_df_csv(ok_path, out_df)
    caches["ikey2mnx"] = ikey2mnx
    caches["xref2mnx"] = xref2mnx
    persist_cache(workdir, caches)

# ----------------------------- Merge final -----------------------------
def merge_final(workdir: str, output_csv: str):
    snap = pd.read_csv(os.path.join(workdir, "snapshot.csv"), dtype=str).fillna("")
    # stage A
    a_dir = os.path.join(workdir, "stageA_name2struct")
    a_parts = [os.path.join(a_dir, fn) for fn in os.listdir(a_dir) if fn.endswith(".ok.csv")]
    A = pd.concat([pd.read_csv(p, dtype=str).fillna("") for p in a_parts], ignore_index=True) if a_parts else pd.DataFrame(columns=["row_uid","InChI","InChIKey"])
    # stage B
    b_dir = os.path.join(workdir, "stageB_ikey2cid")
    b_parts = [os.path.join(b_dir, fn) for fn in os.listdir(b_dir) if fn.endswith(".ok.csv")]
    B = pd.concat([pd.read_csv(p, dtype=str).fillna("") for p in b_parts], ignore_index=True) if b_parts else pd.DataFrame(columns=["row_uid","PubChemID"])
    # stage C
    c_dir = os.path.join(workdir, "stageC_cid2xrefs")
    c_parts = [os.path.join(c_dir, fn) for fn in os.listdir(c_dir) if fn.endswith(".ok.csv")]
    C = pd.concat([pd.read_csv(p, dtype=str).fillna("") for p in c_parts], ignore_index=True) if c_parts else pd.DataFrame(columns=["row_uid","ChebiID","KeggCompoundID"])
    # stage D
    d_dir = os.path.join(workdir, "stageD_struct2mnx")
    d_parts = [os.path.join(d_dir, fn) for fn in os.listdir(d_dir) if fn.endswith(".ok.csv")]
    D = pd.concat([pd.read_csv(p, dtype=str).fillna("") for p in d_parts], ignore_index=True) if d_parts else pd.DataFrame(columns=["row_uid","MetaNetXID"])

    # Left-join by row_uid in A→B→C→D order (column ownership respected)
    out = snap.copy()
    for df_add in (A, B, C, D):
        if not df_add.empty:
            out = out.merge(df_add, on="row_uid", how="left")

    # Clean NaNs -> empty; enforce column order
    for c in ["InChI","InChIKey","ChebiID","PubChemID","KeggCompoundID","MetaNetXID"]:
        if c not in out.columns: out[c] = ""
    out = out[["substrate","InChI","InChIKey","ChebiID","PubChemID","KeggCompoundID","MetaNetXID"]].fillna("")

    ensure_dir(os.path.join(workdir, "final"))
    atomic_write_df_csv(output_csv, out)

# ----------------------------- Orchestrator CLI -----------------------------
def cmd_run(args):
    input_csv = args.input_csv
    workdir   = args.workdir
    batch_sz  = args.batch_size
    stage     = args.stage.upper()
    resume    = args.resume

    # init snapshot + batches if needed
    build_manifest_and_batches(input_csv, workdir, batch_sz)
    manifest = read_json(os.path.join(workdir, "manifest.json"))
    total_batches = manifest["total_batches"]

    batch_files = list_batches(workdir)

    # stage selector
    stage_map = {
        "A": run_stage_A_batch,
        "B": run_stage_B_batch,
        "C": run_stage_C_batch,
        "D": run_stage_D_batch,
    }

    stages = ["A","B","C","D"] if stage == "ALL" else [stage]
    for st in stages:
        print(f"[RUN] Stage {st} | batches={total_batches} | resume={resume} | workdir={workdir}")
        for i, bf in enumerate(batch_files, start=1):
            # skip batch if dependency not satisfied (for B/C/D)
            if st != "A":
                # require stage A for B/D; stage B for C
                if st in ("B","D"):
                    req = os.path.join(workdir, "stageA_name2struct", os.path.basename(bf).replace(".uids",".ok.csv"))
                    if not os.path.exists(req):
                        print(f"[SKIP] Stage {st} batch {i}/{total_batches} waiting for Stage A: {os.path.basename(req)}")
                        continue
                if st in ("C",):
                    req = os.path.join(workdir, "stageB_ikey2cid", os.path.basename(bf).replace(".uids",".ok.csv"))
                    if not os.path.exists(req):
                        print(f"[SKIP] Stage C batch {i}/{total_batches} waiting for Stage B: {os.path.basename(req)}")
                        continue
            try:
                stage_map[st](workdir, bf, resume=resume)
                print(f"[BATCH {i}/{total_batches}] Stage {st} OK")
            except Exception as ex:
                print(f"[ERROR] Stage {st} batch {i}/{total_batches} failed: {ex}", file=sys.stderr)
                # do NOT stop the whole run; continue to next batch
                continue

def cmd_merge(args):
    workdir    = args.workdir
    output_csv = args.output_csv or os.path.join(workdir, "final", "final.csv")
    ensure_dir(os.path.dirname(output_csv))
    merge_final(workdir, output_csv)
    print(f"[DONE] Merged → {output_csv}")

def main():
    ap = argparse.ArgumentParser(description="BRENDA substrates → InChI/InChIKey/PubChem/KEGG/ChEBI/MNX (batched & checkpointed)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # run
    ap_run = sub.add_parser("run", help="Run one stage or all stages with batch checkpoints")
    ap_run.add_argument("stage", choices=["A","B","C","D","all"], help="Which stage to run")
    ap_run.add_argument("input_csv", help="Input CSV with column 'substrate'")
    ap_run.add_argument("--workdir", required=True, help="Working directory to store artifacts")
    ap_run.add_argument("--batch-size", type=int, default=5, help="Rows per batch (default: 5)")
    ap_run.add_argument("--resume", action="store_true", help="Skip already-finished batches")
    ap_run.set_defaults(func=cmd_run)

    # merge
    ap_merge = sub.add_parser("merge", help="Merge snapshot + stages into final CSV")
    ap_merge.add_argument("input_csv", help="Input CSV (for snapshot integrity check / layout)")
    ap_merge.add_argument("--workdir", required=True, help="Working directory used in run")
    ap_merge.add_argument("--output-csv", help="Final merged CSV path (default: workdir/final/final.csv)")
    ap_merge.set_defaults(func=cmd_merge)

    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
