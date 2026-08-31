# MATLAB–Python Bridge Contract

## Purpose

This document defines the boundary between the ECOMAP Python web server and
MATLAB bridge functions. It makes responses machine-checkable and gives later
bridge batches one shared failure vocabulary.

The terms **MUST**, **MUST NOT**, **SHOULD**, and **MAY** are normative. This
contract applies to every MATLAB function exposed to Python through the MATLAB
Engine; internal MATLAB-only functions are outside its scope.

## Function shape

MATLAB bridge functions MUST:

1. use a name of the form `mdpXxx`;
2. live in `scripts/web/matlab/bridge/`;
3. accept only values supported by the MATLAB Engine boundary;
4. build the four-field response envelope below; and
5. return that envelope as a JSON-encoded character vector.

A typical function has this shape:

```matlab
function jsonStr = mdpExample(inputPath)
%MDPEXAMPLE Perform one bridge operation.
    bridge_log('Example', 'starting');
    value = performOperation(inputPath);
    payload = make_ok(value);
    jsonStr = jsonencode(payload);
end
```

Python invokes the function with `nargout=1`, converts the returned MATLAB
character value to `str`, and decodes it with `json.loads`.

Bridge functions SHOULD catch expected domain failures and translate them to a
canonical error code. Unexpected programming errors MAY propagate as MATLAB
Engine errors because they represent defects rather than normal responses.

## Return struct

Before JSON encoding, every response MUST be a MATLAB `struct` with exactly
these four fields:

| Field | MATLAB form | Meaning |
| --- | --- | --- |
| `ok` | logical scalar | `true` on success; `false` on every represented error |
| `error_code` | char vector | `''` on success; one canonical code on failure |
| `error_message` | char vector | `''` on success; one short sentence on failure |
| `result` | any JSON-serialisable value | payload on success; ignored on failure |

Success envelope:

```matlab
payload = struct('ok', true, ...
                 'error_code', '', ...
                 'error_message', '', ...
                 'result', result);
```

Failure envelope:

```matlab
payload = struct('ok', false, ...
                 'error_code', 'err_param_invalid', ...
                 'error_message', 'The model path is required.', ...
                 'result', []);
```

Failure responses SHOULD use `[]` for `result`. Callers MUST ignore that field
when `ok` is `false`. The four fields retain their names and types even when the
result is empty. A bridge MUST NOT return a bare payload, a bare error string,
or an envelope with additional fields.

## JSON serialization

The Engine-visible return value MUST be the JSON string produced by
`jsonencode(payload)`. Python decodes and validates the envelope before
exposing the operation result.

The `result` value MUST be JSON-serialisable. Supported shapes include:

- primitive logical, numeric, and character values;
- cells containing primitives;
- structs whose fields contain primitives;
- one-dimensional arrays; and
- nested struct arrays composed of the supported values above.

Payloads MUST NOT contain function handles, `containers.Map` values, arbitrary
user objects, or other values that `jsonencode` cannot represent reliably.
For tabular output, use an explicit JSON-safe shape such as `columns` plus a
cell array of `rows`.

## Error codes

The following catalogue is canonical and exhaustive for contract version 1.

| Code | Meaning |
| --- | --- |
| `err_init_fail` | Project or state initialisation failed. |
| `err_param_invalid` | A required parameter is missing or has the wrong type or format. |
| `err_model_format` | A model file is unreadable or has the wrong schema. |
| `err_no_biomass` | No biomass reaction can be resolved from the input. |
| `err_no_target` | No target reaction can be resolved from the input. |
| `err_docker_missing` | A required Docker image or service is unavailable. |
| `err_no_proteomics` | The proteomics dataset is empty or unreadable. |
| `err_gurobi_license` | The Gurobi licence check failed; LP/MIP is not solvable. |
| `err_raven_notfound` | RAVEN Toolbox is not on the path when required. |
| `err_sluice_data` | One or more SLUICE inputs are missing or malformed. |
| `err_kcat_merge` | Merging kcat predictions failed. |
| `err_presto_data` | PRESTO input or calibration data is missing or malformed. |
| `err_oom` | Execution ran out of memory; the caller should reduce scope. |
| `err_cancelled` | The user cancelled the run and clean state was preserved. |

For an error not named separately, reuse the closest canonical code. For
example, malformed SBML maps to `err_model_format`; do not invent a convenience
code.

If a genuinely new category is unavoidable, update this document and Python's
`_BRIDGE_ERROR_CODES` in the same commit before returning it. The MATLAB helper
also rejects unknown codes, so both sides enforce the catalogue.

## Logging

Every bridge function MUST write one informational line to standard output.
The line MUST begin with `[mdp<Name>]`, where `<Name>` is the function suffix,
for example: `[mdpRunFseof] completed in 2.4 s`.

Use `bridge_log('RunFseof', 'completed in %.1f s', elapsed)` for the standard
prefix. Include timing when helpful. Do not include secrets, full datasets, or
multiline stack traces.

## Paths

Python MUST pass filesystem paths as absolute `char` values, never relative
paths. Bridge functions MAY normalise separators or canonicalise a path, but
MUST NOT resolve one against MATLAB's current working directory.

An unreadable or schema-invalid model maps to `err_model_format`; other invalid
path parameters map to `err_param_invalid` or the nearest domain-specific code.

## IDs

When a caller needs a stable `projectId`, `modelId`, or `jobId`, return it inside
`result` as:

```matlab
result = struct('id', callerProvidedId);
```

A bridge MUST preserve a caller-provided ID and MUST NOT invent an ID when the
caller did not provide one. IDs are payload data and do not replace envelope
fields.

## Cancellation

Bridge functions SHOULD honour `app.cancellation.Token`-style cancellation at
safe checkpoints during long work. A token-style `false` abort MUST return
`ok=false`, `error_code='err_cancelled'`, and one short `error_message`.

Cancellation MUST preserve clean state and MUST NOT expose partial output as a
successful result. Token plumbing is deferred; these response semantics apply
now.

## Versioning

This document defines bridge contract version 1. Backward-compatible additions
belong inside `result`; the four envelope fields and their types remain fixed.

Changes to envelope fields, JSON boundary behaviour, or the error catalogue
MUST update MATLAB helpers, Python validation, and unit tests together.
Producers and consumers SHOULD land those changes atomically.

## mdpBuildOkoIntervalsFromHomologs

Build OKO+ kcat intervals from UniProt cross-species enzyme retrieval +
UNIKP/CatPred multi-organism prediction. Outputs CSVs in the legacy
`ecoli_kcat_preds.csv` format (6 columns: rxn, uniprot, min, max, mean,
mode) compatible with OKO+ solver.

### Inputs

| Param         | Type        | Description                                   |
|---------------|-------------|-----------------------------------------------|
| ec_model_id   | str         | ecModel identifier (e.g. 'eciML1515')         |
| predictors    | list[str]   | subset of {'UniKP', 'CatPred'}                |
| manager_path  | str         | optional ParameterManager.m path              |

### Output envelope

```python
{
    'ok': bool,
    'error_code': str,
    'error_message': str,
    'result': {
        'predictor_csv_paths': {
            'UniKP': str,    # absolute path
            'CatPred': str,  # absolute path
        },
        'n_candidates_per_enzyme': [
            {'rxn': str, 'nHomologs': int},
            ...
        ],
        'elapsed_seconds': float,
    }
}
```

### Errors

| Code               | When                                       |
|--------------------|--------------------------------------------|
| `err_param_invalid`| predictors not in {'UniKP', 'CatPred'}     |
| `err_model_format` | Cannot resolve model_id                    |
| `err_docker_missing`| Predictor Docker unavailable               |
| `err_internal`     | Any other error                            |
