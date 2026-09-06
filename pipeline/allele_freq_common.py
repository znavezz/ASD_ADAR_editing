"""
Shared allele-frequency helpers — the single source of truth for turning VEP
`Extra` columns into rows of `public.allele_frequencies` (one row per variant that
carries any population frequency).

Used by both:
  * phase2_db.py            — populates the table during a normal pipeline run

Keeping the column map, the value coercion, the variant lookup and the upsert
mutation here guarantees the pipeline and the backfill write byte-identical rows.

Frequencies reach VEP via `--everything` (gnomAD `--af_gnomade`/`--af_gnomadg` +
1000G + `--max_af`) — see phase1_preprocess.py. This table stores the gnomAD
per-population AFs (`gnomade_*` exomes, `gnomadg_*` genomes) plus VEP's `MAX_AF`
(cross-database max over 1000G/ESP/gnomAD — a VEP-computed field, NOT gnomAD-provided)
and `MAX_AF_POPS`. On the GRCh37 cache the gnomAD side is exomes r2.1 only; genome
and MID/REMAINING columns come through NULL until a run provides them.
"""

# gnomAD allele-frequency + MAX_AF columns -> VEP "Extra" key (exact case VEP emits).
FREQ_FLOAT_COLUMNS = {
    # gnomAD exomes
    "gnomade_af": "gnomADe_AF",
    "gnomade_afr_af": "gnomADe_AFR_AF",
    "gnomade_amr_af": "gnomADe_AMR_AF",
    "gnomade_asj_af": "gnomADe_ASJ_AF",
    "gnomade_eas_af": "gnomADe_EAS_AF",
    "gnomade_fin_af": "gnomADe_FIN_AF",
    "gnomade_mid_af": "gnomADe_MID_AF",
    "gnomade_nfe_af": "gnomADe_NFE_AF",
    "gnomade_remaining_af": "gnomADe_REMAINING_AF",
    "gnomade_sas_af": "gnomADe_SAS_AF",
    # gnomAD genomes (genomes carry an extra AMI = Amish population)
    "gnomadg_af": "gnomADg_AF",
    "gnomadg_afr_af": "gnomADg_AFR_AF",
    "gnomadg_ami_af": "gnomADg_AMI_AF",
    "gnomadg_amr_af": "gnomADg_AMR_AF",
    "gnomadg_asj_af": "gnomADg_ASJ_AF",
    "gnomadg_eas_af": "gnomADg_EAS_AF",
    "gnomadg_fin_af": "gnomADg_FIN_AF",
    "gnomadg_mid_af": "gnomADg_MID_AF",
    "gnomadg_nfe_af": "gnomADg_NFE_AF",
    "gnomadg_remaining_af": "gnomADg_REMAINING_AF",
    "gnomadg_sas_af": "gnomADg_SAS_AF",
    # VEP-computed cross-database maximum (1000G/ESP/gnomAD)
    "max_af": "MAX_AF",
}

# text column(s) -> VEP key
FREQ_TEXT_COLUMNS = {
    "max_af_pops": "MAX_AF_POPS",
}

# every non-key column, in table order — also the upsert `update_columns` set.
FREQ_VALUE_COLUMNS = list(FREQ_FLOAT_COLUMNS) + list(FREQ_TEXT_COLUMNS)

_NULL_TOKENS = {"", "-", ".", "nan", "none", "null", "na"}


def _to_float(x):
    if x is None:
        return None
    if isinstance(x, float) and x != x:  # NaN
        return None
    s = str(x).strip()
    if s.lower() in _NULL_TOKENS:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _to_text(x):
    if x is None:
        return None
    if isinstance(x, float) and x != x:
        return None
    s = str(x).strip()
    return None if s.lower() in _NULL_TOKENS else s


def extract_frequencies(row):
    """Map one VEP record (dict or pandas Series keyed by Extra key) to a dict of
    allele_frequencies table columns. Returns None only if EVERY value is absent, so
    no all-NULL row is inserted; any variant with any frequency (incl. a 1000G-only
    MAX_AF) is kept."""
    get = row.get  # dict and pandas Series both support .get
    out = {}
    for col, key in FREQ_FLOAT_COLUMNS.items():
        out[col] = _to_float(get(key))
    for col, key in FREQ_TEXT_COLUMNS.items():
        out[col] = _to_text(get(key))
    return out if any(v is not None for v in out.values()) else None


def variant_key_from_uploaded(uploaded_variation):
    """VEP `Uploaded_variation` is `chr:pos:ref:alt` (phase1 sets the VCF ID to
    exactly this). Returns (chr, pos, ref, alt) or None if it can't be parsed."""
    parts = str(uploaded_variation).split(":")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    try:
        pos = int(pos)
    except ValueError:
        return None
    return (str(chrom), pos, ref, alt)


_VARIANT_LOOKUP_QUERY = """
query VariantLookup($limit: Int!, $offset: Int!) {
  variants(order_by: {id: asc}, limit: $limit, offset: $offset) {
    id
    ref
    alt
    coordinate { chr start }
  }
}
"""


def build_variant_lookup(gql, page=5000, log=None):
    """Return {(chr, start, ref, alt): variant_id} for all variants (hg19)."""
    lut, offset = {}, 0
    while True:
        rows = gql(_VARIANT_LOOKUP_QUERY, {"limit": page, "offset": offset})["variants"]
        if not rows:
            break
        for v in rows:
            co = v.get("coordinate") or {}
            lut[(str(co.get("chr")), co.get("start"), v.get("ref"), v.get("alt"))] = v["id"]
        if log:
            log(f"  variant lookup: {len(lut)} loaded...")
        if len(rows) < page:
            break
        offset += page
    return lut


_UPSERT = (
    """
mutation UpsertAlleleFrequencies($objects: [allele_frequencies_insert_input!]!) {
  insert_allele_frequencies(
    objects: $objects
    on_conflict: {constraint: allele_frequencies_pkey, update_columns: [__COLS__]}
  ) {
    affected_rows
  }
}
""".replace("__COLS__", ", ".join(FREQ_VALUE_COLUMNS))
)


def build_objects(records, lookup, log=None):
    """Turn VEP records (list of dict/Series) into allele_frequencies insert objects,
    keyed by variant_id via `lookup`. Returns (objects, stats)."""
    objects, seen = [], set()
    n_no_match = n_no_freq = 0
    for rec in records:
        key = variant_key_from_uploaded(rec.get("Uploaded_variation"))
        if key is None:
            continue
        vid = lookup.get(key)
        if vid is None:
            n_no_match += 1
            continue
        if vid in seen:
            continue
        vals = extract_frequencies(rec)
        if vals is None:
            n_no_freq += 1
            continue
        seen.add(vid)
        objects.append({"variant_id": vid, **vals})
    stats = {"matched": len(objects), "no_variant_match": n_no_match,
             "no_frequency_data": n_no_freq}
    if log:
        log(f"  allele_frequencies objects: {stats}")
    return objects, stats


def upsert(gql, objects, batch=500, log=None):
    """Batch-upsert allele_frequencies rows; returns total affected_rows."""
    total = 0
    for i in range(0, len(objects), batch):
        chunk = objects[i:i + batch]
        res = gql(_UPSERT, {"objects": chunk})
        total += res["insert_allele_frequencies"]["affected_rows"]
        if log:
            log(f"  upserted {total}/{len(objects)}...")
    return total
