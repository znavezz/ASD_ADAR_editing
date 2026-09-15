"""Stage 2: load the annotated variants into PostgreSQL and attach every external score.

The longest stage - about seven and a half hours for the published run - because each
variant acquires annotation from several sources, and CADD arrives over the network in
batches. It creates the schema from `db/schema.sql`, bulk-loads the staging tables, then
populates the normalised model:

  * **Variants and their effects** - consequence terms collapsed to the project's region
    names, codons, CDS and protein positions, NMD status, variant class.
  * **Scores** - SIFT and PolyPhen from the VEP annotation; CADD from the Ensembl REST
    API (`?CADD=1`), which is why this stage needs network access; AlphaMissense where
    available.
  * **Allele frequencies** - gnomAD exome per-population values plus VEP's MAX_AF, via
    `allele_freq_common.py`. The GRCh37 cache carries no genome frequencies, so every
    `gnomadg_*` column is NULL by construction.
  * **Gene-level annotation** - SFARI gene scores, and GTEx v10 brain expression keyed
    by unversioned ENSG.

Of the 525,687 staged rows, 171 fail preprocessing and are excluded, leaving the
published 329,279 unique variants.

Run through `scripts/run_pipeline.py`. It writes, so the target database is an explicit
argument; see `--env-file`.
"""

# --- imports this stage needs -------------------------------------------------
# Previously inherited from the runner's scope via exec(). A stage that cannot say
# what it imports cannot be read, tested, or run on its own, so each now imports
# for itself; only genuine shared state still arrives from the runner.
import os
import re
import json
import time
import logging
import subprocess
# The class, not the module: this stage calls datetime.utcnow(). It used to inherit
# the name from the runner's scope via exec(), and acquired `import datetime` when the
# stages were given explicit imports - which made every call an AttributeError.
from datetime import datetime
import requests
import pandas as pd
from typing import Optional, Tuple
from joblib import dump, load

df_varicarta = hp.load_vcf(VARICARTA_PATH)
df_varicarta_expanded = hp.parse_info_column(df_varicarta)

df_varicarta_expanded = df_varicarta_expanded.rename(columns={"id": "varicarta_id"})
df_varicarta_expanded = df_varicarta_expanded.rename(columns={"#CHROM": "CHROM"})
df_varicarta_expanded = df_varicarta_expanded.rename(columns={"ID": "variant_id"})
# Strip whitespace from REF/ALT to match preprocess.ipynb
df_varicarta_expanded["REF"] = df_varicarta_expanded["REF"].str.strip()
df_varicarta_expanded["ALT"] = df_varicarta_expanded["ALT"].str.strip()

# Filter to only variants that passed preprocessing (standard alleles, REF!=ALT, ref-validated)
df_clean_keys = pd.read_csv(path_varicarta_standard_positive, sep='\t', dtype={"#CHROM": str})[["#CHROM", "POS", "REF", "ALT"]].rename(columns={"#CHROM": "CHROM"})
df_clean_keys = df_clean_keys.drop_duplicates()

# Ensure CHROM types match (both as string) to avoid int vs str mismatch on merge
df_varicarta_expanded["CHROM"] = df_varicarta_expanded["CHROM"].astype(str)
df_clean_keys["CHROM"] = df_clean_keys["CHROM"].astype(str)

before = len(df_varicarta_expanded)
df_varicarta_expanded = df_varicarta_expanded.merge(
    df_clean_keys, on=["CHROM", "POS", "REF", "ALT"], how="inner"
)
after = len(df_varicarta_expanded)
logging.info(f"Filtered staging variants: {before} -> {after} (removed {before - after} that didn't pass preprocess)")

unique_variants = df_varicarta_expanded.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
logging.info(f"Unique variants: {len(unique_variants)}")


# Check: are there variants with same (CHROM, POS, REF, ALT) but different stop_hg19?
dupes = df_varicarta_expanded.groupby(["CHROM", "POS", "REF", "ALT"])["stop_hg19"].nunique()
conflicting = dupes[dupes > 1]
logging.info(f"Variants with conflicting stop_hg19: {len(conflicting)}")
if len(conflicting) > 0:
    logging.info(conflicting.head(20))
    # Show example rows
    example_key = conflicting.index[0]
    logging.info(f"\nExample: {example_key}")
    mask = (df_varicarta_expanded["CHROM"] == example_key[0]) & \
           (df_varicarta_expanded["POS"] == example_key[1]) & \
           (df_varicarta_expanded["REF"] == example_key[2]) & \
           (df_varicarta_expanded["ALT"] == example_key[3])
    logging.info(df_varicarta_expanded.loc[mask, ["CHROM", "POS", "REF", "ALT", "stop_hg19", "variant_id"]].drop_duplicates())



df_varicarta_annotated = hp.load_vcf(os.path.join(vep_output_dir, vep_output_file_name + ".txt"), "Uploaded_variation")

df_varicarta_annotated_expanded = hp.parse_info_column(df_varicarta_annotated, "Extra")
df_varicarta_annotated_expanded = df_varicarta_annotated_expanded.rename(columns={"#Uploaded_variation": "Uploaded_variation"})

# df_varicarta = your DataFrame

df_varicarta_expanded = hp.normalize_missing(
    df_varicarta_expanded,
    null_tokens={'.'},   # Varicarta missing marker
    strip_whitespace=True
)

# Drop the hg19 validation column before saving
df_varicarta_expanded.to_csv(os.path.join(DB_DIR, "staging_varicarta_variants_raw.csv"), index=False, na_rep="")

# df_vep = your DataFrame

df_varicarta_annotated_expanded = hp.normalize_missing(
    df_varicarta_annotated_expanded,
    null_tokens={'-'},   # VEP missing marker
    strip_whitespace=True
)

df_varicarta_annotated_expanded.to_csv(
    os.path.join(DB_DIR, "staging_vep_annotations_raw.csv"),
    index=False,
    na_rep=""
)




# ── Start Docker Compose ──
project_dir = str(PROJECT_ROOT)

env = os.environ.copy()
env["USER_ID"] = str(os.getuid())
env["GROUP_ID"] = str(os.getgid())

subprocess.run(
    ["docker", "compose", "up", "-d"],
    cwd=project_dir, env=env, check=True
)

logging.info("Waiting for container to start...")
time.sleep(5)

# ── Wait for Postgres & load schema ──
logging.info(f"Waiting for Postgres to be ready (DB: {POSTGRES_DB})...")
for _ in range(30):
    result = subprocess.run(
        ["docker", "exec", POSTGRES_DOCKER_CONTAINER, "pg_isready", "-U", POSTGRES_USER, "-d", POSTGRES_DB],
        capture_output=True
    )
    if result.returncode == 0:
        break
    logging.info(f"Waiting for database '{POSTGRES_DB}' to be ready...")
    time.sleep(2)
else:
    raise RuntimeError(f"Postgres not accepting connections after 60s. Check docker logs.")

# pg_isready succeeds before init scripts finish — wait for the DB to actually exist
logging.info(f"Waiting for database '{POSTGRES_DB}' to be created...")
for _ in range(15):
    result = subprocess.run(
        ["docker", "exec", POSTGRES_DOCKER_CONTAINER, "psql", "-U", POSTGRES_USER, "-d", POSTGRES_DB, "-c", "SELECT 1"],
        capture_output=True
    )
    if result.returncode == 0:
        break
    time.sleep(2)
else:
    raise RuntimeError(f"Database '{POSTGRES_DB}' not created after 30s. Check docker logs.")

logging.info("Postgres is ready. Loading schema...")
# The schema is tracked source, identical for every environment, so it is read from
# the repository rather than from DB_DIR. DB_DIR names where a particular database
# keeps its data, and pointing it at a scratch stack must not hide the DDL.
schema_path = PROJECT_ROOT / "db" / "schema.sql"
with open(schema_path, "r") as f:
    result = subprocess.run(
        ["docker", "exec", "-i", POSTGRES_DOCKER_CONTAINER, "psql", "-U", POSTGRES_USER, "-d", POSTGRES_DB],
        stdin=f, capture_output=True, text=True
    )
    if result.stderr:
        logging.warning(f"Schema load warnings (may be harmless on re-run):\n{result.stderr.strip()}")

logging.info("✓ Schema loaded successfully")



# Import Hasura metadata (run right after applying schema above)


# Tracked source, like the schema beside it. Reading it from DB_DIR meant a run
# pointed at a scratch stack picked up whatever stale copy happened to sit in that
# directory - which is how the bystander_codons relationship went missing and took
# three published numbers with it.
metadata_path = PROJECT_ROOT / "db" / "hasura_metadata.json"
with open(metadata_path) as f:
    metadata = json.load(f)

metadata_url = HASURA_URL.replace("/v1/graphql", "/v1/metadata")
payload = json.dumps({"type": "replace_metadata", "version": 2, "args": metadata})

# Hasura may take a few seconds to start — retry up to 30s
for attempt in range(15):
    result = subprocess.run(
        ["curl", "-s", "-X", "POST", metadata_url,
         "-H", f"x-hasura-admin-secret: {HASURA_ADMIN_SECRET}",
         "-H", "Content-Type: application/json",
         "-d", payload],
        capture_output=True, text=True
    )
    if result.stdout.strip():
        break
    logging.info(f"Waiting for Hasura to be ready... (attempt {attempt + 1})")
    time.sleep(2)

if not result.stdout.strip():
    raise RuntimeError("Hasura did not respond after retries. Check HASURA_URL and container status.")

response = json.loads(result.stdout)
if "message" in response and response["message"] == "success":
    logging.info("Metadata imported successfully")
else:
    logging.warning(f"Hasura metadata response: {json.dumps(response, indent=2)[:500]}")


# ── Load CSV data into Postgres ──
logging.info("Loading staging_varicarta_variants_raw.csv into Postgres...")

copy_varicarta_sql = r"""\copy public.staging_varicarta_variants_raw(chrom, pos, varicarta_vcf_id, ref, alt, qual, filter, info, cadd13_phred, cadd13_raw, cadd_phred, cadd_raw, fathmm_pred, fathmm_score, gerp_rs, lrt_pred, lrt_score, lr_pred, lr_score, mutationassessor_pred, mutationassessor_score, mutationtaster_pred, mutationtaster_score, polyphen2_hdiv_pred, polyphen2_hdiv_score, polyphen2_hvar_pred, polyphen2_hvar_score, radialsvm_pred, radialsvm_score, sift_pred, sift_score, siphy_29way_logodds, vest3_score, aa_change, category, clinvar_20150629, code_change, cytoband, exac03, func, gene_detail, gene_symbol, varicarta_id, inheritance, paper_id, paper_key, phylop100way_vertebrate, phylop46way_placental, pid, protein_change, sample_id, sequencing_study_type, stop_hg19, subject_id, validation, validation_method, validation_reported) FROM '/csv/staging_varicarta_variants_raw.csv' WITH (FORMAT csv, HEADER true);"""

result = subprocess.run(
    ["docker", "exec", "-i", POSTGRES_DOCKER_CONTAINER, "psql", "-v", "ON_ERROR_STOP=1", "-U", POSTGRES_USER, "-d", POSTGRES_DB],
    input=copy_varicarta_sql, capture_output=True, text=True
)
logging.info(f"[psql stdout] {result.stdout}")
if result.stderr:
    logging.info(f"[psql stderr] {result.stderr}")
if result.returncode != 0:
    raise RuntimeError(f"Varicarta CSV import failed with exit code {result.returncode}")
logging.info("✓ staging_varicarta_variants_raw loaded")


logging.info("Loading staging_vep_annotations_raw.csv into Postgres...")

copy_vep_sql = r"""\copy public.staging_vep_annotations_raw(uploaded_variation, location, allele, gene, feature, feature_type, consequence, cdna_position, cds_position, protein_position, amino_acids, codons, existing_variation, extra, impact, strand, variant_class, symbol, symbol_source, hgnc_id, biotype, canonical, exon, hgvsc, gnomade_af, gnomade_afr_af, gnomade_amr_af, gnomade_asj_af, gnomade_eas_af, gnomade_fin_af, gnomade_nfe_af, gnomade_sas_af, max_af, max_af_pops, phastcons100, intron, distance, hgvs_offset, minimised, af, afr_af, amr_af, eas_af, eur_af, sas_af, ensp, uniparc, ccds, swissprot, trembl, gene_pheno, pli_gene_value, domains, hgvsp, clin_sig, pheno, somatic, sift, polyphen, blosum62, am_class, am_genome, am_pathogenicity, am_protein_variant, am_transcript_id, am_uniprot_id, nmd, loftool, pubmed, motif_name, motif_pos, high_inf_pos, motif_score_change, transcription_factors, flags, mirna) FROM '/csv/staging_vep_annotations_raw.csv' WITH (FORMAT csv, HEADER true);"""

result = subprocess.run(
    ["docker", "exec", "-i", POSTGRES_DOCKER_CONTAINER, "psql", "-v", "ON_ERROR_STOP=1", "-U", POSTGRES_USER, "-d", POSTGRES_DB],
    input=copy_vep_sql, capture_output=True, text=True
)
logging.info(f"[psql stdout] {result.stdout}")
if result.stderr:
    logging.info(f"[psql stderr] {result.stderr}")
if result.returncode != 0:
    raise RuntimeError(f"VEP CSV import failed with exit code {result.returncode}")
logging.info("✓ staging_vep_annotations_raw loaded")




UPSERT_ORG_WITH_BUILDS = """
mutation UpsertOrgWithBuild($obj: organisms_insert_input!) {
  insert_organisms(
    objects: [$obj]
    on_conflict: {
      constraint: organisms_taxon_id_key
      update_columns: [abbrev]
    }
  ) {
    returning {
      id
      taxon_id
      scientific_name
      common_name
      abbrev
      genome_builds(where: {build_name: {_eq: "hg19"}}) {
        id
        build_name
        assembly_name
        assembly_accession
        release_date
      }
    }
  }
}
"""

build_obj = {
  "build_name": "hg19",
  "assembly_name": "GRCh37",
  "assembly_accession": "GCF_000001405",
  "release_date": "2009-02-27",
  "notes": "Varicarta staging is hg19; fasta GRCh37 primary assembly",
}

vars1 = {
  "obj": {
    "taxon_id": 9606,
    "scientific_name": "Homo sapiens",
    "common_name": "Human",
    "abbrev": "H. sapiens",
    "genome_builds": {
      "data": [build_obj],
      "on_conflict": {
        "constraint": "genome_builds_organism_id_build_name_key",
        "update_columns": ["assembly_name", "assembly_accession", "release_date", "notes"],
      },
    },
  }
}


org = gql(UPSERT_ORG_WITH_BUILDS, vars1)["insert_organisms"]["returning"][0]
organism_id = org["id"]

# because we filtered by hg19 in returning:
genome_build_id = org["genome_builds"][0]["id"]
logging.info(f"organism_id: {organism_id}")
logging.info(f"genome_build_id: {genome_build_id}")

SET_DEFAULT_BUILD = """
mutation SetDefaultBuild($organism_id: Int!, $genome_build_id: Int!) {
  update_organisms_by_pk(
    pk_columns: {id: $organism_id}
    _set: {default_genome_build_id: $genome_build_id}
  ) {
    id
    default_genome_build_id
  }
}
"""

res2 = gql(SET_DEFAULT_BUILD, {
  "organism_id": organism_id,
  "genome_build_id": genome_build_id
})
logging.info(f"default_genome_build_id: {res2['update_organisms_by_pk']['default_genome_build_id']}")


UPSERT_SOURCE = """
mutation UpsertSource($obj: sources_insert_input!) {
  insert_sources(
    objects: [$obj],
    on_conflict: {
    constraint: sources_name_release_key, 
    update_columns: [name]
    }
  ) {
    returning { id name url description release}
  }
}
"""

# 1) source: Ensembl
ensembl = gql(UPSERT_SOURCE, {"obj": {"name": "Ensembl", "url": "https://www.ensembl.ensembl/index.html", "description": "Ensembl is a public and open project providing access to genomes, annotations, tools and methods. Its goal is to enable genomic science by providing high-quality, integrated and consistent annotation on all cellular genomes within a harmonious, scalable and accessible infrastructure."}})["insert_sources"]["returning"][0]
ensembl_id = ensembl["id"]

# 1) source: HGNC
hgnc = gql(UPSERT_SOURCE, {"obj": {"name": "HGNC", "url": "https://www.genenames.org/", "description": "The HGNC is responsible for approving unique symbols and names for human loci, including protein coding genes, ncRNA genes and pseudogenes, to allow unambiguous scientific communication."}})["insert_sources"]["returning"][0]
hgnc_id = hgnc["id"]

varicarta = gql(UPSERT_SOURCE, {"obj": {"name": "Varicarta", "url": "https://varicarta.msl.ubc.ca/index", "description": "The VariCarta web application and database is an initiative aimed at collecting, reconciling and consistently cataloguing literature-derived genomic variants found in ASD subjects. We put effort into careful curation of the data, standardizing variant reporting, adding comprehensive annotation and identifying overlaps. We hope that the scientific community will find VariCarta to be a useful resource for ASD genomic variants published in peer-reviewed literature."}})["insert_sources"]["returning"][0]
varicarta_id = varicarta["id"]

gtex = gql(UPSERT_SOURCE, {"obj": {"name": "GTEx", "url": "https://www.gtexportal.org/home/", "description": "The Genotype-Tissue Expression (GTEx) Portal is a comprehensive public resource for researchers studying tissue and cell-specific gene expression and regulation across individuals, development, and species, with data from 3 NIH projects."}})["insert_sources"]["returning"][0]
gtex_id = gtex["id"]

uniparc = gql(UPSERT_SOURCE, {"obj": {"name": "UniParc", "url": "https://www.uniprot.org/uniparc/", "description": "The UniProt Archive (UniParc) is a comprehensive and non-redundant database of protein sequences. These sequences are sourced from public sequence databases, and each unique sequence is stored in a UniParc entry with a stable unique identifier (UPI). A UPI is never removed, changed or reassigned to a different sequence. In addition to the protein sequence, a UniParc entry contains cross-references to all source database entries in which the sequence exists or existed, with a date range that shows when the sequence was first and last seen in each source entry. In this way UniParc tracks sequence changes in the source databases and archives the history of all changes."}})["insert_sources"]["returning"][0]
uniparc_id = uniparc["id"]

ccds = gql(UPSERT_SOURCE, {"obj": {"name": "CCDS", "url": "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi", "description": "The Consensus CDS (CCDS) project is a collaborative effort to identify a core set of human and mouse protein coding regions that are consistently annotated and of high quality. The long term goal is to support convergence towards a standard set of gene annotations."}})["insert_sources"]["returning"][0]
ccds_id = ccds["id"]

swissprot = gql(UPSERT_SOURCE, {"obj": {"name": "SwissProt", "url": "https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=*", "description": "Swiss-Prot (from UniProtKB) is a manually annotated and reviewed protein sequence database, which strives to provide a high level of annotation, a minimal level of redundancy and high level of integration with other databases."}})["insert_sources"]["returning"][0]
swissprot_id = swissprot["id"]
swissprot_name = swissprot["name"]

trembl = gql(UPSERT_SOURCE, {"obj": {"name": "TrEMBL", "url": "https://www.uniprot.org/uniprotkb?query=*&facets=reviewed%3Afalse", "description": "TrEMBL (from UniProtKB) is a computer-annotated supplement to the manually annotated and reviewed UniProtKB/Swiss-Prot database that contains all the translations of EMBL nucleotide sequence entries not yet integrated into Swiss-Prot."}})["insert_sources"]["returning"][0]
trembl_id = trembl["id"]

vep = gql(UPSERT_SOURCE, {"obj": {"name": "VEP", "url": "https://www.ensembl.org/info/docs/tools/vep/index.html", "description": "Ensembl VEP predicts the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on gene transcripts and protein sequence, as well as regulatory regions. It reports reference data including gene and variant phenotype associations and population allele frequencies to facilitate variant prioritisation and interpretation."}})["insert_sources"]["returning"][0]
vep_id = vep["id"]

cath_Gene3D_id = gql(UPSERT_SOURCE, {"obj": {"name": "CATH-Gene3D", "url": "https://www.cathdb.info/", "description": "Gene3D is a database of protein domain structure annotations for protein sequences. It uses a library of profile hidden Markov models (HMMs) derived from CATH superfamilies to predict structural domain boundaries and assignments for millions of protein sequences from Ensembl, UniProt, and RefSeq."}})["insert_sources"]["returning"][0]["id"]

HAMAP_id = gql(UPSERT_SOURCE, {"obj": {"name": "HAMAP", "url": "https://hamap.expasy.org/", "description": "HAMAP (High-quality Automated and Manual Annotation of Proteins) is a system for the classification of protein sequences into families and subfamilies and the subsequent annotation of these sequences. HAMAP uses manually curated family profiles and annotation rules to provide high-quality functional annotation to protein sequences from UniProtKB."}})["insert_sources"]["returning"][0]["id"]

panther_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "PANTHER", 
    "url": "https://www.pantherdb.org", 
    "description": "The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a comprehensive system that classifies genes and proteins into families, subfamilies, molecular functions, biological processes and pathways."
}})["insert_sources"]["returning"][0]["id"]

pirsf_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "PIRSF", 
    "url": "https://proteininformationresource.org", 
    "description": "The PIRSF classification system provides a hierarchical classification of whole proteins to reflect their evolutionary relationships of full-length proteins and their domain architecture."
}})["insert_sources"]["returning"][0]["id"]

prosite_patterns_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "PROSITE_patterns", 
    "url": "https://prosite.expasy.org/", 
    "description": "PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns to identify them."
}})["insert_sources"]["returning"][0]["id"]

prosite_profiles_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "PROSITE_profiles", 
    "url": "https://prosite.expasy.org/", 
    "description": "PROSITE profiles are weight matrices (profiles) used for sensitive identification of protein domains and families, complementing the pattern-based approach."
}})["insert_sources"]["returning"][0]["id"]

pfam_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "Pfam", 
    "url": "http://pfam.xfam.org/", 
    "description": "Pfam is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)."
}})["insert_sources"]["returning"][0]["id"]

prints_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "Prints", 
    "url": "https://www.uniprot.org/database/DB-0082", 
    "description": "PRINTS is a compendium of protein fingerprints. A fingerprint is a group of conserved motifs used to characterise a protein family or domain."
}})["insert_sources"]["returning"][0]["id"]

smart_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "SMART", 
    "url": "http://smart.embl-heidelberg.de", 
    "description": "SMART (a Simple Modular Architecture Research Tool) allows the identification and annotation of genetically mobile domains and the analysis of domain architectures. More than 1400 domain families found in signalling, extracellular and chromatin-associated proteins are detectable. These domains are extensively annotated with respect to phyletic distributions, functional class, tertiary structures and functionally important residues. Each domain found in a non-redundant protein database as well as search parameters and taxonomic information are stored in a relational database system. Web user interface to this database allow searches for proteins containing specific combinations of domains in defined taxa. For all the details, please refer to the publications on SMART, listed below."
}})["insert_sources"]["returning"][0]["id"]

superfamily_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "Superfamily", 
    "url": "https://www.ebi.ac.uk/interpro/entry/ssf/#table", 
    "description": "SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes based on a hidden Markov model library representing all proteins of known structure (SCOP)."
}})["insert_sources"]["returning"][0]["id"]

tigrfam_id = gql(UPSERT_SOURCE, {"obj": {
    "name": "TIGRFAM", 
    "url": "https://www.ncbi.nlm.nih.gov/Structure/cdd/docs/tigrfams.html", 
    "description": "TIGRFAMs is a collection of protein families, featuring curated multiple sequence alignments, hidden Markov models (HMMs) and annotation, primarily for prokaryotic proteins."
}})["insert_sources"]["returning"][0]["id"]


# --- Not a Source ---
# Low_complexity_(Seg) is an algorithm output, not a database/source.


logging.info(f"ensembl_id: {ensembl_id}")
logging.info(f"hgnc_id: {hgnc_id}")
logging.info(f"varicarta_id: {varicarta_id}")
logging.info(f"gtex_id: {gtex_id}")
logging.info(f"uniparc_id: {uniparc_id}")
logging.info(f"ccds_id: {ccds_id}")
logging.info(f"swissprot_id: {swissprot_id}")
logging.info(f"swissprot_name: {swissprot_name}")
logging.info(f"trembl_id: {trembl_id}")
logging.info(f"vep_id: {vep_id}")
logging.info(f"cath_Gene3D_id: {cath_Gene3D_id}")
logging.info(f"HAMAP_id: {HAMAP_id}")
logging.info(f"panther_id: {panther_id}")
logging.info(f"pirsf_id: {pirsf_id}")
logging.info(f"prosite_patterns_id: {prosite_patterns_id}")
logging.info(f"prosite_profiles_id: {prosite_profiles_id}")
logging.info(f"pfam_id: {pfam_id}")
logging.info(f"prints_id: {prints_id}")
logging.info(f"smart_id: {smart_id}")
logging.info(f"superfamily_id: {superfamily_id}")
logging.info(f"tigrfam_id: {tigrfam_id}")




FETCH_STAGING_BATCH = """
query FetchStaging($last_id: bigint!, $limit: Int!) {
  staging_varicarta_variants_raw(
    where: {id: {_gt: $last_id}},
    order_by: {id: asc},
    limit: $limit
  ) {
    id
    chrom
    pos
    stop_hg19
    ref
    alt
    varicarta_vcf_id
    varicarta_id
    inheritance
    subject_id
    paper_key
    paper_id
  }
}
"""
MUTATION_INGEST_BATCH = """
mutation IngestBatch($objects: [variants_insert_input!]!) {
  insert_variants(
    objects: $objects
    on_conflict: {
      constraint: variants_coordinate_id_ref_alt_key
      update_columns: [ref_validated]
    }
  ) {
    affected_rows
    returning { id }
  }
}
"""


logging.info(f"Start insering db from staging varicarta table")

genome_build_id = 1  # hg19
READ_BATCH_SIZE = 5000      # rows fetched from staging per GraphQL query
WRITE_BATCH_SIZE = 500     # how many variants per insert mutation (<= READ_BATCH_SIZE)

last_id = 0
total_ok = 0
errors = []

while True:
    # 1) fetch a batch from staging
    batch_data = gql(FETCH_STAGING_BATCH, {
        "last_id": last_id,
        "limit": READ_BATCH_SIZE,
    })
    rows = batch_data["staging_varicarta_variants_raw"]
    if not rows:
        break

    # process rows in chunks of WRITE_BATCH_SIZE
    for i in range(0, len(rows), WRITE_BATCH_SIZE):
        chunk = rows[i : i + WRITE_BATCH_SIZE]
        variant_objects = []
        staging_ids = []

        for row in chunk:
            staging_ids.append(row["id"])
            last_id = row["id"]
            alt = row["alt"]
            if '/' in alt:
                l,r = alt.split('/')
                if l == row["ref"]:
                    alt = r
                elif r == row["ref"]:
                    alt = l

        

            # ---------- subject insert (reusable) ----------
            subject_insert = None
            if row["subject_id"] is not None:
                subject_insert = {
                    "subject_key": str(row["subject_id"]),
                }

            # ---------- variants_subjects ----------
            vs_data = []
            if subject_insert is not None:
                vs_data.append({
                    "subject": {
                        "data": subject_insert,
                        "on_conflict": {
                            "constraint": "subjects_subject_key_key",
                            "update_columns": ["subject_key"],  # no-op
                        },
                    }
                })

            # ---------- variants_papers + papers_subjects ----------
            vp_data = []
            if row["paper_key"] is not None:
                paper_data = {
                    "paper_key": str(row["paper_key"]),
                    "external_paper_id": str(row["paper_id"]) if row["paper_id"] is not None else None,
                }

                # link paper<->subject through join table
                if subject_insert is not None:
                    paper_data["papers_subjects"] = {
                        "data": [
                            {
                                "subject": {
                                    "data": subject_insert,
                                    "on_conflict": {
                                        "constraint": "subjects_subject_key_key",
                                        "update_columns": ["subject_key"],  # no-op
                                    },
                                }
                            }
                        ],
                        "on_conflict": {
                            "constraint": "papers_subjects_pkey",   # (paper_id, subject_id)
                            "update_columns": ["paper_id"],         # no-op
                        },
                    }

                vp_data.append({
                    "paper": {
                        "data": paper_data,
                        "on_conflict": {
                            "constraint": "papers_paper_key_key",
                            "update_columns": ["external_paper_id"],  # will fill later if missing
                        },
                    }
                })


            # ---------- main variant object ----------
            variant_obj = {
                "ref": row["ref"],
                "alt": alt,
                "ref_validated": True,
                "varicarta_vcf_id": int(row["varicarta_vcf_id"]),
                "varicarta_internal_id": row["varicarta_id"],
                "inheritance": row["inheritance"],

                "coordinate": {
                    "data": {
                        "genome_build_id": genome_build_id,
                        "chr": str(row["chrom"]),
                        "start": int(row["pos"]),
                        "end": (
                            int(row["stop_hg19"])
                            if row["stop_hg19"] is not None
                            else int(row["pos"]) + len(row["ref"]) - 1
                        ),
                        "strand": "+"
                    },
                    "on_conflict": {
                        "constraint": "coordinates_genome_build_id_chr_start_end_strand_key",
                        "update_columns": ["strand"],
                    },
                },
                "variants_sources": {
                    "data": [
                        {
                            "source_id": varicarta_id
                        }
                    ],
                    "on_conflict": {
                        "constraint": "variants_sources_pkey",
                        "update_columns": ["variant_id"],
                    },
                },
            }
            if vs_data:
                variant_obj["variants_subjects"] = {
                    "data": vs_data,
                    "on_conflict": {
                        # name from your error message
                        "constraint": "variants_subjects_pkey",
                        # any column is fine as long as it exists; this is a no-op
                        "update_columns": ["variant_id"],
                    },
                }

            if vp_data:
                variant_obj["variants_papers"] = {
                    "data": vp_data,
                    "on_conflict": {
                        # name from your error message
                        "constraint": "variants_papers_pkey",
                        # same trick – no-op update
                        "update_columns": ["variant_id"],
                    },
                }


            variant_objects.append(variant_obj)

        # 2) write this chunk in one mutation
        try:
            result = gql(MUTATION_INGEST_BATCH, {"objects": variant_objects})
            inserted = result["insert_variants"]["returning"]
            total_ok += len(inserted)
        except Exception as e:
            # log batch error (you can also fall back to row-by-row here if you like)
            errors.append({
                "staging_ids": staging_ids,
                "error": str(e),
            })

    logging.info(f"Processed up to staging id={last_id}, ok={total_ok}, errors={len(errors)}")

logging.info(f"DONE. Total OK: {total_ok}, Errors: {len(errors)}")
logging.info(errors)


GET_DOMAINS = """
query GetDomains {
  staging_vep_annotations_raw(where: { domains: { _is_null: false } }) {
    domains
  }
}

"""
rows = gql(GET_DOMAINS, {})["staging_vep_annotations_raw"]

domain_sources = set()

for r in rows:
    domains = r["domains"]
    if not domains:
        continue

    for item in domains.split(","):
        domain_sources.add(item.split(":", 1)[0])

logging.info(sorted(domain_sources))


expression_data_path = os.environ["EXPRESSION_DATA_PATH"]
sample_attributes_path = os.environ["SAMPLE_ATTRIBUTES_PATH"]
subject_attributes_path = os.environ["SUBJECT_ATTRIBUTES_PATH"]    
# Strategy: Read in chunks with lower memory usage
try:
    expression_data_df = pd.read_csv(
        expression_data_path,
        sep='\t',
        skiprows=2,
        compression='gzip',
        low_memory=False,
        dtype={'Description': 'string'},  # Optimize string columns
        engine='c'  # Faster C engine (default, but explicit)
    )
    logging.info(f"✓ Loaded expression data: {expression_data_df.shape[0]:,} genes × {expression_data_df.shape[1]:,} samples")
except Exception as e:
    logging.info(f"⚠️  Fast load failed: {e}")
    logging.info("→ Retrying with python engine (slower but more robust)...")
    expression_data_df = pd.read_csv(
        expression_data_path,
        sep='\t',
        skiprows=2,
        compression='gzip',
        engine='python',
        on_bad_lines='skip'
    )

expression_data_df.set_index("Name", inplace=True)

sample_attributes_df = pd.read_csv(sample_attributes_path, sep='\t', low_memory=False)
sample_attributes_df.set_index("SAMPID", inplace=True)


organs = sample_attributes_df["SMTS"].unique().tolist()
UPSERT_ORGAN = """
mutation UpsertOrgan($obj: organs_insert_input!) {
  insert_organs(
    objects: [$obj],
    on_conflict: {
    constraint: organs_name_key, 
    update_columns: [name]
    }
  ) {
    returning { id name }
  }
}
"""
for organ in organs:
    organ_name = gql(UPSERT_ORGAN, {"obj": {"name": organ}})["insert_organs"]["returning"][0]["name"]
    logging.info(organ_name)


UPSERT_TISSUE_WITH_ORGAN = """
mutation UpsertTissueWithOrgan($obj: tissues_insert_input!) {
  insert_tissues(
    objects: [$obj]
    on_conflict: {
      constraint: tissues_name_organ_id_key
      update_columns: [name, organ_id]
    }
  ) {
    returning {
      id
      name
      organ { id name }
    }
  }
}
"""

tissues = sample_attributes_df["SMTSD"].dropna().unique().tolist()

for raw in tissues:
    s = str(raw).strip()
    if not s or s == ".":
        continue

    parts = s.split("-", 1)
    organ_name = parts[0].strip()
    tissue_name = parts[1].strip() if len(parts) > 1 else organ_name

    variables = {
        "obj": {
            "name": tissue_name,
            "organ": {
                "data": {"name": organ_name},
                "on_conflict": {
                    "constraint": "organs_name_key",
                    "update_columns": ["name"]
                }
            }
        }
    }

    res = gql(UPSERT_TISSUE_WITH_ORGAN, variables)
    row = res["insert_tissues"]["returning"][0]
    logging.info(f'{row["organ"]["name"]} -> {row["name"]}')


# --- tissues lookup
ALL_TISSUES_Q = """
query AllTissues {
  tissues { id name organ { name } }
}
"""
rows = gql(ALL_TISSUES_Q, {})["tissues"]
tissue_id_by_key = {(r["organ"]["name"], r["name"]): r["id"] for r in rows}

def split_organ_tissue(label: str) -> tuple[str, str]:
    """
    GTEx SMTSD examples:
      "Brain - Cortex" -> ("Brain", "Cortex")
      "Pancreas"       -> ("Pancreas", "Pancreas")
    """
    s = str(label).strip()
    parts = s.split("-", 1)
    organ = parts[0].strip()
    tissue = parts[1].strip() if len(parts) > 1 else organ
    return organ, tissue

def tissue_id_from_smtsd(smtsd: str) -> int | None:
    organ, tissue = split_organ_tissue(smtsd)
    return tissue_id_by_key.get((organ, tissue))

# --- sources lookup
GET_SOURCES = """query GetSources { sources { id name } }"""
sources = {r["name"].lower(): r["id"] for r in gql(GET_SOURCES, {})["sources"]}

SOURCE_ENSEMBL = sources["ensembl"]
SOURCE_HGNC    = sources["hgnc"]
SOURCE_GTEX    = sources["gtex"]
SOURCE_UNIPARC    = sources["uniparc"]
SOURCE_SWISSPROT    = sources["swissprot"]
SOURCE_TREMBL    = sources["trembl"]
SOURCE_CCDS    = sources["ccds"]
SOURCE_VEP    = sources["vep"]



def build_gtex_expr_by_ensg_fast(expression_data_df, sample_attributes_df):
    sample_ids = [c for c in expression_data_df.columns if c in sample_attributes_df.index]
    smtsd_by_sample = sample_attributes_df.loc[sample_ids, "SMTSD"].astype(str)

    smtsd_to_tid = {}
    for smtsd in smtsd_by_sample.dropna().unique():
        tid = tissue_id_from_smtsd(smtsd)
        if tid is not None:
            smtsd_to_tid[smtsd] = tid

    expr_samples = expression_data_df[sample_ids]
    ensg = expression_data_df.index.astype(str).str.split(".").str[0]

    parts = []
    for smtsd, tid in smtsd_to_tid.items():
        cols = smtsd_by_sample[smtsd_by_sample == smtsd].index
        if len(cols) == 0:
            continue

        sub = expr_samples.loc[:, cols]
        mean_tpm = sub.mean(axis=1).to_numpy()
        median_tpm = sub.median(axis=1).to_numpy()

        parts.append(pd.DataFrame({
            "ensg": ensg.to_numpy(),
            "tissue_id": tid,
            "source_id": SOURCE_GTEX,
            "mean_tpm": mean_tpm,
            "median_tpm": median_tpm,
        }))

    long_df = pd.concat(parts, ignore_index=True).dropna(subset=["mean_tpm", "median_tpm"], how="all")

    # Convert long_df -> dict[ENSG] -> list[dict] (still some Python, but far less)
    out = {}
    for ens, grp in long_df.groupby("ensg", sort=False):
        out[ens] = grp[["source_id", "tissue_id", "mean_tpm", "median_tpm"]].to_dict("records")
    return out


cache_file = os.path.join(DB_DIR, "gtex_expr_by_ensg_cache.joblib")

# Save to cache
if not os.path.exists(cache_file):
    gtex_expr_by_ensg = build_gtex_expr_by_ensg_fast(expression_data_df, sample_attributes_df)
    dump(gtex_expr_by_ensg, cache_file, compress=3)
    logging.info(f"✓ Saved cache to {cache_file}")
else:
    # Load from cache
    gtex_expr_by_ensg = load(cache_file)
    logging.info(f"✓ Loaded from cache: {len(gtex_expr_by_ensg):,} genes")


FETCH_VEP_DATA = """
query FetchVepGenesDistinct($last_id: bigint!, $limit: Int!) {
  staging_vep_annotations_raw(
    where: {id: {_gt: $last_id}},
    order_by: {id: asc},
    limit: $limit
  ) {
    id
    gene
    symbol
    hgnc_id
    uploaded_variation
    allele
    feature
    feature_type
    consequence
    cdna_position
    cds_position
    protein_position
    amino_acids
    codons
    existing_variation
    strand
    variant_class
    symbol_source
    biotype
    canonical
    exon
    phastcons100
    intron
    distance
    ensp
    uniparc
    ccds
    swissprot
    trembl
    domains
    clin_sig
    pheno
    lovd
    sift
    polyphen
    am_class
    am_genome
    am_pathogenicity
    am_protein_variant
    am_transcript_id
    am_uniprot_id
    nmd
    revel
    condel
    blosum62
    loftool
    mirna
    gene_pheno
    pli_gene_value

  }
}
"""

UPSERT_VARIANTS_FULL_BATCH = """
mutation UpsertVariantsFullBatch($objects: [variants_insert_input!]!) {
  insert_variants(
    objects: $objects
    on_conflict: {
      constraint: variants_coordinate_id_ref_alt_key
      update_columns: [ref_validated class nmd_escaping_variant created_at rs_id ref_codon_id alt_codon_id position_in_codon cdna_position cds_position protein_position exon intron gene_ensg]
    }
  ) {
    affected_rows
  }
}
"""



# Consequence term → region mapping lives in helpers (single source of truth),
# so phase3's bystanders_consequences uses the exact same mapping as
# variants_consequences here. See hp.GENOMIC_CATEGORIES.
genomic_categories = hp.GENOMIC_CATEGORIES

# cds = [
#         "Missense", "Synonymous","Nonsense",
#         "Frameshift", "Start_lost", "Stop_retained_variant"
#     ]


# STEP 1: Create a mapping dictionary from consequence term to region
consequence_to_region = hp.CONSEQUENCE_TO_REGION
logging.info(consequence_to_region)

all_consequences = hp.ALL_CONSEQUENCES
logging.info(all_consequences)


# ═══════════════════════════════════════════════════════════════════
# GENETIC CODE CONSTANTS
# ═══════════════════════════════════════════════════════════════════

GENETIC_CODE = {
    'TTT': {'full_name': 'Phenylalanine', '3_letter': 'Phe', '1_letter': 'F'},
    'TTC': {'full_name': 'Phenylalanine', '3_letter': 'Phe', '1_letter': 'F'},
    'TTA': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'TTG': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'TCT': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'TCC': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'TCA': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'TCG': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'TAT': {'full_name': 'Tyrosine', '3_letter': 'Tyr', '1_letter': 'Y'},
    'TAC': {'full_name': 'Tyrosine', '3_letter': 'Tyr', '1_letter': 'Y'},
    'TAA': {'full_name': 'Termination', '3_letter': 'Ter', '1_letter': '*'},
    'TAG': {'full_name': 'Termination', '3_letter': 'Ter', '1_letter': '*'},
    'TGT': {'full_name': 'Cysteine', '3_letter': 'Cys', '1_letter': 'C'},
    'TGC': {'full_name': 'Cysteine', '3_letter': 'Cys', '1_letter': 'C'},
    'TGA': {'full_name': 'Termination', '3_letter': 'Ter', '1_letter': '*'},
    'TGG': {'full_name': 'Tryptophan', '3_letter': 'Trp', '1_letter': 'W'},
    'CTT': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'CTC': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'CTA': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'CTG': {'full_name': 'Leucine', '3_letter': 'Leu', '1_letter': 'L'},
    'CCT': {'full_name': 'Proline', '3_letter': 'Pro', '1_letter': 'P'},
    'CCC': {'full_name': 'Proline', '3_letter': 'Pro', '1_letter': 'P'},
    'CCA': {'full_name': 'Proline', '3_letter': 'Pro', '1_letter': 'P'},
    'CCG': {'full_name': 'Proline', '3_letter': 'Pro', '1_letter': 'P'},
    'CAT': {'full_name': 'Histidine', '3_letter': 'His', '1_letter': 'H'},
    'CAC': {'full_name': 'Histidine', '3_letter': 'His', '1_letter': 'H'},
    'CAA': {'full_name': 'Glutamine', '3_letter': 'Gln', '1_letter': 'Q'},
    'CAG': {'full_name': 'Glutamine', '3_letter': 'Gln', '1_letter': 'Q'},
    'CGT': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'CGC': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'CGA': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'CGG': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'ATT': {'full_name': 'Isoleucine', '3_letter': 'Ile', '1_letter': 'I'},
    'ATC': {'full_name': 'Isoleucine', '3_letter': 'Ile', '1_letter': 'I'},
    'ATA': {'full_name': 'Isoleucine', '3_letter': 'Ile', '1_letter': 'I'},
    'ATG': {'full_name': 'Methionine', '3_letter': 'Met', '1_letter': 'M'},
    'ACT': {'full_name': 'Threonine', '3_letter': 'Thr', '1_letter': 'T'},
    'ACC': {'full_name': 'Threonine', '3_letter': 'Thr', '1_letter': 'T'},
    'ACA': {'full_name': 'Threonine', '3_letter': 'Thr', '1_letter': 'T'},
    'ACG': {'full_name': 'Threonine', '3_letter': 'Thr', '1_letter': 'T'},
    'AAT': {'full_name': 'Asparagine', '3_letter': 'Asn', '1_letter': 'N'},
    'AAC': {'full_name': 'Asparagine', '3_letter': 'Asn', '1_letter': 'N'},
    'AAA': {'full_name': 'Lysine', '3_letter': 'Lys', '1_letter': 'K'},
    'AAG': {'full_name': 'Lysine', '3_letter': 'Lys', '1_letter': 'K'},
    'AGT': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'AGC': {'full_name': 'Serine', '3_letter': 'Ser', '1_letter': 'S'},
    'AGA': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'AGG': {'full_name': 'Arginine', '3_letter': 'Arg', '1_letter': 'R'},
    'GTT': {'full_name': 'Valine', '3_letter': 'Val', '1_letter': 'V'},
    'GTC': {'full_name': 'Valine', '3_letter': 'Val', '1_letter': 'V'},
    'GTA': {'full_name': 'Valine', '3_letter': 'Val', '1_letter': 'V'},
    'GTG': {'full_name': 'Valine', '3_letter': 'Val', '1_letter': 'V'},
    'GCT': {'full_name': 'Alanine', '3_letter': 'Ala', '1_letter': 'A'},
    'GCC': {'full_name': 'Alanine', '3_letter': 'Ala', '1_letter': 'A'},
    'GCA': {'full_name': 'Alanine', '3_letter': 'Ala', '1_letter': 'A'},
    'GCG': {'full_name': 'Alanine', '3_letter': 'Ala', '1_letter': 'A'},
    'GAT': {'full_name': 'Aspartic acid', '3_letter': 'Asp', '1_letter': 'D'},
    'GAC': {'full_name': 'Aspartic acid', '3_letter': 'Asp', '1_letter': 'D'},
    'GAA': {'full_name': 'Glutamic acid', '3_letter': 'Glu', '1_letter': 'E'},
    'GAG': {'full_name': 'Glutamic acid', '3_letter': 'Glu', '1_letter': 'E'},
    'GGT': {'full_name': 'Glycine', '3_letter': 'Gly', '1_letter': 'G'},
    'GGC': {'full_name': 'Glycine', '3_letter': 'Gly', '1_letter': 'G'},
    'GGA': {'full_name': 'Glycine', '3_letter': 'Gly', '1_letter': 'G'},
    'GGG': {'full_name': 'Glycine', '3_letter': 'Gly', '1_letter': 'G'}
}


def split_id_version(s):
    s = (s or "").strip()
    if not s:
        return None, 0
    if "." in s:
        a, b = s.split(".", 1)
        try:
            return a, int(b)
        except:
            return a, 0
    return s, 0




_SCORE_RE = re.compile(
    r"""
    ^\s*
    (?P<label>[A-Za-z _+-]+?)     # qualitative label
    \s*
    (?:\(\s*(?P<value>-?\d*\.?\d+)\s*\))?  # optional numeric in ()
    \s*$
    """,
    re.VERBOSE,
)

def parse_qualitative_quantitative(
    s: Optional[str],
) -> Tuple[Optional[str], Optional[float]]:
    """
    Parse strings like:
      'tolerated(0.3)'
      'probably_damaging (0.987)'
      'benign'
      None

    Returns:
      (qualitative, quantitative)
    """
    if not s:
        return None, None

    m = _SCORE_RE.match(str(s))
    if not m:
        return s.strip(), None

    label = m.group("label").strip()
    value = m.group("value")

    return label or None, float(value) if value is not None else None


# use the constraint name you showed in Hasura:
PROTEIN_DOMAINS_ON_CONFLICT = {
    "constraint": "protein_domains_source_id_accession_number_key",
    "update_columns": ["source_id"],  
}

DOMAIN_KEY_TO_SOURCE_ID = {
    "Gene3D": cath_Gene3D_id,
    "HAMAP": HAMAP_id,
    "PANTHER": panther_id,
    "PIRSF": pirsf_id,
    "PROSITE_patterns": prosite_patterns_id,
    "PROSITE_profiles": prosite_profiles_id,
    "Pfam": pfam_id,
    "Prints": prints_id,
    "SMART": smart_id,
    "Superfamily": superfamily_id,
    "TIGRFAM": tigrfam_id,
}

EXCLUDED_DOMAIN_KEYS = {
    "Cleavage_site_(Signalp)",
    "Coiled-coils_(Ncoils)",
    "Transmembrane_helices",
    "Low_complexity_(Seg)",
}

def parse_domains_for_insert(domains_str: str):
    """
    Returns LIST of join rows for protein_domains_proteins.data.

    Each element:
      {
        "protein_domain": {
          "data": {"source_id": ..., "accession_number": ...},
          "on_conflict": {...}
        }
      }
    """
    if not domains_str:
        return []

    links = []
    for token in str(domains_str).split(","):
        token = token.strip()
        if not token or ":" not in token:
            continue

        key, accession = token.split(":", 1)
        key = key.strip()
        accession = accession.strip()

        if not key or not accession:
            continue
        if key in EXCLUDED_DOMAIN_KEYS:
            continue

        source_id = DOMAIN_KEY_TO_SOURCE_ID.get(key)
        if source_id is None:
            continue

        links.append({
            "protein_domain": {
                "data": {
                    "source_id": source_id,
                    "accession_number": accession,
                },
                "on_conflict": PROTEIN_DOMAINS_ON_CONFLICT,
            }
        })

    return links


VARIANT_CLASS_MAP = {
    "snv": "SNV",
    "deletion": "Deletion",
    "insertion": "Insertion",
    "substitution": "Substitution",
    "indel": "Indel",
    "sequence_alteration": "Sequence_alteration",
}

def normalize_variant_class(v):
    if v is None:
        return None
    s = str(v).strip()
    # keep exact SNV if it comes already correct
    if s == "SNV":
        return "SNV"
    s_l = s.lower()
    return VARIANT_CLASS_MAP.get(s_l, None)  # or "Sequence_alteration" / "SNV" as fallback


logging.info(f"Start inserting annotated variants data...")


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

BATCH_READ   = 5000
BATCH_INSERT = 500

last_id = 0
total_variants = 0
batch_idx = 0
errors = []
total_ok = 0

t0_global = time.time()

while True:
    batch_idx += 1
    t0_batch = time.time()

    rows = gql(
        FETCH_VEP_DATA,
        {"last_id": last_id, "limit": BATCH_READ}
    )["staging_vep_annotations_raw"]

    if not rows:
        elapsed = time.time() - t0_global
        logging.info(f"\n✅ Done. Total variants processed: {total_variants:,} in {elapsed/60:.1f} min")
        break

    last_id = rows[-1]["id"]

    objects = []


    # process rows in chunks of WRITE_BATCH_SIZE
    for i in range(0, len(rows), BATCH_INSERT):
        chunk = rows[i : i + BATCH_INSERT]
        variant_objects = []
        staging_ids = []

        for row in chunk:
            staging_ids.append(row["id"])
            last_id = row["id"]

            chrom, start, ref, alt = row["uploaded_variation"].split(":")
            alt = row["allele"] if "/" in alt else alt
            start = int(start)
            end = start + len(ref) - 1
            # ---------- Gene Expression ----------
            gene_expression = None
            expressions_data = []
            if row["gene"] is not None:
                seen = set()
                for e in gtex_expr_by_ensg.get(str(row["gene"]), []):
                    tissue_id = e["tissue_id"]
                    if tissue_id is None:
                        continue
                    if tissue_id in seen:
                        continue
                    seen.add(tissue_id)
                    expressions_data.append(
                        {
                            "source_id": SOURCE_GTEX,
                            "tissue_id": tissue_id,
                            "mean_tpm": e["mean_tpm"],
                            "median_tpm": e["median_tpm"],
                        }
                    )
                gene_expression = {
                    "data": expressions_data,
                    "on_conflict": {
                        "constraint": "gene_expressions_tissue_id_gene_ensg_source_id_key",
                        "update_columns": ["mean_tpm", "median_tpm"],
                    },
                }



            # ---------- Variant Scores ----------
            scores_data = []

            # ---------- Variant Feature ----------
            feature_obj = None
            if row["feature"] is not None:
                if row["biotype"] is None:
                    feature_obj = {
                        "data": {
                            "type": row["feature_type"],
                            "identifier": row["feature"].split(".")[0],
                            "version_number": row["feature"].split(".")[1] if "." in row["feature"] else 0,
                        },
                        "on_conflict": {
                            "constraint": "features_type_identifier_version_number_key",
                            "update_columns": ["type"],
                        },
                    }
                else:
                    feature_obj = {
                        "data": {
                            "type": row["feature_type"],
                            "identifier": row["feature"].split(".")[0],
                            "version_number": row["feature"].split(".")[1] if "." in row["feature"] else 0,
                            "biotype": {
                                "data": {
                                    "name": row["biotype"]
                                },
                                "on_conflict": {
                                    "constraint": "biotypes_name_key",
                                    "update_columns": ["name"],
                                },
                            } 
                        },
                        "on_conflict": {
                            "constraint": "features_type_identifier_version_number_key",
                            "update_columns": ["type"],
                        },
                    }

            # ---------- Variant Consequences ----------
            consequences_data = []
            consequences = [c.strip() for c in str(row["consequence"]).split(",")] if row["consequence"] is not None else []
            for consequence in consequences:
                if consequence not in all_consequences:
                    mapped_region = "Other"
                else:
                    mapped_region = consequence_to_region.get(consequence, "Other")
                consequences_data.append({
                    "name": mapped_region,
                })

            # ---------- Codons data ----------
            ref_codon_object = None
            alt_codon_object = None
            position_in_codon = None
            if row["codons"] is not None:
                codons_parts = row["codons"].split("/")
                ref_codon = codons_parts[0].strip() if len(codons_parts) > 1 else None
                alt_codon = codons_parts[1].strip() if len(codons_parts) > 1 else None
                if len(codons_parts) == 2:
                    if ref_codon is not None and len(ref_codon) != 3:
                        ref_codon = None
                    if alt_codon is not None and len(alt_codon) != 3:
                        alt_codon = None
                if ref_codon is not None:
                    position_in_codon = next((i + 1 for i, c in enumerate(ref_codon) if c.isupper()), None)
                    ref_codon_object = {
                        "data": {
                            "is_stop": True,
                            "nt1": ref_codon[0].upper(),
                            "nt2": ref_codon[1].upper(),
                            "nt3": ref_codon[2].upper()

                        } if ref_codon.upper().replace('U', 'T') in ["TAA", "TAG", "TGA"] else {
                            "amino_acid": {
                                "data": {
                                    "full_name": GENETIC_CODE[ref_codon.upper()]["full_name"],
                                    "short_name": GENETIC_CODE[ref_codon.upper()]["3_letter"],
                                    "letter": GENETIC_CODE[ref_codon.upper()]["1_letter"]
                                },
                                "on_conflict": {
                                    "constraint": "amino_acids_full_name_key",
                                    "update_columns": ["full_name"],
                                },
                            },
                            "nt1": ref_codon[0].upper(),
                            "nt2": ref_codon[1].upper(),
                            "nt3": ref_codon[2].upper(),
                            "is_stop": False
                        },
                        "on_conflict": {
                            "constraint": "codons_nt1_nt2_nt3_key",
                            "update_columns": ["is_stop"],
                        },
                    }
                if alt_codon is not None:
                    alt_codon_object = {
                        "data": {
                            "is_stop": True,
                            "nt1": alt_codon[0].upper(),
                            "nt2": alt_codon[1].upper(),
                            "nt3": alt_codon[2].upper()

                        } if alt_codon.upper().replace('U', 'T') in ["TAA", "TAG", "TGA"] else {
                            "amino_acid": {
                                "data": {
                                    "full_name": GENETIC_CODE[alt_codon.upper()]["full_name"],
                                    "short_name": GENETIC_CODE[alt_codon.upper()]["3_letter"],
                                    "letter": GENETIC_CODE[alt_codon.upper()]["1_letter"]
                                },
                                "on_conflict": {
                                    "constraint": "amino_acids_full_name_key",
                                    "update_columns": ["full_name"],
                                },
                            },
                            "nt1": alt_codon[0].upper(),
                            "nt2": alt_codon[1].upper(),
                            "nt3": alt_codon[2].upper(),
                            "is_stop": False
                        },
                        "on_conflict": {
                            "constraint": "codons_nt1_nt2_nt3_key",
                            "update_columns": ["is_stop"],
                        },
                    }

            transcripts_obj = None

            if row.get("feature_type") == "Transcript" and row.get("feature"):
                enst, enst_ver = split_id_version(row["feature"])

                # Optional fields from VEP if you have them (adjust to your columns):
                ensembl_canonical = True if str(row.get("canonical")).upper() in ("YES","1","TRUE") else None
                protein_obj = None
                if row.get("ensp") is not None:
                    ensp, ensp_ver = split_id_version(row["ensp"])
                    
                    domain_links = []
                    if row.get("domains"):
                        domain_links = parse_domains_for_insert(str(row["domains"]))  # LIST

                    protein_obj = {
                        "data": {
                            "ensp": ensp,
                            "proteins_annotations": {
                                "data": [
                                    {
                                        "source_id": SOURCE_ENSEMBL,
                                        "identifier": ensp,
                                        "version_number": ensp_ver,
                                    }
                                ],
                                "on_conflict": {
                                    "constraint": "proteins_annotations_source_id_identifier_version_number_key",
                                    "update_columns": ["identifier", "version_number"],
                                },
                            },

                            # ✅ correct: array rel insert
                            "protein_domains_proteins": {
                                "data": domain_links,
                                "on_conflict": {
                                    "constraint": "protein_domains_proteins_pkey",
                                    "update_columns": ["protein_ensp"],
                                }
                            },
                        },
                        "on_conflict": {
                            "constraint": "proteins_pkey",
                            "update_columns": ["ensp"],
                        },
                    }

                    if row["uniparc"] is not None:
                        uniparc_annotation_data = {
                            "source_id": SOURCE_UNIPARC,
                            "identifier": str(row["uniparc"])
                        }
                        protein_obj["data"]["proteins_annotations"]["data"].append(uniparc_annotation_data)
                    if row["ccds"] is not None:
                        ccds_annotation_data = {
                            "source_id": SOURCE_CCDS,
                            "identifier": str(row["ccds"])
                        }
                        protein_obj["data"]["proteins_annotations"]["data"].append(ccds_annotation_data)
                    if row["swissprot"] is not None:
                        swissprot_annotation_data = {
                            "source_id": SOURCE_SWISSPROT,
                            "identifier": str(row["swissprot"])
                        }
                        protein_obj["data"]["proteins_annotations"]["data"].append(swissprot_annotation_data)
                    if row["trembl"] is not None:
                        trembl_annotation_data = []
                        for trembl_id in str(row["trembl"]).split(","):
                            trembl_annotation_data.append({
                                "source_id": SOURCE_TREMBL,
                                "identifier": trembl_id.strip()
                            })
                        protein_obj["data"]["proteins_annotations"]["data"] += trembl_annotation_data

                transcript_data = {
                    # If you add columns to transcripts (recommended):
                    "enst": enst,

                    # Existing columns:
                    "ensembl_canonical": True if str(ensembl_canonical).upper() in ("YES","1","TRUE") else None,

                    "protein": protein_obj if protein_obj is not None else None,


                    # If you want to link protein/sequence later, you can add nested objects here
                    # protein: {...}
                    # sequence: {...}

                    "transcripts_annotations": {
                        "data": [
                            {
                                "source_id": SOURCE_ENSEMBL,
                                "identifier": enst,
                                "version_number": enst_ver,
                            }
                        ],
                        "on_conflict": {
                            # use the FIXED constraint name after you create it
                            "constraint": "transcripts_annotations_source_id_identifier_version_number_key",
                            "update_columns": ["identifier", "version_number"],
                        },
                    },
                }

                transcripts_obj = {
                    "data": [transcript_data],
                    "on_conflict": {
                        # if you added transcripts_enst_version_key:
                        "constraint": "transcripts_pkey",
                        "update_columns": ["ensembl_canonical"],
                    },
                }

            # ---------- Variant scores (JOIN rows) ----------
            vq_rows = []
            vqt_rows = []

            def add_score(*, kind: str, name: str, score_type: str, value):
                """
                kind: 'qual' or 'quan'
                """
                if value is None:
                    return

                if kind == "qual":
                    vq_rows.append({
                        "qualitative_score": {
                            "data": {
                                "source_id": SOURCE_VEP,
                                "name": name,
                                "score_type": score_type,
                                "value": str(value),  # usually enum/text
                            },
                            "on_conflict": {
                                "constraint": "qualitative_scores_source_id_score_type_name_value_key",
                                "update_columns": ["name"],
                            },
                        }
                    })
                elif kind == "quan":
                    vqt_rows.append({
                        "quantitative_score": {
                            "data": {
                                "source_id": SOURCE_VEP,
                                "name": name,
                                "score_type": score_type,
                                "value": float(value),
                            },
                            "on_conflict": {
                                "constraint": "quantitative_scores_source_id_score_type_name_value_key",
                                "update_columns": ["name"],
                            },
                        }
                    })
                else:
                    raise ValueError(f"Unknown kind: {kind}")


            if row["sift"] is not None:
                ql, qn = parse_qualitative_quantitative(str(row["sift"]))
                add_score(kind="qual", name="SIFT", score_type="Pathogenicity", value=ql)
                add_score(kind="quan", name="SIFT", score_type="Pathogenicity", value=qn)

            if row["polyphen"] is not None:
                ql, qn = parse_qualitative_quantitative(str(row["polyphen"]))
                add_score(kind="qual", name="PolyPhen", score_type="Pathogenicity", value=ql)
                add_score(kind="quan", name="PolyPhen", score_type="Pathogenicity", value=qn)
            if row["phastcons100"] is not None and len(ref) == 1 and len(alt) == 1:
                add_score(kind="quan", name="PhastCons100",   score_type="Preservation",   value=row.get("phastcons100").split(",")[0])
            if row["blosum62"] is not None:
                add_score(kind="quan", name="BLOSUM62",       score_type="Pathogenicity",  value=row.get("blosum62"))
            if row["am_class"] is not None:
                add_score(kind="qual", name="AlphaMissense",  score_type="Pathogenicity",  value=row.get("am_class"))
            if row["am_pathogenicity"] is not None:
                add_score(kind="quan", name="AlphaMissense",  score_type="Pathogenicity",  value=row.get("am_pathogenicity"))
            if row["gene_pheno"] is not None:
                add_score(kind="quan", name="GENE_PHENO",     score_type="Other",          value=row.get("gene_pheno"))
            if row["pli_gene_value"] is not None:
                add_score(kind="quan", name="pLI",            score_type="Pathogenicity",  value=row.get("pli_gene_value"))
            if row["loftool"] is not None:
                add_score(kind="quan", name="LoFtool",        score_type="Pathogenicity",  value=row.get("loftool"))



            # ---------- main variant object ----------


            variant_obj = {
                "ref": ref,
                "alt": alt,
                "ref_validated": True,
                "created_at": datetime.utcnow().isoformat(),
                "alt_codon": alt_codon_object if alt_codon_object is not None else None,  
                "ref_codon": ref_codon_object if ref_codon_object is not None else None,


                "exon": str(row["exon"]) if row["exon"] is not None else None,
                "intron": str(row["intron"]) if row["intron"] is not None else None,
                "nmd_escaping_variant": True if str(row["nmd"]) == "NMD_escaping_variant" else False,

                "coordinate": {
                    "data": {
                        "genome_build_id": 1,
                        "chr": str(chrom),
                        "start": start,
                        "end": end,
                        "strand": "+"
                    },
                    "on_conflict": {
                        "constraint": "coordinates_genome_build_id_chr_start_end_strand_key",
                        "update_columns": ["strand"],
                    },
                },

                "gene": {
                    "data": {
                        "symbol": str(row["symbol"]) if row["symbol"] is not None else None,
                        "ensg": str(row["gene"]) if row["gene"] is not None else None,
                        "strand": "+" if str(row["strand"]) == "1" else "-",
                        "genes_annotations": {
                            "data": [
                                {
                                    "source_id": SOURCE_ENSEMBL,
                                    "identifier": str(row["gene"]) if row["gene"] is not None else None,
                                },
                                {
                                    "source_id": SOURCE_HGNC,
                                    "identifier": str(row["hgnc_id"]) if row["hgnc_id"] is not None else None
                                }
                            ] if row["hgnc_id"] is not None else [
                                {
                                    "source_id": SOURCE_ENSEMBL,
                                    "identifier": str(row["gene"]),
                                }
                            ],
                            "on_conflict": {
                                "constraint": "genes_annotations_source_id_identifier_version_number_key",
                                "update_columns": ["identifier"],
                            },
                        },
                        "gene_expressions": gene_expression,
                        "transcripts": transcripts_obj, 
                    },
                    "on_conflict": {
                        "constraint": "genes_pkey",
                        "update_columns": ["symbol"],
                    }
                } if row["gene"] is not None else None,
                "position_in_codon": position_in_codon,
                "cdna_position": str(row["cdna_position"]) if row["cdna_position"] is not None else None,
                "cds_position": str(row["cds_position"]) if row["cds_position"] is not None else None,
                "protein_position": str(row["protein_position"]) if row["protein_position"] is not None else None,
                "class": normalize_variant_class(row.get("variant_class")),


            }
            if feature_obj is not None:
                variant_obj["variants_features"] = {
                    "data": {
                        "feature": feature_obj,
                    },
                    "on_conflict": {
                        "constraint": "variants_features_pkey",
                        "update_columns": ["variant_id"],
                    },
                }

            if consequences_data:
                variants_consequences_data = []
                for c in consequences_data:
                    variants_consequences_data.append({
                        # IMPORTANT: this field name must match your object-relationship name
                        # in variants_consequences -> consequences
                        "consequence": {
                            "data": c,
                            "on_conflict": {
                                "constraint": "consequences_name_key",
                                "update_columns": ["name"],
                            },
                        }
                    })

                variant_obj["variants_consequences"] = {
                    "data": variants_consequences_data,
                    # IMPORTANT: this constraint must be a UNIQUE key on the join table,
                    # typically (variant_id, consequence_id), NOT the PK id.
                    "on_conflict": {
                        "constraint": "variants_consequences_pkey",
                        "update_columns": ["consequence_id"],
                    },
                }

            if vq_rows:
                variant_obj["variants_qualitative_scores"] = {
                    "data": vq_rows,
                    "on_conflict": {
                        # IMPORTANT: this should be UNIQUE on (variant_id, qualitative_score_id), not the PK id
                        "constraint": "variants_qualitative_scores_pkey",
                        "update_columns": ["qualitative_score_id"],
                    },
                }

            if vqt_rows:
                variant_obj["variants_quantitative_scores"] = {
                    "data": vqt_rows,
                    "on_conflict": {
                        "constraint": "variants_quantitative_scores_pkey",
                        "update_columns": ["quantitative_score_id"],
                    },
                }
            
            variant_obj["variants_sources"] = {
                "data": [
                    {
                        "source_id": SOURCE_VEP,
                    }
                ],
                "on_conflict": {
                    "constraint": "variants_sources_pkey",
                    "update_columns": ["source_id"],
                },
            }
            if row["existing_variation"] is not None:
                variant_obj["rs_id"] = str(row["existing_variation"]).strip()



            variant_objects.append(variant_obj)

        # 2) write this chunk in one mutation
        try:
            result = gql(UPSERT_VARIANTS_FULL_BATCH, {"objects": variant_objects})
            inserted = result["insert_variants"]
            total_ok += inserted.get("affected_rows", 0)

        except Exception as e:
            # log batch error (you can also fall back to row-by-row here if you like)
            errors.append({
                "staging_ids": staging_ids,
                "error": str(e),
            })

    logging.info(f"Processed up to staging id={last_id}, ok={total_ok}, errors={len(errors)}")

logging.info(f"DONE. Total OK: {total_ok}, Errors: {len(errors)}")

logging.info(f"Total errors: {len(errors)}")
for err in errors[:20]:  # show first 20
    logging.info(f"  staging_ids: {err['staging_ids'][:3]}...  error: {err['error'][:300]}")

VARIANTS_COUNT = """
query VariantsCount {
  variants_aggregate {
    aggregate {
      count
      }
    }
  }
"""

res = gql(VARIANTS_COUNT, {})
logging.info(f"Total variants in DB: {res['variants_aggregate']['aggregate']['count']}")



# ── 1. Staging table counts & ID gaps ──
STAGING_VEP_INFO = """
query StagingVepInfo {
  staging_vep_annotations_raw_aggregate {
    aggregate {
      count
      max { id }
      min { id }
    }
  }
  variants_aggregate {
    aggregate { count }
  }
  coordinates_aggregate {
    aggregate { count }
  }
}
"""
res = gql(STAGING_VEP_INFO, {})
staging = res["staging_vep_annotations_raw_aggregate"]["aggregate"]
variants_count = res["variants_aggregate"]["aggregate"]["count"]
coords_count = res["coordinates_aggregate"]["aggregate"]["count"]

logging.info(f"Staging VEP records count: {staging['count']}")
logging.info(f"Staging VEP ID range: {staging['min']['id']} – {staging['max']['id']}")
logging.info(f"ID gaps (max - min + 1 - count): {staging['max']['id'] - staging['min']['id'] + 1 - staging['count']}")
logging.info(f"Variants in DB: {variants_count}")
logging.info(f"Coordinates in DB: {coords_count}")
logging.info(f"Staging − Variants = {staging['count'] - variants_count} (duplicate uploaded_variations)")

# ── 2. NULL uploaded_variation / allele ──
STAGING_INVALID = """
query StagingInvalid {
  staging_vep_annotations_raw_aggregate(
    where: {
        allele: {_is_null: true}
    }
  ) {
    aggregate { count }
  }
}
"""
res = gql(STAGING_INVALID, {})
logging.info(f"\nVEP annotations with NULL allele: {res['staging_vep_annotations_raw_aggregate']['aggregate']['count']}")

# ── 3. NULL / empty consequence ──
NO_CONSEQUENCES = """
query NoConsequences {
  staging_vep_annotations_raw_aggregate(
    where: {
      _or: [
        {consequence: {_is_null: true}},
        {consequence: {_eq: ""}}
      ]
    }
  ) {
    aggregate { count }
  }
}
"""
res = gql(NO_CONSEQUENCES, {})
logging.info(f"VEP annotations with no consequence: {res['staging_vep_annotations_raw_aggregate']['aggregate']['count']}")



# Step 1: Query VEP annotations with clinical significance data
GET_VEP_CLINSIG = """
query GetVepClinSig($last_id: bigint!, $limit: Int!) {
    staging_vep_annotations_raw(
        where: {
            id: {_gt: $last_id}
            clin_sig: {_is_null: false}
        }
        order_by: {id: asc}
        limit: $limit
    ) {
        id
        uploaded_variation
        allele
        clin_sig
    }
}
"""

# Step 2: Mutation to upsert clinical significance scores
UPSERT_CLINSIG_SCORES = """
mutation UpsertClinSigScores($objects: [variants_qualitative_scores_insert_input!]!) {
    insert_variants_qualitative_scores(
        objects: $objects
        on_conflict: {
            constraint: variants_qualitative_scores_pkey
            update_columns: [qualitative_score_id]
        }
    ) {
        affected_rows
    }
}
"""

logging.info("✅ ClinVar processing queries defined")

# Step 3: Process VEP annotations and build score links
logging.info("Processing ClinVar clinical significance data from VEP annotations...")

# Add a query to find variant ID by coordinates
GET_VARIANT_ID = """
query GetVariantId($chr: chr!, $start: Int!, $ref: String!, $alt: String!) {
  variants(
    where: {
      ref: {_eq: $ref}
      alt: {_eq: $alt}
      coordinate: {
        chr: {_eq: $chr}
        start: {_eq: $start}
        genome_build_id: {_eq: 1}
      }
    }
    limit: 1
  ) {
    id
  }
}
"""

# Simplified mutation - just link variant_id to score
UPSERT_CLINSIG_SCORES_SIMPLE = """
mutation UpsertClinSigScores($objects: [variants_qualitative_scores_insert_input!]!) {
    insert_variants_qualitative_scores(
        objects: $objects
        on_conflict: {
            constraint: variants_qualitative_scores_pkey
            update_columns: [qualitative_score_id]
        }
    ) {
        affected_rows
    }
}
"""

BATCH_READ = 5000
BATCH_INSERT = 500

# 🔄 RESUMABILITY: Change this to resume from a specific ID if interrupted
START_FROM_ID = 0  # Set to last processed ID to resume

last_id = START_FROM_ID
total_processed = 0
total_inserted = 0
total_skipped = 0
errors = []

try:
    while True:
        # Fetch batch of VEP annotations with clin_sig
        try:
            rows = gql(GET_VEP_CLINSIG, {
                "last_id": last_id,
                "limit": BATCH_READ
            })["staging_vep_annotations_raw"]
        except Exception as e:
            logging.info(f"⚠️  Error fetching batch at id={last_id}: {str(e)[:150]}")
            errors.append({"type": "fetch", "last_id": last_id, "error": str(e)})
            break
        
        if not rows:
            break
        
        last_id = rows[-1]["id"]
        
        # Process in smaller insert batches
        for i in range(0, len(rows), BATCH_INSERT):
            chunk = rows[i : i + BATCH_INSERT]
            score_objects = []
            
            for row in chunk:
                total_processed += 1
                
                # Parse uploaded_variation to get variant coordinates
                try:
                    chrom, start, ref, alt = row["uploaded_variation"].split(":")
                    # Use allele from VEP (handles multi-allelic correctly)
                    alt = row["allele"] if "/" in alt else alt
                    start = int(start)
                except (ValueError, AttributeError):
                    continue
                
                # Find variant ID by coordinates
                try:
                    variant_result = gql(GET_VARIANT_ID, {
                        "chr": str(chrom),
                        "start": start,
                        "ref": ref,
                        "alt": alt
                    })
                    if not variant_result["variants"]:
                        total_skipped += 1
                        continue
                    
                    variant_id = variant_result["variants"][0]["id"]
                except Exception as e:
                    total_skipped += 1
                    continue
                
                # Parse clinical significance (can be multiple values separated by ,)
                clin_sig_raw = str(row["clin_sig"]).strip()
                if not clin_sig_raw or clin_sig_raw.lower() in ['-', 'nan', 'none', '']:
                    continue
                
                # Split multiple ClinVar classifications (e.g., "pathogenic,likely_pathogenic")
                clin_sig_values = clin_sig_raw.split(',')
                
                for clin_sig_value in clin_sig_values:
                    clin_sig_value = clin_sig_value.strip()
                    if not clin_sig_value:
                        continue
                    
                    # Build the simplified insert object - just link existing variant to score
                    score_objects.append({
                        "variant_id": variant_id,
                        "qualitative_score": {
                            "data": {
                                "source_id": SOURCE_VEP,
                                "name": "ClinVar",
                                "score_type": "Clinical_significance",
                                "value": clin_sig_value
                            },
                            "on_conflict": {
                                "constraint": "qualitative_scores_source_id_score_type_name_value_key",
                                "update_columns": ["name"]
                            }
                        }
                    })
            
            # Insert batch
            if score_objects:
                try:
                    result = gql(UPSERT_CLINSIG_SCORES_SIMPLE, {"objects": score_objects})
                    inserted = result["insert_variants_qualitative_scores"]["affected_rows"]
                    total_inserted += inserted
                except Exception as e:
                    errors.append({
                        "type": "insert",
                        "batch_start_id": chunk[0]["id"],
                        "error": str(e)
                    })
                    logging.info(f"⚠️  Error in batch starting at VEP id {chunk[0]['id']}: {str(e)[:100]}")
        
        logging.info(f"Processed up to VEP id={last_id}, processed={total_processed:,}, inserted={total_inserted:,}, skipped={total_skipped:,}")

except KeyboardInterrupt:
    logging.info(f"\n⚠️  Interrupted by user at VEP id={last_id}")
except Exception as e:
    logging.info(f"\n⚠️  Unexpected error: {str(e)}")
    errors.append({"type": "unexpected", "last_id": last_id, "error": str(e)})

logging.info(f"\n{'='*60}")
logging.info("CLINVAR PROCESSING SUMMARY")
logging.info(f"{'='*60}")
logging.info(f"Total VEP rows processed: {total_processed:,}")
logging.info(f"Total score links inserted: {total_inserted:,}")
logging.info(f"Total variants not found (skipped): {total_skipped:,}")
logging.info(f"Last processed VEP id: {last_id}")
logging.info(f"Errors: {len(errors)}")
logging.info(f"{'='*60}")

if errors:
    logging.info("\n⚠️  Errors encountered:")
    for err in errors[:5]:
        logging.info(f"  {err.get('type', 'unknown')}: {err.get('error', 'N/A')[:150]}")
    if len(errors) > 5:
        logging.info(f"  ... and {len(errors) - 5} more errors")
    logging.info(f"\n💡 To resume from last processed ID, set START_FROM_ID = {last_id}")

# Verify ClinVar data was added successfully
GET_CLINVAR_STATS = """
query GetClinVarStats {
  qualitative_scores_aggregate(
    where: {
      name: {_eq: "ClinVar"}
      score_type: {_eq: "Clinical_significance"}
    }
  ) {
    aggregate {
      count
    }
    nodes {
      value
      variants_qualitative_scores_aggregate {
        aggregate { count }
      }
    }
  }
}
"""

clinvar_stats = gql(GET_CLINVAR_STATS, {})
total_clinvar_scores = clinvar_stats["qualitative_scores_aggregate"]["aggregate"]["count"]
logging.info(f"\n✅ Total unique ClinVar clinical significance values: {total_clinvar_scores}")

if total_clinvar_scores > 0:
    logging.info("\nBreakdown by classification:")
    for node in clinvar_stats["qualitative_scores_aggregate"]["nodes"]:
        variant_count = node["variants_qualitative_scores_aggregate"]["aggregate"]["count"]
        logging.info(f"  • {node['value']}: {variant_count:,} variants")



UPSERT_SOURCE = """
mutation UpsertSource($obj: sources_insert_input!) {
  insert_sources(
    objects: [$obj],
    on_conflict: {
    constraint: sources_name_release_key, 
    update_columns: [name]
    }
  ) {
    returning { id name url description release}
  }
}
"""

# 1) source: Ensembl
sfari = gql(UPSERT_SOURCE, {"obj": {"name": "SFARI Gene", "url": "https://gene.sfari.org/", "description": "SFARI Gene utilizes a systems biology approach, linking information on autism candidate genes within its original “Human Gene” module to corresponding data from a diverse array of supplementary data modules. ASD risk genes are then scored using a set of annotation rules developed in consultation with an external advisory board and classified into specific categories based on the evidence supporting their link to autism."}})["insert_sources"]["returning"][0]
sfari_id = sfari["id"]

logging.info(f"sfari_id: {sfari_id}")

path_sfari_genes = resolve(os.environ["SFARI_GENES"])

FETCH_GENES_ENSG = """
query FetchGenesEnsg($last: String!, $limit: Int!) {
  genes(
    where: { ensg: { _gt: $last } }
    order_by: { ensg: asc }
    limit: $limit
  ) {
    ensg
  }
}
"""

UPSERT_QUANT_SCORE = """
mutation UpsertQuantScore($score: quantitative_scores_insert_input!) {
  insert_quantitative_scores(
    objects: [$score],
    on_conflict: {
      constraint: quantitative_scores_source_id_score_type_name_value_key
      update_columns: [name, score_type, value, source_id]
    }
  ) {
    returning { id }
  }
}
"""


UPSERT_QUAL_SCORE = """
mutation UpsertQualScore($score: qualitative_scores_insert_input!) {
  insert_qualitative_scores(
    objects: [$score],
    on_conflict: {
      constraint: qualitative_scores_source_id_score_type_name_value_key
      update_columns: [name, score_type, value, source_id]
    }
  ) {
    returning { id }
  }
}
"""

LINK_GENE_QUANT = """
mutation LinkGeneQuant($link: genes_quantitative_scores_insert_input!) {
  insert_genes_quantitative_scores(
    objects: [$link],
    on_conflict: {
      constraint: genes_quantitative_scores_pkey
      update_columns: [gene_ensg, quantitative_score_id]
    }
  ) {
    affected_rows
  }
}
"""

LINK_GENE_QUAL = """
mutation LinkGeneQual($link: genes_qualitative_scores_insert_input!) {
  insert_genes_qualitative_scores(
    objects: [$link],
    on_conflict: {
      constraint: genes_qualitative_scores_pkey
      update_columns: [gene_ensg, qualitative_score_id]
    }
  ) {
    affected_rows
  }
}
"""

def fetch_all_gene_ensg(limit=5000):
    existing = set()
    last = ""
    while True:
        res = gql(FETCH_GENES_ENSG, {"last": last, "limit": limit})
        rows = res["genes"]
        if not rows:
            break
        for g in rows:
            if g["ensg"]:
                existing.add(g["ensg"])
        last = rows[-1]["ensg"]
        logging.info(f"Loaded genes ensg: {len(existing):,} (last={last})")
    return existing


def safe_str(x):
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.lower() in ("nan", "na", "none"):
        return None
    return s

def safe_float(x):
    s = safe_str(x)
    if s is None:
        return None
    try:
        return float(s)
    except Exception:
        return None

FETCH_GENES_BY_SYMBOLS = """
query FetchGenesBySymbols($symbols: [String!]!) {
  genes(where: { symbol: { _in: $symbols } }) {
    ensg
    symbol
  }
}
"""

def load_sfari_gene_scores(sfari_csv_path: str, sfari_source_id: int, *, progress_every=2000):
    # 0) load existing genes once
    existing_genes = fetch_all_gene_ensg()
    logging.info(f"✅ Existing genes loaded: {len(existing_genes):,}")

    sfari_df = pd.read_csv(sfari_csv_path, dtype=str, low_memory=False)

    required = ["gene-symbol", "gene-score", "syndromic", "ensembl-id", "number-of-reports"]
    missing = [c for c in required if c not in sfari_df.columns]
    if missing:
        raise ValueError(f"Missing columns in SFARI CSV: {missing}")

    # 0b) Build symbol->ensg map for rows missing ensembl-id
    no_ensg_symbols = [
        safe_str(s) for s in sfari_df.loc[
            sfari_df["ensembl-id"].apply(lambda x: safe_str(x) is None), "gene-symbol"
        ] if safe_str(s)
    ]
    symbol_to_ensg = {}
    if no_ensg_symbols:
        found = gql(FETCH_GENES_BY_SYMBOLS, {"symbols": no_ensg_symbols})["genes"]
        symbol_to_ensg = {g["symbol"]: g["ensg"] for g in found}
        logging.info(f"✅ Symbol fallback: {len(symbol_to_ensg)}/{len(no_ensg_symbols)} resolved")

    t0 = time.time()
    ok = 0
    resolved_by_symbol = 0
    skipped_no_ensg = 0
    skipped_missing_gene = 0
    skipped_no_values = 0
    errs = []

    n = len(sfari_df)

    for idx, r in sfari_df.iterrows():
        if (idx + 1) % progress_every == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            logging.info(
                f"[{idx+1:,}/{n:,}] ok={ok:,} "
                f"resolved_by_symbol={resolved_by_symbol:,} "
                f"skip_missing_gene={skipped_missing_gene:,} "
                f"skip_no_values={skipped_no_values:,} "
                f"errors={len(errs):,} "
                f"({rate:,.1f} rows/s)"
            )

        ensg = safe_str(r["ensembl-id"])
        if not ensg:
            # Fallback: resolve by gene symbol
            symbol = safe_str(r["gene-symbol"])
            ensg = symbol_to_ensg.get(symbol) if symbol else None
            if ensg:
                resolved_by_symbol += 1
            else:
                skipped_no_ensg += 1
                continue

        # Skip if gene doesn't exist in DB
        if ensg not in existing_genes:
            skipped_missing_gene += 1
            continue

        gene_score = safe_float(r["gene-score"])          # quantitative
        syndromic  = safe_str(r["syndromic"])             # qualitative
        n_reports  = safe_float(r["number-of-reports"])   # quantitative

        if gene_score is None and syndromic is None and n_reports is None:
            skipped_no_values += 1
            continue

        try:
            # --- quantitative: SFARI Gene Score ---
            if gene_score is not None:
                res = gql(
                    UPSERT_QUANT_SCORE,
                    {"score": {
                        "source_id": sfari_source_id,
                        "name": "SFARI Gene Score",
                        "score_type": "SFARI",
                        "value": float(gene_score),
                    }},
                )
                qid = res["insert_quantitative_scores"]["returning"][0]["id"]
                gql(LINK_GENE_QUANT, {"link": {"gene_ensg": ensg, "quantitative_score_id": qid}})

            # --- qualitative: SFARI Syndromic ---
            if syndromic is not None:
                res = gql(
                    UPSERT_QUAL_SCORE,
                    {"score": {
                        "source_id": sfari_source_id,
                        "name": "SFARI Syndromic",
                        "score_type": "SFARI",
                        "value": syndromic,
                    }},
                )
                qsid = res["insert_qualitative_scores"]["returning"][0]["id"]
                gql(LINK_GENE_QUAL, {"link": {"gene_ensg": ensg, "qualitative_score_id": qsid}})

            # --- quantitative: SFARI #Reports ---
            if n_reports is not None:
                res = gql(
                    UPSERT_QUANT_SCORE,
                    {"score": {
                        "source_id": sfari_source_id,
                        "name": "SFARI #Reports",
                        "score_type": "SFARI",
                        "value": float(n_reports),
                    }},
                )
                qid2 = res["insert_quantitative_scores"]["returning"][0]["id"]
                gql(LINK_GENE_QUANT, {"link": {"gene_ensg": ensg, "quantitative_score_id": qid2}})

            ok += 1

        except Exception as e:
            errs.append({"ensg": ensg, "error": str(e)})
            if len(errs) <= 10:
                logging.info(f"ERR {ensg} {e}")

    elapsed = time.time() - t0
    logging.info("\n✅ DONE SFARI load")
    logging.info(f"Processed rows: {n:,} in {elapsed/60:.1f} min")
    logging.info(f"ok: {ok:,}")
    logging.info(f"resolved_by_symbol: {resolved_by_symbol:,}")
    logging.info(f"skipped_no_ensg: {skipped_no_ensg:,}")
    logging.info(f"skipped_missing_gene: {skipped_missing_gene:,}")
    logging.info(f"skipped_no_values: {skipped_no_values:,}")
    logging.info(f"errors: {len(errs):,}")

    return errs


errs = load_sfari_gene_scores(path_sfari_genes, sfari_id)
if errs:
    logging.info("Some errors occurred:")
    for e in errs:
        logging.info(e)


Q = """
query VariantsWithVsWithoutSfariGeneScore($sfari_score_name: String!) {
  total: variants_aggregate { aggregate { count } }
  with_sfari: variants_aggregate(
    where: { gene: { genes_quantitative_scores: { quantitative_score: { name: { _eq: $sfari_score_name } } } } }
  ) { aggregate { count } }
  no_gene: variants_aggregate(where: { gene_ensg: { _is_null: true } }) { aggregate { count } }
}
"""

res = gql(Q, {"sfari_score_name": "SFARI Gene Score"})

total = res["total"]["aggregate"]["count"]
with_sfari = res["with_sfari"]["aggregate"]["count"]
no_gene = res["no_gene"]["aggregate"]["count"]

without_sfari = total - with_sfari  # includes null-gene variants too
without_sfari_but_has_gene = total - with_sfari - no_gene

logging.info(f"total: {total}")
logging.info(f"with_sfari: {with_sfari}")
logging.info(f"without_sfari (incl no_gene): {without_sfari}")
logging.info(f"no_gene: {no_gene}")
logging.info(f"without_sfari (but has gene): {without_sfari_but_has_gene}")



res = gql(Q, {"sfari_score_name": "SFARI Gene Score"})
without_sfari_but_has_gene = total - with_sfari - no_gene

without_sfari = total - with_sfari  # includes null-gene variants too

total = res["total"]["aggregate"]["count"]

with_sfari = res["with_sfari"]["aggregate"]["count"]
no_gene = res["no_gene"]["aggregate"]["count"]


# ═══════════════════════════════════════════════════════════════
# CADD score population (ported from cadd.ipynb)
# Queries VEP HGVS endpoint for every SNV variant and upserts
# CADD_raw + CADD_phred into variants_quantitative_scores.
# ═══════════════════════════════════════════════════════════════

logging.info("=" * 60)
logging.info("CADD SCORE POPULATION STARTED")
logging.info("=" * 60)

FETCH_VARIANTS_FOR_CADD = """
query FetchVariantsForCADD($offset: Int!, $limit: Int!) {
    variants(
        offset: $offset
        limit: $limit
        order_by: {id: asc}
        where: {
            class: {_eq: "SNV"}
        }
    ) {
        id
        ref
        alt
        coordinate {
            chr
            start
        }
    }
}
"""

UPSERT_CADD_SCORES = """
mutation UpsertCADDScores($objects: [variants_quantitative_scores_insert_input!]!) {
    insert_variants_quantitative_scores(
        objects: $objects
        on_conflict: {
            constraint: variants_quantitative_scores_pkey
            update_columns: [quantitative_score_id]
        }
    ) {
        affected_rows
    }
}
"""

CADD_VEP_SERVER = os.environ["GRCH37_VEP_SERVER"]
CADD_VEP_EXT = os.environ["VEP_HGVS_EXT"]
CADD_VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

logging.info(f"VEP Server: {CADD_VEP_SERVER}")
logging.info(f"VEP Endpoint: {CADD_VEP_EXT}")

def _extract_cadd(vep_result):
    """Return (cadd_raw, cadd_phred) from a single VEP result object, or (None, None)."""
    for key in ('transcript_consequences', 'regulatory_feature_consequences',
                'intergenic_consequences', 'motif_feature_consequences'):
        for c in vep_result.get(key, []):
            cadd_raw = c.get('cadd_raw')
            if cadd_raw is not None:
                return cadd_raw, c.get('cadd_phred')
    return None, None

def _fetch_hgvs_individually(variants_batch):
    """Query VEP one variant at a time — used when a batch gets HTTP 400."""
    results = {v["variant_id"]: {"cadd_raw": None, "cadd_phred": None} for v in variants_batch}
    for v in variants_batch:
        hgvs = v["hgvs_notation"]
        for attempt in range(5):
            try:
                r = requests.post(
                    CADD_VEP_SERVER + CADD_VEP_EXT,
                    headers=CADD_VEP_HEADERS,
                    data=json.dumps({"hgvs_notations": [hgvs]}),
                    timeout=60,
                )
                if r.status_code == 429:
                    time.sleep(2 ** attempt)
                    continue
                if r.status_code == 400:
                    logging.warning(f"  Individual 400 for {hgvs} — skipping")
                    break
                r.raise_for_status()
                cadd_raw, cadd_phred = _extract_cadd(r.json()[0])
                results[v["variant_id"]] = {"cadd_raw": cadd_raw, "cadd_phred": cadd_phred}
                time.sleep(0.1)
                break
            except Exception as e:
                if attempt == 4:
                    logging.warning(f"  Individual query failed for {hgvs}: {str(e)[:60]}")
                else:
                    time.sleep(2 ** attempt)
    return results

def fetch_vep_cadd_batch(variants_batch):
    """Fetch CADD scores for a batch using VEP HGVS endpoint. Falls back per-variant on HTTP 400."""
    hgvs_list = [v["hgvs_notation"] for v in variants_batch]
    id_by_hgvs = {v["hgvs_notation"]: v["variant_id"] for v in variants_batch}
    payload = json.dumps({"hgvs_notations": hgvs_list})

    max_retries = 5
    for attempt in range(max_retries):
        try:
            r = requests.post(CADD_VEP_SERVER + CADD_VEP_EXT, headers=CADD_VEP_HEADERS, data=payload, timeout=120)

            if r.status_code == 429:
                wait_time = 2 ** attempt
                logging.warning(f"Rate limited, waiting {wait_time}s (attempt {attempt+1}/{max_retries})")
                time.sleep(wait_time)
                continue

            if r.status_code == 400:
                logging.warning("Batch got HTTP 400 — falling back to individual queries")
                return _fetch_hgvs_individually(variants_batch)

            r.raise_for_status()

            results = {v["variant_id"]: {"cadd_raw": None, "cadd_phred": None} for v in variants_batch}
            for vep_result in r.json():
                input_str = vep_result.get("input", "")
                variant_id = id_by_hgvs.get(input_str)
                if variant_id is None:
                    continue
                cadd_raw, cadd_phred = _extract_cadd(vep_result)
                results[variant_id] = {"cadd_raw": cadd_raw, "cadd_phred": cadd_phred}

            time.sleep(0.5)
            return results

        except Exception as e:
            if attempt == max_retries - 1:
                logging.error(f"Batch failed after {max_retries} retries: {str(e)[:100]}")
                return {v["variant_id"]: {"cadd_raw": None, "cadd_phred": None} for v in variants_batch}
            else:
                logging.warning(f"Retry {attempt+1}/{max_retries} after error: {str(e)[:50]}")
                time.sleep(2 ** attempt)

    return {}

CADD_BATCH_SIZE = 200
CADD_INSERT_BATCH_SIZE = 1000
CADD_QUERY_LIMIT = 10000

logging.info(f"Batch size (VEP): {CADD_BATCH_SIZE}")
logging.info(f"Insert batch size: {CADD_INSERT_BATCH_SIZE}")
logging.info(f"Query limit: {CADD_QUERY_LIMIT}")
logging.info("=" * 60)

cadd_offset = 0
cadd_total_processed = 0
cadd_total_inserted = 0
cadd_failed_batches = []

try:
    while True:
        logging.info(f"{'─' * 60}")
        logging.info(f"Processing offset {cadd_offset} (limit {CADD_QUERY_LIMIT})")

        try:
            variants_data = gql(FETCH_VARIANTS_FOR_CADD, {"offset": cadd_offset, "limit": CADD_QUERY_LIMIT})
            cadd_variants = variants_data["variants"]
        except Exception as e:
            logging.error(f"Failed to fetch variants at offset {cadd_offset}: {e}")
            cadd_failed_batches.append({"offset": cadd_offset, "error": str(e), "type": "fetch"})
            break

        if not cadd_variants:
            logging.info("No more variants to process")
            break

        # Build HGVS notations — DB stores forward-strand alleles, no revcomp needed
        variants_with_hgvs = []
        for v in cadd_variants:
            chrom = v["coordinate"]["chr"]
            start = v["coordinate"]["start"]
            ref = v["ref"]
            alt = v["alt"]
            hgvs_notation = f"{chrom}:g.{start}{ref}>{alt}"
            variants_with_hgvs.append({
                "variant_id": v["id"],
                "hgvs_notation": hgvs_notation,
            })

        logging.info(f"Built {len(variants_with_hgvs)} HGVS notations")

        all_cadd_results = {}
        num_batches = (len(variants_with_hgvs) + CADD_BATCH_SIZE - 1) // CADD_BATCH_SIZE
        logging.info(f"Fetching CADD for {num_batches} VEP batches...")

        for i in range(0, len(variants_with_hgvs), CADD_BATCH_SIZE):
            batch = variants_with_hgvs[i:i + CADD_BATCH_SIZE]
            batch_num = (i // CADD_BATCH_SIZE) + 1
            logging.info(f"  VEP batch {batch_num}/{num_batches} ({len(batch)} variants)...")
            batch_results = fetch_vep_cadd_batch(batch)
            all_cadd_results.update(batch_results)

        score_inserts = []
        cadd_raw_count = 0
        cadd_phred_count = 0

        for variant_id, scores in all_cadd_results.items():
            if scores.get("cadd_raw") is not None:
                score_inserts.append({
                    "variant_id": variant_id,
                    "quantitative_score": {
                        "data": {
                            "source_id": SOURCE_VEP,
                            "name": "CADD_raw",
                            "score_type": "Pathogenicity",
                            "value": float(scores["cadd_raw"])
                        },
                        "on_conflict": {
                            "constraint": "quantitative_scores_source_id_score_type_name_value_key",
                            "update_columns": ["name"]
                        }
                    }
                })
                cadd_raw_count += 1

            if scores.get("cadd_phred") is not None:
                score_inserts.append({
                    "variant_id": variant_id,
                    "quantitative_score": {
                        "data": {
                            "source_id": SOURCE_VEP,
                            "name": "CADD_phred",
                            "score_type": "Pathogenicity",
                            "value": float(scores["cadd_phred"])
                        },
                        "on_conflict": {
                            "constraint": "quantitative_scores_source_id_score_type_name_value_key",
                            "update_columns": ["name"]
                        }
                    }
                })
                cadd_phred_count += 1

        logging.info(f"Prepared {len(score_inserts)} score inserts ({cadd_raw_count} raw, {cadd_phred_count} phred)")

        insert_success_count = 0
        insert_fail_count = 0

        for i in range(0, len(score_inserts), CADD_INSERT_BATCH_SIZE):
            chunk = score_inserts[i:i + CADD_INSERT_BATCH_SIZE]
            try:
                result = gql(UPSERT_CADD_SCORES, {"objects": chunk})
                affected = result["insert_variants_quantitative_scores"]["affected_rows"]
                cadd_total_inserted += affected
                insert_success_count += 1
                logging.info(f"  Insert batch {insert_success_count}: ✓ {affected} scores")
            except Exception as e:
                insert_fail_count += 1
                logging.error(f"  Insert batch failed: {str(e)[:100]}")
                cadd_failed_batches.append({
                    "offset": cadd_offset,
                    "insert_batch_start": i,
                    "error": str(e),
                    "type": "insert"
                })

        cadd_total_processed += len(cadd_variants)
        cadd_offset += CADD_QUERY_LIMIT

        logging.info(f"✅ Progress: {cadd_total_processed:,} processed, {cadd_total_inserted:,} inserted")

except KeyboardInterrupt:
    logging.warning("⚠️  CADD process interrupted by user")
except Exception as e:
    logging.error(f"⚠️  Unexpected error during CADD population: {e}", exc_info=True)

logging.info(f"{'=' * 60}")
logging.info("CADD SCORE POPULATION COMPLETE")
logging.info(f"{'=' * 60}")
logging.info(f"Total variants processed: {cadd_total_processed:,}")
logging.info(f"Total CADD scores inserted: {cadd_total_inserted:,}")
logging.info(f"Failed batches: {len(cadd_failed_batches)}")

if cadd_failed_batches:
    logging.warning("⚠️  FAILED BATCHES:")
    for i, batch in enumerate(cadd_failed_batches, 1):
        logging.warning(f"{i}. Type: {batch['type']}, Offset: {batch.get('offset')}, Error: {batch['error'][:100]}")

logging.info(f"{'=' * 60}")


# ── Step 1: Fetch the 123 SNV missense variants without SIFT from DB ──

FETCH_MISSING_SIFT = """
query FetchMissingSift {
  variants(where: {
    class: {_eq: "SNV"},
    cds_position: {_is_null: false},
    variants_consequences: {consequence: {name: {_eq: "Missense"}}},
    _not: {variants_quantitative_scores: {quantitative_score: {name: {_eq: "SIFT"}}}}
  }) {
    id
    ref
    alt
    coordinate {
      chr
      start
    }
    gene {
      symbol
    }
  }
}
"""

missing_sift_variants = gql(FETCH_MISSING_SIFT, {})["variants"]
logging.info(f"Found {len(missing_sift_variants)} SNV missense variants without SIFT score")

# ── Step 2: Query VEP API for SIFT scores ──

VEP_SERVER = "https://grch37.rest.ensembl.org"
VEP_HGVS_EXT = "/vep/human/hgvs"
VEP_SIFT_PARAMS = "?SIFT=b&pick=1"  # b = both prediction and score
VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
VEP_BATCH_SIZE = 200

def build_hgvs(chrom, pos, ref, alt):
    return f"{chrom}:g.{pos}{ref}>{alt}"

sift_results = {}  # variant_id -> {sift_prediction, sift_score}

for batch_start in range(0, len(missing_sift_variants), VEP_BATCH_SIZE):
    batch = missing_sift_variants[batch_start:batch_start + VEP_BATCH_SIZE]

    hgvs_notations = []
    id_by_hgvs = {}
    for v in batch:
        hgvs = build_hgvs(v["coordinate"]["chr"], v["coordinate"]["start"], v["ref"], v["alt"])
        hgvs_notations.append(hgvs)
        id_by_hgvs[hgvs] = v["id"]

    payload = json.dumps({"hgvs_notations": hgvs_notations})

    for attempt in range(5):
        try:
            r = requests.post(
                VEP_SERVER + VEP_HGVS_EXT + VEP_SIFT_PARAMS,
                headers=VEP_HEADERS,
                data=payload,
                timeout=120,
            )
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            break
        except Exception as e:
            if attempt == 4:
                logging.info(f"Batch starting at {batch_start} failed: {e}")
                continue
            time.sleep(2 ** attempt)

    vep_responses = r.json()

    for vep_result in vep_responses:
        input_hgvs = vep_result.get("input") or vep_result.get("id", "")
        variant_id = id_by_hgvs.get(input_hgvs)
        if variant_id is None:
            continue

        sift_prediction = None
        sift_score = None
        for key in ("transcript_consequences",):
            for cons in vep_result.get(key, []):
                if cons.get("sift_prediction") is not None:
                    sift_prediction = cons["sift_prediction"]
                    sift_score = cons.get("sift_score")
                    break
            if sift_prediction is not None:
                break

        sift_results[variant_id] = {
            "sift_prediction": sift_prediction,
            "sift_score": sift_score,
        }

    logging.info(f"  Processed {min(batch_start + VEP_BATCH_SIZE, len(missing_sift_variants))}/{len(missing_sift_variants)}")
    time.sleep(0.5)

# ── Step 3: Summarize results ──
has_sift = sum(1 for v in sift_results.values() if v["sift_score"] is not None)
no_sift = sum(1 for v in sift_results.values() if v["sift_score"] is None)

logging.info(f"\n=== SIFT Validation Results ===")
logging.info(f"Total variants checked: {len(sift_results)}")
logging.info(f"VEP returned SIFT:      {has_sift}")
logging.info(f"VEP confirmed no SIFT:  {no_sift}")

if has_sift > 0:
    logging.info(f"\n⚠️  {has_sift} variants DO have SIFT scores from VEP — these were missed during ingestion!")
    logging.info("\nVariants with SIFT scores from VEP:")
    for vid, scores in sift_results.items():
        if scores["sift_score"] is not None:
            v = next(x for x in missing_sift_variants if x["id"] == vid)
            logging.info(f"  id={vid} gene={v['gene']['symbol']} {v['coordinate']['chr']}:{v['coordinate']['start']} "
                  f"{v['ref']}>{v['alt']}  SIFT={scores['sift_prediction']}({scores['sift_score']})")
else:
    logging.info("\n✅ Confirmed: VEP also returns no SIFT for all 123 variants. No data was missed.")


# ═══════════════════════════════════════════════════════════════
# Phase 1: Baseline counts — variants table by class
# ═══════════════════════════════════════════════════════════════
data = gql("""
{
  total: variants_aggregate { aggregate { count } }
  snv: variants_aggregate(where: {class: {_eq: "SNV"}}) { aggregate { count } }
  deletion: variants_aggregate(where: {class: {_eq: "Deletion"}}) { aggregate { count } }
  insertion: variants_aggregate(where: {class: {_eq: "Insertion"}}) { aggregate { count } }
  substitution: variants_aggregate(where: {class: {_eq: "Substitution"}}) { aggregate { count } }
  indel: variants_aggregate(where: {class: {_eq: "Indel"}}) { aggregate { count } }
  seq_alt: variants_aggregate(where: {class: {_eq: "Sequence_alteration"}}) { aggregate { count } }
  null_class: variants_aggregate(where: {class: {_is_null: true}}) { aggregate { count } }
}
""", {})

logging.info("═══ Phase 1: Variants table — class distribution ═══")
total = data["total"]["aggregate"]["count"]
logging.info(f"\nTotal variants:          {total}")
logging.info(f"  SNV:                   {data['snv']['aggregate']['count']}")
logging.info(f"  Deletion:              {data['deletion']['aggregate']['count']}")
logging.info(f"  Insertion:             {data['insertion']['aggregate']['count']}")
logging.info(f"  Substitution:          {data['substitution']['aggregate']['count']}")
logging.info(f"  Indel:                 {data['indel']['aggregate']['count']}")
logging.info(f"  Sequence_alteration:   {data['seq_alt']['aggregate']['count']}")
null_count = data['null_class']['aggregate']['count']
logging.info(f"  NULL (unclassified):   {null_count}")

snv_count = data['snv']['aggregate']['count']
accounted = sum(data[k]['aggregate']['count'] for k in ['snv','deletion','insertion','substitution','indel','seq_alt','null_class'])
logging.info(f"\n  Accounted total:       {accounted}")
if accounted != total:
    logging.info(f"  ⚠️  Unaccounted:       {total - accounted}")

# ═══════════════════════════════════════════════════════════════
# Phase 2: Staging VEP — count by variant_class
# ═══════════════════════════════════════════════════════════════
staging = gql("""
{
  total: staging_vep_annotations_raw_aggregate { aggregate { count } }
  snv: staging_vep_annotations_raw_aggregate(where: {variant_class: {_eq: "SNV"}}) { aggregate { count } }
  null_class: staging_vep_annotations_raw_aggregate(where: {variant_class: {_is_null: true}}) { aggregate { count } }
}
""", {})

logging.info("\n═══ Phase 2: Staging VEP — row counts ═══")
logging.info(f"Total staging rows:      {staging['total']['aggregate']['count']}")
logging.info(f"  SNV rows:              {staging['snv']['aggregate']['count']}")
logging.info(f"  NULL class rows:       {staging['null_class']['aggregate']['count']}")

# ═══════════════════════════════════════════════════════════════
# Phase 3: True unique SNV count in staging (Python-side dedup)
# ═══════════════════════════════════════════════════════════════
logging.info("\n═══ Phase 3: True unique SNVs in staging (fetching...) ═══")
all_snv_ids = set()
offset = 0
batch_size = 50000
while True:
    batch = gql(f"""
    {{
      staging_vep_annotations_raw(
        where: {{variant_class: {{_eq: "SNV"}}}}
        limit: {batch_size}
        offset: {offset}
        order_by: {{id: asc}}
      ) {{
        uploaded_variation
      }}
    }}
    """, {})
    rows = batch["staging_vep_annotations_raw"]
    if not rows:
        break
    for r in rows:
        all_snv_ids.add(r["uploaded_variation"])
    offset += len(rows)
    logging.info(f"  Fetched {offset} rows, {len(all_snv_ids)} unique so far...")

logging.info(f"\nTrue unique SNV variants in staging: {len(all_snv_ids)}")
logging.info(f"SNV variants in variants table:      {snv_count}")
logging.info(f"Difference:                          {len(all_snv_ids) - snv_count}")

if null_count > 0:
    logging.info(f"\n⚠️  {null_count} variants have NULL class — these may be SNVs that weren't classified by VEP ingestion")
    logging.info(f"   NULL + SNV = {null_count + snv_count}")

logging.info("\n═══ Diagnosis ═══")
if snv_count + null_count >= len(all_snv_ids) * 0.95:
    logging.info("→ Most likely cause: Many variants have class=NULL because VEP ingestion didn't reach them.")
    logging.info("  The NULL-class variants are likely the missing SNVs. Re-running VEP ingestion should fix this.")
elif len(all_snv_ids) < 200000:
    logging.info("→ Staging has fewer unique SNVs than expected. The 282K count was inflated by duplicate uploaded_variation values.")
else:
    logging.info("→ Unexpected gap. Check VEP ingestion logs for errors or skipped batches.")



# ═══════════════════════════════════════════════════════════════════════════════
# Population allele frequencies  →  public.allele_frequencies  (one row per variant)
# ═══════════════════════════════════════════════════════════════════════════════
# VEP emits gnomAD (exome+genome) AFs + VEP's own MAX_AF in the `Extra` column (via
# --everything, see phase1_preprocess.py); those keys were expanded into columns of
# df_varicarta_annotated_expanded at the top of this phase. Here we map each VEP row
# to its variant_id and upsert one allele_frequencies row per variant that has any
# frequency value. Extraction/lookup/upsert live in allele_freq_common.py.
# Idempotent (on_conflict upserts every column).
import allele_freq_common as afc

logging.info("Populating allele_frequencies table from VEP annotations...")
af_lookup = afc.build_variant_lookup(gql, log=logging.info)
logging.info(f"  variant lookup built: {len(af_lookup)} variants")

af_objects, af_stats = afc.build_objects(
    df_varicarta_annotated_expanded.to_dict("records"), af_lookup, log=logging.info
)
logging.info(f"  allele_frequencies: {af_stats}")

af_inserted = afc.upsert(gql, af_objects, log=logging.info)
logging.info(f"✓ allele_frequencies populated: {af_inserted} rows upserted")
