--
-- PostgreSQL database dump
--

\restrict 7zFEJ1s4ZqFLoV3ZeieUXPeMX1jMav1Tar3WQmEEIllsslEDziTJk7xSKmfZDGt

-- Dumped from database version 16.11 (Debian 16.11-1.pgdg13+1)
-- Dumped by pg_dump version 16.11 (Debian 16.11-1.pgdg13+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: SCHEMA public; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON SCHEMA public IS '';


--
-- Name: pgcrypto; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pgcrypto WITH SCHEMA public;


--
-- Name: EXTENSION pgcrypto; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION pgcrypto IS 'cryptographic functions';


--
-- Name: amino_acid_full_name; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.amino_acid_full_name AS ENUM (
    'Alanine',
    'Arginine',
    'Asparagine',
    'Aspartic acid',
    'Cysteine',
    'Glutamic acid',
    'Glutamine',
    'Glycine',
    'Histidine',
    'Isoleucine',
    'Leucine',
    'Lysine',
    'Methionine',
    'Phenylalanine',
    'Proline',
    'Serine',
    'Threonine',
    'Tryptophan',
    'Tyrosine',
    'Valine'
);


--
-- Name: amino_acid_letter; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.amino_acid_letter AS ENUM (
    'A',
    'R',
    'N',
    'D',
    'C',
    'E',
    'Q',
    'G',
    'H',
    'I',
    'L',
    'K',
    'M',
    'F',
    'P',
    'S',
    'T',
    'W',
    'Y',
    'V'
);


--
-- Name: amino_acid_short_name; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.amino_acid_short_name AS ENUM (
    'Ala',
    'Arg',
    'Asn',
    'Asp',
    'Cys',
    'Glu',
    'Gln',
    'Gly',
    'His',
    'Ile',
    'Leu',
    'Lys',
    'Met',
    'Phe',
    'Pro',
    'Ser',
    'Thr',
    'Trp',
    'Tyr',
    'Val'
);


--
-- Name: chr; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.chr AS ENUM (
    '1',
    '2',
    '3',
    '4',
    '5',
    '6',
    '7',
    '8',
    '9',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18',
    '19',
    '20',
    '21',
    '22',
    'X',
    'Y',
    'MT'
);


--
-- Name: edit_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.edit_type AS ENUM (
    'Fix',
    'Rescue',
    'Improve'
);


--
-- Name: feature_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.feature_type AS ENUM (
    'Transcript',
    'RegulatoryFeature',
    'MotifFeature'
);


--
-- Name: genome_version; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.genome_version AS ENUM (
    'hg19',
    'hg38',
    'mm10'
);


--
-- Name: guide_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.guide_type AS ENUM (
    'pre_mRNA',
    'mature_mRNA'
);


--
-- Name: nucleotide; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.nucleotide AS ENUM (
    'A',
    'C',
    'G',
    'T',
    'U',
    'I'
);


--
-- Name: score_value_kind; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.score_value_kind AS ENUM (
    'quantitative',
    'qualitative',
    'both'
);


--
-- Name: seq_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.seq_type AS ENUM (
    'Protein',
    'DNA',
    'RNA'
);


--
-- Name: strand; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.strand AS ENUM (
    '+',
    '-'
);


--
-- Name: variant_class; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.variant_class AS ENUM (
    'SNV',
    'Deletion',
    'Insertion',
    'Substitution',
    'Indel',
    'Sequence_alteration'
);


--
-- Name: is_ga_variant(integer); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.is_ga_variant(variant_id_param integer) RETURNS boolean
    LANGUAGE plpgsql STABLE
    AS $$
DECLARE
    v_ref character varying;
    v_alt character varying;
    v_gene_strand public.strand;
BEGIN
    -- Get variant ref, alt, and gene strand
    SELECT v.ref, v.alt, g.strand
    INTO v_ref, v_alt, v_gene_strand
    FROM public.variants v
    LEFT JOIN public.genes g ON g.ensg = v.gene_ensg
    WHERE v.id = variant_id_param;
    
    -- Check if it's a G>A variant
    -- Case 1: G>A on any strand
    IF v_ref = 'G' AND v_alt = 'A' THEN
        RETURN TRUE;
    END IF;
    
    -- Case 2: C>T on minus strand gene (equivalent to G>A on coding strand)
    IF v_ref = 'C' AND v_alt = 'T' AND v_gene_strand = '-' THEN
        RETURN TRUE;
    END IF;
    
    -- Not a G>A variant
    RETURN FALSE;
END;
$$;


--
-- Name: is_not_ga_variant(integer); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.is_not_ga_variant(editing_context_id_param integer) RETURNS boolean
    LANGUAGE plpgsql STABLE
    AS $$
DECLARE
    v_ref character varying;
    v_alt character varying;
    v_gene_strand public.strand;
BEGIN
    -- Get variant ref, alt, and gene strand via editing_context
    SELECT v.ref, v.alt, g.strand
    INTO v_ref, v_alt, v_gene_strand
    FROM public.editing_contexts ec
    JOIN public.variants v ON v.id = ec.variant_id
    LEFT JOIN public.genes g ON g.ensg = v.gene_ensg
    WHERE ec.id = editing_context_id_param;
    
    -- Check if it's a G>A variant on the coding strand
    -- Case 1: G>A on plus strand or no strand info
    IF v_ref = 'G' AND v_alt = 'A' AND (v_gene_strand = '+' OR v_gene_strand IS NULL) THEN
        RETURN FALSE;
    END IF;
    
    -- Case 2: C>T on minus strand (equivalent to G>A on coding strand)
    IF v_ref = 'C' AND v_alt = 'T' AND v_gene_strand = '-' THEN
        RETURN FALSE;
    END IF;
    
    -- Not a G>A variant on coding strand
    RETURN TRUE;
END;
$$;


--
-- Name: set_has_adjacent_u(); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.set_has_adjacent_u() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
BEGIN
  NEW.has_adjacent_u :=
    COALESCE(NEW.has_u_before, false) OR COALESCE(NEW.has_u_after, false);
  RETURN NEW;
END;
$$;


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: amino_acids; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.amino_acids (
    id integer NOT NULL,
    full_name public.amino_acid_full_name NOT NULL,
    short_name public.amino_acid_short_name NOT NULL,
    letter public.amino_acid_letter NOT NULL
);


--
-- Name: amino_acids_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.amino_acids ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.amino_acids_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: biotypes; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.biotypes (
    id integer NOT NULL,
    name character varying NOT NULL
);


--
-- Name: biotypes_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.biotypes ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.biotypes_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: bystanders; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.bystanders (
    id integer NOT NULL,
    cds_position integer,
    am_pathogenicity double precision,
    am_class character varying,
    sift_score double precision,
    sift_prediction character varying,
    ref_codon_id integer,
    alt_codon_id integer,
    cadd_raw double precision,
    cadd_phred double precision,
    has_u_before boolean NOT NULL,
    has_u_after boolean NOT NULL,
    has_adjacent_u boolean GENERATED ALWAYS AS ((COALESCE(has_u_before, false) OR COALESCE(has_u_after, false))) STORED,
    guide_id integer NOT NULL,
    genomic_position integer,
    in_cds boolean NOT NULL,
    variant_id integer NOT NULL
);


--
-- Name: bystanders_consequences; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.bystanders_consequences (
    bystander_id integer NOT NULL,
    consequence_id integer NOT NULL
);


--
-- Name: bystanders_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.bystanders ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.bystanders_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: codons; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.codons (
    id integer NOT NULL,
    amino_acid_id integer,
    is_stop boolean NOT NULL,
    nt1 public.nucleotide NOT NULL,
    nt2 public.nucleotide NOT NULL,
    nt3 public.nucleotide NOT NULL,
    CONSTRAINT codons_stop_consistency CHECK (((is_stop AND (amino_acid_id IS NULL)) OR ((NOT is_stop) AND (amino_acid_id IS NOT NULL))))
);


--
-- Name: codons_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.codons_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: codons_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.codons_id_seq OWNED BY public.codons.id;


--
-- Name: consequences; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.consequences (
    id integer NOT NULL,
    name character varying NOT NULL
);


--
-- Name: consequences_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.consequences ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.consequences_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: coordinates; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.coordinates (
    id integer NOT NULL,
    chr public.chr NOT NULL,
    start integer NOT NULL,
    "end" integer NOT NULL,
    strand public.strand,
    genome_build_id integer NOT NULL
);


--
-- Name: coordinates_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.coordinates ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.coordinates_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: exons; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.exons (
    id integer NOT NULL,
    exon_position_in_transcript integer,
    sequence_id integer,
    transcript_enst character varying NOT NULL
);


--
-- Name: features; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.features (
    id integer NOT NULL,
    type public.feature_type NOT NULL,
    identifier character varying NOT NULL,
    version_number integer DEFAULT 0 NOT NULL,
    biotype_id integer
);


--
-- Name: features_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.features ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.features_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: gene_expressions; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_expressions (
    id integer NOT NULL,
    tissue_id integer NOT NULL,
    median_tpm double precision NOT NULL,
    source_id integer NOT NULL,
    mean_tpm double precision NOT NULL,
    gene_ensg character varying NOT NULL
);


--
-- Name: gene_expressions_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.gene_expressions ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.gene_expressions_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: genes; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genes (
    symbol character varying,
    ensg character varying NOT NULL,
    strand public.strand NOT NULL
);


--
-- Name: genes_annotations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genes_annotations (
    id integer NOT NULL,
    source_id integer NOT NULL,
    version_number integer DEFAULT 0 NOT NULL,
    identifier text NOT NULL,
    gene_ensg character varying NOT NULL
);


--
-- Name: genes_annotations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.genes_annotations ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.genes_annotations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: genes_qualitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genes_qualitative_scores (
    gene_ensg character varying NOT NULL,
    qualitative_score_id integer NOT NULL
);


--
-- Name: genes_quantitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genes_quantitative_scores (
    gene_ensg character varying NOT NULL,
    quantitative_score_id integer NOT NULL
);


--
-- Name: genome_builds; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genome_builds (
    id integer NOT NULL,
    organism_id integer NOT NULL,
    build_name text NOT NULL,
    assembly_name text NOT NULL,
    assembly_accession text,
    release_date date,
    notes text
);


--
-- Name: genome_builds_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.genome_builds ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.genome_builds_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: guides; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.guides (
    id integer NOT NULL,
    hits_85 integer,
    hits_90 integer,
    hits_95 integer,
    hits_100 integer,
    sequence text NOT NULL,
    type public.guide_type NOT NULL
);


--
-- Name: guides_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.guides ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.guides_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: introns; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.introns (
    intron_position_in_transcript integer,
    sequence_id integer,
    id integer NOT NULL,
    transcript_enst character varying NOT NULL
);


--
-- Name: introns_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.introns_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: introns_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.introns_id_seq OWNED BY public.introns.id;


--
-- Name: organisms; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.organisms (
    id integer NOT NULL,
    taxon_id integer NOT NULL,
    scientific_name text NOT NULL,
    common_name text NOT NULL,
    abbrev text,
    default_genome_build_id integer
);


--
-- Name: organisms_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.organisms ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.organisms_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: organs; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.organs (
    name character varying NOT NULL,
    id integer NOT NULL
);


--
-- Name: organs_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.organs_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: organs_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.organs_id_seq OWNED BY public.organs.id;


--
-- Name: papers; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.papers (
    paper_key character varying NOT NULL,
    id integer NOT NULL,
    external_paper_id character varying NOT NULL
);


--
-- Name: papers_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.papers_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: papers_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.papers_id_seq OWNED BY public.papers.id;


--
-- Name: papers_id_seq1; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.papers ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.papers_id_seq1
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: papers_subjects; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.papers_subjects (
    subject_id integer NOT NULL,
    paper_id integer NOT NULL
);


--
-- Name: populations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.populations (
    id integer NOT NULL,
    code character varying NOT NULL,
    description character varying
);


--
-- Name: populations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.populations ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.populations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: protein_domains; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.protein_domains (
    id integer NOT NULL,
    source_id integer,
    accession_number character varying NOT NULL
);


--
-- Name: protein_domains_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.protein_domains ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.protein_domains_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: protein_domains_proteins; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.protein_domains_proteins (
    protein_domain_id integer NOT NULL,
    protein_ensp character varying NOT NULL
);


--
-- Name: proteins; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.proteins (
    ensp character varying NOT NULL
);


--
-- Name: proteins_annotations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.proteins_annotations (
    source_id integer NOT NULL,
    version_number integer DEFAULT 0 NOT NULL,
    identifier character varying NOT NULL,
    protein_ensp character varying NOT NULL,
    id bigint NOT NULL
);


--
-- Name: proteins_annotations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.proteins_annotations ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.proteins_annotations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: qualitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.qualitative_scores (
    id integer NOT NULL,
    name text NOT NULL,
    score_type text NOT NULL,
    value text NOT NULL,
    source_id integer,
    population_id integer
);


--
-- Name: qualitative_scores_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.qualitative_scores ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.qualitative_scores_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: quantitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.quantitative_scores (
    id integer NOT NULL,
    source_id integer,
    population_id integer,
    name text NOT NULL,
    score_type text NOT NULL,
    value double precision NOT NULL
);


--
-- Name: quantitative_scores_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.quantitative_scores ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.quantitative_scores_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: sources; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.sources (
    id integer NOT NULL,
    name character varying NOT NULL,
    url text,
    description text,
    release integer DEFAULT 0 NOT NULL
);


--
-- Name: sources_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.sources ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.sources_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: staging_varicarta_variants_raw; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.staging_varicarta_variants_raw (
    id bigint NOT NULL,
    chrom text,
    pos text,
    varicarta_vcf_id text,
    ref text,
    alt text,
    qual text,
    filter text,
    info text,
    cadd13_phred text,
    cadd13_raw text,
    cadd_phred text,
    cadd_raw text,
    fathmm_pred text,
    fathmm_score text,
    gerp_rs text,
    lrt_pred text,
    lrt_score text,
    lr_pred text,
    lr_score text,
    mutationassessor_pred text,
    mutationassessor_score text,
    mutationtaster_pred text,
    mutationtaster_score text,
    polyphen2_hdiv_pred text,
    polyphen2_hdiv_score text,
    polyphen2_hvar_pred text,
    polyphen2_hvar_score text,
    radialsvm_pred text,
    radialsvm_score text,
    sift_pred text,
    sift_score text,
    siphy_29way_logodds text,
    vest3_score text,
    aa_change text,
    category text,
    clinvar_20150629 text,
    code_change text,
    cytoband text,
    exac03 text,
    func text,
    gene_detail text,
    gene_symbol text,
    varicarta_id text,
    inheritance text,
    paper_id text,
    paper_key text,
    phylop100way_vertebrate text,
    phylop46way_placental text,
    pid text,
    protein_change text,
    sample_id text,
    sequencing_study_type text,
    stop_hg19 text,
    subject_id text,
    validation text,
    validation_method text,
    validation_reported text
);


--
-- Name: staging_varicarta_variants_raw_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.staging_varicarta_variants_raw ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.staging_varicarta_variants_raw_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: staging_vep_annotations_raw; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.staging_vep_annotations_raw (
    id bigint NOT NULL,
    uploaded_variation text,
    location text,
    allele text,
    gene text,
    feature text,
    feature_type text,
    consequence text,
    cdna_position text,
    cds_position text,
    protein_position text,
    amino_acids text,
    codons text,
    existing_variation text,
    extra text,
    impact text,
    strand text,
    variant_class text,
    symbol text,
    symbol_source text,
    hgnc_id text,
    biotype text,
    canonical text,
    exon text,
    hgvsc text,
    gnomade_af text,
    gnomade_afr_af text,
    gnomade_amr_af text,
    gnomade_asj_af text,
    gnomade_eas_af text,
    gnomade_fin_af text,
    gnomade_nfe_af text,
    gnomade_sas_af text,
    max_af text,
    max_af_pops text,
    phastcons100 text,
    intron text,
    distance text,
    hgvs_offset text,
    minimised text,
    af text,
    afr_af text,
    amr_af text,
    eas_af text,
    eur_af text,
    sas_af text,
    enformer_sad text,
    enformer_sar text,
    ensp text,
    uniparc text,
    ccds text,
    swissprot text,
    trembl text,
    gene_pheno text,
    pli_gene_value text,
    domains text,
    hgvsp text,
    clin_sig text,
    pheno text,
    lovd text,
    somatic text,
    sift text,
    polyphen text,
    blosum62 text,
    condel text,
    revel text,
    am_class text,
    am_genome text,
    am_pathogenicity text,
    am_protein_variant text,
    am_transcript_id text,
    am_uniprot_id text,
    nmd text,
    loftool text,
    pubmed text,
    motif_name text,
    motif_pos text,
    high_inf_pos text,
    motif_score_change text,
    transcription_factors text,
    flags text,
    flaglrg text,
    mirna text
);


--
-- Name: staging_vep_annotations_raw_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.staging_vep_annotations_raw ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.staging_vep_annotations_raw_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: subjects; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.subjects (
    subject_key text NOT NULL,
    id integer NOT NULL
);


--
-- Name: subjects_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.subjects ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.subjects_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: tissues; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.tissues (
    id integer NOT NULL,
    name character varying NOT NULL,
    organ_id integer NOT NULL
);


--
-- Name: tissues_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.tissues ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.tissues_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: transcripts; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.transcripts (
    sequence_id integer,
    ensembl_canonical boolean,
    mane_select boolean,
    enst character varying NOT NULL,
    gene_ensg character varying,
    protein_ensp character varying
);


--
-- Name: transcripts_annotations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.transcripts_annotations (
    id integer NOT NULL,
    source_id integer NOT NULL,
    version_number integer DEFAULT 0 NOT NULL,
    identifier character varying NOT NULL,
    transcript_enst character varying NOT NULL
);


--
-- Name: transcripts_annotations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.transcripts_annotations ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.transcripts_annotations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: variants; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants (
    id integer NOT NULL,
    coordinate_id integer NOT NULL,
    class public.variant_class,
    nmd_escaping_variant boolean,
    created_at timestamp without time zone,
    ref character varying NOT NULL,
    alt character varying NOT NULL,
    ref_validated boolean DEFAULT false NOT NULL,
    varicarta_vcf_id integer,
    varicarta_internal_id integer,
    inheritance character varying,
    rs_id character varying,
    ref_codon_id integer,
    alt_codon_id integer,
    position_in_codon integer,
    cdna_position character varying,
    cds_position character varying,
    protein_position character varying,
    exon character varying,
    intron character varying,
    gene_ensg character varying,
    CONSTRAINT check_codon_position CHECK (((position_in_codon IS NULL) OR (position_in_codon = ANY (ARRAY[1, 2, 3])))),
    CONSTRAINT variants_ref_validated_chk CHECK ((ref_validated = true))
);


--
-- Name: variants_consequences; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_consequences (
    consequence_id integer NOT NULL,
    variant_id integer NOT NULL
);


--
-- Name: variants_features; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_features (
    variant_id integer NOT NULL,
    feature_id integer NOT NULL
);


--
-- Name: variants_guides; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_guides (
    variant_id integer NOT NULL,
    guide_id integer NOT NULL,
    variant_idx integer NOT NULL,
    edit_type public.edit_type,
    rescue_sift double precision
);


--
-- Name: COLUMN variants_guides.variant_idx; Type: COMMENT; Schema: public; Owner: -
--

COMMENT ON COLUMN public.variants_guides.variant_idx IS '1-based position of the variant within the guide sequence.
Typically 21 for n=20 extraction (21st nucleotide), but varies near boundaries.';


--
-- Name: variants_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

ALTER TABLE public.variants ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.variants_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);


--
-- Name: variants_papers; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_papers (
    variant_id integer NOT NULL,
    paper_id integer NOT NULL
);


--
-- Name: variants_qualitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_qualitative_scores (
    variant_id integer NOT NULL,
    qualitative_score_id integer NOT NULL
);


--
-- Name: variants_quantitative_scores; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_quantitative_scores (
    variant_id integer NOT NULL,
    quantitative_score_id integer NOT NULL
);


--
-- Name: variants_sources; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_sources (
    variant_id integer NOT NULL,
    source_id integer NOT NULL
);


--
-- Name: variants_subjects; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variants_subjects (
    variant_id integer NOT NULL,
    subject_id integer NOT NULL
);


--
-- Name: codons id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.codons ALTER COLUMN id SET DEFAULT nextval('public.codons_id_seq'::regclass);


--
-- Name: introns id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.introns ALTER COLUMN id SET DEFAULT nextval('public.introns_id_seq'::regclass);


--
-- Name: organs id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organs ALTER COLUMN id SET DEFAULT nextval('public.organs_id_seq'::regclass);


--
-- Name: amino_acids amino_acids_full_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amino_acids
    ADD CONSTRAINT amino_acids_full_name_key UNIQUE (full_name);


--
-- Name: amino_acids amino_acids_letter_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amino_acids
    ADD CONSTRAINT amino_acids_letter_key UNIQUE (letter);


--
-- Name: amino_acids amino_acids_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amino_acids
    ADD CONSTRAINT amino_acids_pkey PRIMARY KEY (id);


--
-- Name: amino_acids amino_acids_short_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amino_acids
    ADD CONSTRAINT amino_acids_short_name_key UNIQUE (short_name);


--
-- Name: biotypes biotypes_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.biotypes
    ADD CONSTRAINT biotypes_name_key UNIQUE (name);


--
-- Name: biotypes biotypes_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.biotypes
    ADD CONSTRAINT biotypes_pkey PRIMARY KEY (id);


--
-- Name: bystanders_consequences bystanders_consequences_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders_consequences
    ADD CONSTRAINT bystanders_consequences_pkey PRIMARY KEY (bystander_id, consequence_id);


--
-- Name: bystanders bystanders_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders
    ADD CONSTRAINT bystanders_pkey PRIMARY KEY (id);


--
-- Name: codons codons_nt1_nt2_nt3_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.codons
    ADD CONSTRAINT codons_nt1_nt2_nt3_key UNIQUE (nt1, nt2, nt3);


--
-- Name: codons codons_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.codons
    ADD CONSTRAINT codons_pkey PRIMARY KEY (id);


--
-- Name: consequences consequences_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.consequences
    ADD CONSTRAINT consequences_name_key UNIQUE (name);


--
-- Name: consequences consequences_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.consequences
    ADD CONSTRAINT consequences_pkey PRIMARY KEY (id);


--
-- Name: coordinates coordinates_genome_build_id_chr_start_end_strand_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.coordinates
    ADD CONSTRAINT coordinates_genome_build_id_chr_start_end_strand_key UNIQUE (genome_build_id, chr, start, "end", strand);


--
-- Name: coordinates coordinates_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.coordinates
    ADD CONSTRAINT coordinates_pkey PRIMARY KEY (id);


--
-- Name: exons exons_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.exons
    ADD CONSTRAINT exons_pkey PRIMARY KEY (id);


--
-- Name: features features_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_pkey PRIMARY KEY (id);


--
-- Name: features features_type_identifier_version_number_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_type_identifier_version_number_key UNIQUE (type, identifier, version_number);


--
-- Name: gene_expressions gene_expressions_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_expressions
    ADD CONSTRAINT gene_expressions_pkey PRIMARY KEY (id);


--
-- Name: gene_expressions gene_expressions_tissue_id_gene_ensg_source_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_expressions
    ADD CONSTRAINT gene_expressions_tissue_id_gene_ensg_source_id_key UNIQUE (tissue_id, gene_ensg, source_id);


--
-- Name: genes_annotations genes_annotations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_annotations
    ADD CONSTRAINT genes_annotations_pkey PRIMARY KEY (id);


--
-- Name: genes_annotations genes_annotations_source_id_identifier_version_number_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_annotations
    ADD CONSTRAINT genes_annotations_source_id_identifier_version_number_key UNIQUE (source_id, identifier, version_number);


--
-- Name: genes genes_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_pkey PRIMARY KEY (ensg);


--
-- Name: genes_qualitative_scores genes_qualitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_qualitative_scores
    ADD CONSTRAINT genes_qualitative_scores_pkey PRIMARY KEY (gene_ensg, qualitative_score_id);


--
-- Name: genes_quantitative_scores genes_quantitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_quantitative_scores
    ADD CONSTRAINT genes_quantitative_scores_pkey PRIMARY KEY (gene_ensg, quantitative_score_id);


--
-- Name: genome_builds genome_builds_organism_id_build_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_builds
    ADD CONSTRAINT genome_builds_organism_id_build_name_key UNIQUE (organism_id, build_name);


--
-- Name: genome_builds genome_builds_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_builds
    ADD CONSTRAINT genome_builds_pkey PRIMARY KEY (id);


--
-- Name: guides guides_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.guides
    ADD CONSTRAINT guides_pkey PRIMARY KEY (id);


--
-- Name: guides guides_sequence_type_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.guides
    ADD CONSTRAINT guides_sequence_type_key UNIQUE (sequence, type);


--
-- Name: introns introns_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.introns
    ADD CONSTRAINT introns_id_key UNIQUE (id);


--
-- Name: introns introns_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.introns
    ADD CONSTRAINT introns_pkey PRIMARY KEY (id);


--
-- Name: organisms organisms_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms
    ADD CONSTRAINT organisms_pkey PRIMARY KEY (id);


--
-- Name: organisms organisms_taxon_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms
    ADD CONSTRAINT organisms_taxon_id_key UNIQUE (taxon_id);


--
-- Name: organs organs_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organs
    ADD CONSTRAINT organs_id_key UNIQUE (id);


--
-- Name: organs organs_name_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organs
    ADD CONSTRAINT organs_name_key UNIQUE (name);


--
-- Name: organs organs_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organs
    ADD CONSTRAINT organs_pkey PRIMARY KEY (id);


--
-- Name: papers papers_external_paper_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers
    ADD CONSTRAINT papers_external_paper_id_key UNIQUE (external_paper_id);


--
-- Name: papers papers_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers
    ADD CONSTRAINT papers_id_key UNIQUE (id);


--
-- Name: papers papers_paper_key_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers
    ADD CONSTRAINT papers_paper_key_key UNIQUE (paper_key);


--
-- Name: papers papers_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers
    ADD CONSTRAINT papers_pkey PRIMARY KEY (id);


--
-- Name: papers_subjects papers_subjects_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers_subjects
    ADD CONSTRAINT papers_subjects_pkey PRIMARY KEY (subject_id, paper_id);


--
-- Name: populations populations_code_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.populations
    ADD CONSTRAINT populations_code_key UNIQUE (code);


--
-- Name: populations populations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.populations
    ADD CONSTRAINT populations_pkey PRIMARY KEY (id);


--
-- Name: protein_domains protein_domains_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains
    ADD CONSTRAINT protein_domains_pkey PRIMARY KEY (id);


--
-- Name: protein_domains_proteins protein_domains_proteins_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains_proteins
    ADD CONSTRAINT protein_domains_proteins_pkey PRIMARY KEY (protein_ensp, protein_domain_id);


--
-- Name: protein_domains protein_domains_source_id_accession_number_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains
    ADD CONSTRAINT protein_domains_source_id_accession_number_key UNIQUE (source_id, accession_number);


--
-- Name: proteins_annotations proteins_annotations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proteins_annotations
    ADD CONSTRAINT proteins_annotations_pkey PRIMARY KEY (id);


--
-- Name: proteins_annotations proteins_annotations_source_id_identifier_version_number_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proteins_annotations
    ADD CONSTRAINT proteins_annotations_source_id_identifier_version_number_key UNIQUE (source_id, identifier, version_number);


--
-- Name: proteins proteins_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proteins
    ADD CONSTRAINT proteins_pkey PRIMARY KEY (ensp);


--
-- Name: qualitative_scores qualitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.qualitative_scores
    ADD CONSTRAINT qualitative_scores_pkey PRIMARY KEY (id);


--
-- Name: qualitative_scores qualitative_scores_source_id_score_type_name_value_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.qualitative_scores
    ADD CONSTRAINT qualitative_scores_source_id_score_type_name_value_key UNIQUE (source_id, score_type, name, value);


--
-- Name: quantitative_scores quantitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.quantitative_scores
    ADD CONSTRAINT quantitative_scores_pkey PRIMARY KEY (id);


--
-- Name: quantitative_scores quantitative_scores_source_id_score_type_name_value_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.quantitative_scores
    ADD CONSTRAINT quantitative_scores_source_id_score_type_name_value_key UNIQUE (source_id, score_type, name, value);


--
-- Name: sources sources_name_release_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.sources
    ADD CONSTRAINT sources_name_release_key UNIQUE (name, release);


--
-- Name: sources sources_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.sources
    ADD CONSTRAINT sources_pkey PRIMARY KEY (id);


--
-- Name: staging_varicarta_variants_raw staging_varicarta_variants_raw_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.staging_varicarta_variants_raw
    ADD CONSTRAINT staging_varicarta_variants_raw_pkey PRIMARY KEY (id);


--
-- Name: staging_vep_annotations_raw staging_vep_annotations_raw_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.staging_vep_annotations_raw
    ADD CONSTRAINT staging_vep_annotations_raw_pkey PRIMARY KEY (id);


--
-- Name: subjects subjects_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subjects
    ADD CONSTRAINT subjects_pkey PRIMARY KEY (id);


--
-- Name: subjects subjects_subject_key_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subjects
    ADD CONSTRAINT subjects_subject_key_key UNIQUE (subject_key);


--
-- Name: tissues tissues_name_organ_id_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tissues
    ADD CONSTRAINT tissues_name_organ_id_key UNIQUE (name, organ_id);


--
-- Name: tissues tissues_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tissues
    ADD CONSTRAINT tissues_pkey PRIMARY KEY (id);


--
-- Name: transcripts_annotations transcripts_annotations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts_annotations
    ADD CONSTRAINT transcripts_annotations_pkey PRIMARY KEY (id);


--
-- Name: transcripts_annotations transcripts_annotations_source_id_identifier_version_number_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts_annotations
    ADD CONSTRAINT transcripts_annotations_source_id_identifier_version_number_key UNIQUE (source_id, identifier, version_number);


--
-- Name: transcripts transcripts_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts
    ADD CONSTRAINT transcripts_pkey PRIMARY KEY (enst);


--
-- Name: variants_consequences variants_consequences_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_consequences
    ADD CONSTRAINT variants_consequences_pkey PRIMARY KEY (consequence_id, variant_id);


--
-- Name: variants variants_coordinate_id_ref_alt_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_coordinate_id_ref_alt_key UNIQUE (coordinate_id, ref, alt);


--
-- Name: variants_features variants_features_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_features
    ADD CONSTRAINT variants_features_pkey PRIMARY KEY (variant_id, feature_id);


--
-- Name: variants_guides variants_guides_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_guides
    ADD CONSTRAINT variants_guides_pkey PRIMARY KEY (variant_id, guide_id);


--
-- Name: variants_papers variants_papers_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_papers
    ADD CONSTRAINT variants_papers_pkey PRIMARY KEY (variant_id, paper_id);


--
-- Name: variants variants_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_pkey PRIMARY KEY (id);


--
-- Name: variants_qualitative_scores variants_qualitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_qualitative_scores
    ADD CONSTRAINT variants_qualitative_scores_pkey PRIMARY KEY (variant_id, qualitative_score_id);


--
-- Name: variants_quantitative_scores variants_quantitative_scores_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_quantitative_scores
    ADD CONSTRAINT variants_quantitative_scores_pkey PRIMARY KEY (variant_id, quantitative_score_id);


--
-- Name: variants_sources variants_sources_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_sources
    ADD CONSTRAINT variants_sources_pkey PRIMARY KEY (variant_id, source_id);


--
-- Name: variants_subjects variants_subjects_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_subjects
    ADD CONSTRAINT variants_subjects_pkey PRIMARY KEY (variant_id, subject_id);


--
-- Name: amino_acid_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX amino_acid_id_idx ON public.codons USING btree (amino_acid_id);


--
-- Name: bystanders_guide_genomic_pos_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX bystanders_guide_genomic_pos_idx ON public.bystanders USING btree (guide_id, genomic_position);


--
-- Name: bystanders_guide_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX bystanders_guide_id_idx ON public.bystanders USING btree (guide_id);


--
-- Name: bystanders_in_cds_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX bystanders_in_cds_idx ON public.bystanders USING btree (in_cds);


--
-- Name: bystanders_variant_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX bystanders_variant_id_idx ON public.bystanders USING btree (variant_id);


--
-- Name: common_name_key; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX common_name_key ON public.organisms USING btree (common_name);


--
-- Name: coordinates_genome_build_chr_start_end_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX coordinates_genome_build_chr_start_end_idx ON public.coordinates USING btree (genome_build_id, chr, start, "end");


--
-- Name: coordinates_genome_build_chr_start_end_strand_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX coordinates_genome_build_chr_start_end_strand_idx ON public.coordinates USING btree (genome_build_id, chr, start, "end", strand);


--
-- Name: coordinates_genome_build_chr_start_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX coordinates_genome_build_chr_start_idx ON public.coordinates USING btree (genome_build_id, chr, start);


--
-- Name: features_type_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX features_type_idx ON public.features USING btree (type);


--
-- Name: gene_expressions_gene_ensg_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_expressions_gene_ensg_idx ON public.gene_expressions USING btree (gene_ensg);


--
-- Name: gene_expressions_tissue_id_gene_ensg_key; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_expressions_tissue_id_gene_ensg_key ON public.gene_expressions USING btree (tissue_id, gene_ensg);


--
-- Name: gene_expressions_tissue_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_expressions_tissue_id_idx ON public.gene_expressions USING btree (tissue_id);


--
-- Name: genes_annotations_gene_ensg_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX genes_annotations_gene_ensg_idx ON public.genes_annotations USING btree (gene_ensg);


--
-- Name: genes_annotations_source_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX genes_annotations_source_id_idx ON public.genes_annotations USING btree (source_id);


--
-- Name: genes_symbol; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX genes_symbol ON public.genes USING btree (symbol);


--
-- Name: nt1_nt2_nt3_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX nt1_nt2_nt3_idx ON public.codons USING btree (nt1, nt2, nt3);


--
-- Name: organs_name_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX organs_name_idx ON public.organs USING btree (name);


--
-- Name: papers_subjects_paper_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX papers_subjects_paper_id_idx ON public.papers_subjects USING btree (paper_id);


--
-- Name: papers_subjects_subject_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX papers_subjects_subject_id_idx ON public.papers_subjects USING btree (subject_id);


--
-- Name: papers_subjects_subject_id_paper_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX papers_subjects_subject_id_paper_id_idx ON public.papers_subjects USING btree (subject_id, paper_id);


--
-- Name: protein_domains_proteins_protein_domain_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX protein_domains_proteins_protein_domain_id_idx ON public.protein_domains_proteins USING btree (protein_domain_id);


--
-- Name: protein_domains_proteins_protein_ensp_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX protein_domains_proteins_protein_ensp_idx ON public.protein_domains_proteins USING btree (protein_ensp);


--
-- Name: scientific_name_key; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX scientific_name_key ON public.organisms USING btree (scientific_name);


--
-- Name: sequence_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX sequence_id_idx ON public.exons USING btree (sequence_id);


--
-- Name: source_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX source_id ON public.proteins_annotations USING btree (source_id);


--
-- Name: tissues_organ_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX tissues_organ_id_idx ON public.tissues USING btree (organ_id);


--
-- Name: transcripts_annotations_transcript_enst_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX transcripts_annotations_transcript_enst_idx ON public.transcripts_annotations USING btree (transcript_enst);


--
-- Name: transcripts_gene_ensg_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX transcripts_gene_ensg_idx ON public.transcripts USING btree (gene_ensg);


--
-- Name: transcripts_protein_ensp_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX transcripts_protein_ensp_idx ON public.transcripts USING btree (protein_ensp);


--
-- Name: transcripts_sequence_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX transcripts_sequence_id_idx ON public.transcripts USING btree (sequence_id);


--
-- Name: variants_class_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_class_idx ON public.variants USING btree (class);


--
-- Name: variants_consequences_consequence_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_consequences_consequence_id_idx ON public.variants_consequences USING btree (consequence_id);


--
-- Name: variants_coordinate_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_coordinate_id_idx ON public.variants USING btree (coordinate_id);


--
-- Name: variants_coordinate_id_ref_alt_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX variants_coordinate_id_ref_alt_idx ON public.variants USING btree (coordinate_id, ref, alt);


--
-- Name: variants_features_feature_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_features_feature_id_idx ON public.variants_features USING btree (feature_id);


--
-- Name: variants_features_variant_id_feature_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX variants_features_variant_id_feature_id_idx ON public.variants_features USING btree (variant_id, feature_id);


--
-- Name: variants_gene_ensg_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_gene_ensg_idx ON public.variants USING btree (gene_ensg);


--
-- Name: variants_guides_guide_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_guides_guide_id_idx ON public.variants_guides USING btree (guide_id);


--
-- Name: variants_papers_paper_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_papers_paper_id_idx ON public.variants_papers USING btree (paper_id);


--
-- Name: variants_sources_source_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_sources_source_id_idx ON public.variants_sources USING btree (source_id);


--
-- Name: variants_subjects_subject_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variants_subjects_subject_id_idx ON public.variants_subjects USING btree (subject_id);


--
-- Name: variants_subjects_variant_id_subject_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX variants_subjects_variant_id_subject_id_idx ON public.variants_subjects USING btree (variant_id, subject_id);


--
-- Name: bystanders trg_set_has_adjacent_u; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER trg_set_has_adjacent_u BEFORE INSERT OR UPDATE OF has_u_before, has_u_after ON public.bystanders FOR EACH ROW EXECUTE FUNCTION public.set_has_adjacent_u();


--
-- Name: bystanders bystanders_alt_codon_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders
    ADD CONSTRAINT bystanders_alt_codon_id_fkey FOREIGN KEY (alt_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: bystanders_consequences bystanders_consequences_bystander_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders_consequences
    ADD CONSTRAINT bystanders_consequences_bystander_id_fkey FOREIGN KEY (bystander_id) REFERENCES public.bystanders(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: bystanders_consequences bystanders_consequences_consequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders_consequences
    ADD CONSTRAINT bystanders_consequences_consequence_id_fkey FOREIGN KEY (consequence_id) REFERENCES public.consequences(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: bystanders bystanders_ref_codon_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders
    ADD CONSTRAINT bystanders_ref_codon_id_fkey FOREIGN KEY (ref_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: bystanders bystanders_variant_id_guide_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.bystanders
    ADD CONSTRAINT bystanders_variant_id_guide_id_fkey FOREIGN KEY (variant_id, guide_id) REFERENCES public.variants_guides(variant_id, guide_id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: codons codons_amino_acid_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.codons
    ADD CONSTRAINT codons_amino_acid_id_fkey FOREIGN KEY (amino_acid_id) REFERENCES public.amino_acids(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: coordinates coordinates_genome_build_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.coordinates
    ADD CONSTRAINT coordinates_genome_build_id_fkey FOREIGN KEY (genome_build_id) REFERENCES public.genome_builds(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: exons exons_transcript_enst_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.exons
    ADD CONSTRAINT exons_transcript_enst_fkey FOREIGN KEY (transcript_enst) REFERENCES public.transcripts(enst) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: features features_biotype_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_biotype_id_fkey FOREIGN KEY (biotype_id) REFERENCES public.biotypes(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: gene_expressions gene_expressions_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_expressions
    ADD CONSTRAINT gene_expressions_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: gene_expressions gene_expressions_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_expressions
    ADD CONSTRAINT gene_expressions_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: gene_expressions gene_expressions_tissue_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_expressions
    ADD CONSTRAINT gene_expressions_tissue_id_fkey FOREIGN KEY (tissue_id) REFERENCES public.tissues(id);


--
-- Name: genes_annotations genes_annotations_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_annotations
    ADD CONSTRAINT genes_annotations_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: genes_annotations genes_annotations_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_annotations
    ADD CONSTRAINT genes_annotations_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: genes_qualitative_scores genes_qualitative_scores_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_qualitative_scores
    ADD CONSTRAINT genes_qualitative_scores_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: genes_qualitative_scores genes_qualitative_scores_qualitative_score_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_qualitative_scores
    ADD CONSTRAINT genes_qualitative_scores_qualitative_score_id_fkey FOREIGN KEY (qualitative_score_id) REFERENCES public.qualitative_scores(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: genes_quantitative_scores genes_quantitative_scores_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_quantitative_scores
    ADD CONSTRAINT genes_quantitative_scores_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: genes_quantitative_scores genes_quantitative_scores_quantitative_score_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genes_quantitative_scores
    ADD CONSTRAINT genes_quantitative_scores_quantitative_score_id_fkey FOREIGN KEY (quantitative_score_id) REFERENCES public.quantitative_scores(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: genome_builds genome_builds_organism_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_builds
    ADD CONSTRAINT genome_builds_organism_id_fkey FOREIGN KEY (organism_id) REFERENCES public.organisms(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: introns introns_transcript_enst_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.introns
    ADD CONSTRAINT introns_transcript_enst_fkey FOREIGN KEY (transcript_enst) REFERENCES public.transcripts(enst) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: organisms organisms_default_genome_build_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms
    ADD CONSTRAINT organisms_default_genome_build_id_fkey FOREIGN KEY (default_genome_build_id) REFERENCES public.genome_builds(id) ON UPDATE RESTRICT ON DELETE SET NULL;


--
-- Name: papers_subjects papers_subjects_paper_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers_subjects
    ADD CONSTRAINT papers_subjects_paper_id_fkey FOREIGN KEY (paper_id) REFERENCES public.papers(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: papers_subjects papers_subjects_subject_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.papers_subjects
    ADD CONSTRAINT papers_subjects_subject_id_fkey FOREIGN KEY (subject_id) REFERENCES public.subjects(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: protein_domains_proteins protein_domains_proteins_protein_domain_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains_proteins
    ADD CONSTRAINT protein_domains_proteins_protein_domain_id_fkey FOREIGN KEY (protein_domain_id) REFERENCES public.protein_domains(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: protein_domains_proteins protein_domains_proteins_protein_ensp_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains_proteins
    ADD CONSTRAINT protein_domains_proteins_protein_ensp_fkey FOREIGN KEY (protein_ensp) REFERENCES public.proteins(ensp) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: protein_domains protein_domains_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_domains
    ADD CONSTRAINT protein_domains_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id);


--
-- Name: proteins_annotations proteins_annotations_protein_ensp_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proteins_annotations
    ADD CONSTRAINT proteins_annotations_protein_ensp_fkey FOREIGN KEY (protein_ensp) REFERENCES public.proteins(ensp) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: proteins_annotations proteins_annotations_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proteins_annotations
    ADD CONSTRAINT proteins_annotations_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: qualitative_scores qualitative_scores_population_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.qualitative_scores
    ADD CONSTRAINT qualitative_scores_population_id_fkey FOREIGN KEY (population_id) REFERENCES public.populations(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: qualitative_scores qualitative_scores_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.qualitative_scores
    ADD CONSTRAINT qualitative_scores_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: quantitative_scores quantitative_scores_population_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.quantitative_scores
    ADD CONSTRAINT quantitative_scores_population_id_fkey FOREIGN KEY (population_id) REFERENCES public.populations(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: quantitative_scores quantitative_scores_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.quantitative_scores
    ADD CONSTRAINT quantitative_scores_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: tissues tissues_organ_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tissues
    ADD CONSTRAINT tissues_organ_id_fkey FOREIGN KEY (organ_id) REFERENCES public.organs(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: transcripts_annotations transcripts_annotations_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts_annotations
    ADD CONSTRAINT transcripts_annotations_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: transcripts_annotations transcripts_annotations_transcript_enst_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts_annotations
    ADD CONSTRAINT transcripts_annotations_transcript_enst_fkey FOREIGN KEY (transcript_enst) REFERENCES public.transcripts(enst) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: transcripts transcripts_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts
    ADD CONSTRAINT transcripts_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: transcripts transcripts_protein_ensp_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcripts
    ADD CONSTRAINT transcripts_protein_ensp_fkey FOREIGN KEY (protein_ensp) REFERENCES public.proteins(ensp) ON UPDATE SET NULL ON DELETE SET NULL;


--
-- Name: variants variants_alt_codon_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_alt_codon_id_fkey FOREIGN KEY (alt_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: variants_consequences variants_consequences_consequence_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_consequences
    ADD CONSTRAINT variants_consequences_consequence_id_fkey FOREIGN KEY (consequence_id) REFERENCES public.consequences(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_consequences variants_consequences_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_consequences
    ADD CONSTRAINT variants_consequences_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants variants_coordinate_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_coordinate_id_fkey FOREIGN KEY (coordinate_id) REFERENCES public.coordinates(id);


--
-- Name: variants_features variants_features_feature_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_features
    ADD CONSTRAINT variants_features_feature_id_fkey FOREIGN KEY (feature_id) REFERENCES public.features(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_features variants_features_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_features
    ADD CONSTRAINT variants_features_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants variants_gene_ensg_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_gene_ensg_fkey FOREIGN KEY (gene_ensg) REFERENCES public.genes(ensg) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: variants_guides variants_guides_guide_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_guides
    ADD CONSTRAINT variants_guides_guide_id_fkey FOREIGN KEY (guide_id) REFERENCES public.guides(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_guides variants_guides_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_guides
    ADD CONSTRAINT variants_guides_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_papers variants_papers_paper_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_papers
    ADD CONSTRAINT variants_papers_paper_id_fkey FOREIGN KEY (paper_id) REFERENCES public.papers(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_papers variants_papers_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_papers
    ADD CONSTRAINT variants_papers_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_qualitative_scores variants_qualitative_scores_qualitative_score_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_qualitative_scores
    ADD CONSTRAINT variants_qualitative_scores_qualitative_score_id_fkey FOREIGN KEY (qualitative_score_id) REFERENCES public.qualitative_scores(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_qualitative_scores variants_qualitative_scores_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_qualitative_scores
    ADD CONSTRAINT variants_qualitative_scores_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_quantitative_scores variants_quantitative_scores_quantitative_score_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_quantitative_scores
    ADD CONSTRAINT variants_quantitative_scores_quantitative_score_id_fkey FOREIGN KEY (quantitative_score_id) REFERENCES public.quantitative_scores(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_quantitative_scores variants_quantitative_scores_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_quantitative_scores
    ADD CONSTRAINT variants_quantitative_scores_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants variants_ref_codon_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants
    ADD CONSTRAINT variants_ref_codon_id_fkey FOREIGN KEY (ref_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


--
-- Name: variants_sources variants_sources_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_sources
    ADD CONSTRAINT variants_sources_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.sources(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_sources variants_sources_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_sources
    ADD CONSTRAINT variants_sources_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_subjects variants_subjects_subject_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_subjects
    ADD CONSTRAINT variants_subjects_subject_id_fkey FOREIGN KEY (subject_id) REFERENCES public.subjects(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- Name: variants_subjects variants_subjects_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variants_subjects
    ADD CONSTRAINT variants_subjects_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--

\unrestrict 7zFEJ1s4ZqFLoV3ZeieUXPeMX1jMav1Tar3WQmEEIllsslEDziTJk7xSKmfZDGt


--
-- Neighbor drill ("Improve" edit_type): per-variant codon-edit options + their scores.
-- Tables: neighbor_edits, neighbor_edit_scores, neighbor_edits_consequences.
-- Mirrors the bystanders / *_consequences conventions (identity id, RESTRICT codon
-- FKs, CASCADE variant/consequence FKs).
--

CREATE TABLE public.neighbor_edits (
    id integer NOT NULL,
    variant_id integer NOT NULL,
    pre_codon_id integer,
    post_codon_id integer,
    edited_positions integer[] NOT NULL,
    edited_cds_positions integer[],
    restores_reference boolean DEFAULT false NOT NULL,
    is_best boolean DEFAULT false NOT NULL,
    improves boolean DEFAULT false NOT NULL
);

ALTER TABLE public.neighbor_edits ALTER COLUMN id ADD GENERATED BY DEFAULT AS IDENTITY (
    SEQUENCE NAME public.neighbor_edits_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1
);

ALTER TABLE ONLY public.neighbor_edits
    ADD CONSTRAINT neighbor_edits_pkey PRIMARY KEY (id);

ALTER TABLE ONLY public.neighbor_edits
    ADD CONSTRAINT neighbor_edits_variant_id_edited_positions_key UNIQUE (variant_id, edited_positions);

CREATE INDEX neighbor_edits_variant_id_idx ON public.neighbor_edits USING btree (variant_id);

ALTER TABLE ONLY public.neighbor_edits
    ADD CONSTRAINT neighbor_edits_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY public.neighbor_edits
    ADD CONSTRAINT neighbor_edits_pre_codon_id_fkey FOREIGN KEY (pre_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;

ALTER TABLE ONLY public.neighbor_edits
    ADD CONSTRAINT neighbor_edits_post_codon_id_fkey FOREIGN KEY (post_codon_id) REFERENCES public.codons(id) ON UPDATE RESTRICT ON DELETE RESTRICT;


CREATE TABLE public.neighbor_edit_scores (
    neighbor_edit_id integer NOT NULL,
    source character varying NOT NULL,
    value double precision,
    label character varying,
    improvement double precision
);

ALTER TABLE ONLY public.neighbor_edit_scores
    ADD CONSTRAINT neighbor_edit_scores_pkey PRIMARY KEY (neighbor_edit_id, source);

ALTER TABLE ONLY public.neighbor_edit_scores
    ADD CONSTRAINT neighbor_edit_scores_neighbor_edit_id_fkey FOREIGN KEY (neighbor_edit_id) REFERENCES public.neighbor_edits(id) ON UPDATE CASCADE ON DELETE CASCADE;


CREATE TABLE public.neighbor_edits_consequences (
    neighbor_edit_id integer NOT NULL,
    consequence_id integer NOT NULL
);

ALTER TABLE ONLY public.neighbor_edits_consequences
    ADD CONSTRAINT neighbor_edits_consequences_pkey PRIMARY KEY (neighbor_edit_id, consequence_id);

ALTER TABLE ONLY public.neighbor_edits_consequences
    ADD CONSTRAINT neighbor_edits_consequences_neighbor_edit_id_fkey FOREIGN KEY (neighbor_edit_id) REFERENCES public.neighbor_edits(id) ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY public.neighbor_edits_consequences
    ADD CONSTRAINT neighbor_edits_consequences_consequence_id_fkey FOREIGN KEY (consequence_id) REFERENCES public.consequences(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- bystander_codon view: flags whether a bystander sits INSIDE the variant's own codon, so a
-- "deleterious bystander outside the codon" (e.g. CADD >= 25 & in_variant_codon = false) can be
-- expressed as a plain Hasura/GraphQL filter instead of per-row arithmetic. The check is in
-- genomic space (pre_mRNA bystanders carry genomic_position) and strand-aware; it assumes a
-- contiguous codon (split codons, ~0.7%, are not handled).
--

CREATE VIEW public.bystander_codon AS
 SELECT b.id,
    b.variant_id,
    b.guide_id,
    b.in_cds,
    b.genomic_position,
    b.cds_position,
    b.cadd_phred,
    b.cadd_raw,
    b.sift_score,
    b.sift_prediction,
    CASE
        WHEN (ge.strand = '+') THEN (b.genomic_position BETWEEN ((co.start - v.position_in_codon) + 1) AND ((co.start - v.position_in_codon) + 3))
        WHEN (ge.strand = '-') THEN (b.genomic_position BETWEEN ((co.start + v.position_in_codon) - 3) AND ((co.start + v.position_in_codon) - 1))
        ELSE NULL::boolean
    END AS in_variant_codon
   FROM public.bystanders b
     JOIN public.variants v ON ((v.id = b.variant_id))
     JOIN public.coordinates co ON ((co.id = v.coordinate_id))
     JOIN public.genes ge ON ((ge.ensg = v.gene_ensg));



--
-- Name: allele_frequencies; Type: TABLE; Schema: public; Owner: -
--
-- One row per variant that carries any population allele frequency, parsed from VEP
-- `Extra` columns (via --everything -> gnomAD --af_gnomade/--af_gnomadg + 1000G +
-- --max_af; see pipeline/phase1_preprocess.py). `gnomade_*` = gnomAD exomes, `gnomadg_*` =
-- gnomAD genomes (per-population + combined). `max_af`/`max_af_pops` are VEP's
-- cross-database maximum over 1000G/ESP/gnomAD (VEP-computed, not gnomAD-provided) and
-- so a row can exist on a 1000G-only frequency with all gnomAD columns NULL. Populated
-- in pipeline/phase2_db.py, via pipeline/allele_freq_common.py.
--

CREATE TABLE public.allele_frequencies (
    variant_id integer NOT NULL,
    gnomade_af double precision,
    gnomade_afr_af double precision,
    gnomade_amr_af double precision,
    gnomade_asj_af double precision,
    gnomade_eas_af double precision,
    gnomade_fin_af double precision,
    gnomade_mid_af double precision,
    gnomade_nfe_af double precision,
    gnomade_remaining_af double precision,
    gnomade_sas_af double precision,
    gnomadg_af double precision,
    gnomadg_afr_af double precision,
    gnomadg_ami_af double precision,
    gnomadg_amr_af double precision,
    gnomadg_asj_af double precision,
    gnomadg_eas_af double precision,
    gnomadg_fin_af double precision,
    gnomadg_mid_af double precision,
    gnomadg_nfe_af double precision,
    gnomadg_remaining_af double precision,
    gnomadg_sas_af double precision,
    max_af double precision,
    max_af_pops text
);

--
-- Name: allele_frequencies allele_frequencies_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.allele_frequencies
    ADD CONSTRAINT allele_frequencies_pkey PRIMARY KEY (variant_id);

--
-- Name: allele_frequencies allele_frequencies_variant_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.allele_frequencies
    ADD CONSTRAINT allele_frequencies_variant_id_fkey FOREIGN KEY (variant_id) REFERENCES public.variants(id) ON UPDATE CASCADE ON DELETE CASCADE;
