"""Default parameters for soybean primer design, matching NCBI Primer-BLAST defaults."""

# Primer3 design parameters
DEFAULT_PRIMER_OPT_SIZE = 22
DEFAULT_PRIMER_MIN_SIZE = 20
DEFAULT_PRIMER_MAX_SIZE = 27
DEFAULT_PRIMER_OPT_TM = 60.0
DEFAULT_PRIMER_MIN_TM = 57.0
DEFAULT_PRIMER_MAX_TM = 63.0
DEFAULT_PRIMER_MAX_DIFF_TM = 3.0
DEFAULT_PRIMER_PRODUCT_MIN = 100
DEFAULT_PRIMER_PRODUCT_MAX = 400
DEFAULT_NUM_RETURN = 5
DEFAULT_PRIMER_GC_MIN = 35.0
DEFAULT_PRIMER_GC_MAX = 65.0
DEFAULT_MAX_POLY_X = 5
DEFAULT_GC_CLAMP = 0
DEFAULT_MAX_END_STABILITY = 9.0
DEFAULT_MAX_END_GC = 5

# Biochemical conditions (matching NCBI Primer-BLAST)
DEFAULT_MONO_CATIONS = 50.0       # mM KCl
DEFAULT_DIVA_CATIONS = 1.5        # mM MgCl2
DEFAULT_CON_DNTPS = 0.6           # mM dNTPs
DEFAULT_SALT_CORRECTION = 1       # SantaLucia 1998
DEFAULT_TM_METHOD = 1             # SantaLucia 1998

# Specificity checking
DEFAULT_MAX_OFF_TARGET_MISMATCH = 1
DEFAULT_MIN_3PRIME_MATCH = 7
DEFAULT_MAX_AMPLICON_SIZE = 4000
DEFAULT_BLAST_EVALUE = 1000
DEFAULT_BLAST_WORD_SIZE = 7
DEFAULT_BLAST_MAX_TARGET_SEQS = 2000

# BLAST paths
DEFAULT_BLASTN_PATH = "C:/Users/leo/blast/bin/blastn.exe"
DEFAULT_MAKEBLASTDB_PATH = "C:/Users/leo/blast/bin/makeblastdb.exe"

# Default file paths (relative to project root)
DEFAULT_GENOME_FA = "Gmax_508_v4.0.fa"
DEFAULT_GFF3 = "Gmax_508_Wm82.a4.v1.gene.gff3"
DEFAULT_TRANSCRIPT_FA = "Gmax_508_Wm82.a4.v1.transcript_primaryTranscriptOnly.fa"
DEFAULT_ANNOTATION_TSV = "Soybean_Arabidopsis_Complete_Annotation.tsv"

# Off-target sensitivity presets
OFF_TARGET_PRESETS = {
    "strict":   {"max_mismatches": 0, "min_3prime": 8,  "min_coverage": 0.85},
    "standard": {"max_mismatches": 1, "min_3prime": 7,  "min_coverage": 0.75},
    "sensitive":{"max_mismatches": 2, "min_3prime": 6,  "min_coverage": 0.60},
    "relaxed":  {"max_mismatches": 3, "min_3prime": 5,  "min_coverage": 0.50},
}
