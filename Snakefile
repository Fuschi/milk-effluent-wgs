# ─────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────
import pandas as pd

# Read metadata table
samples = pd.read_csv("config/sample_metadata.tsv", sep="\t")

# Extract sample info
SAMPLES = samples["sample_id"].tolist()
BIOMES = samples["biome"].tolist()
UNIQUE_BIOMES = sorted(set(BIOMES))
SAMPLE_TO_BIOME = dict(zip(SAMPLES, BIOMES))

# ─────────────────────────────────────────────────────────────
# RULE INCLUDES
# ─────────────────────────────────────────────────────────────
include: "rules/cleaning.smk"

# ─────────────────────────────────────────────────────────────
# FINAL TARGETS
# ─────────────────────────────────────────────────────────────
rule all:
    input:
        expand("snakestream/reads_clean/{sample}_R{pe}_clean.fastq.gz", sample=SAMPLES, pe=["1", "2"])
