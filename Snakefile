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
BIOME_TO_SAMPLE={"Milk":[],"Sewage":[])
for i in SAMPLE_TO_BIOME.keys():
	BIOME_TO_SAMPLE[SAMPLE_TO_BIOME[i]]=BIOME_TO_SAMPLE[SAMPLE_TO_BIOME[i]]+[i]
#General resources
def default_resources():
    return dict(
        requeue=0,
        trigger=0,
    )
# ─────────────────────────────────────────────────────────────
# RULE INCLUDES
# ─────────────────────────────────────────────────────────────
include: "rules/cleaning.smk"
include: "rules/stats_reads.smk"

# ─────────────────────────────────────────────────────────────
# FINAL TARGETS
# ─────────────────────────────────────────────────────────────
rule all:
    input:
        expand("snakestream/reads_clean/{sample}_R{pe}_clean.fastq.gz", sample=SAMPLES, pe=["1", "2"]),
        expand("snakestream/qc/trim/{sample}_R{pe}_trim_fastqc.html", sample=SAMPLES, pe=["1", "2"]),
        expand("snakestream/qc/raw/{sample}_R{pe}_fastqc.html", sample=SAMPLES, pe=["1", "2"]),
        "snakestream/stats/seqkit_raw_reads.tsv",
        "snakestream/stats/seqkit_trimmed_reads.tsv",
        "snakestream/stats/seqkit_cleaned_reads.tsv"
    resources:
        mem_mb=1000,
        qos="normal",
        time="00:05:00",
        **default_resources(),
