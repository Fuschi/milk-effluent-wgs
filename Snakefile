import pandas as pd

# ─────────────────────────────────────────────────────────────
#                         Configuration
# ─────────────────────────────────────────────────────────────

# Read the sample configuration file with sample_id and biome (e.g. milk, sewage)
samples = pd.read_csv("config/sample_metadata.tsv", sep="\t")

# List of sample IDs
SAMPLES = samples["sample_id"].tolist()
# Corresponding list of biome labels
BIOMES = samples["biome"].tolist()
# Unique biomes (e.g. ["milk", "sewage"])
UNIQUE_BIOMES = sorted(set(BIOMES))

# Mapping sample → biome (used for grouping)
SAMPLE_TO_BIOME = dict(zip(SAMPLES, BIOMES))

# ─────────────────────────────────────────────────────────────
#                            Rules
# ─────────────────────────────────────────────────────────────

# Declares the final output files expected from the workflow
# ─────────────────────────────────────────────────────────────
rule all:
    input:
        expand("snakestream/reads_clean/{sample}_R{pe}_clean.fastq.gz", sample=SAMPLES, pe=["1", "2"])

# Download Bos taurus host genome from UCSC
# ─────────────────────────────────────────────────────────────
rule bostaurus_reference:
    output:
        fastagz = protected("snakestream/host/bosTau9.fa.gz"),
        fasta = temp("snakestream/host/bosTau9.fa"),
        index_files = protected(expand("snakestream/host/bosTau9.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ]))
    threads: 4
    conda:
        "hostremoval"
    shell:
        """
        wget -O {output.fastagz} \
             https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.fa.gz

        gunzip -c {output.fastagz} > {output.fasta}

        bowtie2-build --threads {threads} --seed 42 {output.fasta} snakestream/host/bosTau9
        """

# Filter reads by minimum length using reformat.sh (from BBTools)
# ─────────────────────────────────────────────────────────────
rule filter_length:
    input:
        r1 = "snakestream/reads_raw/{sample}_R1.fastq.gz",
        r2 = "snakestream/reads_raw/{sample}_R2.fastq.gz"
    output:
        r1_filt = temp("snakestream/reads_clean/{sample}_R1_lenfilt.fastq.gz"),
        r2_filt = temp("snakestream/reads_clean/{sample}_R2_lenfilt.fastq.gz")
    conda:
        "bbmap"
    shell:
        """
        reformat.sh in1={input.r1} in2={input.r2} \
                    out1={output.r1_filt} out2={output.r2_filt} \
                    minlength=50
        """

# Remove host reads using Bowtie2 and Samtools
# ─────────────────────────────────────────────────────────────
rule remove_host:
    input:
        r1 = "snakestream/reads_clean/{sample}_R1_lenfilt.fastq.gz",
        r2 = "snakestream/reads_clean/{sample}_R2_lenfilt.fastq.gz",
        index_check = "snakestream/host/bosTau9.1.bt2"
    output:
        r1_clean = protected("snakestream/reads_clean/{sample}_R1_clean.fastq.gz"),
        r2_clean = protected("snakestream/reads_clean/{sample}_R2_clean.fastq.gz"),
        bam =      temp("snakestream/reads_clean/{sample}_hostaligned.bam"),
        unmapped = temp("snakestream/reads_clean/{sample}_unmapped.bam"),
        sorted =   temp("snakestream/reads_clean/{sample}_unmapped_sorted.bam")
    threads: 4
    params:
        index = "snakestream/host/bosTau9"
    conda:
        "hostremoval"
    shell:
        """
        bowtie2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} -S /dev/stdout |
        samtools view -bS - > {output.bam}

        samtools view -b -f 12 -F 256 {output.bam} > {output.unmapped}
        samtools sort -n -m 5G -@ {threads} {output.unmapped} -o {output.sorted}

        samtools fastq -@ {threads} -1 {output.r1_clean} -2 {output.r2_clean} \
            -0 /dev/null -s /dev/null -n {output.sorted}
        """

