# Cellsnp-lite

## 0. Set the environment and installation
```bash
conda create -n cellsnp-env -c bioconda -c conda-forge cellsnp-lite samtools bcftools htslib zlib liblzma -y
conda activate cellsnp-env

cellsnp-lite -V
```


## 1. Basic commends

Example.1
```bash
cellsnp-lite -s your_bam_file.bam -b your_barcodes.tsv -O result_output -p 10 --minMAF 0.05 -minCOUNT 100 --gzip -f hg38_genome.fa --genotype
```


## 2. Important parameters
```bash
-f, --refseq FILE
```