# cellsnp-lite

## 0. Set the environment and installation
```bash
conda create -n cellsnp-env -c bioconda -c conda-forge cellsnp-lite samtools bcftools htslib zlib liblzma -y
conda activate cellsnp-env

cellsnp-lite -V
```


## 1. Basic commends

```bash
cellsnp-lite -s your_bam_file.bam -b your_barcodes.tsv -O result_output -p 10 --minMAF 0.05 -minCOUNT 20 --gzip -f hg38_genome.fa --genotype
```

## 2. Important parameters
```bash
-f, --refseq FILE
-s, --samFile STR
-S, --samFileList FILE
-R, --regionsVCF FILE
--genotype
--minCOUNT 
--minMAF
```

**-f** or **--refseq** could provide the normal real genome for reference. Because most of the time we analyzed the abnormal genomes or genomes from diseases.

**-R** or **--regionsVC** normally we use default because we want to the all of the information.

**--genotype** if use, do the genotyping in addition to counting

**--minCOUNT** is 20 in default, it is reasonable.

**--minMAF** is 0.00 in default, we could try 0.1(when you **UNDERSTAND** this!).