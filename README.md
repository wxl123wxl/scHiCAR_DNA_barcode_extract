# sci-hicar

snakemake --latency-wait 50 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &


## This pipeline was used for generating pseudo-bulk fastq files by extracting DNA reads from the same clusters or cell types that identidied by scHiCAR RNA library.

To run this pipeline, you need to install the following software:
- **seqtk**: [install](https://github.com/lh3/seqtk)
- **pigz**: [install](https://zlib.net/pigz/)
- **nextflow**: [install](https://www.nextflow.io/docs/latest/install.html)
  
#### 1. trim specific sequence in the 5 end of read1 and read2
```bash

# bash code

cutadapt -j {threads} -e 0.15 --action=retain --discard-untrimmed  -g 'NNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG;min_overlap=30' -G 'NNNNNNNNNNNNNNGTCATAGCTGTTTCCTGTA;min_overlap=19' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
```
#### 2. Add DNA barcodes to the Seurat object generated from RNA library
```r

# R code

library(Seurat)

rna_seurat_object<-readRDS("rna_seurat.rds") # read personalized analysis results by Seurat based on filtered matrix from 1_RNA directory

dna_barcode<-read.table("total_RNA_DNA_barcode.txt",sep='\t',header=F,row.names=1) #1st column is RNA barcode and 2nd column is DNA barcode

dna_bd<-dna_barcode$V2

names(dna_bd)<-rownames(dna_barcode)
