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
#### 2. extract the three barcodes at line of read name after "@"
```r

script/raw_fq_update.py
```
