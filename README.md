## This pipeline was used for extract the cell barcode, trim adaptor sequence from DNA fastq file and write the three barcode at line of read name following "@"

snakemake --latency-wait 50 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &


  
#### 1. trim specific sequence in the 5 end of read1 and read2
```bash

# bash code

cutadapt -j {threads} -e 0.15 --action=retain --discard-untrimmed  -g 'NNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG;min_overlap=30' -G 'NNNNNNNNNNNNNNGTCATAGCTGTTTCCTGTA;min_overlap=19' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
```
#### 2. extract the three barcodes and write the three barcode at line of read name following "@"
```

script/raw_fq_update.py
```
#### 3. extract total barcodes list, count each barcode
```

export TMPDIR=/work/xw171/tmp
awk  '{{if(NR%4==1) print substr($0,2,18)}}' {input} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output}
```
#### 4. compare the extracted barcodes with the whitelist
```

script/barcode_hash_v2_ME.py
```
#### 5. correct the barcode which barcode only have one mismatch
```

script/fq_barcode_correction_R1_ME.py
```

#### 5. compress fastq file
```

pigz -p {threads} {input}
```
#### 6. remove the ME sequence in the fastq file
```

cutadapt -Z -j {threads} -e 0.2  -g file:ME_index\
        -o 05_ME_fq/{wildcards.sample}_cutM13R_cut{{name}}_L001_R1_001.fastq.gz \
        -p 05_ME_fq/{wildcards.sample}_cutM13R_cut{{name}}_L001_R2_001.fastq.gz \
        {input[0]}  {input[1]}
```
