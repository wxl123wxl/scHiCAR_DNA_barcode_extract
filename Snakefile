shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")
import collections
configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON'])) ##  fastq for each lane
SAMPLES = sorted(FILES.keys())

#chromsizes = config['chromsizes']
cool_bin       = config['cool_bin']
TARGETS = []
genome = config['genome']
tmp='/work/xw171/tmp'
barcode = json.load(open('barcode.json')) # json file to include barcode file location.

genome_version = 'hs'
smooth_window  = 150
shiftsize      = -75
pval_thresh    = 0.01


def load_group(grouptxt): # group per sample key
    groupid = collections.defaultdict(list)
    with open(grouptxt,'r') as f:
        for line in f:
            line = line.strip()
            if not line =='':
                group, value = line.split("\t")
                groupid[group].append(value)
    return groupid

groupid = load_group('group.txt')



# frag_path = config['frag_path']

TARGETS.extend(expand("05_ME_fq/{sample}_cutM13R_cutME_L001_R2_001.fastq.gz",sample = SAMPLES ))
#TARGETS.extend(expand("06_barcode_fq/{sample}_cutM13R_cutME_L001_mouse_muscle_day7_25_R1_001.fastq.gz",sample = SAMPLES ))

localrules: targetfiles
rule targetfiles:
    input: TARGETS


rule cut_R2_tso:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
        #r3 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R1_001.fastq",
        #r4 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R2_001.fastq"
    output:
        r1 = "00_raw_fq_update/{sample}_cut_adapt3_L001_R1_001.fastq.gz",
        r2 = "00_raw_fq_update/{sample}_cut_adapt3_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        """
        cutadapt -j {threads} -e 0.15 --action=retain --discard-untrimmed  -g 'NNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG;min_overlap=30' -G 'NNNNNNNNNNNNNNGTCATAGCTGTTTCCTGTA;min_overlap=19' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
        """


## run for each fastq pairs individually. 

rule raw_fq_trim:
    input:
        r1 = "00_raw_fq_update/{sample}_cut_adapt3_L001_R1_001.fastq.gz",
        r2 = "00_raw_fq_update/{sample}_cut_adapt3_L001_R2_001.fastq.gz"
    output: 
        r1 = "01_raw_fq_update/{sample}_index_L001_R1_001.fastq",
        r2 = "01_raw_fq_update/{sample}_index_L001_R2_001.fastq"
        #r1 = temp("01_raw_fq_update/{sample}_index_L001_R1_001.fastq"),
        #r2 = temp("01_raw_fq_update/{sample}_index_L001_R2_001.fastq")
    script:
        "script/raw_fq_update.py"

rule barcode_QC: ## extract total barcodes list
    input: "01_raw_fq_update/{sample}_index_L001_R1_001.fastq"
    output: "02_barcode_info/{sample}_raw_barcode_count.txt"
    threads: 8 
    shell:
        """
        export TMPDIR=/work/xw171/tmp
        awk  '{{if(NR%4==1) print substr($0,2,18)}}' {input} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output} 
        """

rule find_right_barcodes: # mismatch correction for each barcode
    input: "02_barcode_info/{sample}_raw_barcode_count.txt"
    output: sum = "02_barcode_info/{sample}.barcode_final_summary",
            map = "02_barcode_info/{sample}.barcode_final_map",
            log = "02_barcode_info/{sample}.barcode_log"
    script:
        "script/barcode_hash_v2_ME.py"


rule read_barcode_correction:
    input :
        "01_raw_fq_update/{sample}_index_L001_R1_001.fastq",
        "01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
        "02_barcode_info/{sample}.barcode_final_map"
    output :
        "03_corrected_fq/{sample}_L001_R1_001.fastq",
        "03_corrected_fq/{sample}_L001_R2_001.fastq"
    log: "00_log/{sample}_L001_R1_corrected.log"
    script:
        "script/fq_barcode_correction_R1_ME.py"

#rule read2_barcode_correction:
#    input :
#        "01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
#        "02_barcode_info/{sample}.barcode_final_map"
#    output :
#        "03_corrected_fq/{sample}_L001_R2_001.fastq"
#    log:"00_log/{sample}_L001_R2_corrected.log"   
#    script:
#        "script/fq_barcode_correction_R1_ME.py"

rule r1_zip:
    input  : "03_corrected_fq/{sample}_L001_R1_001.fastq"
    output : "03_corrected_fq/{sample}_L001_R1_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule r2_zip:
    input  : "03_corrected_fq/{sample}_L001_R2_001.fastq"
    output : "03_corrected_fq/{sample}_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"


rule ME:
    input:
        "03_corrected_fq/{sample}_L001_R1_001.fastq.gz",
        "03_corrected_fq/{sample}_L001_R2_001.fastq.gz"
    output:
        #"03_ME_fq/{sample}_L001_R1_001.fastq.gz",
        #"03_ME_fq/{sample}_L001_R2_001.fastq.gz"
        "05_ME_fq/{sample}_cutM13R_cutME_L001_R1_001.fastq.gz",
        "05_ME_fq/{sample}_cutM13R_cutME_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        "cutadapt -Z -j {threads} -e 0.2  -g file:ME_index\
        -o 05_ME_fq/{wildcards.sample}_cutM13R_cut{{name}}_L001_R1_001.fastq.gz \
        -p 05_ME_fq/{wildcards.sample}_cutM13R_cut{{name}}_L001_R2_001.fastq.gz \
        {input[0]}  {input[1]}"


rule split_barcode:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
        #"05_ME_fq/{sample}_cutM13R_cutME_L001_R1_001.fastq.gz",
        #"05_ME_fq/{sample}_cutM13R_cutME_L001_R2_001.fastq.gz"
    output:
        #"03_ME_fq/{sample}_L001_R1_001.fastq.gz",
        #"03_ME_fq/{sample}_L001_R2_001.fastq.gz"
        "06_barcode_fq/{sample}_cutM13R_cutME_L001_mouse_muscle_day7_25_R1_001.fastq.gz",
        "06_barcode_fq/{sample}_cutM13R_cutME_L001_mouse_muscle_day7_25_R2_001.fastq.gz"
    threads: 11
    shell:
        "python barcode_split.py --forward_fastq {input[0]} --reverse_fastq {input[1]} --barcode_file mouse_muscle_barcod_DNAchange.txt --output_dir 06_barcode_fq"

