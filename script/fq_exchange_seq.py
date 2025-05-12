#!/hpc/home/xw171/miniconda3/envs/ArchRenv/bin/python
import gzip
import argparse

#exchange sequence in fq.gz
parser = argparse.ArgumentParser()
parser.add_argument("--fq_file", help="fastq file", type=str)
args = parser.parse_args()
fq = args.fq_file
fq_file = fq.split("/",)[-1]
fq_file = fq_file.split(".gz",)[0]
exchange_fq = "exchange" + fq_file 


def update_fastq(fq, exchange_fq): ## process two files

    f_fq = gzip.open(fq, 'rt')
    f_exchange_fq = open(exchange_fq, 'w')

    while True:
        
        cur_fq_name = f_fq.readline().strip()[1:]
        cur_fq_read = f_fq.readline().strip()
        cur_fq_plus = f_fq.readline().strip()
        cur_fq_qual = f_fq.readline().strip()
    
        
        position1 = cur_fq_read[:100]  ## 6bp tn5 index sequence
        position2 = cur_fq_read[100:110]  ## 6bp tn5 index sequence
        position3 = cur_fq_read[110:]  ## 6bp tn5 index sequence
        
        cur_fq_read = position1 + position3 + position2 #+ tn5_index_i5+tn5_index_i7)
        
    
        f_exchange_fq.write(cur_fq_name+"\n")
        f_exchange_fq.write(cur_fq_read+"\n")
        f_exchange_fq.write(cur_fq_plus+"\n")
        f_exchange_fq.write(cur_fq_qual+"\n")
        
    f_fq.close()
    f_exchange_fq.close()
    
update_fastq(fq, exchange_fq)
os.system('pigz -p 11 '+exchange_fq)
