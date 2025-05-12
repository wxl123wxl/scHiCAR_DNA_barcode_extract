#"mouse_muscle_barcode.txt"的第5列为”sample_barcode”，"mouse_muscle_barcode.txt"的第4列为对应的”sample_name”。fastq.gz文件的“name”行以“@”开始，紧接着6个字母为“sample_barcode”。写一个脚本，将fastq.gz 按照”sample_barcode”分为多个文件，新文件以对应的”sample_name”命名。
import os
import gzip
import argparse
# Extract library_name from file path
def extract_library_name(file_path):
    base_name = os.path.basename(file_path)
    library_name = base_name.split('_R1_')[0]  # Extract the part before "_R1_001.fastq.gz"
    return library_name

# Function to read fastq file in chunks of 4 lines (for paired reads)
def read_fastq_pair(forward_file, reverse_file):
    with gzip.open(forward_file, 'rt') as f_fw, gzip.open(reverse_file, 'rt') as f_rv:
        while True:
            # Read forward read
            name_fw = f_fw.readline().strip()
            if not name_fw:
                break
            sequence_fw = f_fw.readline().strip()
            plus_fw = f_fw.readline().strip()
            quality_fw = f_fw.readline().strip()

            # Read reverse read
            name_rv = f_rv.readline().strip()
            sequence_rv = f_rv.readline().strip()
            plus_rv = f_rv.readline().strip()
            quality_rv = f_rv.readline().strip()

            yield (name_fw, sequence_fw, plus_fw, quality_fw), (name_rv, sequence_rv, plus_rv, quality_rv)

# Function to write forward and reverse reads to corresponding files
def write_to_fastq(file_dict, sample_name, library_name, forward_read, reverse_read, output_dir):
    if sample_name not in file_dict:
        # Create output file paths
        forward_path = os.path.join(output_dir, f"{library_name}_{sample_name}_R1_001.fastq.gz")
        reverse_path = os.path.join(output_dir, f"{library_name}_{sample_name}_R2_001.fastq.gz")

        # Open two files for forward and reverse reads with library_name in the filename (gzipped output)
        file_dict[sample_name] = {
            'forward': gzip.open(forward_path, 'wt'),
            'reverse': gzip.open(reverse_path, 'wt')
        }
    
    # Write forward read
    file_dict[sample_name]['forward'].write("\n".join(forward_read) + "\n")
    
    # Write reverse read
    file_dict[sample_name]['reverse'].write("\n".join(reverse_read) + "\n")

# Load the barcode to sample_name mapping from mouse_muscle_barcode.txt
def load_barcode_mapping(barcode_file):
    name_to_barcodes = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 5:
                sample_barcode = columns[4]
                sample_name = columns[3]
                if sample_name not in name_to_barcodes:
                    name_to_barcodes[sample_name] = []
                name_to_barcodes[sample_name].append(sample_barcode)
    return name_to_barcodes

# Main function to split fastq files by barcode
def split_fastq_by_barcode(forward_fastq, reverse_fastq, barcode_file, output_dir):
    # Load barcode mapping
    name_to_barcodes = load_barcode_mapping(barcode_file)

    # Initialize a dictionary to keep track of open file handles
    paired_files = {}

    # Extract library_name from the forward fastq file name
    library_name = extract_library_name(forward_fastq)

    # Process both forward and reverse fastq files simultaneously
    for forward_read, reverse_read in read_fastq_pair(forward_fastq, reverse_fastq):
        sample_barcode = forward_read[0][1:7]  # Extract first 6 characters after '@' from forward read name
        sample_name = None
        # Find the corresponding sample_name based on sample_barcode
        for name, barcodes in name_to_barcodes.items():
            if sample_barcode in barcodes:
                sample_name = name
                break
        if sample_name:
            # Write forward and reverse read to the corresponding files
            write_to_fastq(paired_files, sample_name, library_name, forward_read, reverse_read, output_dir)

    # Close all open files
    for file_dict in paired_files.values():
        file_dict['forward'].close()
        file_dict['reverse'].close()

    print("Splitting completed.")

# Command-line argument parser
parser = argparse.ArgumentParser(description="python split_fastq.py --forward_fastq path/to/library_name_R1_001.fastq.gz --reverse_fastq path/to/library_name_R2_001.fastq.gz --barcode_file path/to/mouse_muscle_barcode.txt")
parser.add_argument("--forward_fastq", required=True, help="Path to the forward fastq.gz file (e.g. library_name_R1_001.fastq.gz)")
parser.add_argument("--reverse_fastq", required=True, help="Path to the reverse fastq.gz file (e.g. library_name_R2_001.fastq.gz)")
parser.add_argument("--barcode_file", required=True, help="Path to the barcode file (e.g. mouse_muscle_barcode.txt)")
parser.add_argument("--output_dir", required=True, help="Directory to save the output fastq.gz files")

args = parser.parse_args()

# Call the main function to split the fastq files by barcode
split_fastq_by_barcode(args.forward_fastq, args.reverse_fastq, args.barcode_file, args.output_dir)
