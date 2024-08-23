import argparse
import random
from Bio import SeqIO
import time
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Cut sequences into single-end or paired-end short reads with specific overlap.")
    parser.add_argument("input_fasta", help="Input FASTA file containing the long sequence.")
    parser.add_argument("output_prefix", help="Output FASTQ file stem for the extracted sequences.")
    parser.add_argument("-l", "--read_length", type=int, default=100, help="Length of the reads (default: 100bp).")
    parser.add_argument("-c", "--overlap", type=int, default=50, help="Overlap length between reads (default: 50bp).")
    parser.add_argument("-s", "--insize", type=int, default=150, help="Single-end if equal to read_length; Or insert size for paired-end reads (default: 150bp).")
    return parser.parse_args()

def write_read(read_seq, identifier, output_file):
    mode = 'at' if output_file.endswith('.gz') else 'a'
    with gzip.open(output_file, mode) if output_file.endswith('.gz') else open(output_file, mode) as f:
        f.write(f"@{identifier}\n")
        f.write(f"{read_seq}\n")
        f.write("+\n")
        f.write(f"{'-' * len(read_seq)}\n")

def extract_reads(id, sequence, read_length, overlap, insize, output_prefix):
    step = insize - overlap
    #timestamp = str(int(time.time()))
    
    for i in range(0, len(sequence) - read_length + 1, step):
        read_1 = sequence[i:i + read_length]
        id_1 = f"{id}:{i}/1"
        write_read(read_1, id_1, f"{output_prefix}_1.fq")
        
        if insize > read_length:
            paired_offset = i + insize - read_length
            if paired_offset + read_length <= len(sequence):
                read_2 = sequence[paired_offset:paired_offset + read_length]
                id_2 = f"{id}:{i}/2"
                write_read(read_2, id_2, f"{output_prefix}_2.fq")

def main():
    args = parse_args()

    # Load the long sequence from the FASTA file
    gzip_open = args.input_fasta.endswith('.gz')
    with (gzip.open(args.input_fasta, "rt") if gzip_open else open(args.input_fasta, "r")) as fasta_file:
        for fasta_record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(fasta_record.seq)
            identifier = fasta_record.id
            # Extract and write reads
            extract_reads(identifier, sequence, args.read_length, args.overlap, args.insize, args.output_prefix)

    print(f"Extracted reads from {args.input_fasta} with length {args.read_length}bp and overlap {args.overlap}bp.")

if __name__ == "__main__":
    main()
