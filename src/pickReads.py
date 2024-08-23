import argparse
import random
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Extract short sequences from a long FASTA sequence and output in FASTQ format based on insert size.")
    parser.add_argument("input_fasta", help="Input FASTA file containing the long sequence.")
    parser.add_argument("output_prefix", help="Output FASTQ file stem for the extracted sequences with '_1.fq' and '_2.fq' suffixes for paired-end.")
    parser.add_argument("-r", "--read_length", type=int, default=100, help="Length of the reads (default: 100bp).")
    parser.add_argument("-i", "--insert_size", type=int, default=300, help="Insert size for paired-end reads (default: 300bp).")
    parser.add_argument("-b", "--barcode", default="NA", help="Assign a barcode tag to each record.")
    parser.add_argument("-n", "--num_reads", type=int, default=1, help="Number of reads to extract (default: 1).")
    return parser.parse_args()

def write_read(read_seq, identifier, output_file):
    with open(output_file, "a") as f:
        f.write(f"@{identifier}\n")
        f.write(f"{read_seq}\n")
        f.write("+\n")
        f.write(f"{'-' * len(read_seq)}\n")

def extract_and_write_reads(barcode, fasta_record, read_length, insize, num_reads, insert_size, output):
    sequence = str(fasta_record.seq)
    for _ in range(num_reads):
        start = random.randint(0, len(sequence) - read_length)
        read = sequence[start:start + read_length]
        pos = ":" + str(start)
        bc = "#" + barcode
        identifier = fasta_record.id + pos + bc + ("/1" if read_length <= insert_size else "")

        write_read(read, identifier, f"{output}_{'1' if read_length <= insert_size else ''}.fq")

        if read_length < insert_size:
            # For paired-end, write the second read
            start_2 = start + insize - read_length + random.randint(-10,10)
            read_2 = sequence[start_2:start_2 + read_length]
            pos_2 = ":" + str(start_2)
            write_read(read_2, fasta_record.id + pos_2 + bc + "/2", f"{output}_2.fq")

def main():
    args = parse_args()

    # Load the long sequence from the FASTA file
    fasta_record = next(SeqIO.parse(args.input_fasta, "fasta"), None)
    if not fasta_record:
        raise ValueError("No sequence found in the input FASTA file.")

    # Determine read length based on insert size
    read_length = args.read_length if args.read_length <= args.insert_size else args.insert_size

    # Extract and write reads
    extract_and_write_reads(args.barcode, fasta_record, read_length, args.insert_size, args.num_reads, args.insert_size, args.output_prefix)

    print(f"Extracted {args.num_reads} reads of length {read_length}bp from {args.input_fasta}.")

if __name__ == "__main__":
    main()