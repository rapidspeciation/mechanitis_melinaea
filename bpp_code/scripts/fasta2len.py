import argparse
from Bio import SeqIO

# Set up argument parser
parser = argparse.ArgumentParser(description="Print the name and length of each sequence in a FASTA file.")
parser.add_argument("fasta_file", help="Path to the input FASTA file")

# Parse the arguments
args = parser.parse_args()

# Parse the FASTA file and print the name and length of each sequence
for record in SeqIO.parse(args.fasta_file, "fasta"):
    sequence_name = record.id
    sequence_length = len(record.seq)
    print(f"{sequence_name}\t1\t{sequence_length}")
