import argparse
import pysam
import time

time_a=time.perf_counter()

# Define the command line arguments
parser = argparse.ArgumentParser(description='Calculate the number of called positions in sliding windows.')
parser.add_argument('--chromosome', type=str, help='The name of the chromosome to process.')
parser.add_argument('--chrom_size', type=int, help='The size of the chromosome.')
parser.add_argument('--vcf_file', type=str, help='The name of the VCF file to process.')
parser.add_argument('--window-size', type=int, default=1000, help='The size of the sliding window.')
parser.add_argument('--step-size', type=int, default=500, help='The step size between consecutive windows.')

# Parse the command line arguments
args = parser.parse_args()

# Open the VCF file using PySam
vcf = pysam.VariantFile(args.vcf_file)
nsamples = len(vcf.header.samples)

# print header
print(f"scaffold\tstart\tend\tpropCalled")
# Iterate over the sliding windows
for start_pos in range(0, args.chrom_size - args.window_size + 1, args.step_size):
    end_pos = start_pos + args.window_size - 1

    # Count the number of called positions in the window
    num_called = 0
    prop_called = 0
    for rec in vcf.fetch(args.chromosome, start_pos, end_pos):
        num_called += len(rec.samples)/nsamples
        prop_called = num_called / args.window_size
    # Print the results for this window
    print(f"{args.chromosome}\t{start_pos}\t{end_pos}\t{prop_called}")

# Close the VCF file
vcf.close()

# Print the total run time
time_b=time.perf_counter()
runtime=time_b-time_a
# print(f"Total run time: {runtime:.2f} seconds")
