import sys
import string

def read_breakpoints(breakpoints_file):
    breakpoints = {}
    with open(breakpoints_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            breaks = parts[1].split(",")
            scaffold = parts[0]
            breakpoints[scaffold] = [int(bp) for bp in breaks[0:]] + [20000000]
    return breakpoints

def split_bed_file(bed_file, breakpoints):
    with open(bed_file, 'r') as f:
        lines = f.readlines()

        # Write BED lines to appropriate output files based on breakpoints
        # for scaffold in breakpoints:
        output_files = {}
        bk_index = 0
        suffix = chr(ord('`')+1+bk_index)
        output_files[bk_index] = open(f'{out_dir}{chromCode}_{suffix}.list', 'w')

        for line in lines:
            parts = line.strip().split(".")
            chrom, start, end = parts[:3]
            start, end = int(start), int(end)

            if chrom in breakpoints:
                if end <= breakpoints[chrom][bk_index]:
                    output_files[bk_index].write(line)

                else:
                    bk_index += 1
                    suffix = chr(ord('`')+1+bk_index)
                    output_files[bk_index] = open(f'{out_dir}{chromCode}_{suffix}.list', 'w')
                    output_files[bk_index].write(line)

            # else:
            #     # If scaffold not in breakpoints, write to a default file
            #     output_files["default"].write(line)

        # Close output files
        for file in output_files.values():
            file.close()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <bed_file> <breakpoints_file>")
        sys.exit(1)

    bed_file = sys.argv[1]
    breakpoints_file = sys.argv[2]
    out_dir = sys.argv[3]
    chromCode = sys.argv[4]
    
    breakpoints = read_breakpoints(breakpoints_file)
    # print(breakpoints)
    split_bed_file(bed_file, breakpoints)
