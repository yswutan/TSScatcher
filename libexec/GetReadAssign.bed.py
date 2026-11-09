import pysam
import argparse
from tqdm import tqdm

def extract_max_zw_reads(input_bam, output_bed, tag_name, threshold):
    # Initialize dictionaries to store the maximum scores and corresponding reads
    max_score = {}
    max_read = {}

    # Open the input BAM file using pysam
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        total_reads = in_bam.mapped + in_bam.unmapped  # Estimate the total number of reads

        for read in tqdm(in_bam, total=total_reads, desc="Processing reads", unit="read"):
            read_name = read.query_name
            tags = dict(read.get_tags())

            # Check if the specified tag exists
            if tag_name not in tags:
                continue

            score = tags[tag_name]

            # Filter out scores less than or equal to the threshold
            if score <= threshold:
                continue

            # Update the max_score and max_read dictionaries
            if read_name in max_score:
                if score > max_score[read_name]:
                    max_score[read_name] = score
                    max_read[read_name] = read
            else:
                max_score[read_name] = score
                max_read[read_name] = read

    # Write the output to a BED file
    with open(output_bed, "w") as bedfile:
        for read_name, read in max_read.items():
            chrom = read.reference_name
            start = read.reference_start
            end = start + 1  # End position is one base after the start
            strand = '+' if not read.is_reverse else '-'  # Determine strand
            bedfile.write(f"{chrom}\t{start}\t{end}\t{read_name}\t{max_score[read_name]}\t{strand}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BAM file by ZW tag, extract max score reads, and create a BED file with start sites")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bed", help="Output BED file")
    parser.add_argument("tag_name", help="Tag name to filter and score reads")
    parser.add_argument("threshold", type=float, help="Score threshold for filtering reads")

    args = parser.parse_args()
    extract_max_zw_reads(args.input_bam, args.output_bed, args.tag_name, args.threshold)