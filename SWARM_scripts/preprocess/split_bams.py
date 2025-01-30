import pysam
import argparse

def split_bam_file(input_bam_path, output_prefix, num_chunks):
    
    if input_bam_path[-4:] == ".sam":
        print("INPUT IS SAM",flush=True)
        input_bam = pysam.AlignmentFile(input_bam_path, "r")
        total_reads = sum(1 for _ in input_bam)
        input_bam.close()
        input_bam = pysam.AlignmentFile(input_bam_path, "r")

    elif input_bam_path[-4:] == ".bam":
        print("INPUT IS NOT SAM",flush=True)
        input_bam = pysam.AlignmentFile(input_bam_path, "rb")
        total_reads = input_bam.mapped
    else:
        raise ("Input path extension must end in .bam or .sam")

    #total_reads = input_bam.mapped
    #total_reads = sum(1 for _ in input_bam)
    print("ORIGINAL READS =",total_reads)
    reads_per_chunk = total_reads // num_chunks

    current_chunk = 0
    current_chunk_reads = 0
    for read in input_bam:
        if current_chunk_reads == 0:
            if input_bam_path[-4:] == ".sam":
                print("INPUT IS SAM",flush=True)
                output_bam = pysam.AlignmentFile(f"{output_prefix}_{current_chunk + 1}.sam", "w", header=input_bam.header)
            else:
                print("INPUT IS NOT SAM",flush=True)
                output_bam = pysam.AlignmentFile(f"{output_prefix}_{current_chunk + 1}.bam", "wb", header=input_bam.header)
        output_bam.write(read)
        current_chunk_reads += 1

        # Move to the next chunk if the current chunk is full
        if current_chunk_reads >= reads_per_chunk and current_chunk < num_chunks - 1:
            print("STARTING NEW CHUNK")
            output_bam.close()
            current_chunk += 1
            current_chunk_reads = 0

    if current_chunk_reads > 0:
        output_bam.close()

    input_bam.close()
    print("INPUT CLOSED")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a BAM file into multiple chunks.")

    # Add command-line arguments
    parser.add_argument("--input_bam" ,"-i" ,help="Path to the input BAM file.")
    parser.add_argument("--output_prefix" ,"-o" ,help="Path for the output BAM files.")
    parser.add_argument("--num_chunks" ,"-n" ,type=int, help="Number of chunks to create.")

    # Parse command-line arguments
    args = parser.parse_args()

    # Split the BAM file into chunks
    split_bam_file(args.input_bam, args.output_prefix, args.num_chunks)

