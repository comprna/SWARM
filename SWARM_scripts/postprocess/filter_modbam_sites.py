import pysam
import array


def get_read_to_ref_mapping(read):
    read_to_ref = {}
    ref_pos = read.reference_start  # Start of the alignment on the reference
    read_pos = 0  # Position in the read

    if read.cigartuples is None:  # unmapped or no CIGAR
        return read_to_ref

    # Iterate through the CIGAR tuples (operation, length)
    for cig_op, cig_len in read.cigartuples:
        if cig_op == 0:  # Alignment Match (M)
            for i in range(cig_len):
                read_to_ref[read_pos] = ref_pos
                ref_pos += 1
                read_pos += 1
        elif cig_op == 1:  # Insertion to reference
            read_pos += cig_len
        elif cig_op == 2:  # Deletion from reference
            ref_pos += cig_len
        elif cig_op == 3:  # Skipped region (N) → intron in spliced alignment
            ref_pos += cig_len
        elif cig_op == 4 or cig_op == 5:  # Soft/hard clipping
            read_pos += cig_len

    return read_to_ref



def load_bed(bed_file):
    """Load BED into a dict of chrom -> set(positions)."""
    positions = set()
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            # BED is usually 0-based, half-open. Convert to 1-based if needed.

            for pos in range(start + 1, end + 1):  # now 1-based inclusive
                positions.add(f"{chrom}_{pos}")
    return positions


def filter_mods_by_bed(in_bam, out_bam, bed_file):

    keep_positions = load_bed(bed_file)
    print(keep_positions)
    infile = pysam.AlignmentFile(in_bam, "rb")
    outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)

    for read in infile:
        if read.has_tag('MM') and read.has_tag('ML'):

            # Extract the MM and ML tags
            mm_tag = read.get_tag('MM')
            ml_tag = read.get_tag('ML')

            ml_probs = list(ml_tag)
            read_to_ref = get_read_to_ref_mapping(read)
            contig = infile.get_reference_name(read.reference_id)
            read_id = read.query_name
            mod_dct = read.modified_bases

            try:
                key, item = next(iter(mod_dct.items()))
                mm_position, probability = next(iter(item))
            except:
                print(read_id, "skipped")
                continue

            mm_new_values = []
            # print(mod_dct)
            for key, item in mod_dct.items():
                for mm_position, probability in sorted(item,key = lambda x:x[0]):
                    if mm_position in read_to_ref:
                        ref_pos = read_to_ref[mm_position]
                        #print(f"{contig}_{ref_pos + 1}  {probability}")
                        # print(f"{contig}_{ref_pos + 1}")
                        if f"{contig}_{ref_pos + 1}" in keep_positions:
                            print("validated_pos",probability)
                            mm_new_values.append(probability)
                        else:
                            mm_new_values.append(0)
                    else:
                        mm_new_values.append(0)


        #print(new_ml)
        if len(mm_new_values) > 0:
            read.set_tag("ML", mm_new_values)

        outfile.write(read)
        # exit()
    infile.close()
    outfile.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python filter_mods.py in.bam out.bam regions.bed \n")
        sys.exit(1)

    in_bam, out_bam, bed_file = sys.argv[1:]
    filter_mods_by_bed(in_bam, out_bam, bed_file)
