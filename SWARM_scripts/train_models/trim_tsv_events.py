"""
@author: Stefan Prodic

Output limited number of event align tsv signals per 9mer
"""

import os
import copy
import time
import argparse
from collections import Counter


BUFF_LIM = 50
def write_trimmed_nanopolish(file_path, lines_5mers, ninemer_seq,cbuf,buffer_lines, process_index = ""):
    output_path = file_path + "_trimmed.tsv"

    # if cbuf == 1000:
        # print("ITS 1000")

    if cbuf >= BUFF_LIM:
        buffer_lines.append(lines_5mers)
        # if cbuf == 1000:
        #     print("OK 1000")
        #     print(len(buffer_lines))
            
        cbuf = 0
        
        with open(output_path, "a+") as trimmed_file:
            for lines_9mers in buffer_lines:
                for lines_5mer in lines_9mers:
                    
                    trimmed_file.write(lines_5mer)
                #trimmed_file.write("END\n")
                #print(len(lines_5merss), lines_5merss[0])
            #sys.exit()
        buffer_lines = []
    else:
        cbuf +=1
        buffer_lines.append(copy.deepcopy(lines_5mers))
    #print(cbuf,"cbf")
    return cbuf, buffer_lines



def parse_nanopolish(nanopolish_file):
    """
    Main function for the preprocessing.
    Iterates through the nanopolish file seeking 5 consecutive 5-mers.
    Outputs the vector containing currents, distances from expected currents for all bases, dwelling time, and qscores.
    """

    print("starting to parse nanopolish", nanopolish_file)
    start_time_np = time.time()
    counter_9mers = Counter()
    OUTPUT_file = PATH_OUT


    process_index = ""


    with open(nanopolish_file) as nf:
        prev_pos = -9999
        prev_read, prev_contig, current_seq = "", "", ""

        with open(OUTPUT_file + "_trimmed.tsv", "a+") as trimmed_file:
            trimmed_file.write(nf.readline())
        
        consecutive = 0
        sample_lst = []
        cb,lb = 0,[]
        c = 0
        for line in nf:
            c += 1
            if c % 1000000 == 0:
                print(c, time.time() - start_time_np, flush = True)
            lst = line.split("\t")
            contig, pos, kmer, read_name = lst[:4]
            try:
                pos = int(pos)
            except:
                print("invalid line")
                continue
            model_kmer = lst[9]
            samples = lst[-1]
            if read_name != prev_read:
                new_read = True
            else:
                new_read = False

            if prev_pos != pos:  # if new kmer
                if model_kmer == kmer:  # if model kmer is reference

                    if pos != prev_pos + 1 or new_read:  # if the consecutive 5-mer chain is broken or new read

                        if consecutive >= 5:  # process the previously built 9-mer if built

                            sample_nested_lst.append(sample_lst)
                            lines_that_yield_9mer.append(lines_that_yield_5mer)
                            
                            if counter_9mers[current_seq] < LIMIT_OUT:
                                counter_9mers[current_seq] += 1
                                cb, lb = write_trimmed_nanopolish(OUTPUT_file,lines_that_yield_9mer, current_seq,cb, lb, process_index)



                        #### reset 9-mer building as chain is broken
                        prev_read = read_name
                        prev_contig = contig
                        consecutive = 1
                        prev_pos = pos
                        current_seq = kmer
                        sample_nested_lst = []
                        lines_that_yield_9mer = []
                        sample_lst = [float(i) for i in samples.split(",")]
                        lines_that_yield_5mer = line

                    elif consecutive >= 5:
                        ### process the previously built 9-mer
                        sample_nested_lst.append(sample_lst)
                        lines_that_yield_9mer.append(lines_that_yield_5mer)

                        if counter_9mers[current_seq] < LIMIT_OUT:
                            counter_9mers[current_seq] += 1
                            cb, lb = write_trimmed_nanopolish(OUTPUT_file, lines_that_yield_9mer, current_seq,cb,lb, process_index)
                        #### execute vector operations

                        #### slide the 9-mer start
                        del sample_nested_lst[0]
                        del lines_that_yield_9mer[0]
                        current_seq = current_seq[1:] + kmer[4]
                        sample_lst = [float(i) for i in samples.split(",")]
                        lines_that_yield_5mer = line
                        prev_pos = pos
                        prev_read = read_name
                        prev_contig = contig
                        consecutive += 1

                    else:  # if building 6-9mer and new kmer seen
                        sample_nested_lst.append(sample_lst)  ### append currents for previous vector
                        lines_that_yield_9mer.append(lines_that_yield_5mer)

                        sample_lst = [float(i) for i in samples.split(",")]
                        lines_that_yield_5mer = line
                        prev_pos = pos
                        prev_read = read_name
                        prev_contig = contig
                        consecutive += 1
                        current_seq += kmer[4]

            elif model_kmer == kmer:  ### if another event for the already seen kmer position
                #    if consecutive >= 5 and current_seq[4] == "A":
                #           count_dwell_error["{}_{}".format(contig,pos)] +=1
                sample_lst += [float(i) for i in samples.split(",")]
                lines_that_yield_5mer+= line

    if consecutive >= 5:  # process the last 9-mer if built
        sample_nested_lst.append(sample_lst)
        lines_that_yield_9mer.append(lines_that_yield_5mer)

        if counter_9mers[current_seq] < LIMIT_OUT:
            counter_9mers[current_seq] += 1
            cb, lb = write_trimmed_nanopolish(OUTPUT_file, lines_that_yield_9mer, current_seq, 1000, lb, process_index)
    if cb > 0:
        
        with open(OUTPUT_file + "_trimmed.tsv", "a+") as trimmed_file:
            print("HERE")
            for lines_9mers in lb:
                for lines_5mer in lines_9mers:
                    trimmed_file.write(lines_5mer)
        #### execute vector operations

    print(counter_9mers.most_common(5), nanopolish_file, flush=True)
    print("Nanopolish parsed in", (time.time() - start_time_np) / 60, "minutes", flush=True)

    print(nanopolish_file, "done", flush=True)

    return counter_9mers


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='trim_nanopolish_v1.0', description=
    """ 
    This script takes nanopolish/f5c eventalign and trims it for input number of signals per 9mer suitable for training. 
    
    """, usage='python trim_eventalign-tsv.v1.0.py -i <nanopolish_file.txt> -o <out_dir> --limit-out 500  \nversion: %(prog)s')

    OPTIONAL = parser._action_groups.pop()

    REQUIRED = parser.add_argument_group('required arguments')

    REQUIRED.add_argument("-i", "--input_nanopolish",
                          help="Nanopolish file. Run nanopolish with the following flags: " \
                               " nanopolish eventalign --reads <in.fasta>" \
                               "--bam <in.bam> --genome <genome.fa> --print-read-names" \
                               "--scale-events --samples  > <out.txt>",
                          metavar='\b',
                          required=True)


    REQUIRED.add_argument("-o", "--out_path",
                          help="Output path",
                          metavar='\b',
                          required=True)

    OPTIONAL.add_argument("--limit_out",
                          help="Input the number of outputs per 9mer",
                          type = int,
                          default = 500)

    OPTIONAL.add_argument('-v', '--version',
                          action='version',
                          version='%(prog)s')

    OPTIONAL.add_argument("--out_counter",
                          help="Output counts of all 9mers for each split nanopolish file",
                          action="store_true")


    parser._action_groups.append(OPTIONAL)

    ARGS = parser.parse_args()
    # required arg
    NANOPOLISH_PATH = ARGS.input_nanopolish
    PATH_OUT = ARGS.out_path
    #optional args
    LIMIT_OUT = ARGS.limit_out
    OUT_COUNTER = ARGS.out_counter

    NUCLEOTIDES = ["A", "C", "G", "T"]
    start_time_total = time.time()

    final_counter = parse_nanopolish(NANOPOLISH_PATH)

    total_signals = sum(final_counter.values())
    print("Total of ", total_signals, " preprocessed")

    if OUT_COUNTER:
        with open(PATH_OUT + ".counts", "w+") as out_counts:
            out_counts.write("\n".join(["\t".join([str(x) for x in y]) for y in sorted(final_counter.items())]))

    print("All done in", (time.time() - start_time_total) / 60, "minutes")
