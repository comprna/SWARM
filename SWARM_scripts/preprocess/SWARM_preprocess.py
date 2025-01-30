"""
@author: Stefan Prodic


Same as 2.033 but outputs qscore_readBase for target base at the end of the signal's key (for downstream multi-feature model2). 0_D if deleted

This code takes a nanopolish file and a bam file to generate the input for the
4DDQ model (current,distances from A/C/G/T bases, dwelling time, q-score of the reference base)
Current smoothing: pedestrian interpolation and compression (fixed current reversing!)
q-score = avg 9mer (only non-indel bases) if read aligned with a gap at that position
Vector length needs to be divisible by 9 to evenly output all q-scores within a 9mer
"""

import os
import sys
import time
import math
import pysam
import random
import argparse
import numpy as np
import pandas as pd
import _pickle as cPickle
from random import shuffle
from collections import Counter
from multiprocessing import Pool


def split_nanopolish(file, n_split, temp_path, file_len=None):
    """
    Splits the nanopolish file in desired number of chunks at 'breakpoints'
    where a consecutive chain of 5-mer positions is broken
    """
    print("splitting nanopolish")
    start_time_np = time.time()
    if file_len == None:
        with open(file) as f:
            f.readline()
            file_len = 0
            for line in f:
                file_len += 1
    temp_path += os.path.basename(file)  # output split file as [outputdir]/tmp/[nanopolish]_[index]
    index, i, previous_pos = 0, 0, 0
    step = file_len // n_split
    limit = step
    output_file = open(temp_path + "_" + str(index), "w+")
    outputs_list = [temp_path + "_" + str(index)]
    with open(file) as np_file:
        np_file.readline()
        for line in np_file:
            i += 1
            read, current_pos, = line.split("\t")[:2]
            try:
                current_pos = int(current_pos)
            except ValueError:
                continue

            if i > limit and (current_pos > previous_pos + 1 or current_pos < previous_pos):
                limit += step
                index += 1
                output_file.close()
                output_file = open(temp_path + "_" + str(index), "w+")
                outputs_list.append(temp_path + "_" + str(index))

            output_file.write(line)
            previous_pos = current_pos

    output_file.close()
    print("splitting done in ", (time.time() - start_time_np) / 60, "minutes")

    return outputs_list


def get_qscores(indexed_reads, read_name, contig, temp_path, HEADER_BAM, suffix=""):
    """
    Finds the read from the pysam read-indexed alignment
    Writes and indexes a new bam file with only the desired read
    Returns a dictionary storing qscores and bases at reference positions
    """

    try:
        read_object = next(indexed_reads.find(read_name))  # get the read pysam object from the bam
    except:
        print(read_name + "  skipped, not present in the provided bam.", flush=True)
        return {}

    read_bam_name = temp_path + read_name + "_" + suffix + ".bam"
    with pysam.Samfile(read_bam_name, "wb", header=HEADER_BAM) as read_bam_file:  # write a bam containing the read
        read_bam_file.write(read_object)
    pysam.index(read_bam_name)  # index the read bam
    read_alignment = pysam.AlignmentFile(read_bam_name, "rb")  # load the indexed read
    output_dict = {}
    for pileupcol in read_alignment.pileup(contig, min_base_quality=0,
                                           max_depth=9999999):  # go through reference positions
        ref_position = pileupcol.pos
        for pileupread in pileupcol.pileups:
            if pileupread.query_position:  # if not an indel store basecall qscore and called base
                output_dict[ref_position] = (pileupread.alignment.query_qualities[pileupread.query_position],
                                             pileupread.alignment.query_sequence[pileupread.query_position])
    try:
        os.remove(read_bam_name)
        os.remove(read_bam_name + ".bai")
        return output_dict

    except FileNotFoundError:
        print("tmp bam or bai not deleted", read_bam_name)
        return output_dict


def get_ref_qscore(ref_base, read_base, qscore):
    """
    Computes the qscore of the reference base
    """
    if ref_base == read_base:
        return qscore
    else:
        try:
            p = 10 ** ((-1 * qscore) / 10)  # get probability of nonref base from qscore
            return -10 * math.log10(1 - p)  # return qscore of ref base
        except:
            print("Invalid qscore", qscore, flush=True)
            return 0


def interpolate_signal(original_vector, target_scale):
    """
        Pedestrian smoothing devised by Eduardo Eyras
        """
    # get original vector length
    original_len = len(original_vector)
    interp_op = original_len - 1
    interp_target = target_scale - 1
    # get sliding window size
    min_step_size = interp_target // interp_op

    # get increased sliding window needed to capture remaining values
    max_step_size = min_step_size + 1
    remainder = interp_target % interp_op

    # get a list of all window sizes needed to compress to desired scale
    window_values = [max_step_size] * remainder + [min_step_size] * (interp_op - remainder)

    shuffle(window_values)  # shuffle the list for randomly distributed smoothness
    index = 0
    output_vector = []

    for index, window_size in enumerate(window_values):  # get average for each window to compress original values
        low_value = original_vector[index]
        up_value = original_vector[index + 1]
        increment = (up_value - low_value) / (window_size)
        add = 0
        for _ in range(window_size):
            output_vector.append(low_value + add)
            add += increment
    output_vector.append(original_vector[-1])
    return output_vector


def smooth_event(original_vector, target_scale):
    """
    Pedestrian smoothing devised by Eduardo Eyras
    """

    # original_vector.reverse()  # reverse to account for 3'->5' direction of nanopore

    if len(original_vector) < target_scale:
        new_event = interpolate_signal(original_vector, target_scale)
        return new_event
    else:
        # get original vector length
        original_len = len(original_vector)

        # get sliding window size
        min_step_size = original_len // target_scale

        # get increased sliding window needed to capture remaining values
        max_step_size = min_step_size + 1
        remainder = original_len % target_scale

        # get a list of all window sizes needed to compress to desired scale
        window_values = [max_step_size] * remainder + [min_step_size] * (target_scale - remainder)

        random.shuffle(window_values)  # shuffle the list for randomly distributed smoothness
        index = 0
        output_vector = []

        for window_size in window_values:  # get average for each window to compress original values
            output_vector.append(sum(original_vector[index: index + window_size]) / window_size)
            # output_vector.append(np.median(original_vector[index: index + window_size]))
            index += window_size

        return output_vector


def get_4DDQ_list(smooth_event_2D, nine_mer, original_lengths, expected_dict, qscore_dict, start_pos, length=20):
    """
    outputs the vector used for the SWARM model 1
    [[current, distances from 4 bases, dwelling time, qscore(ref nucleotide)] x length x 5]
    """

    reference_9mers = [nine_mer[:4] + nucleotide + nine_mer[5:] for nucleotide in NUCLEOTIDES]
    out_lst = []
    qscore_lst = []
    qscore_non_0 = []
    for i in range(9):
        ref_base = nine_mer[i]
        try:
            guppy_qscore, guppy_base = qscore_dict[start_pos + i]

        except KeyError:  # if indel within 9mer in this read (or read skipped)
            ref_qscore = 0

        else:
            ref_qscore = get_ref_qscore(ref_base, guppy_base, guppy_qscore)  # get the reference q score if mismatch
            qscore_non_0.append(ref_qscore)
        # if ref_qscore is None:    # skip if q score cannot be obtained
        #     ref_qscore = 0

        qscore_lst.append(ref_qscore)
    if len(qscore_non_0) == 0:
        return None
    avg_q = sum(qscore_non_0) / len(qscore_non_0)
    qscore_lst = [y if y != 0 else avg_q for y in qscore_lst]
    qscore_step = (5 * length) // 9  # calculate the number of reps for the same qscore to fit it into the vector size


    for kmer_i in range(5):  # build the vector for the CNN input
        kmer_signals = smooth_event_2D[kmer_i]
        for current_index, current in enumerate(kmer_signals):
            inner_lst = [current]
            for ref_9mer in reference_9mers:
                ref_5mer = ref_9mer[kmer_i:kmer_i + 5]
                expected_signal, expected_std = expected_dict[ref_5mer]
                inner_lst.append(expected_signal - current)
            inner_lst.append(original_lengths[kmer_i])
            inner_lst.append(qscore_lst[(kmer_i * length + current_index) // qscore_step])
            out_lst.append(inner_lst)

    return np.array([[round(y, 3) for y in x] for x in out_lst])


def handle_ninemer(contig, pos, kmer_seq, read_id, sample_nested_lst,
                   model_kmer_dict, counter_dict, qscore_dict, length=20, nucleotide_mode="ALL4",
                   OUTPUT_file=None, outs_dict=None):
    if LIMIT_OUT:  ### don't preprocess if reached limit of outputed signals for this 9mer in this process
        if counter_dict[kmer_seq] == LIMIT_OUT:
            return None

    ### handle outputs for a single base or multiple bases
    if nucleotide_mode != "ALL4":
        if kmer_seq[4] != nucleotide_mode:
            return None
        else:

            start_9mer_pos = pos - 4
            length_lst = [len(events) for events in sample_nested_lst]  # store dwelling time before smoothing
            combined_vectors = get_4DDQ_list(
                [smooth_event(events, length) for events in sample_nested_lst],
                kmer_seq,
                length_lst,
                model_kmer_dict,
                qscore_dict,
                start_9mer_pos,
                length)
            if combined_vectors is not None:
                if pos in qscore_dict:
                    q_base = f"{qscore_dict[pos][0]}_{qscore_dict[pos][1]}"
                else:
                    q_base = "0_D"
                with open(OUTPUT_file + ".pickle",
                          "ab") as out_file:  # output vectors (position of the start of the 9mer)
                    cPickle.dump({"_".join([contig, str(start_9mer_pos), kmer_seq, read_id,q_base]): combined_vectors},
                                 out_file)
                counter_dict[kmer_seq] += 1
    else:

        start_9mer_pos = pos - 4
        length_lst = [len(event) for event in sample_nested_lst]
        combined_vectors = get_4DDQ_list(
            [smooth_event(events, length) for events in sample_nested_lst],
            kmer_seq,
            length_lst,
            model_kmer_dict,
            qscore_dict,
            start_9mer_pos,
            length)
        if combined_vectors is not None:
            # print("_".join([contig, str(start_9mer_pos), kmer_seq, read_id]))
            # print(combined_vectors)
            if pos in qscore_dict:
                    q_base = f"{qscore_dict[pos][0]}_{qscore_dict[pos][1]}"
            else:
                    q_base = "0_D"
            with open(outs_dict[kmer_seq[4]], "ab") as out_file:
                cPickle.dump({"_".join([contig, str(start_9mer_pos), kmer_seq, read_id,q_base]): combined_vectors}, out_file)
            counter_dict[kmer_seq] += 1
    return combined_vectors


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

    if INPUT_BASE == "ALL4":
        outs_dict = {nucleotide: "{}_{}.pickle".format(OUTPUT_file, nucleotide) for nucleotide in ["A", "G", "T", "C"]}
    else:
        outs_dict = None
        OUTPUT_file += "_" + INPUT_BASE

    if CPU_NUMBER > 1:
        process_index = str(nanopolish_file).split("_")[-1]
    else:
        process_index = ""

    BAM_ALIGNMENT = pysam.AlignmentFile(PATH_BAM, "rb")
    HEADER_BAM = BAM_ALIGNMENT.header.copy()
    READ_INDEXED = pysam.IndexedReads(BAM_ALIGNMENT, multiple_iterators=True)
    READ_INDEXED.build()

    TARGETS_DICT = {}
    if ARGS.positions:
        print("DOING positions")
        with open(ARGS.positions) as pos_bed:
            pos_bed.readline()
            for line in pos_bed:
                chr, pos = line.strip().split("\t")[:2]
                TARGETS_DICT[f"{chr}_{pos}"] = 0




    # print("BAM INDEXED BY READS")
    with open(nanopolish_file) as nf:
        prev_pos = -9999
        prev_read, prev_contig, current_seq = "", "", ""
        nf.readline()
        consecutive = 0
        sample_lst = []
        c = 0
        for line in nf:
            c += 1
            if c % 1000 == 0:
                print(c, time.time() - start_time_np, flush=True)
            lst = line.split("\t")
            contig, pos, kmer, read_name = lst[:4]
            pos = int(pos)
            model_kmer = lst[9]
            samples = lst[-1]
            skip_read = False
            if read_name != prev_read:
                new_read = True
                if prev_read:
                    prev_qscore_dict = qscore_dict
                    qscore_dict = get_qscores(READ_INDEXED, read_name, contig, TEMP_PATH, HEADER_BAM, process_index)
                else:
                    qscore_dict = get_qscores(READ_INDEXED, read_name, contig, TEMP_PATH, HEADER_BAM, process_index)

                # check if read is not in the bam (if qscore dict is empty)
                # if so, preprocess last 9mer of the previous read and skip preprocessing this read
                if not qscore_dict:
                    skip_read = True
                else:
                    skip_read = False
            else:
                new_read = False
                if skip_read:
                    continue

            if prev_pos != pos:  # if new kmer
                if True:  # if model kmer is reference

                    if pos != prev_pos + 1 or new_read:  # if the consecutive 5-mer chain is broken or new read

                        if consecutive >= 5:  # process the previously built 9-mer if built

                            sample_nested_lst.append(sample_lst)

                            if new_read:  #### execute vector operations  with previous read information

                                if not TARGETS_DICT or f"{prev_contig}_{prev_pos}" in TARGETS_DICT:

                                    combined_vector = handle_ninemer(prev_contig,
                                                                     prev_pos,
                                                                     current_seq,
                                                                     prev_read,
                                                                     sample_nested_lst,
                                                                     MODEL_KMER_DICT,
                                                                     counter_9mers,
                                                                     prev_qscore_dict,
                                                                     length=VECTOR_LEN,
                                                                     nucleotide_mode=INPUT_BASE,
                                                                     OUTPUT_file=OUTPUT_file,
                                                                     outs_dict=outs_dict)
                                    if combined_vector is not None and TARGETS_DICT:
                                        TARGETS_DICT[f"{prev_contig}_{prev_pos}"]+=1
                                # print(combined_vector, "new read", current_seq, sep="\n")
                            else:  #### execute vector operations  with current read information

                                if not TARGETS_DICT or f"{prev_contig}_{prev_pos}" in TARGETS_DICT:
                                    combined_vector = handle_ninemer(prev_contig,
                                                                     prev_pos,
                                                                     current_seq,
                                                                     prev_read,
                                                                     sample_nested_lst,
                                                                     MODEL_KMER_DICT,
                                                                     counter_9mers,
                                                                     qscore_dict,
                                                                     length=VECTOR_LEN,
                                                                     nucleotide_mode=INPUT_BASE,
                                                                     OUTPUT_file=OUTPUT_file,
                                                                     outs_dict=outs_dict)
                                    if combined_vector is not None and TARGETS_DICT:
                                        TARGETS_DICT[f"{prev_contig}_{prev_pos}"]+=1
                            # print(combined_vector, "broken chain", current_seq , sep = "\n")

                        #### reset 9-mer building as chain is broken
                        prev_read = read_name
                        prev_contig = contig
                        consecutive = 1
                        prev_pos = pos
                        current_seq = kmer
                        sample_nested_lst = []

                        if NOISY:
                            sample_lst = [float(i) + random.randint(-NOISE_LEVEL, NOISE_LEVEL) for i in
                                          samples.split(",")]
                        else:
                            sample_lst = [float(i) for i in samples.split(",")]

                        sample_lst.reverse()

                    elif consecutive >= 5:
                        ### process the previously built 9-mer
                        sample_nested_lst.append(sample_lst)
                        #### execute vector operations

                        if not TARGETS_DICT or f"{prev_contig}_{prev_pos}" in TARGETS_DICT:

                            combined_vector = handle_ninemer(prev_contig,
                                                             prev_pos,
                                                             current_seq,
                                                             prev_read,
                                                             sample_nested_lst,
                                                             MODEL_KMER_DICT,
                                                             counter_9mers,
                                                             qscore_dict,
                                                             length=VECTOR_LEN,
                                                             nucleotide_mode=INPUT_BASE,
                                                             OUTPUT_file=OUTPUT_file,
                                                             outs_dict=outs_dict)
                            if combined_vector is not None and TARGETS_DICT:
                                TARGETS_DICT[f"{prev_contig}_{prev_pos}"] += 1

                        # print(current_seq, "c is 5 new 9mer", "\n".join([str(x) for x in combined_vector]), prev_pos, sep="\n")
                        # sys.exit()
                        #### slide the 9-mer start
                        del sample_nested_lst[0]
                        current_seq = current_seq[1:] + kmer[4]
                        if NOISY:
                            sample_lst = [float(i) + random.randint(-NOISE_LEVEL, NOISE_LEVEL) for i in
                                          samples.split(",")]
                        else:
                            sample_lst = [float(i) for i in samples.split(",")]
                        sample_lst.reverse()

                        prev_pos = pos
                        prev_read = read_name
                        prev_contig = contig
                        consecutive += 1

                    else:  # if building 6-9mer and new kmer seen
                        sample_nested_lst.append(sample_lst)  ### append currents for previous vector
                        if NOISY:
                            sample_lst = [float(i) + random.randint(-NOISE_LEVEL, NOISE_LEVEL) for i in
                                          samples.split(",")]
                        else:
                            sample_lst = [float(i) for i in samples.split(",")]
                        sample_lst.reverse()
                        prev_pos = pos
                        prev_read = read_name
                        prev_contig = contig
                        consecutive += 1
                        current_seq += kmer[4]

            elif True:  ### if another event for the already seen kmer position
                #    if consecutive >= 5 and current_seq[4] == "A":
                #           count_dwell_error["{}_{}".format(contig,pos)] +=1
                if NOISY:
                    rev_lst = [float(i) + random.randint(-NOISE_LEVEL, NOISE_LEVEL) for i in samples.split(",")]
                else:
                    rev_lst = [float(i) for i in samples.split(",")]
                rev_lst.reverse()
                sample_lst += rev_lst

    if consecutive >= 5:  # process the last 9-mer if built
        sample_nested_lst.append(sample_lst)

        #### execute vector operations
        if not TARGETS_DICT or f"{prev_contig}_{prev_pos}" in TARGETS_DICT:
            combined_vector = handle_ninemer(prev_contig,
                                             prev_pos,
                                             current_seq,
                                             prev_read,
                                             sample_nested_lst,
                                             MODEL_KMER_DICT,
                                             counter_9mers,
                                             qscore_dict,
                                             length=VECTOR_LEN,
                                             nucleotide_mode=INPUT_BASE,
                                             OUTPUT_file=OUTPUT_file, outs_dict=outs_dict)

            if combined_vector is not None and TARGETS_DICT:
                TARGETS_DICT[f"{prev_contig}_{prev_pos}"] += 1
        # print(combined_vector, "c is 5 last", current_seq, sep="\n")

    print("Nanopolish parsed in", (time.time() - start_time_np) / 60, "minutes", flush=True)

    print(nanopolish_file, "done", flush=True)

    return (counter_9mers,TARGETS_DICT)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='SWARM_preprocess_4D_D_Qscore_SP_v1.0', description=
    """ 
    This script takes nanopolish and bam files to extract 
    signals from 9-mers centered around desired nucleotides.
    It creates a file/s that will be the input of the predict_model_1.py 
    using the 7 features (current, 4 distances, dwelling time, and qscore).
    """, usage='python SWARM_preprocess_4D_D_Qscore_SP_v1.0.py -i <nanopolish_file.txt> ' \
               '-m <kmer_model> -o <out_dir> -b <bam_file> -n 10 --input_base A    \nversion: %(prog)s')

    OPTIONAL = parser._action_groups.pop()

    REQUIRED = parser.add_argument_group('required arguments')

    REQUIRED.add_argument("-i", "--input_nanopolish",
                          help="Nanopolish file. Run nanopolish with the following flags: " \
                               " nanopolish eventalign --reads <in.fasta>" \
                               "--bam <in.bam> --genome <genome.fa> --print-read-names" \
                               "--scale-events --samples  > <out.txt>",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("-m", "--kmer_model",
                          help="File containing all the expected signal k-mer means",
                          metavar='\b',
                          required=True)


    REQUIRED.add_argument("-b", "--bam",
                          help="Alignment file in bam format",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("-o", "--out_path",
                          help="Output path",
                          metavar='\b',
                          required=True)

    OPTIONAL.add_argument('-v', '--version',
                          action='version',
                          version='%(prog)s')

    OPTIONAL.add_argument("-n", "--cpu",
                          help='Number of cores to use',
                          type=int,
                          default=1
                          )
    OPTIONAL.add_argument("-l", "--vlen",
                          help='Number of smooth current measurements per 5-mer',
                          type=int,
                          default=20
                          )

    OPTIONAL.add_argument("--input_base",
                          help="Input base as desired centre of the 9mers. Default is ALL4 to preprocess all 4 bases",
                          metavar='\b',
                          default="ALL4")

    OPTIONAL.add_argument("-w", "--wc_l",
                          help="Input the length of the nanopolish file",
                          metavar='\b',
                          default=None)

    OPTIONAL.add_argument("--limit_out",
                          help="Input the number of outputs per 9mer within chunks of the split nanopolish file",
                          metavar='\b',
                          default=None)

    OPTIONAL.add_argument("--out_counter",
                          help="Output counts of all 9mers for each split nanopolish file",
                          action="store_true")

    OPTIONAL.add_argument("--current_noise",
                          help="Introduce noise to current",
                          action="store_true")

    OPTIONAL.add_argument("--noise_level",
                          help="Level of noise",
                          default=0,
                          type=int)

    OPTIONAL.add_argument("--positions",
                          help="Target positions in .bed format (with column names)",
                          default=None)

    parser._action_groups.append(OPTIONAL)

    ARGS = parser.parse_args()

    # required arg
    NANOPOLISH_PATH = ARGS.input_nanopolish
    MODEL_KMER_PATH = ARGS.kmer_model
    PATH_OUT = ARGS.out_path
    PATH_BAM = ARGS.bam
    # optional args
    CPU_NUMBER = ARGS.cpu
    VECTOR_LEN = ARGS.vlen
    INPUT_BASE = ARGS.input_base
    FILE_LEN = ARGS.wc_l
    LIMIT_OUT = ARGS.limit_out
    OUT_COUNTER = ARGS.out_counter
    NOISY = ARGS.current_noise
    NOISE_LEVEL = ARGS.noise_level

    NUCLEOTIDES = ["A", "C", "G", "T"]

    start_time_total = time.time()
    TEMP_PATH = os.path.dirname(PATH_OUT) + "/tmp/"

    if not os.path.exists(TEMP_PATH):
        os.mkdir(TEMP_PATH)

    model_kmer = pd.read_csv(MODEL_KMER_PATH, sep=',')
    # MODEL_KMER_DICT = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    MODEL_KMER_DICT = model_kmer.set_index('model_kmer').apply(
        lambda x: [x['model_mean'], x['model_stdv']], axis=1).to_dict()

    if LIMIT_OUT:
        if not LIMIT_OUT.isdigit():
            raise Exception(" --limit_out needs to be an integer")

    if VECTOR_LEN % 9 != 0:
        raise Exception(" -l needs to be a multiple of 9 to fit qscores")



    if CPU_NUMBER == 1:
        parse_nanopolish(NANOPOLISH_PATH)

    else:
        temp_files = split_nanopolish(NANOPOLISH_PATH, CPU_NUMBER, TEMP_PATH, FILE_LEN)
        print(temp_files, CPU_NUMBER)
        p = Pool(CPU_NUMBER)
        time.sleep(1)
        outs = p.map(parse_nanopolish, temp_files)
        counters = [x[0] for x in outs]
        targets = [x[1] for x in outs]

        if OUT_COUNTER:
            final_counter = Counter()
            for counter in counters:
                final_counter.update(counter)
            with open(PATH_OUT + ".counts", "w+") as out_counts:
                out_counts.write("\n".join(["\t".join([str(x) for x in y]) for y in sorted(final_counter.items())]))


        if ARGS.positions:
            total_targets = Counter()
            for counter in targets:
                total_targets.update(counter)

            with open(ARGS.positions) as in_bed:
                with open(PATH_OUT + ".targets.bed","w+") as out_bed:
                    out_bed.write(in_bed.readline().strip() + "\tsignals\n")
                    for line in in_bed:
                        chr,pos = line.strip().split("\t")[:2]
                        out_bed.write(line.strip()  + "\t" +str(total_targets[f"{chr}_{pos}"]) + "\n")

        time.sleep(1)
        p.close()

        if CPU_NUMBER > 1:
            for tmp in temp_files:
                os.remove(tmp)

    print("All done in", (time.time() - start_time_total) / 60, "minutes")
