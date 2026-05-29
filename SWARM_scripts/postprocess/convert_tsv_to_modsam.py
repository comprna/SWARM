import os
import sys
import time
import math
import pysam
import argparse




parser = argparse.ArgumentParser(prog='convert_model1_to_modsam v1.0')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--m1_input",
                      help="path to model1 prediction tsv output",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-b", "--bam",
                      help="path to input bam file used for f5C",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument("--model2",
                          help="Model2 predictions",
                          default=None)

OPTIONAL.add_argument("--cutoff",
                          help="Model2 cutoff",
                          default=0,
                          type=float)

OPTIONAL.add_argument("--min_m2_stoich",
                          help="Model2 cutoff",
                          default=0.1,
                          type=float)

OPTIONAL.add_argument('-v', '--version',
                      action='version',
                      version='%(prog)s')

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required args

file_out = ARGS.file_out



def get_MM_tag(MMs):
    point = 0
    MM = "N+N"
    for pos in MMs:
        delta = int(pos) - point
        MM+=f",{delta}"
        point+= 1 + delta
    return MM + ";"


model2_set = set()
if ARGS.model2:
    with open(ARGS.model2) as inf:
        inf.readline()
        for line in inf:
            line_lst = line.strip().split("\t")
            contig,pos = line_lst[:2]
            stoich,m2_p = list(map(float,line_lst[-2:]))
            if stoich <= ARGS.min_m2_stoich:
                m2_p = 0
            if m2_p > ARGS.cutoff:
                model2_set.add(f"{contig}_{pos}")


prev_read_name, MMs, MLs = "",[],[]
BAM_ALIGNMENT = pysam.AlignmentFile(ARGS.bam, "rb")
HEADER_BAM = BAM_ALIGNMENT.header.copy()
READ_INDEXED = pysam.IndexedReads(BAM_ALIGNMENT, multiple_iterators=True)
READ_INDEXED.build()
with open(ARGS.m1_input) as inf:
    with pysam.AlignmentFile(file_out + ".mod.sam", "wh", header=HEADER_BAM) as modsam_file:
        for line in inf:
            key, m1_p, label = line.strip().split("\t")
            m1_p = float(m1_p)
            key_lst = key.split("_")

            # last = ID_lst[-1]
            # qscore = ""
            # base = ""
            # ### Get qscore and base if they are in the key
            # if len(last) == 1 and not last.isnumeric():
            #     qscore = ID_lst.pop()
            #     base = ID_lst.pop()
            # if ID_lst[-1].isnumeric():  # handle within read position if --modsam was used
            #     current_ID = '_'.join(ID_lst[:-2])
            # else:
            #     current_ID = '_'.join(ID_lst[:-1])

            last = key_lst[-1]
            if len(last) == 1 and not last.isnumeric():
                qscore = key_lst.pop()
                base = key_lst.pop()
            if not key_lst[-1].isnumeric():
                raise ("Preprocessing format not supported. ReadID must be followed by integer. Exiting")
            kmer, readName, readPos = key_lst[-3:]
            site_index = "_".join(key_lst[:-3])
            if readName == prev_read_name:
                if int(readPos) > 0 or len(MMs) == 0:
                    try:
                        if not ARGS.model2 or site_index in model2_set:
                            col=math.floor(m1_p*255)
                        else:
                            col = 0
                    except Exception as e:
                        print("skipped Nan probability")
                    else:
                        MMs.append(readPos)
                        MLs.append(col)
            elif MMs:
                input_read.set_tag("MM", get_MM_tag(MMs), value_type="Z")
                input_read.set_tag("ML", MLs)

                modsam_file.write(input_read)
                try:
                    input_read = next(READ_INDEXED.find(readName))  # get the read pysam object from the bam
                except:
                    raise(readName + " not present in the provided bam. Exiting.")
                prev_read_name = input_read.query_name
                MMs.clear()
                MLs.clear()
            else:
                try:
                    input_read = next(READ_INDEXED.find(readName))  # get the read pysam object from the bam
                except:
                    raise(readName + " not present in the provided bam. Exiting.")
                prev_read_name = input_read.query_name
                MMs.clear()
                MLs.clear()

BAM_ALIGNMENT.close()





