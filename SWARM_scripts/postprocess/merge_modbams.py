import os
import sys
import math
import time
import pysam
import argparse
import numpy as np


MAXPROB = 0.9999
MAXSCALED = -np.log10(1-MAXPROB)

start = time.time()
def scale_p(x):
    # Apply a logarithmic transformation to probabilities
    if x == 255:
        return 255
    p=x/255
    scaled_values = -np.log10(1-p )
    # Scale the transformed values to the range [0, 255]
    scaled_values = math.ceil((scaled_values /MAXSCALED) * 255)

    return scaled_values



parser = argparse.ArgumentParser()
OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

OPTIONAL.add_argument("-i", "--input_bam",
                      metavar='\b',
                      required=False)

REQUIRED.add_argument("-o", "--out_path",
                      help="Output path",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-m', '--merge',
                      default=None,
                      required=False)

OPTIONAL.add_argument('-t', '--threshold',
                      default=0.9,
                      type=float,
                      required=False)

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

mods_dict = {"pU":["T","17802"], "m5C":["C","c"],"m6A":["A","a"]}
input_files,input_mods = [],[]

if ARGS.merge:
    with open(ARGS.merge) as f:
        for line in f:
            path,mod = line.strip().split(",")
            input_files.append(path)
            input_mods.append(mod)
else:
    input_mods = ["pU"]
    input_files = [ARGS.input_bam]


file_out = ARGS.out_path
original_bam = input_files[0]
THRESHOLD = ARGS.threshold
bam_objects_lst = []
creads = 0
for i,SAM_PATH in enumerate(input_files):
    BAM_ALIGNMENT = pysam.AlignmentFile(SAM_PATH, "rb")
    HEADER_BAM = BAM_ALIGNMENT.header.copy()
    READ_INDEXED = pysam.IndexedReads(BAM_ALIGNMENT)
    READ_INDEXED.build()
    bam_objects_lst.append(READ_INDEXED)

with pysam.AlignmentFile(original_bam,"rb") as input_bamfile:
    with pysam.AlignmentFile(file_out, "wb", template=input_bamfile) as merged_modbam_file:
        for input_read in input_bamfile:
            read_name = input_read.query_name
            merged_MM = ""
            merged_ML = []
            for i,mod in enumerate(input_mods):
                try:
                    read_object = next(bam_objects_lst[i].find(read_name))  # get the read pysam object from the bam
                except:
                    print(read_name + "  skipped, not present in the provided bam for",mod, flush=True)
                    pass
                else:
                    MM = read_object.get_tag("MM")
                    if len(MM[4:]) == 0:
                        print(read_name, "MM tag empty in bam for", mod, flush=True)
                    else:
                        merged_MM = f"{merged_MM}N+{mods_dict[mod][1]}?,{MM[4:]}"
                        ML=read_object.get_tag("ML")
                        merged_ML+=[x if x > THRESHOLD*255 else 0 for x in ML]

            if merged_MM:
                input_read.set_tag("MM", merged_MM, value_type="Z")
                input_read.set_tag("ML", merged_ML)
                merged_modbam_file.write(input_read)
                pass
                    # print(ML)
            else:
                print(read_name)
            creads+=1
            print(creads, time.time() - start)




