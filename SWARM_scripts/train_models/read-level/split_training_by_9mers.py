import sys
import argparse
import numpy as np
import _pickle as cPickle
from collections import Counter

"""
@author: Stefan Prodic
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='split_training_by_9mers', description=
    """ 
    """, usage='python split_training_by_9mers.py -i <SM11_preprocessed.p> --counts <SM23_IVT-i1c1_i2c1.counts> -l 500 -o <output_path> ')

    OPTIONAL = parser._action_groups.pop()

    REQUIRED = parser.add_argument_group('required arguments')


    REQUIRED.add_argument("-i", "--input",
                          help="Preprocessed pickle for the desired modification",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("--counts",
                          help="signals per 9mer for the desired modification",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("--limit",
                          help="Maximal number of signals per 9mer",
                          type=int,
                          required=True)

    REQUIRED.add_argument("-o", "--out_path",
                          help="Output path",
                          metavar='\b',
                          required=True)

    OPTIONAL.add_argument("--train_percent",
                          help="Input the proportion of training data",
                          type=float,
                          default = 0.6)

    OPTIONAL.add_argument("--validate_percent",
                          help="Input the proportion of validation data",
                          type=float,
                          default=0.2)


    parser._action_groups.append(OPTIONAL)

    ARGS = parser.parse_args()

    # required arg
    INPUT_PATH = ARGS.input
    COUNTS_PATH = ARGS.counts
    PATH_OUT = ARGS.out_path
    LIMIT = ARGS.limit
    # optional args
    TRAIN_PERCENT = ARGS.train_percent
    VALIDATE_PERCENT = ARGS.validate_percent


    limits_dict = {}
    with open(COUNTS_PATH) as f:
        for line in f:
            ninemer, count = line.strip().split("\t")
            count = min( int(count), LIMIT)
            limits_dict[ninemer] = [count * TRAIN_PERCENT, count * (TRAIN_PERCENT + VALIDATE_PERCENT)]


    with open(INPUT_PATH , "rb") as pf:

        mod_c_new = Counter()
        while True:
            try:
                dct = cPickle.load(pf)
                key,item = next(iter(dct.items()))
            except:
                break

            ninemer = key.split("_")[-2]
            #print(ninemer)
            train_limit, validate_limit = limits_dict[ninemer]
            mod_c_new[ninemer] += 1
            seen_count = mod_c_new[ninemer]
            if seen_count <= train_limit:

                with open(PATH_OUT + "_training_data.p", "ab") as of:
                    cPickle.dump(item,of)

            elif seen_count <= validate_limit:
                with open(PATH_OUT + "_validation_data.p", "ab") as of:
                    cPickle.dump(item,of)

            elif seen_count <= LIMIT:
                with open(PATH_OUT + "_testing_data.p", "ab") as of:
                    cPickle.dump(item,of)
