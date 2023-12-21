#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Stefan Prodic
"""

import sys
import pickle
import argparse
import numpy as np
import _pickle as cPickle

def load_all_data(path):
    elements = []
    with open(path, 'rb') as signal_in:
        while True:
            try:
                elements.append(cPickle.load(signal_in))
                
            except Exception as e:
                print(e)
                print(path)
                break
    return elements


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='assemble_data', description=
    """ 
    """, usage='python assemble_data.py --input_positive </path/m6A> ' \
               '--input_negative </path/nonmod> -o <out_dir>')

    OPTIONAL = parser._action_groups.pop()

    REQUIRED = parser.add_argument_group('required arguments')


    REQUIRED.add_argument("--input_positive",
                          help="Split data prefix for the desired modification",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("--input_negative",
                          help="Split data prefix for the unmodified or negative modification/s",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("-o", "--out_path",
                          help="Output path",
                          metavar='\b',
                          required=True)

    OPTIONAL.add_argument("--positive_label",
                          help="Label for the positive modification",
                          default=1,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--negative_label",
                          help="Labels for the negatives",
                          default=0,
                          type = int,
                          required=False)




    parser._action_groups.append(OPTIONAL)

    ARGS = parser.parse_args()

    # required arg
    INPUT_POSITIVE = ARGS.input_positive
    INPUT_NEGATIVE = ARGS.input_negative
    PATH_OUT = ARGS.out_path
    # optional args
    POSITIVE_LABEL = ARGS.positive_label
    NEGATIVE_LABEL = ARGS.negative_label

    for data_set in ["training", "validation", "testing"]:
        positive_data = load_all_data("{}_{}_data.p".format(INPUT_POSITIVE, data_set))
        positive_labels = [POSITIVE_LABEL for _ in range(len(positive_data))]

        negative_data = load_all_data("{}_{}_data.p".format(INPUT_NEGATIVE, data_set))
        negative_labels = [NEGATIVE_LABEL for _ in range(len(negative_data))]
        print("positive data/labels",len(positive_data), len(positive_labels),"negative data/labels",len(negative_data),len(negative_labels))

        data = np.array(positive_data + negative_data)
        labels = np.array(positive_labels + negative_labels)
        with open("{}_final_{}_data.p".format(PATH_OUT, data_set), "ab") as of:
            pickle.dump(data, of)

        with open("{}_final_{}_labels.p".format(PATH_OUT, data_set), "ab") as of:
            pickle.dump(labels, of)


