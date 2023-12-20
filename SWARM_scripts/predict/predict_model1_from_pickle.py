#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: stefan prodic
"""
import os
import sys
import time
import argparse
import subprocess
import numpy as np
import _pickle as cPickle
from tensorflow.keras import Input
from tensorflow.keras.models import Model
from network_27082022 import readwise_prediction_network



parser = argparse.ArgumentParser(prog='SWARM_predict_model1 v0.1', description=
""" 
This script takes an ID+signal file generated using WARM_preprocess.py and predict modification status \
per read and per 9-mer.

""", usage='python SWARM_preprocess.py -i <path_to_signlas+IDs_file> ' \
           '-m <path_to_DL_model> -l <label> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--signals_input",
                      help="path to the ID+signal file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to trainned model 1 DL model",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-l", "--label",
                      help="label of the condition of the sample, e.g. WT_rep1",
                      required=True)

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument("-r", "--resume",
                      action='store_true',
                      help="Continue predictions from previous file",
                      default=False,
                      required=False)

OPTIONAL.add_argument("--vector_len",
                          help="Size of the smooth window",
                          type=int,
                          default=36,
                          required=False)

OPTIONAL.add_argument("-f", "--features",
                          help="Number of model features",
                          type=int,
                          default=7,
                          required=False)

OPTIONAL.add_argument("--model_labels",
                          help="Number of model labels",
                          type=int,
                          default=2,
                          required=False)

OPTIONAL.add_argument("--limit",
                          help="Number of predictions",
                          type=int,
                          default=-1,
                          required=False)

OPTIONAL.add_argument("--batch_size",
                          help="Number of predictions loaded at a time on the GPU/CPU",
                          type=int,
                          default=20000,
                          required=False)

OPTIONAL.add_argument('-v', '--version',
                      action='version',
                      version='%(prog)s')

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
signals_input = ARGS.signals_input
DL_model = ARGS.DL_model
label = ARGS.label
file_out = ARGS.file_out
#optional args
vector_len=ARGS.vector_len
number_of_features = ARGS.features
number_of_labels = ARGS.model_labels
LIMIT = ARGS.limit
BATCH_SIZE =ARGS.batch_size
resume = ARGS.resume

# load the trainned model
inputs = Input(shape=(1,int(vector_len)*5, number_of_features))
output = readwise_prediction_network(inputs, number_of_labels)
model = Model(inputs=inputs, outputs=output)
model.load_weights(DL_model)

if LIMIT == -1:
    LIMIT = float("inf")

if resume is not True:
    if os.path.isfile(file_out):
        print('WARNING: read level prediction file already exists, please delete it or change the output name')
        sys.exit()

if resume is True:
    try:
        result = subprocess.run(["wc", "-l", file_out], stdout=subprocess.PIPE, text=True)
        total_lines = int(result.stdout.split(' ')[0])
        print('Previous number of predictions ', total_lines)
    except:
        print('Previous file not found: ' + file_out)
        total_lines = 0
else:
    total_lines = 0

### load the stored data
counter = 0
IDs = []
signals = []
start_loading = time.time()
with open(signals_input, 'rb') as signal_in:
    with open(file_out, 'a') as f_out:
        while True:
                counter += 1
                if counter > LIMIT:
                    if IDs:
                        predictions = model.predict(np.array(signals), batch_size = len(signals))
                        f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                               enumerate(predictions)])+"\n")

                    print('All signals have been processed', counter)
                    break

                if counter <= total_lines:
                    try:
                        seen_signal = cPickle.load(signal_in)
                    except:
                        break
                    continue
                try:
                    #key,signal = next(iter(cPickle.load(signal_in).items()))
                    pickle_dict = cPickle.load(signal_in)
                except Exception as ex:
                    print(ex)

                    if IDs:
                        #signal = np.array(signals)
                        predictions = model.predict(np.array(signals), batch_size = len(signals))
                        f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                               enumerate(predictions)])+"\n")

                    print('All signals have been processed', counter)
                    break

                else:
                    key, signal = next(iter(pickle_dict.items()))
                    IDs.append(key)
                    #signals.append(np.expand_dims(np.array(list(pickle_dict.values())),axis=1))
                    signals.append(np.expand_dims(signal, axis=0))
                    # to avoid loading everything predict every 20k singnals
                    if counter % BATCH_SIZE == 0:
                        if len(signals) > 0:
                            print("loading done in", time.time() - start_loading)
                            prediction_start = time.time()
                            predictions = model.predict(np.array(signals), batch_size = len(signals))
                            print("predicting done in", time.time() - prediction_start)
                            writing_start = time.time()
                            f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                   enumerate(predictions)])+"\n")
                            print("writing done in", time.time() - writing_start)
                            IDs, signals = [],[]
                            start_loading = time.time()
                        else:
                            continue
