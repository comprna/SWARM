#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: stefan
"""
import os
import sys
import time
import struct
import argparse
import subprocess
import numpy as np
from tensorflow.keras import Input
from tensorflow.keras.models import Model
from network_27082022 import readwise_prediction_network



parser = argparse.ArgumentParser(prog='predict_model1_parallel v1.0', description=
""" 
This script takes an ID+signal file generated using SWARM_preprocess_* and predicts RNA modification status \
per read and per 9-mer.

""", usage='python predict_model1_parallel -i <path_to_signlas+IDs_file> ' \
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
                          help="Limit to a given number of predictions",
                          type=int,
                          default=-1,
                          required=False)

OPTIONAL.add_argument("--batch_size",
                          help="Number of predictions loaded at a time on the GPU/CPU",
                          type=int,
                          default=4096,
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
BATCH_SIZE = ARGS.batch_size


def read_signal(file):
    char_len = struct.unpack('i', file.read(4))[0]
    contigName = file.read(char_len).decode('utf-8')

    ninemerPos = str(struct.unpack('i', file.read(4))[0])

    ninemer = file.read(9).decode('utf-8')

    char_len = struct.unpack('i', file.read(4))[0]
    readName = file.read(char_len).decode('utf-8')
    readPos = str(struct.unpack('i', file.read(4))[0])
    key = "_".join([str(x) for x in [contigName,ninemerPos,ninemer,readName,readPos]])
    #print(char_len, name)
    # Read the 2D double array (float values) from the file
    array_data = file.read(4 * 180 * 7)
    return key, np.frombuffer(array_data, dtype=np.float32).reshape((1,180, 7))




# load the trainned model
inputs = Input(shape=(1,int(vector_len)*5, number_of_features))
output = readwise_prediction_network(inputs, number_of_labels)
model = Model(inputs=inputs, outputs=output)
model.load_weights(DL_model)

if LIMIT == -1:
    LIMIT = float("inf")


### load the stored data
counter = 0
IDs = []
signals = []
start_loading = time.time()
PREP_I = 0
in_dir = os.path.dirname(signals_input)

waited = 0
WAIT_LIMIT = 600
WAIT_INTERVAL=0.01
DONE=False
with open(file_out, 'a') as f_out:
    while True:
        files = os.listdir(in_dir)
        if waited > WAIT_LIMIT or DONE:
            if DONE:
                break
            else: # look for the last prep file and predict
                curr_file = signals_input + "_" + str(PREP_I)
                if curr_file.split("/")[-1] in files:
                    with open(curr_file, 'rb') as signal_in:
                        while True:
                            if counter == LIMIT:
                                if IDs:
                                    predictions = model.predict(np.array(signals), batch_size=len(signals))
                                    f_out.write(
                                        "\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                   enumerate(predictions)]) + "\n")
                                    IDs, signals = [], []
                                print('All signals have been processed', counter)
                                DONE = True
                                break
                            try:
                                key, signal = read_signal(signal_in)
                            except Exception as ex:
                                print(ex)
                                if IDs:
                                    predictions = model.predict(np.array(signals), batch_size=len(signals))
                                    f_out.write(
                                        "\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                   enumerate(predictions)]) + "\n")
                                    IDs, signals = [], []
                                print('All signals have been processed', counter)
                                DONE=True
                                break
                            else:
                                IDs.append(key)
                                signals.append(signal)
                                counter += 1
                                if counter % BATCH_SIZE == 0:
                                    if len(signals) > 0:
                                        print("loading done in", time.time() - start_loading)
                                        prediction_start = time.time()
                                        predictions = model.predict(np.array(signals), batch_size=len(signals))
                                        print("predicting done in", time.time() - prediction_start)
                                        writing_start = time.time()
                                        f_out.write("\n".join(
                                            [f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                             enumerate(predictions)]) + "\n")
                                        print("writing done in", time.time() - writing_start)
                                        IDs, signals = [], []
                                        start_loading = time.time()
                                    else:
                                        continue
                    os.remove(curr_file)
                    DONE = True
                    break

        next_file = signals_input + "_" + str(PREP_I + 1)
        if next_file.split("/")[-1] in files:
            curr_file = signals_input + "_" + str(PREP_I)
            PREP_I+=1
            print("FOUND! waited:", waited, len(files), files)
            with open(curr_file, 'rb') as signal_in:
                    while True:
                            if counter == LIMIT:
                                if IDs:
                                    predictions = model.predict(np.array(signals),batch_size=len(signals))
                                    f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                           enumerate(predictions)])+"\n")
                                    IDs, signals = [],[]
                                print('All signals have been processed', counter)
                                DONE = True
                                break

                            try:
                                key, signal = read_signal(signal_in)
                            except Exception as ex:

                                print(ex)

                                if IDs:
                                    #signal = np.array(signals)
                                    predictions = model.predict(np.array(signals),batch_size=len(signals))
                                    f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                           enumerate(predictions)])+"\n")
                                    IDs, signals = [],[]
                                print('All signals have been processed', counter)
                                break

                            else:
                                IDs.append(key)
                                signals.append(signal)
                                counter+=1
                                if counter % BATCH_SIZE == 0:
                                    if len(signals) > 0:
                                        print("loading done in", time.time() - start_loading)
                                        prediction_start = time.time()
                                        predictions = model.predict(np.array(signals),batch_size=len(signals))
                                        print("predicting done in", time.time() - prediction_start)
                                        writing_start = time.time()
                                        f_out.write("\n".join([f"{IDs[index]}\t{prediction[1]}\t{label}" for index, prediction in
                                                               enumerate(predictions)])+"\n")
                                        print("writing done in", time.time() - writing_start)
                                        IDs, signals = [],[]
                                        start_loading = time.time()
                                    else:
                                        continue
            os.remove(curr_file)
            waited = 0
        else:

            done_file = signals_input + "_DONE"
            if done_file.split("/")[-1] in files:
                waited = WAIT_LIMIT
                os.remove(done_file)
            time.sleep(WAIT_INTERVAL)
            waited+=WAIT_INTERVAL

