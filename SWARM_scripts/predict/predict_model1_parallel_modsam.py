#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: stefan
"""
import os
import sys
import time
import math
import pysam
import struct
import argparse
import subprocess
import numpy as np
from tensorflow.keras import Input
from tensorflow.keras.models import Model
from network_27082022 import readwise_prediction_network



parser = argparse.ArgumentParser(prog='predict_model1_SWARM_bin_parallel v1.0', description=
""" 
This script takes an ID+signal file generated using SWARM_preprocess_* and predicts RNA modification status \
per read and per 9-mer.

""", usage='python predict_model1_SWARM_bin_parallel -i <path_to_signlas+IDs_file> ' \
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
OPTIONAL.add_argument("--sam",
                          help="Input sam for making modsam",
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

OPTIONAL.add_argument("--arch",
                          help="Model1 arch from Mini/Mid/Large. Mini is default",
                          default="Mini",
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
SAM_PATH = ARGS.sam


def read_signal(file):
    char_len = struct.unpack('i', file.read(4))[0]
    contigName = file.read(char_len).decode('utf-8')

    ninemerPos = str(struct.unpack('i', file.read(4))[0])

    ninemer = file.read(9).decode('utf-8')

    char_len = struct.unpack('i', file.read(4))[0]
    readName = file.read(char_len).decode('utf-8')
    readPos = str(struct.unpack('i', file.read(4))[0])
    #key = "_".join([str(x) for x in [contigName,ninemerPos,ninemer,readName,readPos]])
    qscore = str(struct.unpack('i', file.read(4))[0])
    base=file.read(1).decode('utf-8')
    key = [contigName,ninemerPos,ninemer,readName,readPos,qscore,base]
    #print(char_len, name)
    # Read the 2D double array (float values) from the file
    array_data = file.read(4 * 180 * 7)
    return key, np.frombuffer(array_data, dtype=np.float32).reshape(180, 7)


def get_MM_tag(MMs):
    point = 0
    MM = "N+N"
    for pos in MMs:
        delta = int(pos) - point
        MM+=f",{delta}"
        point+= 1 + delta
    return MM + ";"

def write_outputs(IDs, predictions, f_out, modsam_file, input_read_iterator, input_read, MMs, MLs ):
    input_read_name = input_read.query_name

    for index,prediction in enumerate(predictions):
        ID = IDs[index]
        f_out.write(f"{'_'.join(ID)}\t{prediction[0]}\t{label}\n")
        try:
            contigName,ninemerPos,ninemer,readName,readPos,qscore,base = ID
        except:
            print("Error splitting",ID)
            sys.exit()
        if readName == input_read_name:

            if int(readPos) != 0 or len(MMs) == 0:
                try:
                    col=math.floor(prediction[0]*255)
                except Exception as e:
                    print("skipped Nan probability")
                else:
                    MMs.append(readPos)
                    MLs.append(col)
            # else:
                # print(readPos, prediction[0])
        elif MMs:
            input_read.set_tag("MM",get_MM_tag(MMs),value_type="Z")
            input_read.set_tag("ML", MLs)

            modsam_file.write(input_read)
            try:
                input_read = next(input_read_iterator)
            except StopIteration:
                raise ("Last processed read is not in the input sam. Exiting.")
            input_read_name = input_read.query_name
            # print("CLEARING after writing")
            MMs.clear()
            MLs.clear()
        else:
            try:
                input_read = next(input_read_iterator)
            except StopIteration:
                raise ("Last processed read is not in the input sam. Exiting.")

            input_read_name = input_read.query_name
            MMs.clear()
            MLs.clear()

    return input_read



# load the trainned model
#inputs = Input(shape=(1,int(vector_len)*5, number_of_features))
#output = readwise_prediction_network(inputs, number_of_labels)
#model = Model(inputs=inputs, outputs=output)
#model.load_weights(DL_model)


# load the trainned model
inputs = Input(shape=(int(vector_len)*5, number_of_features))
if ARGS.arch == "Large":
            from network_2132024 import readwise_prediction_network
            outputs = readwise_prediction_network(inputs, 1)

elif ARGS.arch == "Mid":
            from DL_models import build_Jasper
            outputs = build_Jasper(inputs, Deep=True)

elif ARGS.arch == "Mini":
            from network_21122023 import build_jasper_model
            outputs = build_jasper_model(inputs)
else:
            raise ("--arch  must be either Mini or Mid or Large")

model = Model(inputs=inputs, outputs=outputs)

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



with pysam.AlignmentFile(SAM_PATH,"r") as input_samfile:
    with pysam.AlignmentFile(file_out +".mod.sam","wh",header=input_samfile.header) as modsam_file:
        input_read_iterator = iter(input_samfile)
        input_read = next(input_read_iterator)
        MMs,MLs = [],[]
        with open(file_out +".pred.tsv", 'a') as f_out:
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
                                            input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                       input_read_iterator, input_read, MMs, MLs)
                                            IDs, signals = [], []
                                        print('All signals have been processed', counter, f"{len(MMs)} MMs")
                                        DONE = True
                                        break
                                    try:
                                        key, signal = read_signal(signal_in)
                                    except Exception as ex:
                                        print(ex)
                                        if IDs:
                                            predictions = model.predict(np.array(signals), batch_size=len(signals))
                                            input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                       input_read_iterator, input_read, MMs, MLs)
                                            IDs, signals = [], []
                                        print('All signals have been processed', counter, f"{len(MMs)} MMs")
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
                                                input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                           input_read_iterator, input_read, MMs, MLs)
                                                print("writing done in", time.time() - writing_start)
                                                IDs, signals = [], []
                                                start_loading = time.time()
                                            else:
                                                continue
                            os.remove(curr_file)
                            if len(MMs) > 0:
                                # print("MMs", MMs)
                                input_read.set_tag("MM", get_MM_tag(MMs), value_type="Z")
                                input_read.set_tag("ML",MLs)
                                modsam_file.write(input_read)
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
                                            input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                       input_read_iterator, input_read, MMs, MLs)
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
                                            input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                       input_read_iterator, input_read, MMs, MLs)
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
                                                input_read = write_outputs(IDs, predictions, f_out, modsam_file,
                                                                           input_read_iterator, input_read, MMs, MLs)
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




