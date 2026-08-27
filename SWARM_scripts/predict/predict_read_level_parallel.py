import os
import sys
import time
import random
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

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument("--arch",
                          help="Model1 arch from Mini/Mid/Large. Mini is default",
                          default="Mini",
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
                          help="Limit to a given number of predictions",
                          type=int,
                          default=-1,
                          required=False)

OPTIONAL.add_argument("--batch_size",
                          help="Number of predictions loaded at a time on the GPU/CPU",
                          type=int,
                          default=4096,
                          required=False)

OPTIONAL.add_argument("-n","--nworkers",
                          help="Number of preprocessing workers",
                          type=int,
                          default=1,
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


class PreprocessWorker():
    def __init__(self, path ):
        self.path = path
        self.prep_i = 0
        self.current_file = f"{self.path}_{self.prep_i}".split("/")[-1]
        self.current_path = f"{self.path}_{self.prep_i}"
        self.next_file = f"{self.path}_{self.prep_i +1}".split("/")[-1]
        self.done_file = f"{self.path}_DONE".split("/")[-1]
        self.done_path = f"{self.path}_DONE"

    def update_path(self):
        self.prep_i+=1
        self.current_file = f"{self.path}_{self.prep_i}".split("/")[-1]
        self.current_path = f"{self.path}_{self.prep_i}"
        self.next_file = f"{self.path}_{self.prep_i +1}".split("/")[-1]


def read_signal(file):
    char_len = struct.unpack('i', file.read(4))[0]
    contigName = file.read(char_len).decode('utf-8')

    ninemerPos = str(struct.unpack('i', file.read(4))[0])

    ninemer = file.read(9).decode('utf-8')

    char_len = struct.unpack('i', file.read(4))[0]
    readName = file.read(char_len).decode('utf-8')
    readPos = str(struct.unpack('i', file.read(4))[0])
    key = "_".join([str(x) for x in [contigName,ninemerPos,ninemer,readName,readPos]])
    qscore = str(struct.unpack('i', file.read(4))[0])
    base=file.read(1).decode('utf-8')
    key = "_".join([str(x) for x in [contigName,ninemerPos,ninemer,readName,readPos,qscore,base]])
    #print(char_len, name)
    # Read the 2D double array (float values) from the file
    array_data = file.read(4 * 180 * 7)
    return key, np.frombuffer(array_data, dtype=np.float32).reshape(180, 7)




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
LIMIT_REACHED = False

workers_list = []
for workerN in range(ARGS.nworkers):
    workers_list.append(PreprocessWorker(f"{signals_input}_{workerN}"))

with open(file_out, 'a') as f_out:
    while True:
        files = os.listdir(in_dir)
        if waited > WAIT_LIMIT:
            print(f"Unexpected long wait for the next worker preprocessed chunk ({waited} seconds). Terminating all workers.", file=sys.stderr)
            break

        elif DONE:
            print(f"All done")
            break


        found_worker = None
        for worker in workers_list:
            if worker.next_file in files:
                found_worker = worker
                break

        random.shuffle(workers_list)
        if found_worker is not None:

            curr_file = found_worker.current_file
            curr_path = found_worker.current_path
            found_worker.update_path()
            print("FOUND! waited:", waited, len(files), files)
            with open(curr_path, 'rb') as signal_in:
                    while True:
                            if counter == LIMIT:
                                if IDs:
                                    predictions = model.predict(np.array(signals),batch_size=len(signals))
                                    f_out.write("\n".join([f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
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
                                    f_out.write("\n".join([f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
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
                                        f_out.write("\n".join([f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
                                                               enumerate(predictions)])+"\n")
                                        print("writing done in", time.time() - writing_start)
                                        IDs, signals = [],[]
                                        start_loading = time.time()
                                    else:
                                        continue
            os.remove(curr_path)
            waited = 0
        else:
            updated_workers = []
            for worker in workers_list:
                if DONE == True:
                    break
                if worker.done_file in files:
                    waited = 0
                    os.remove(worker.done_path)
                    curr_file = worker.current_file
                    curr_path = worker.current_path
                    if curr_file in files:          #### predict on last file from done worker
                        with open(curr_path, 'rb') as signal_in:
                            while True:
                                if counter == LIMIT:
                                    if IDs:
                                        predictions = model.predict(np.array(signals), batch_size=len(signals))
                                        f_out.write(
                                            "\n".join(
                                                [f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
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
                                            "\n".join(
                                                [f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
                                                 enumerate(predictions)]) + "\n")
                                        IDs, signals = [], []
                                    print('All signals have been processed', counter)
                                    # DONE = True
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
                                                [f"{IDs[index]}\t{prediction[0]}\t{label}" for index, prediction in
                                                 enumerate(predictions)]) + "\n")
                                            print("writing done in", time.time() - writing_start)
                                            IDs, signals = [], []
                                            start_loading = time.time()
                                        else:
                                            continue
                        os.remove(curr_path)
                        # DONE = True
                        # break
                    ### remove worker from list
                    # workers_list.remove(worker)
                else:
                    updated_workers.append(worker)
            if len(updated_workers) == 0:
                DONE = True
            workers_list = updated_workers
            time.sleep(WAIT_INTERVAL)
            waited+=WAIT_INTERVAL




