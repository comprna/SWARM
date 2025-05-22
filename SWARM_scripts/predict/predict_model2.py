#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(prog='predict_model2 v0.1', usage='python predict_model2.py -i <path_to_predictions_model_1> ' \
           '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input",
                      help="path to read-level prediction file from CHEUI_predict_model1.py",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to pretrainned DL model 2",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-c", "--cutoff",
                      help="model 2 probability cutoff for printing sites",
                      metavar='\b',
                      default='-1'
                      )

REQUIRED.add_argument("-d", "--double_cutoff",
                      help="Model 1 probability cutoffs used to calculate the stoichiometry",
                      metavar='\b',
                      default='0.1,0.5',
                      )

REQUIRED.add_argument("-n", "--min_reads",
                      help="Minimun number of reads in a site to include in the analysis, ",
                      metavar='\b',
                      default=20
                      )

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument('--arch',
                      default="Mini",
                      required=False)

REQUIRED.add_argument('--percentile',
                      default=None,
                      required=False)

OPTIONAL.add_argument('-v', '--version',
                      action='version',
                      version='%(prog)s')

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
input_df = ARGS.input
DL_model = ARGS.DL_model
cutoff = float(ARGS.cutoff)
double_cutoff = ARGS.double_cutoff
lower_cutoff = float(double_cutoff.split(',')[0])
upper_cutoff = float(double_cutoff.split(',')[1])
file_out_path = ARGS.file_out
min_reads = int(ARGS.min_reads)

import os
import sys
import time
import random
import numpy as np
from math import floor
import tensorflow as tf
import _pickle as cPickle
from tensorflow.keras import Input
from tensorflow.keras.models import Model

random.seed(42)
BOOTSTRAP_NUM = 21
BATCH_SIZE = 4096


label_dct = {
    "1":"pU_RNA002",
    "2":"pU_RNA004",
    "3":"m5C_RNA002" ,
    "4":"m5C_RNA004",
    "5" :"m6A_RNA002",
    "6" :"m6A_RNA004"
}

MODELS_PATH = "../SWARM_models/"

def convert_p_to_vector_SWARM(probs):

    prob_dist = [0]*100
    for probability_m1 in probs:
        if probability_m1 == 1:
            prob_dist[-1]+=1
        else:
            try:
                index = floor(probability_m1*100)   # get the bin index for this probability
            except Exception as e:
                print(e,"p=",probability_m1)
            else:
                prob_dist[index] +=1

    return (prob_dist)

class convert_p_to_vector_Percentile():
    def __init__(self, path):
        self.path = path
        self.percentiles = [0]
        with open(self.path) as f:
            for line in f:
                self.percentiles.append(float(line.strip()))
        self.percentiles.append(1)

    def __call__(self, probs):
        probs_valid = [x for x in probs if self.is_valid_m1(x)]
        counts, bin_ranges = np.histogram(probs_valid, bins = self.percentiles)
        return counts

    def is_valid_m1(self,value):
        try:
            if 1.0 >= value >= 0.0:
                return True
            return False
        except:
            return False


def convert_p_to_vector_CHEUI(probs):
    prob_dist = [0]*99
    for probability_m1 in probs:
        if probability_m1 < 0.01:
            continue
        elif probability_m1 == 1:
            prob_dist[-1] += 1
        else:
            index = floor(probability_m1*100) -1   # get the bin index for this probability
            prob_dist[index] +=1
    return (prob_dist)

def predict_SWARM(probs):
    prob_vectors = np.array(probs).reshape(len(probs), 100, 1)
    # prob_vectors = np.expand_dims(prob_vectors, axis=0)
    return model.predict(prob_vectors,batch_size = len(prob_vectors))

def predict_CHEUI(probs):

    return model.predict(probs)

def biggerThan100(prob_list):

    return [convert_p_to_vector(random.sample(prob_list, 100)) for _ in range(BOOTSTRAP_NUM)]

def predict_average_SWARM(keys_100, lr_probs):

    for index, key in enumerate(keys_100):
        p_mod = np.mean([1-x[0] for x in lr_probs[index * BOOTSTRAP_NUM: index*BOOTSTRAP_NUM  +BOOTSTRAP_NUM]])
        if p_mod > cutoff and keys_100[index].split("\t")[-1] != "None":
            print(keys_100[index] + '\t' + str(p_mod), file=file_out)

def predict_average_CHEUI(keys_100, lr_probs):

    for index, key in enumerate(keys_100):
        p_mod = np.mean(lr_probs[index * BOOTSTRAP_NUM: index*BOOTSTRAP_NUM  +BOOTSTRAP_NUM])
        if p_mod > cutoff and keys_100[index].split("\t")[-1] != "None":
            print(keys_100[index] + '\t' + str(p_mod), file=file_out)

def get_pMod_SWARM(pmod):
    return 1-pmod[0]

def get_pMod_CHEUI(pmod):
    return pmod[0]

inputs = Input(shape=(100, 1))



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

if ARGS.percentile:
    convert_p_to_vector = convert_p_to_vector_Percentile(ARGS.percentile)
else:
    convert_p_to_vector = convert_p_to_vector_SWARM
predict_function = predict_SWARM
predict_average = predict_average_SWARM
get_pMod = get_pMod_SWARM
predictions_site = []
counter = 0
counter_predictions = 0

if os.path.isfile(file_out_path):
    print('WARNING: site level prediction file already exists, please delete it or change the output name')
    sys.exit()

# code the last one
with open(file_out_path, 'w') as file_out:
    print('contig' + '\t' + 'position' + '\t' + 'site' + '\t' + 'coverage' + \
          '\t' + 'stoichiometry' + '\t' + 'probability', file=file_out)
    prev_ID = ""
    predictions_site, probs_lst, probs_100_lst, keys, keys_100 = [],[],[],[],[]
    start_time = time.time()
    model_name = ""
    with open(input_df, 'r') as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            if not model_name:
                model_key = line[-1]
                try:
                    model_name = label_dct[model_key]
                except KeyError:
                    raise ("Rerun model1 using SWARM_read_level.py ; Model info cannot be inferred from the label (3rd column).\nLabel should be integer 1-6\n"
                           "{1:pU_RNA002, 2:pU_RNA004, 3:m5C_RNA002 , 4:m5C_RNA004, 5:m6A_RNA002, 6 :m6A_RNA004}")
                else:
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    RNAmod,KIT = model_name.split("_")
                    model2_path = os.path.join(script_dir,MODELS_PATH + f"Model2/{KIT}/{RNAmod}/Model_100_epoch_relu.h5")
                    model.load_weights(model2_path)
            else:
                if line[-1] != model_key:
                    raise ("Exiting. Different model1 keys in 3rd column. Make sure that same model1 was used on all signals.")

            ID_lst = line[0].split('_')
            last=ID_lst[-1]
            qscore=""
            base=""
            ### Get qscore and base if they are in the key
            if len(last) == 1 and not last.isnumeric():
                qscore=ID_lst.pop()
                base=ID_lst.pop()
            if ID_lst[-1].isnumeric(): #handle within read position if --modsam was used
                current_ID = '_'.join(ID_lst[:-2])
            else:
                current_ID = '_'.join(ID_lst[:-1])

            if current_ID != prev_ID:
                if len(predictions_site) < min_reads:
                    # do not make predictions for sites with low coverage
                    if line[1] == 'Prediction':
                        continue
                else:

                    try:     #compute stoichiometry using double cutoff
                        mod = [i for i in predictions_site if i > upper_cutoff]
                        no_mod = [i for i in predictions_site if i < lower_cutoff]
                        stoichiometry = len(mod) / (len(mod) + len(no_mod))
                    except:
                        print('Warning: Stoichiometry cannot compute, try less stringent double cutoff values')
                        stoichiometry = 'None'
                    coverage = len(predictions_site)
                    ID_colums = ['_'.join(prev_ID.split('_')[:-2])] + prev_ID.split('_')[-2:]

                    if coverage > 150 and stoichiometry > 0.1:
                        # handle predictions on larger coverages (bootstrap)

                        probs_100_lst.extend(biggerThan100(predictions_site))
                        # store info for the output file

                        keys_100.append(ID_colums[0] + '\t' + ID_colums[1] + '\t' + ID_colums[2] + '\t' + \
                                        str(coverage) + '\t' + str(stoichiometry))
                        if len(probs_100_lst) >= BATCH_SIZE:
                            # predict on a batch
                            lr_probs = predict_function(probs_100_lst)
                            predict_average(keys_100,lr_probs)
                            keys_100,probs_100_lst = [],[]

                    else:
                        # when coverage is under 150 or stoichiometry is low
                        probs_lst.append(convert_p_to_vector(predictions_site))
                        keys.append(ID_colums[0] + '\t' + ID_colums[1] + '\t' + ID_colums[2] + '\t' + \
                                        str(coverage) + '\t' + str(stoichiometry))

                        if len(keys) == BATCH_SIZE:
                            lr_probs = predict_function(probs_lst)

                            for index,prediction in enumerate(lr_probs):
                                p_mod = get_pMod(prediction)
                                if p_mod > cutoff and keys[index].split("\t")[-1] != "None":
                                    print(keys[index] + '\t' + str(p_mod), file=file_out)
                            probs_lst, keys = [],[]


                predictions_site = [float(line[1])]
                prev_ID = current_ID

            else:
                predictions_site.append(float(line[1]))


    if predictions_site:   # handle the last site
        try:  # compute stoichiometry using double cutoff
            mod = [i for i in predictions_site if i > upper_cutoff]
            no_mod = [i for i in predictions_site if i < lower_cutoff]
            stoichiometry = len(mod) / (len(mod) + len(no_mod))
        except:
            print('Warning: Stoichiometry cannot compute, try less stringent double cutoff values')
            stoichiometry = 'None'
        coverage = len(predictions_site)
        ID_colums = ['_'.join(prev_ID.split('_')[:-2])] + prev_ID.split('_')[-2:]

        if coverage > 150 and stoichiometry > 0.1:
            # handle predictions on larger coverages (bootstrap)
            probs_100_lst.extend(biggerThan100(predictions_site))
            # store info for the output file
            keys_100.append(ID_colums[0] + '\t' + ID_colums[1] + '\t' + ID_colums[2] + '\t' + \
                            str(coverage) + '\t' + str(stoichiometry))
            if len(probs_100_lst) >= BATCH_SIZE:
                # predict on a batch
                lr_probs = predict_function(probs_100_lst)
                predict_average(keys_100, lr_probs)
                keys_100, probs_100_lst = [], []

        else:
            # when coverage is under 100 or stoichiometry is low
            probs_lst.append(convert_p_to_vector(predictions_site))
            keys.append(ID_colums[0] + '\t' + ID_colums[1] + '\t' + ID_colums[2] + '\t' + \
                        str(coverage) + '\t' + str(stoichiometry))

            if len(keys) == BATCH_SIZE:
                lr_probs = predict_function(probs_lst)

                for index, prediction in enumerate(lr_probs):
                    p_mod = get_pMod(prediction)
                    if p_mod > cutoff and keys[index].split("\t")[-1] != "None":
                        print(keys[index] + '\t' + str(p_mod), file=file_out)
                probs_lst, keys = [], []


    if len(keys) > 0:
        lr_probs = predict_function(probs_lst)
        for index, prediction in enumerate(lr_probs):
            p_mod = get_pMod(prediction)
            if p_mod > cutoff and keys[index].split("\t")[-1] != "None":
                print(keys[index] + '\t' + str(p_mod), file=file_out)
        probs_lst, keys = [], []

    if len(keys_100) > 0:
        lr_probs = predict_function(probs_100_lst)
        predict_average(keys_100, lr_probs)

    print("All done in",(time.time() - start_time) / 60, "minutes!")

