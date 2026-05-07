import argparse

parser = argparse.ArgumentParser(prog='CHEUI_predict_model2 v0.1', description=
""" 
This script takes a per-read prediction file generated using CHEUI_predict_mode1.py \
and generate RNA modification predictions per-site

""", usage='python CHEUI_predict_model2.py -i <path_to_predictions_model_1> ' \
           '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')




REQUIRED.add_argument("-NM", "--NonModified",
                      help="path to non-modified IVT predictions from model1",
                      metavar='\b',
                      required=True)
REQUIRED.add_argument("-M", "--Modified",
                      help="path to modified IVT predictions from model1",
                      metavar='\b',
                      required=True)
REQUIRED.add_argument('--center',
                      default="A",
                      required=True)
REQUIRED.add_argument("-o", "--dir_out",
                      help="Path to the output directory",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument("-L", "--Loops",
                      help="number of loops to perform over each kmer",
                      metavar='\b',
                      default=40)
                      
OPTIONAL.add_argument("--kmer",
                      help="Unique to use only Unique centered kmers, or All to use all 9mers",
                      metavar='\b',
                      default="All",
                      required=False)

OPTIONAL.add_argument("--features",
                      help="Features used for model 2, can be Standard, Mismatch or Bases",
                      default="Standard",
                      required=False)

OPTIONAL.add_argument("--percentile",
                      help="Path to model1 percentiles",
                      default=None,
                      required=False)

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
Input_nonmod = ARGS.NonModified
Input_mod = ARGS.Modified
CenterNuc = ARGS.center
OutDir = ARGS.dir_out
Loops = ARGS.Loops
Unique = ARGS.kmer
FEATURES =  ARGS.features

nbCPU=48



import sys
import pickle
import random
import pandas as pd
from math import floor
from random import choices
from random import uniform
import multiprocessing as mp
from random import choices
import numpy as np



def load_predictions(path):
    '''
    load the predictions and probabilities
    '''
    df = pd.read_csv(path, sep='\t', names=['site', 'predictions', 'label'])
    split_name = ['_'.join(i.split('_')[-1]) for i in df['site']]
    df['base'] = split_name
    split_name = ['_'.join(i.split('_')[-2]) for i in df['site']]
    df['Qscore'] = split_name
    split_name = ['_'.join(i.split('_')[:-3]) for i in df['site']]
    df['name'] = split_name
    split_kmer = [(i.split('_')[-1]) for i in df['name']]
    df['kmer'] = split_kmer
    return df

def base_to_index(base):
    index=[4]
    if base == 'A':
        index = [0]    
    if base == 'C':
        index = [1]        
    if base == 'G':
        index = [2]        
    if base == 'T':
        index = [3]
    return(index)        


def convert_p_to_vector_Mismatch(pairs):    
    histograms = [0,0]*100
    for m1p, missmatch in pairs:
        if m1p == 1:
            histograms[99]+=1
            histograms[-1]+=int(missmatch != CenterNuc)
        else:
            try:
                index = floor(m1p*100)   # get the bin index for this probability
            except Exception as e:
                print(e,"p=",m1p)
            else:
                histograms[index] +=1
                histograms[100+index] +=int(missmatch != CenterNuc)
    a= np.array(histograms)
    return np.swapaxes(a.reshape(2,100),0,1)    

def convert_p_to_vector_Bases(pairs):
    prob_dist = [0]*100
    basecount = [0]*5
    for m1, nuc in pairs:
        if m1 == 1:
            prob_dist[-1]+=1
            basecount[base_to_index(nuc)[0]]+=1
        else:
            try:
                index = floor(m1*100)   # get the bin index for this probability
            except Exception as e:
                print(e,"p=",m1)
            else:
                prob_dist[index] +=1
                basecount[base_to_index(nuc)[0]]+=1
    for x in basecount:
        prob_dist.append(x)
    return (prob_dist)

def convert_p_to_vector_Standard(pairs):
    prob_dist = [0]*100
    for m1, nuc in pairs:
        if m1 == 1:
            prob_dist[-1]+=1
        else:
            try:
                index = floor(m1*100)   # get the bin index for this probability
            except Exception as e:
                print(e,"p=",m1)
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

    def __call__(self, pairs):
        probs_valid = [x[0] for x in pairs if self.is_valid_m1(x[0])]
        counts, bin_ranges = np.histogram(probs_valid, bins = self.percentiles)
        return counts

    def is_valid_m1(self,value):
        try:
            if 1.0 >= value >= 0.0:
                return True
            return False
        except:
            return False

if FEATURES == "Standard":
    convert_p_to_vector = convert_p_to_vector_Standard
elif FEATURES == "Mismatch":
    convert_p_to_vector = convert_p_to_vector_Mismatch
elif FEATURES == "Bases":
    convert_p_to_vector = convert_p_to_vector_Bases
elif FEATURES == "Percentile":
    if not ARGS.percentile:
        raise("--percentile must be provided")
    convert_p_to_vector = convert_p_to_vector_Percentile(ARGS.percentile)
else:
    raise ("--features must be either Standard or Mismatch or Bases")



no_modify_p = Input_nonmod
modify_p = Input_mod
out_dir = OutDir

modify_df = load_predictions(modify_p)
no_modify_df = load_predictions(no_modify_p)

unique_ids_mod = modify_df['kmer'].unique().tolist()
unique_ids_nomod = no_modify_df['kmer'].unique().tolist()
unique_ids = [value for value in unique_ids_mod if value in unique_ids_nomod]

unique_ids_keep = []
for ids in unique_ids:
    if ids.count(CenterNuc) == 1:
        unique_ids_keep.append(ids)

# prepare the dataset
X_test = []
y_test = []

# create probabilities for N of reads and stoichiometries
weights = []
N_reads_stoichiometry = []

for i in range(90):
    N_reads_stoichiometry.append(i + 10)

# create probabilities for N of reads and stoichiometries
N_reads_cov = []

for i in range(990):
    N_reads_cov.append(i + 10)

weights.reverse()

if Unique == "Unique":
    ids_to_use = unique_ids_keep
else:
    ids_to_use = unique_ids

rounds = int(Loops)
TARGET= (rounds*len(ids_to_use)) //2
counterWT = 0
MAXITER = rounds
iters = 0

def DataWTSimple(id_):
    successattempt=0
    trials=0
    while successattempt < 1:
        trials +=1
        if trials > 50:
            print ("giving up on ", id_) 
            mix=[uniform(0,0.1) for x in range(92)]+ [uniform(0.2,0.8) for x in range(8)] 
            match= [CenterNuc]*100  
            pairs=[[x,y] for x,y in zip(mix,match)]        
            break
        n_reads = choices(N_reads_cov)[0]
        stoichiometry = choices(N_reads_stoichiometry)[0]
        Noise = choices(range(8))[0]
        modify_site_df = modify_df[modify_df['kmer'] == id_]
        non_modify_site_df = no_modify_df[no_modify_df['kmer'] == id_]
        number_read_noise = int((Noise / 100) * n_reads)
        number_read_unmodify = n_reads - number_read_noise
        try:
            modify_site_df_noise = modify_site_df.sample(n=number_read_noise)
            non_modify_site_df_random = non_modify_site_df.sample(n=number_read_unmodify)
            mix = modify_site_df_noise['predictions'].tolist() + non_modify_site_df_random['predictions'].tolist()
            match = modify_site_df_noise['base'].tolist() + non_modify_site_df_random['base'].tolist()
            pairs=[[x,y] for x,y in zip(mix,match)]
            successattempt = 1
        except:
            continue
    return(convert_p_to_vector(pairs))    

def DataWTAdded(id_):
    successattempt=0
    trials=0
    while successattempt < 1:
        trials +=1
        if trials > 50:
            print ("giving up on ", id_) 
            mix=[uniform(0,0.1) for x in range(92)]+ [uniform(0.2,0.8) for x in range(8)]                        
            match= [CenterNuc]*100
            pairs=[[x,y] for x,y in zip(mix,match)]
            break
        n_reads = choices(N_reads_cov)[0]
        stoichiometry = choices(N_reads_stoichiometry)[0]
        Noise = choices(range(8))[0]
        modify_site_df = modify_df[modify_df['kmer'] == id_]
        non_modify_site_df = no_modify_df[no_modify_df['kmer'] == id_]
        number_read_noise = int((Noise / 200) * n_reads)
        number_read_unmodify = n_reads - 2*number_read_noise
        try:
            modify_site_df_noise = modify_site_df.sample(n=number_read_noise)
            non_modify_site_df_random = non_modify_site_df.sample(n=number_read_unmodify)
            mix = modify_site_df_noise['predictions'].tolist() + non_modify_site_df_random['predictions'].tolist() + [uniform(0.2,0.8) for x in range(number_read_noise)]
            match = modify_site_df_noise['base'].tolist() + non_modify_site_df_random['base'].tolist() + choices(["A","C","G","T"], k=number_read_noise)
            pairs=[[x,y] for x,y in zip(mix,match)]
            successattempt = 1
        except:
            continue
    return(convert_p_to_vector(pairs))    
    
    
def DataKO(id_):
    successattempt=0
    trials=0
    while successattempt < 1:
        trials +=1
        if trials > 50:
            print ("giving up on ", id_) 
            mix=[uniform(0,0.1) for x in range(50)]+ [uniform(0.5,1) for x in range(50)]  
            match= [CenterNuc]*100
            pairs=[[x,y] for x,y in zip(mix,match)]
            break
        n_reads = choices(N_reads_cov)[0]
        stoichiometry = choices(N_reads_stoichiometry)[0]
        Noise = choices(range(8))[0]
        modify_site_df = modify_df[modify_df['kmer'] == id_]
        non_modify_site_df = no_modify_df[no_modify_df['kmer'] == id_]
        number_read_modify = int((stoichiometry / 100) * n_reads)
        number_read_unmodify = n_reads - number_read_modify
        try:
            modify_site_df_random = modify_site_df.sample(n=number_read_modify)
            non_modify_site_df_random = non_modify_site_df.sample(n=number_read_unmodify)
            mix = modify_site_df_random['predictions'].tolist() + non_modify_site_df_random['predictions'].tolist()
            match = modify_site_df_random['base'].tolist() + non_modify_site_df_random['base'].tolist()
            pairs = [[x,y] for x,y in zip(mix,match)]
            successattempt = 1
        except Exception as e:
            print("KO attempt failed",e,sep="\n")
            continue
    return(convert_p_to_vector(pairs))    


        

# do KO
for i in range(rounds):
    pool=mp.Pool(nbCPU)
    results=pool.map(DataKO,[id_ for id_ in ids_to_use])
    pool.close()
    pool.join()
    with open(out_dir + 'XKO_train.p', "ab") as pickle_file_ko:
        pickle.dump(results, pickle_file_ko)

print ("Done with KO")

# Do WT Added
for i in range(rounds):
    pool=mp.Pool(nbCPU)
    results=pool.map(DataWTAdded,[id_ for id_ in ids_to_use])
    pool.close()
    pool.join()
    with open(out_dir + 'XWT_added_train.p', "ab") as pickle_file_wt:
        pickle.dump(results, pickle_file_wt)

print ("Done with WT Added")

# Do WT Simple
for i in range(rounds):
    pool=mp.Pool(nbCPU)
    results=pool.map(DataWTSimple,[id_ for id_ in ids_to_use])
    pool.close()
    pool.join()
    with open(out_dir + 'XWT_simple_train.p', "ab") as pickle_file_wt:
        pickle.dump(results, pickle_file_wt)
    
print ("Done with WT Simple")
    

