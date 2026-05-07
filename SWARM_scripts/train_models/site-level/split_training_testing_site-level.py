#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:13:47 2022

@author: admin
"""

import _pickle as cPickle

from sklearn.model_selection import train_test_split
import pickle
import sys
def load_all_data(path):
    '''
    '''
    ### load the stored data 
    counter = 0
    elements = []
    
    with open(path, 'rb') as signal_in:
        while True:
            try:        
                counter +=1
                item = cPickle.load(signal_in)
                elements.extend(item)
                
             # to avoid loading everything predict every 10k singnals
            except:
                break
    return elements

input_file = sys.argv[1]
type_meth=sys.argv[2]
out_dir=sys.argv[3]

label_dct = {"KO": 1, "WT": 0, "WT_added": 0, "WT_simple": 0 }

if type_meth not in label_dct:
    raise (type_meth ," not in conditions, choose from:",label_dct.keys())

type_label= label_dct[type_meth]

data = load_all_data(input_file)
data_labels=len(data)*[type_label]

print(type_meth, type_label, "signals=",len(data))

X_train, X_rem, y_train, y_rem = train_test_split(data,data_labels, train_size=0.6)
del(data,data_labels)

test_size = 0.5
X_valid, X_test, y_valid, y_test = train_test_split(X_rem,y_rem, test_size=0.5)

# training
output1=str(out_dir)+str(type_meth)+'_training_data.p'
output2=str(out_dir)+str(type_meth)+'_training_label.p'
f = open(output1,"wb")
pickle.dump(X_train,f)
f.close()
del(X_train)
f = open(output2,"wb")
pickle.dump(y_train,f)
f.close()
del(y_train)
## validation
output1=str(out_dir)+str(type_meth)+'_validation_data.p'
output2=str(out_dir)+str(type_meth)+'_validation_label.p'
f = open(output1,"wb")
pickle.dump(X_valid,f)
f.close()
del(X_valid)
f = open(output2,"wb")
pickle.dump(y_valid,f)
f.close()
del(y_valid)
## testing
output1=str(out_dir)+str(type_meth)+'_testing_data.p'
output2=str(out_dir)+str(type_meth)+'_testing_label.p'
f = open(output1,"wb")
pickle.dump(X_test,f)
f.close()
del(X_test)
f = open(output2,"wb")
pickle.dump(y_test,f)
f.close()
del(y_test)
