#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:41:22 2023

@author: admin
"""

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from tensorflow.keras.layers import Conv1D, BatchNormalization, ReLU, Dropout, Add, Input, GlobalAveragePooling1D, Dense
def site_level_predictions(inputs):
        x=inputs
# Define the model
    
        x=Conv1D(32, kernel_size=3, activation='relu')(x)
        x=MaxPooling1D(pool_size=2)(x)
        x=Conv1D(64, kernel_size=3, activation='relu')(x)
        x=MaxPooling1D(pool_size=2)(x)
        x=Conv1D(128, kernel_size=3, activation='relu')(x)
        x=MaxPooling1D(pool_size=2)(x)
        x=Flatten()(x)
        x=Dense(128, activation='relu')(x)
        x=Dense(1, activation='sigmoid')(x)
        return x
    


def jasper_sub_block(x, filters, kernel_size, dilation_rate, dropout_rate):
    residual = x
    x = Conv1D(filters=filters, kernel_size=kernel_size, dilation_rate=dilation_rate, padding='same')(x)
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = Dropout(rate=dropout_rate)(x)
    
    x = Conv1D(filters=filters, kernel_size=kernel_size, dilation_rate=dilation_rate, padding='same')(x)
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = Dropout(rate=dropout_rate)(x)
    
    # Residual connection
    residual = Conv1D(filters=filters, kernel_size=1, padding='same')(residual)
    residual = BatchNormalization()(residual)
    
    x = Add()([residual, x])
    return x

def jasper_block(x, num_sub_blocks, filters, kernel_size, dilation_rate, dropout_rate):
    for _ in range(num_sub_blocks):
        x = jasper_sub_block(x, filters, kernel_size, dilation_rate, dropout_rate)
    return x

def build_jasper_model(inputs):
    x=inputs
    
    # First layers for initial processing
    x = Conv1D(filters=32, kernel_size=11, strides=2, padding='same')(x)
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = Dropout(0.2)(x)
    
    x = jasper_block(x, num_sub_blocks=1, filters=32, kernel_size=11, dilation_rate=1, dropout_rate=0.2)
    x = jasper_block(x, num_sub_blocks=1, filters=64, kernel_size=11, dilation_rate=2, dropout_rate=0.2)
    
    x = GlobalAveragePooling1D()(x)
    x = Dense(1, activation='sigmoid')(x)  # Adjust output neurons for your binary classification
    
  
    return x


