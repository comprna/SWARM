#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 10:53:46 2024

@author: admin
"""
import sys
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, DepthwiseConv1D, Add, Dense, Activation, BatchNormalization, MaxPooling1D, Dropout
from tensorflow.keras.regularizers import l2

#afunc = 'swish'
afunc = 'relu'
def conv_layer(x, f, k, s, p):
    return Conv1D(filters=f, kernel_size=k, strides=s, padding=p, kernel_regularizer=l2(5e-4))(x)

def dconv_layer(x, k, s, m, p):
    return DepthwiseConv1D(kernel_size=k, strides=s, depth_multiplier=m, padding=p, kernel_regularizer=l2(5e-4))(x)

def identity_block(x, m):
   
    c = int(x.shape[-1])
    x_shortcut = x
   
    x = conv_layer(x, f=c*m, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=3, s=1, m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)

    x = Add()([x, x_shortcut])
    x = Activation(afunc)(x)

    return x

def residual_conv_block(x, m, s):
   
    c = int(x.shape[-1])
    x_shortcut = x

    x = conv_layer(x, f=c*m, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=3, s=s, m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c*m, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)

    x_shortcut = conv_layer(x_shortcut, f=c*m, k=3, s=s, p='same')
    x_shortcut = BatchNormalization(axis=-1)(x_shortcut)

    x = Add()([x, x_shortcut])
    x = Activation(afunc)(x)

    return x

def combine_conv_block(x, m):

    c = int(x.shape[-1])
   
    x = conv_layer(x, f=c, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=3, s=1, m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c*m, k=1, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)
    
    return x
    
def feature_extractor(x, f,k):
    
    x = conv_layer(x, f=f, k=k, s=1, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)
    x = MaxPooling1D(pool_size=3, strides=2, padding='same')(x)
    
    x = residual_conv_block(x, m=2, s=2)
    x = identity_block(x, m=2)
    x = identity_block(x, m=2)

    x = residual_conv_block(x, m=2, s=2)
    x = identity_block(x, m=2)
    x = identity_block(x, m=2)
    
    return x
    
def readwise_prediction_backbone(x):
    
    x3 = feature_extractor(x, f=48, k=3)
    x7 = feature_extractor(x, f=48, k=7)
    x11 = feature_extractor(x, f=48, k=11)
    
    x = tf.keras.layers.Concatenate()([x3, x7, x11])
    x = combine_conv_block(x, 0.5)

    x = residual_conv_block(x, m=1, s=2)
    x = identity_block(x, m=1)
    x = identity_block(x, m=1)

    x = residual_conv_block(x, m=1, s=2)
    x = identity_block(x, m=1)
    x = identity_block(x, m=1)

    x = residual_conv_block(x, m=1, s=2)
    x = identity_block(x, m=1)
    x = identity_block(x, m=1)
    
    return x

def readwise_prediction_head(x, output_size):

    x = tf.keras.layers.Flatten()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='swish')(x)
    x = Dropout(0.4)(x)
    x = Dense(512, activation='swish')(x)
    x = Dropout(0.2)(x)
    x = Dense(output_size, activation='sigmoid')(x)
    
    return x

def readwise_prediction_network(inputs, output_size):

    x = readwise_prediction_backbone(inputs)
    x = readwise_prediction_head(x, output_size)
    
    return x
    
if __name__ == '__main__':

    input_size = (100, 2)
    output_size = 2
    
    inputs = Input(shape=input_size)
    outputs = readwise_prediction_network(inputs, output_size)

    model = Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])    

    model.summary()
    sys.stdout.flush()
