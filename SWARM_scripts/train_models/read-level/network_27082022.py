import sys
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, DepthwiseConv2D, Add, Dense, Conv2D, Activation, BatchNormalization, MaxPooling2D, Dropout
from tensorflow.keras.regularizers import l2

#afunc = 'swish'
afunc = 'relu'
def conv_layer(x, f, k, s, p):
    return Conv2D(filters=f, kernel_size=k, strides=s, padding=p, kernel_regularizer=l2(5e-4))(x)

def dconv_layer(x, k, s, m, p):
    return DepthwiseConv2D(kernel_size=k, strides=s, depth_multiplier=m, padding=p, kernel_regularizer=l2(5e-4))(x)

def identity_block(x, m):
   
    c = int(x.shape[-1])
    x_shortcut = x
   
    x = conv_layer(x, f=c*m, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=(1,3), s=(1,1), m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)

    x = Add()([x, x_shortcut])
    x = Activation(afunc)(x)

    return x

def residual_conv_block(x, m, s):
   
    c = int(x.shape[-1])
    x_shortcut = x

    x = conv_layer(x, f=c*m, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=(1,3), s=(s,s), m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c*m, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)

    x_shortcut = conv_layer(x_shortcut, f=c*m, k=(1,3), s=(s,s), p='same')
    x_shortcut = BatchNormalization(axis=-1)(x_shortcut)

    x = Add()([x, x_shortcut])
    x = Activation(afunc)(x)

    return x

def combine_conv_block(x, m):

    c = int(x.shape[-1])
   
    x = conv_layer(x, f=c, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = dconv_layer(x, k=(1,3), s=(1,1), m=2, p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)

    x = conv_layer(x, f=c*m, k=(1,1), s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)
    
    return x
    
def feature_extractor(x, f, k):
    
    x = conv_layer(x, f=f, k=k, s=(1,1), p='same')
    x = BatchNormalization(axis=-1)(x)
    x = Activation(afunc)(x)
    x = MaxPooling2D(pool_size=(1, 3), strides=(1, 2), padding='same')(x)
    
    x = residual_conv_block(x, m=2, s=2)
    x = identity_block(x, m=2)
    x = identity_block(x, m=2)

    x = residual_conv_block(x, m=2, s=2)
    x = identity_block(x, m=2)
    x = identity_block(x, m=2)
    
    return x
    
def readwise_prediction_backbone(x):
    
    x3 = feature_extractor(x, f=48, k=(1,3))
    x7 = feature_extractor(x, f=48, k=(1,7))
    x11 = feature_extractor(x, f=48, k=(1,11))
    
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
    x = Dense(output_size, activation='softmax')(x)
    
    return x

def readwise_prediction_network(inputs, output_size):

    x = readwise_prediction_backbone(inputs)
    x = readwise_prediction_head(x, output_size)
    
    return x
    
if __name__ == '__main__':

    input_size = (1,100,2)
    output_size = 3
    
    inputs = Input(shape=input_size)
    outputs = readwise_prediction_network(inputs, output_size)

    model = Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])    

    model.summary()
    sys.stdout.flush()
    