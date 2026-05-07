import os
import pandas
import argparse
import lr_types
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from tensorflow.keras.layers import Input
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import Callback
from tensorflow.keras.callbacks import CSVLogger
from scipy.stats import spearmanr,pearsonr,kstest
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from collections import Counter


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='trainer', description=
    """ 
    """, usage='python trainer.py -i <input_path/prefix> -o <out_dir>')

    OPTIONAL = parser._action_groups.pop()

    REQUIRED = parser.add_argument_group('required arguments')

    REQUIRED.add_argument("-i", "--input",
                          help="Prefix for assembled data",
                          metavar='\b',
                          required=True)

    REQUIRED.add_argument("-o", "--out_path",
                          help="Output path",
                          metavar='\b',
                          required=True)

    OPTIONAL.add_argument("--vector_len",
                          help="Length of the vectors, default is 36",
                          default=36,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--features",
                          help="Number of vectors, default is 7",
                          default=7,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--labels",
                          help="Number of labels, default is 1",
                          default=1,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--fine_tune",
                          help="Start from pretrained model1 path",
                          default=None,
                          required=False)

    OPTIONAL.add_argument("--arch",
                          help="Mini/Mid/Large network",
                          default="Mini",
                          required=False)

    OPTIONAL.add_argument("--weights",
                          help="Type of weights. 0 for positive/negative balance;"
                               "\n1 for subsample balance (validated/synthetic, positive/negative, single/multiple modified 9mer;"
                               "\n2 for subsample balance (validated/synthetic, positive/negative) ",
                          default="0",
                          required=False)


    parser._action_groups.append(OPTIONAL)

    ARGS = parser.parse_args()

    # required arg
    INPUT_PATH = ARGS.input
    PATH_OUT = ARGS.out_path
    # optional args
    VLEN = ARGS.vector_len
    FEATURES = ARGS.features
    LABELS = ARGS.labels
    
    path_signals = INPUT_PATH + '_final_training_data.p'
    path_labels = INPUT_PATH + '_final_training_labels.p'
    path_weights = INPUT_PATH + f"_final_training_weights{ARGS.weights}.p"

    vpath_signals = INPUT_PATH + '_final_validation_data.p'
    vpath_labels = INPUT_PATH + '_final_validation_labels.p'
    vpath_weights = INPUT_PATH + f"_final_validation_weights{ARGS.weights}.p"

    signals = np.load(path_signals, allow_pickle=True)
    print('signal loading done')
    labels = np.load(path_labels, allow_pickle=True)
    vsignals = np.load(vpath_signals, allow_pickle=True)
    print('vsignal loading done')
    vlabels = np.load(vpath_labels, allow_pickle=True)

    weights = np.load(path_weights, allow_pickle=True)
    vweights = np.load(vpath_weights, allow_pickle=True)

    signals, labels, weights = shuffle(signals, labels,weights)
    print('shuffle done')

    vsignals, vlabels, vweights = shuffle(vsignals, vlabels,vweights)

    labels_one_hot = tf.keras.utils.to_categorical(labels, 2)
    vlabels_one_hot = tf.keras.utils.to_categorical(vlabels, 2)
    signals = np.array(signals, dtype='float32')
    vsignals = np.array(vsignals, dtype='float32')
    
    print(Counter([str(x) for x in labels_one_hot]))
    print(Counter([str(x) for x in vlabels_one_hot]))

    cfg = {
        'input_size': ( int(VLEN * 5), FEATURES),
        'output_size': 1,
        'max_epoch': 100,
        'batch_size': 1024,
        'lr': 0,  # 0: cyclic, 1: cosine
        'lr_warmup': 2,
        'patience': 50,
    }

    strategy = tf.distribute.MirroredStrategy()
    print("Number of devices: {}".format(strategy.num_replicas_in_sync))
    with strategy.scope():
        inputs = Input(shape=cfg['input_size'])

        if ARGS.arch == "Large":
            from network_2132024 import readwise_prediction_network
            outputs = readwise_prediction_network(inputs, cfg['output_size'])
        
        elif ARGS.arch == "Mid":
            from DL_models import build_Jasper
            outputs = build_Jasper(inputs, Deep=True)
        
        elif ARGS.arch == "Mini":
            from network_21122023 import build_jasper_model
            outputs = build_jasper_model(inputs)
        else:
            raise ("--arch  must be either Mini or Mid or Large")

        model = Model(inputs=inputs, outputs=outputs)

        if ARGS.fine_tune:
            model.load_weights(ARGS.fine_tune) 
            # Freeze batch normalization layers for fine-tuning
            for layer in model.layers:
                if isinstance(layer, tf.keras.layers.BatchNormalization):
                    layer.trainable = False

        model.compile(optimizer='adam',
                      loss='binary_crossentropy',
                      weighted_metrics=['accuracy'])
    model.summary()

    # callbacks

    early_stopping = EarlyStopping(monitor='val_loss', patience=cfg['patience'])
    
    if not os.path.exists(PATH_OUT+"/models/"):
            os.mkdir(PATH_OUT+"/models/")
    

        
    model_path = PATH_OUT + '/models/Model_100_epoch_relu.h5'
    checkpoint = ModelCheckpoint(model_path, save_best_only=True, monitor='val_loss', mode='min')
    csv_logger = CSVLogger(PATH_OUT + "/Model_train_val_epoch.log.csv")

    if cfg['lr'] == 0:

        batches_per_epoch = signals.shape[0] // cfg['batch_size']  # important

        if ARGS.fine_tune:
            lr_callback = lr_types.CyclicLR(base_lr=0.0000001,
                                        max_lr=0.00001,
                                        step_size=batches_per_epoch,
                                        mode='exp_range',
                                        gamma=0.999995)
        else:
            lr_callback = lr_types.CyclicLR(base_lr=0.00001,
                                        max_lr=0.001,
                                        step_size=batches_per_epoch,
                                        mode='exp_range',
                                        gamma=0.999995)
    else:

        batches_per_epoch = signals.shape[0] // cfg['batch_size']  # important
        total_steps = cfg['max_epoch'] * batches_per_epoch
        warmup_steps = cfg['lr_warmup'] * batches_per_epoch
        lr_callback = lr_types.WarmUpCosineDecayScheduler(learning_rate_base=0.001,
                                                          learning_rate_min=0.00001,
                                                          total_steps=total_steps,
                                                          warmup_learning_rate=0.0,
                                                          warmup_steps=warmup_steps,
                                                          hold_base_rate_steps=0)

    
    
    callback_lst=[early_stopping, checkpoint, lr_callback, csv_logger]
    
    #fit the model
    history = model.fit(signals, labels,
                        batch_size=cfg['batch_size'], epochs=cfg['max_epoch'],sample_weight=weights,
                        validation_data=(vsignals, vlabels,vweights), validation_batch_size=1024,
                        callbacks= callback_lst)
    
    print(history.history.keys())

    # summarize history for loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')

    plt.savefig(PATH_OUT + "/Model_100_epoch_relu.png")
    plt.show()

