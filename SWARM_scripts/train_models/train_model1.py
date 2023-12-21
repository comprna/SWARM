import argparse
import lr_types
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from tensorflow.keras.layers import Input
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import  CSVLogger
from network_27082022 import readwise_prediction_network
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='trainer', description=
    """ 
    """, usage='python train_model1.py -i <input_path/prefix> -o <out_dir>')

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
                          help="Number of labels, default is 2",
                          default=2,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--epochs",
                          help="Number epochs to train, default is 100",
                          default=100,
                          type=int,
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
    EPOCHS = ARGS.epochs

    path_signals = INPUT_PATH + '_final_training_data.p'
    path_labels = INPUT_PATH + '_final_training_labels.p'
    vpath_signals = INPUT_PATH + '_final_validation_data.p'
    vpath_labels = INPUT_PATH + '_final_validation_labels.p'    
    
    signal = np.load(path_signals,allow_pickle=True)
    print('signal loading done')
    label = np.load(path_labels,allow_pickle=True)
    vsignal = np.load(vpath_signals,allow_pickle=True)
    print( 'vsignal loading done')
    vlabel = np.load(vpath_labels,allow_pickle=True)
    signals, labels = shuffle(signal, label)
    print('shuffle done')
    del(signal,label)
    vsignals, vlabels = shuffle(vsignal, vlabel)  
    del(vsignal,vlabel)
    labels_one_hot = tf.keras.utils.to_categorical(labels, LABELS)
    vlabels_one_hot = tf.keras.utils.to_categorical(vlabels, LABELS)
    signals = np.array(signals, dtype='float32')
    vsignals = np.array(vsignals, dtype='float32')
    signals=np.expand_dims(signals, axis=1)
    vsignals=np.expand_dims(vsignals, axis=1)
    cfg = {
        'input_size' : (1,int(VLEN * 5), FEATURES),
        'output_size' : LABELS,
        'max_epoch' : EPOCHS,
        'batch_size' : 1024,
        'lr' : 0, # 0: cyclic, 1: cosine
        'lr_warmup' : 2,
        'patience' : 50,
    }
    
    strategy = tf.distribute.MirroredStrategy()
    print("Number of devices: {}".format(strategy.num_replicas_in_sync))
    with strategy.scope():
        inputs = Input(shape=cfg['input_size'])
        outputs = readwise_prediction_network(inputs, cfg['output_size'])    
        model = Model(inputs=inputs, outputs=outputs)    
        model.compile(optimizer='adam',
                      loss='categorical_crossentropy',
                      metrics=['accuracy'])    
    model.summary()
        
        
    #callbacks
    model_path = PATH_OUT + f'Model_{EPOCHS}_epoch_relu.h5'
    checkpoint = ModelCheckpoint(model_path, save_best_only=True, monitor='val_loss', mode='min')
    csv_logger = CSVLogger( PATH_OUT + f'Model_{EPOCHS}_epoch_relu.log')
    early_stopping = EarlyStopping(monitor='val_loss', patience=cfg['patience'])
    
    if cfg['lr'] == 0:
        
        batches_per_epoch = signals.shape[0] // cfg['batch_size']   #important
        lr_callback = lr_types.CyclicLR(base_lr=0.00001,
                                            max_lr=0.001,
                                            step_size=batches_per_epoch,
                                            mode='exp_range',
                                            gamma=0.999995)
    else:
        
        batches_per_epoch = signals.shape[0] // cfg['batch_size'] #important
        total_steps = cfg['max_epoch']*batches_per_epoch
        warmup_steps = cfg['lr_warmup']*batches_per_epoch
        lr_callback = lr_types.WarmUpCosineDecayScheduler(learning_rate_base=0.001,
                                                    learning_rate_min=0.00001,
                                                    total_steps=total_steps,
                                                    warmup_learning_rate=0.0,
                                                    warmup_steps=warmup_steps,
                                                    hold_base_rate_steps=0) 

    # fit the model
    history = model.fit(signals, labels_one_hot,
                        batch_size=cfg['batch_size'], epochs=cfg['max_epoch'],
                        validation_data=(vsignals,vlabels_one_hot), validation_batch_size=1024,
                        callbacks=[early_stopping, checkpoint, lr_callback,csv_logger])   
    print(history.history.keys())
    
    # summarize history for loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.show()
    plt.savefig(PATH_OUT + f"Model_{EPOCHS}_epoch_relu.png")
    
