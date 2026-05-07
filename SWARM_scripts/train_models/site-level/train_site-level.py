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

class LogCustomTestingMetrics(Callback):

    def __init__(self,test_path,out_path):
        super().__init__()
        self.test_path = test_path
        self.load_test_data(self.test_path)
        self.out_path = out_path

        if not os.path.exists(out_path+"/graphs/"):
            os.mkdir(out_path+"/graphs/")
        self.graphs_path = self.out_path+"/graphs/"

        if not os.path.exists(out_path+"/models/"):
            os.mkdir(out_path+"/models/")
        self.models_path = self.out_path + "/models/"

        self.log_path = self.out_path + "/Model_custom_metrics.log.tsv"

        with open(self.log_path, "w+") as logf:
            logf.write("epoch\tks_d\todds10\todds5\todds1\tavg_delta_stoich\tpearsonr\tspearmanr\n")

    def load_test_data(self,path):
        self.test_df = pandas.read_pickle(path).sample(frac=1)

    def on_epoch_end(self, epoch, logs=None):
        #save model
        self.model.save(self.models_path + f"Model_epoch_{epoch}.h5", overwrite=True)

        #print("Starting epoch testing")
        ##predict with most recent model
        #signals = np.array(self.test_df["signal"].values.tolist())
        #print(signals.shape, len(signals))
        #self.test_df["y_pred"] =self.model.predict(signals,batch_size=1024)


        ## get stoichiometry\predictions for each genomic site
        ## if validated range of TTT etc., take site with max predicted WT stoichiometry

        #WT_stoichs,IVT_stoichs,Validated_stoichs = [],[],[]
        #WT_p, IVT_p = [],[]

        #for range_index, df_range_sub in self.test_df.groupby("Range_index"):
        #    maxWT_stoich, IVT_stoich = -1,-1
        #    maxWT_ps, maxIVT_ps = [],[]
        #    print(len(df_range_sub),"before Nan removal")
        #    df_range_sub.dropna(subset=["y_pred"])
        #    print(len(df_range_sub),"after Nan removal")
        #    for gen_index, df_gen_sub in df_range_sub.groupby("Gen_index"):
        #        df_sub_WT = df_gen_sub[df_gen_sub["IVT-WT"] == "WT"]
        #        WT_stoich = len(df_sub_WT[df_sub_WT["y_pred"]>0.5]) / len(df_sub_WT)
        #        if WT_stoich > maxWT_stoich:
        #            maxWT_stoich = WT_stoich
        #            df_sub_IVT = df_gen_sub[df_gen_sub["IVT-WT"] == "IVT"]
        #            IVT_stoich = len(df_sub_IVT[df_sub_IVT["y_pred"] > 0.5]) / len(df_sub_IVT)
        #            maxWT_ps =df_sub_WT["y_pred"].values.tolist()
        #            maxIVT_ps = df_sub_IVT["y_pred"].values.tolist()

        #    WT_stoichs.append(maxWT_stoich)
        #    IVT_stoichs.append(IVT_stoich)
        #    WT_p.extend(maxWT_ps)
        #    IVT_p.extend(maxIVT_ps)
        #    Validated_stoichs.append(next(iter(df_range_sub["stoichiometry"])))
        #print(len(WT_stoichs),WT_stoichs[0])
        #print(len(IVT_stoichs),IVT_stoichs[0])
        #WT_stoichs = np.array(WT_stoichs)
        #IVT_stoichs = np.array(IVT_stoichs)
        #print(len(WT_stoichs),WT_stoichs[0])
        #print(len(IVT_stoichs),IVT_stoichs[0])
        #Validated_stoichs = np.array(Validated_stoichs)

        # get custom metrics

        #ks_d, odds10, odds5, odds1= self.evaluate_m1_separation(WT_p,IVT_p)
        #avg_delta = self.evaluate_stoich_separation(WT_stoichs,IVT_stoichs,Validated_stoichs,epoch)
        #pearson_r,spearman_r = self.evaluate_correlation(WT_stoichs,Validated_stoichs,epoch)

        #write epoch testing metrics in log file
        #with open(self.log_path, "a") as logf:
        #    logf.write("\t".join([str(x) for x in [epoch, ks_d, odds10, odds5, odds1,avg_delta, pearson_r,spearman_r]]) +"\n")
        
        #save model
        #self.model.save(self.models_path + f"Model_epoch_{epoch}.h5", overwrite=True)

    def evaluate_m1_separation(self,WT_p:list, IVT_p:list):
        try:
            ks = kstest(WT_p, IVT_p)
            ks_d = round(ks[0],4)
        except Exception as e:
            print("KS test failed from ",len(WT_p),len(IVT_p), "WT and IVT probs. Raised:\n", Exception)
            ks_d = "Nan"           

        odds10 = self.evaluate_topk_odds_ratio(WT_p,IVT_p,10)
        odds5 = self.evaluate_topk_odds_ratio(WT_p, IVT_p, 5)
        odds1 = self.evaluate_topk_odds_ratio(WT_p, IVT_p, 1)
        return ks_d, odds10, odds5, odds1

    def evaluate_topk_odds_ratio(self,WT_p, IVT_p, k):
        percentile= np.percentile(WT_p + IVT_p, int(100-k))
        odds_wt = len([x for x in WT_p if x >= percentile]) / len(WT_p)
        odds_IVT = len([x for x in IVT_p if x >= percentile]) / len(IVT_p)
        if odds_IVT > 0:
            return odds_wt/odds_IVT
        else:
            print(int(100-k), "percentile returned 0 IVT at threshold",percentile)
            return "inf"

    def evaluate_stoich_separation(self,WT_stoichs, IVT_stoichs,Validated_stoichs,epoch):
        delta = WT_stoichs - IVT_stoichs
        avg_delta = round(np.mean(delta),3)
        plt.scatter(Validated_stoichs,delta)
        plt.ylabel("WT-IVT stoich",size=15)
        plt.xlabel("Validated_stoich",size=15)
        plt.title(f"Epoch {epoch}; mean delta={avg_delta}",size=15)
        plt.ylim(-1.1,1.1)
        plt.xlim(-0.1,1.1)
        plt.tight_layout()
        plt.savefig(self.graphs_path +f"WT-IVT_separation_epoch{epoch}.png",dpi=100)
        plt.close()
        return avg_delta

    def evaluate_correlation(self,WT_stoichs, Validated_stoichs,epoch):
        try:
            spearman_r = round(spearmanr(WT_stoichs,Validated_stoichs)[0], 3)
            pearson_r = round(pearsonr(WT_stoichs,Validated_stoichs)[0], 3)
        except Exception as e:
            print (f"Correlation failed with exception in epoch {epoch}:\n",Exception,"\n\n")
            spearman_r, pearson_r = 0,0

        plt.plot([0,1],[0,1],c="red",linewidth=1)
        plt.scatter(Validated_stoichs, WT_stoichs)
        plt.ylabel("SWARM stoich", size=15)
        plt.xlabel("Validated stoich", size=15)
        plt.ylim(-0.1,1.1)
        plt.xlim(-0.1,1.1)
        plt.title(f"Epoch {epoch}; spearmanr={spearman_r}; pearsonr={pearson_r}", size=15)
        plt.tight_layout()
        plt.savefig(self.graphs_path + f"Validated_correlation_epoch{epoch}.png", dpi=100)
        plt.close()
        return pearson_r,spearman_r


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
                          help="Length of the vectors, default is 105",
                          default=105,
                          type=int,
                          required=False)

    OPTIONAL.add_argument("--features",
                          help="Number of vectors, default is 1",
                          default=1,
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
                          default=None,
                          required=False)

    OPTIONAL.add_argument("--test_data",
                          help="Path to file with testing signals and labels. Run create_model1_testing_data.py",
                          default=None,
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

    vpath_signals = INPUT_PATH + '_final_validation_data.p'
    vpath_labels = INPUT_PATH + '_final_validation_labels.p'

    signals = np.load(path_signals, allow_pickle=True)
    print('signal loading done')
    labels = np.load(path_labels, allow_pickle=True)
    vsignals = np.load(vpath_signals, allow_pickle=True)
    print('vsignal loading done')
    vlabels = np.load(vpath_labels, allow_pickle=True)


    signals, labels = shuffle(signals, labels)
    print('shuffle done')

    vsignals, vlabels = shuffle(vsignals, vlabels)

    labels_one_hot = tf.keras.utils.to_categorical(labels, 2)
    vlabels_one_hot = tf.keras.utils.to_categorical(vlabels, 2)
    signals = np.array(signals, dtype='float32')
    vsignals = np.array(vsignals, dtype='float32')
    
    signals=signals.reshape( signals.shape[0],VLEN, FEATURES)
    vsignals=vsignals.reshape( vsignals.shape[0],VLEN, FEATURES)
    print(Counter([str(x) for x in labels_one_hot]))
    print(Counter([str(x) for x in vlabels_one_hot]))

    cfg = {
        'input_size': ( VLEN, FEATURES),
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
            model.load_weights(ARGS.DL_model)

        model.compile(optimizer='adam',
                      loss='binary_crossentropy',
                      metrics=['accuracy'])
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

    
    
    if ARGS.test_data:
        CustomMetrics = LogCustomTestingMetrics(ARGS.test_data, ARGS.out_path)
        callback_lst=[early_stopping, checkpoint, lr_callback, csv_logger, CustomMetrics]
    else:
        callback_lst=[early_stopping, checkpoint, lr_callback, csv_logger]
    
    #fit the model
    history = model.fit(signals, labels,
                        batch_size=cfg['batch_size'], epochs=cfg['max_epoch'],
                        validation_data=(vsignals, vlabels), validation_batch_size=1024,
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

