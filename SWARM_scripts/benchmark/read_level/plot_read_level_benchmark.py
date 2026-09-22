from sklearn.metrics import precision_recall_curve, auc
import random
import matplotlib.pyplot as plt

# example key: cc6m_2244_t7_ecorv_2075_GAGTCAAAG_b654b5be-3f1e-43ae-a5d2-022f1fb32f70_12_C
def get_prob_lst(path):
        out_dct = {}
        with open(path) as f:
                for line in f:
                        key, p, l = line.strip().split("\t")
                        key_lst = key.split("_")
                        chr = "_".join(key_lst[:-6])
                        pos = key_lst[-6]
                        pos = int(pos) + 5
                        index = f"{chr}_{pos}"
                        if index not in out_dct:
                                out_dct[index] = [float(p)]
                        elif len(out_dct[index]) < 500:
                                out_dct[index].append(float(p))
        return [p for key,item in out_dct.items() for p in item]

target_mod_lst = ["m6A","m5C","pU"]
sample_mod_lst =  ["m6A","m5C","pU","NM"]

for target_mod in target_mod_lst:
        positive_probs = get_prob_lst(f"outputs/IVT-m1-m2.sample-{target_mod}.target-{target_mod}.pred.tsv")
        unmodified_base = target_mod[-1]
        for negative_mod in sample_mod_lst:
                if negative_mod == target_mod:
                        continue
                negative_probs = get_prob_lst(f"outputs/IVT-m1-m2.sample-{negative_mod}.target-{target_mod}.pred.tsv")
                min_len = min(len(positive_probs), len(negative_probs))
                ypred = random.sample(positive_probs,min_len) + random.sample(negative_probs, min_len)
                ytrue = [1] * min_len + [0] * min_len

                precision, recall, thresholds = precision_recall_curve(ytrue, ypred)
                pr_auc = auc(recall, precision)  # Area under the PR curve
                plt.plot(recall, precision, label = f"SWARM AUC = {round(pr_auc,3)}")
                plt.xlabel('Recall',fontsize=20)
                plt.ylabel('Precision',fontsize=20)
                plt.title(f'RNA002 {target_mod.replace("pU","Ψ")} vs {unmodified_base} ({negative_mod.replace("pU","Ψ")} IVT)',fontsize=20)
                plt.ylim(0.45,1.05)
                plt.legend()
                plt.savefig(f"IVT-m1-m2.sample-{negative_mod}.target-{target_mod}.prc.png",dpi=300)
                plt.close()




