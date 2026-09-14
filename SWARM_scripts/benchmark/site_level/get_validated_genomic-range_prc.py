import sys
import time
import pickle

st = time.time()
### get prediction rate on WT

validated_ranges = sys.argv[1]
predicted_path = sys.argv[2]
out_path = sys.argv[3]

validated_range_dict = {}
with open(validated_ranges) as f:
    f.readline()
    for line in f:
        try:
            chr, coordinate,sequence,plus_chain,minus_chain = line.strip().split(" ")
        except:
            print(line)
            exit()
        if "chr" not in chr:
            chr="chr"+chr
        coordinate = int(coordinate)
        for offset in range(-1*int(minus_chain),int(plus_chain)+1):
            validated_range_dict[f"{chr}_{coordinate + offset}"] =f"{chr}_{coordinate}"
        # validated_set = pickle.load(p)
#validated_set = set()

p_lst, all_p_lst = [],[]
validated_tested = 0

preds_dct_val = {}
preds_dct_all = {}

with open(predicted_path) as f:
    f.readline()
    for line in f:
        line_lst = line.strip().split("\t")
        stoich, prob = line_lst[-2:]
        contig,start,end = line_lst[:3]
        if "chr" not in contig:
            contig = "chr" + contig
        site_index = f"{contig}_{end}"
        
        try:
            if float(stoich) > 0.1:
                pm2 = float(prob)
            else:
                pm2=0
        except:
            print(stoich, prob, "skipped")
            continue
         


        if site_index in validated_range_dict:
            site_index = validated_range_dict[site_index] # overwrite to BID-seq coordinate of TTT range

            if site_index in preds_dct_val:
                    preds_dct_val[site_index] = max(preds_dct_val[site_index], pm2)
            else:
                preds_dct_val[site_index] = pm2

        if site_index in preds_dct_all:
            preds_dct_all[site_index] = max(preds_dct_all[site_index], pm2)
        else:
            preds_dct_all[site_index] = pm2

all_p_lst = [item for key,item in preds_dct_all.items()]
p_lst = [item for key,item in preds_dct_val.items()]

# p_lst = [0.1,0.11,0.2,0.4,0.8,0.89,0.9001,1,1,0.91995,0.9891,0.999994]
p_lst = sorted(p_lst)
all_p_lst = sorted(all_p_lst)
#print(p_lst)
thresholds, fpr_lst = [],[]

# creates range of thresholds 0.9, 0.901, 0.902 ... 0.989, 0.990, 0.9901 ... 0.9990 ... 0.99999999989 ... 1
thresholds = [0, 1/10**7, 1/10**6, 1/10**5,1/10**4]+[x/1000 for x in range(1,900)]
base_n_lst = [0.9]
for n1 in range(1,17):
    base_n  = str(base_n_lst[-1])
    for n2 in range(0,9):
        base_n_lst.append(float(base_n + str(n2)))
        for n3 in range(0,10):
            thresholds.append(float(base_n+str(n2)+str(n3)))
    base_n_lst.append(float(base_n + str(9)))
thresholds.append(1)

index_p, index_t = 0,0

out_vals = []
while index_p < len(p_lst) and index_t < len(thresholds):
    prob = p_lst[index_p]
    thresh = thresholds[index_t]
    if prob < thresh:
        index_p+=1
    else:
        index_t+=1
        out_vals.append([thresh,1- index_p/len(p_lst),len(p_lst) - index_p])
        # print([thresh,1- index_p/len(p_lst),prob])

while index_t < len(thresholds):
    thresh = thresholds[index_t]
    out_vals.append([thresh, 0,0])
    # print([thresh, 0])
    index_t+=1




index_p, index_t = 0,0

out_vals_all = []
while index_p < len(all_p_lst) and index_t < len(thresholds):
    prob = all_p_lst[index_p]
    thresh = thresholds[index_t]
    if prob < thresh:
        index_p+=1
    else:
        index_t+=1
        out_vals_all.append([thresh,1- index_p/len(all_p_lst),len(all_p_lst) - index_p])
        # print([thresh,1- index_p/len(p_lst),prob])

while index_t < len(thresholds):
    thresh = thresholds[index_t]
    out_vals_all.append([thresh, 0,0])
    # print([thresh, 0])
    index_t+=1

with open(out_path,"w+") as of:
    of.write("site_threshold\tvalidated_called\tvalidated_rate\tall_called\tall_rate\tvalidated_precision\n")
    for index,row in enumerate(out_vals):
        thresh,validated_freq,validated_pred = row
        thresh,all_freq, all_pred = out_vals_all[index]
        if all_pred > 0:
            called_validated_percent = validated_pred / all_pred
        else:
            called_validated_percent = 0

        of.write("\t".join([str(x) for x in [thresh,validated_pred, validated_freq,all_pred,all_freq,called_validated_percent]]) + "\n")

# print(len(thresholds))

print("all done in ", time.time() -st, "seconds")

