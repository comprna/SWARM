import time
import pickle
import argparse

parser = argparse.ArgumentParser()
OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input_bed",
                      help="Input bed file, lifted to genome with columns:\n"
                           "chr pos0 pos1    .... coverage stoichiometry site_confidence",
                      required=True)

REQUIRED.add_argument("-v", "--validated",
                      help="Input bed file with a set of validated sites in bed3 format",
                      required=True)

REQUIRED.add_argument("-o", "--output",
                      help="Input bed file with validated sites:\n"
                           "chr pos0 pos1 ",
                      required=True)

OPTIONAL.add_argument("-d", "--discard",
                      help='Path to set with sites to discard from analysis',
                      default=None
                      )

OPTIONAL.add_argument("-a", "--aggregate",
                      help='Method for handling aggregated transcripts at same genomic location.\n'
                           '<max/avg> ; Default is max',
                      default= "max"
                      )

OPTIONAL.add_argument("-m", "--metric",
                      help='Metric for selecting sites.'
                           '<model2/rate> ; Default is model2',
                      default= "model2"
                      )

OPTIONAL.add_argument("--min_stoich",
                      help='Minimum stoichiometry for selecting sites. (Between 0 and 1)',
                      default= -1.0,
                      type=float
                      )

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()
st = time.time()

validated = ARGS.validated
m2_path = ARGS.input_bed
out_path = ARGS.output

if ARGS.metric not in ["model2","rate"]:
    raise ("--metric must be either model2 or rate")

if ARGS.aggregate not in ["max","avg"]:
    raise ("--aggregate must be either max or avg")

if ARGS.discard:
    with open(ARGS.discard,"rb") as p:
        discard_set = pickle.load(p)
else:
    discard_set=set()

validated_set = set()
with open(validated) as f:
    for line in f:
        chr,pos0,pos = line.strip().split("\t")[:3]
        validated_set.add(f"{chr}_{pos}")


p_lst, all_p_lst = [],[]
validated_tested = 0

preds_dct_val = {}
preds_dct_all = {}


#parse site level info to and store site thesholds
with open(m2_path) as f:
    f.readline()
    for line in f:
        line_lst = line.strip().split("\t")
        coverage,stoich, prob = line_lst[-3:]
        contig,start,end = line_lst[:3]
        if "chr" not in contig:
            contig = "chr" + contig
        site_index = f"{contig}_{end}"

        if ARGS.metric == "rate":
            prob = stoich

        try:
            if float(stoich) > ARGS.min_stoich:
                pm2 = float(prob)
            else:
                pm2=0
        except:
            print(stoich, prob, "skipped")
            continue

        if site_index in discard_set:
            continue

        if ARGS.aggregate == "max":
            if site_index in preds_dct_all:
                preds_dct_all[site_index] = max(preds_dct_all[site_index], pm2) # store maximum probability if multiple transcripts
            else:
                preds_dct_all[site_index] = pm2 # store maximum probability if multiple transcripts
            if site_index in validated_set:
                if site_index in preds_dct_val:    # store maximum probability if multiple transcripts
                        preds_dct_val[site_index] = max(preds_dct_val[site_index], pm2)
                else:
                    preds_dct_val[site_index] = pm2

        if ARGS.aggregate == "avg":
            coverage = float(coverage)
            if site_index in preds_dct_all:
                predsum, covsum = preds_dct_all[site_index]
                predsum += pm2*coverage
                covsum+=coverage
                preds_dct_all[site_index] = [predsum, covsum] # store avg weighted probability if multiple transcripts
            else:
                preds_dct_all[site_index] = [pm2*coverage, coverage]

            if site_index in validated_set:
                if site_index in preds_dct_val:
                    predsum, covsum = preds_dct_val[site_index]
                    predsum += pm2 * coverage
                    covsum += coverage
                    preds_dct_val[site_index] = [predsum, covsum]  # store avg weighted probability if multiple transcripts
                else:
                    preds_dct_val[site_index] = [pm2*coverage, coverage]

if ARGS.aggregate == "max":
    all_p_lst = [item for key,item in preds_dct_all.items()]
    p_lst = [item for key,item in preds_dct_val.items()]

if ARGS.aggregate == "avg":
    all_p_lst = [item[0] / item[1] for key, item in preds_dct_all.items()]
    p_lst = [item[0] / item[1] for key, item in preds_dct_val.items()]


p_lst = sorted(p_lst)
all_p_lst = sorted(all_p_lst)

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


#two pointer to get n predicted at thresholds for validated sites
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

while index_t < len(thresholds):
    thresh = thresholds[index_t]
    out_vals.append([thresh, 0,0])
    index_t+=1

#two pointer to get n predicted at thresholds for all sites
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

while index_t < len(thresholds):
    thresh = thresholds[index_t]
    out_vals_all.append([thresh, 0,0])
    index_t+=1

with open(out_path,"w+") as of:
    of.write("site_threshold\tvalidated_called\tvalidated_rate\tall_called\tall_rate\tvalidated_precision\n")
    for index,row in enumerate(out_vals):
        thresh,validated_freq,validated_pred = row
        thresh,all_freq, all_pred = out_vals_all[index]
        if all_pred>0:
            called_validated_percent = validated_pred / all_pred
        else:
            called_validated_percent=0
        of.write("\t".join([str(x) for x in [thresh, validated_pred, validated_freq,all_pred,all_freq,called_validated_percent]]) + "\n")


print("all done in ", time.time() -st, "seconds")

