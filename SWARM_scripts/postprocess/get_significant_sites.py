import pandas
import sys
import argparse
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--output",
                      help="Output path",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument('-m', '--modification',
                      choices = ["m6A", "pU", "m5C"],
                      required=True)

REQUIRED.add_argument('--kit',
                      choices = ["RNA002","RNA004"],
                      required=True)

OPTIONAL.add_argument('--context',
                      default="writer",
                      choices=["writer","all"],
                      required=False)

OPTIONAL.add_argument('-t', '--threshold',
                      default=-1,
                      type=float,
                      required=False)

parser._action_groups.append(OPTIONAL)

args = parser.parse_args()

df_SWARM = pandas.read_csv(args.input,sep="\t")

df_cutoffs = pandas.read_csv(f"{script_dir}/site-level_cutoffs.tsv",sep="\t")

def is_PUS7(seq):
    if seq[2] != "T":
        return False
    if seq[3] == "T":
        return False
    if seq[4] != "T":
        return False
    if seq[5] != "A":
        return False
    if seq[6] not in ["A","G"]:
        return False
    # print("PUS7")
    return True


def is_TRUB1(seq):
    if seq[2] != "G":
        return False
    if seq[3] != "T":
        return False
    if seq[4] != "T":
        return False
    if seq[5] != "C":
         return False
    if seq[7] != "A":
        return False
    #if seq[8] != "C":
    #    return False
    # print("TRUB1")
    return True


def is_DRACH(seq):
    if seq[2] not in ["A","G","T"]:
        return False
    if seq[3] not in ["G","A"]:
        return False
    if seq[4] != "A":
        return False
    if seq[5] !="C":
        return False
    if seq[6] not in ["A","C","T"]:
        return False

    return True

def is_NSUN6(kmer):
    if kmer[4]!= "C":
       return False
    if kmer[5] != "T":
        return False
    elif kmer[6] != "C":
        return False
    elif kmer[7] != "C":
       return False
    elif kmer[8] != "A":
         return False
    return True

def is_motif(seq):
    return is_PUS7(seq) or is_TRUB1(seq) or is_DRACH(seq) or is_NSUN6(seq)

if args.context == "writer":
    if "motif" in df_SWARM.columns:
        df_SWARM_sig= df_SWARM[df_SWARM["motif"].apply(is_motif) == True]
    else:
        df_SWARM_sig = df_SWARM[df_SWARM["site"].apply(is_motif) == True]
        df_SWARM_sig["start"] = df_SWARM_sig["position"] +4
        df_SWARM_sig["end"] = df_SWARM_sig["position"] +5
        df_SWARM_sig = df_SWARM_sig[["contig","start","end","site","coverage","stoichiometry","probability"]]
        #df_SWARM_sig.columns = ["contig","start","end","motif","coverage","stoichiometry","probability"]

if args.threshold == -1:
    cutoff = df_cutoffs[(df_cutoffs["Kit"] == args.kit) & (df_cutoffs["Modification"] == args.modification) & (df_cutoffs["Context"] == args.context)]["Cutoff"].values.tolist()[0]
else:
    cutoff = args.threshold

print("CUTOFF=",cutoff)

df_SWARM_sig = df_SWARM_sig[(df_SWARM_sig["stoichiometry"] > 0.1 ) & (df_SWARM_sig["probability"] > cutoff )]


df_SWARM_sig.to_csv(args.output,sep="\t",index=False,header=False)


