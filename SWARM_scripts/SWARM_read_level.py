import os
import sys
import time
import random
import argparse
import subprocess
import multiprocessing



## predefined paths

MODELS_PATH = "../SWARM_models/"
SCRIPTS_PATH = "../SWARM_scripts/"

MODEL_KMER_cpp = MODELS_PATH + "kmer_model/IVT_model_c++.csv"
MODEL_KMER_py = MODELS_PATH + "kmer_model/IVT_model.pickle"

PREPROCESS_CPP = SCRIPTS_PATH + "preprocess/SWARM_preprocess"
PREPROCESS_py = SCRIPTS_PATH + "preprocess/SWARM_preprocess.py"
regPredict = SCRIPTS_PATH + "predict/predict_model1_parallel.py"
modsamPredict = SCRIPTS_PATH + "predict/predict_model1_parallel_modsam.py"
onlyPredict = SCRIPTS_PATH + "predict/predict_model1_from_pickle.py"

mod_dct = {
    "pU": "T",
    "m5C": "C",
    "ac4C": "C",
    "m6A": "A"
}

def run_script(cmd,error_event):
    try:
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        # Print stdout and stderr from the script
        print(f"stdout from {cmd[0]}\n")
        print(process.stdout)
        sys.stdout.flush()

        print(f"stderr from {cmd[0]}\n")
        print(process.stderr)
        sys.stdout.flush()

    except subprocess.CalledProcessError as e:
        # The attribute to access stderr is e.stderr
        print(f"Error in {cmd[0]}: {e}\n\n")
        print(f"stderr:\n{e.stderr}\n")
        sys.stdout.flush()

        # Raise the error again to propagate it
        error_event.set()
        raise e
    

def main():
    parser = argparse.ArgumentParser(description="Detect RNA modifications at each parsed read/base\n")
    parser.add_argument("-m", "--RNAmod", required=False, help="\nTarget RNA modification\n")
    parser.add_argument("-o", "--out", required=True, help="\nPrefix to the output file. Outputs <prefix>.pred.tsv\n")

    parser.add_argument("--mode", required=False, default="parallel", help="\nparallel / preprocess / predict \n")
    parser.add_argument("-l", "--label", required=False, default="1", help="\nLabel for condition such as WT_rep1 \n")

    parser.add_argument("-s", "--sam", required=False, help="\nPath to the input sam event align\n")
    parser.add_argument("-f", "--fasta", required=False, help="\nPath to the input fasta reference genome\n")
    parser.add_argument("-r", "--raw", required=False, help="\nPath to the input signals in blow5 format\n")
    parser.add_argument("--temp", required=False, help="\nDirectory for temp files\n")
    parser.add_argument("--modsam", help="\nMakes sam file (MM/ML tags) as <prefix>.mod.sam\n", action="store_true")

    parser.add_argument("-n", "--nanopolish", required=False, help="\nPath to the input nanopolish event align\n")
    parser.add_argument("-b", "--bam", required=False, help="\nPath to the input bam file\n")
    parser.add_argument("--base", required=False, help="\nTarget base to preprocess\n")
    parser.add_argument("-t", "--threads", required=False, default=1,type=int, help="\nNumber of threads for preprocessing\n")

    parser.add_argument("-p", "--pickle", required=False, help="\nPath to the input pickle file (model1 output)\n")
    parser.add_argument("--out_counter",help="Output counts of all 9mers for each split nanopolish file",action="store_true")

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    if args.mode == "parallel": # preprocess C++ and predict python (default option)
        if not(args.sam and args.fasta and args.raw and args.RNAmod):
            raise ("Missing args for default parallel prediction.\n --sam --fasta --raw --RNAmod    must all be provided!")

        else:
            RNAmod = args.RNAmod
            try:
                BASE = mod_dct[RNAmod]
            except KeyError:
                raise ("--RNAmod must be one of pU / m6A / m5C / ac4C ;\ncase sensitive!")

            MODEL1_PATH = os.path.join(script_dir,MODELS_PATH + f"Model1/{RNAmod}/Model_100_epoch_relu.h5")

            if args.temp: # if temporary directory provided (ideally on SSD)
                TMP = args.temp
                if TMP[-1] == "/":
                    TMP = TMP[:-1]

                TEMP = f"{TMP}/t{hash(args.out)}{random.randint(0, 10000)}_tmp"

            else: # put temp files in provided output dir
                TEMP = f"{os.path.dirname(args.out)}/t{hash(args.out)}{random.randint(0, 10000)}_tmp"

            if args.modsam: # run prediction script which makes regular+modsam output
                SCRIPT_py = os.path.join(script_dir,modsamPredict)
                print("running modsam prediction")

                args_script_py = ["python3",os.path.join(script_dir,SCRIPT_py),
                              "-i", TEMP,
                              "-o", args.out,
                              "-m", MODEL1_PATH,
                              "--sam", args.sam,
                              "-l", args.label]

            else: # run prediction script which makes just regular tsv output
                SCRIPT_py = os.path.join(script_dir,regPredict)
                print("running regular prediction")
                args_script_py = ["python3",os.path.join(script_dir,SCRIPT_py),
                              "-i", TEMP,
                              "-o", args.out,
                              "-m", MODEL1_PATH,
                              "-l", args.label]

            # Define arguments for each script
            args_script_cpp = [os.path.join(script_dir,PREPROCESS_CPP),
                               "--sam", args.sam,
                               "--raw", args.raw,
                               "--fasta", args.fasta,
                               "-m", os.path.join(script_dir,MODEL_KMER_cpp),
                               "-o", TEMP,
                               "--base", BASE]

            error_event = multiprocessing.Event()

            process_script_cpp = multiprocessing.Process(target=run_script, args=[args_script_cpp,error_event])
            process_script_py = multiprocessing.Process(target=run_script, args=[args_script_py,error_event])
            # Start both processes
            process_script_py.start()
            process_script_cpp.start()

            while process_script_py.is_alive() or process_script_cpp.is_alive():
                # Check if there was an error in either process
                if error_event.is_set():
                    print("Terminating due to an error in one of the processes.")
                    # Terminate the other process
                    if process_script_py.is_alive():
                        process_script_py.terminate()
                        process_script_py.join()
                    if process_script_cpp.is_alive():
                        process_script_cpp.terminate()
                        process_script_cpp.join()
                    exit(1)

                # Add a delay to avoid excessive checking
                time.sleep(20)
            process_script_py.join()
            process_script_cpp.join()
            print("Both processes have completed successfully.")



    elif args.mode == "preprocess":
        if not(args.nanopolish and args.bam):
            raise ("Missing args for preprocess mode. --nanopolish and --bam    must both be provided.")
        else:
            args_script_py = ["python3", os.path.join(script_dir,PREPROCESS_py),
                              "-i", args.nanopolish,
                              "-o", args.out,
                              "-m", os.path.join(script_dir,MODEL_KMER_py),
                              "-b", args.bam,
                              "-n", args.threads]
            if args.out_couter:
                args_script_py.append("--out_counter")
            if args.base:
                args_script_py+= ["--input_base",args.base]
            subprocess.run(args_script_py)

    elif args.mode == "predict":
        if not(args.pickle and args.RNAmod):
            raise ("Missing args for preprocess mode. --pickle --RNAmod   must be provided.")
        else:
            RNAmod = args.RNAmod

            try:
                BASE = mod_dct[RNAmod]
            except KeyError:
                raise ("--RNAmod must be one of pU / m6A / m5C / ac4C ;\ncase sensitive!")

            MODEL1_PATH = os.path.join(script_dir,MODELS_PATH + f"Model1/{RNAmod}/Model_100_epoch_relu.h5")
            args_script_py = ["python3", os.path.join(script_dir,onlyPredict),
                              "-i", args.pickle,
                              "-o", args.out,
                              "-m", MODEL1_PATH,
                              "-l",args.label]


            subprocess.run(args_script_py)
    else:
        raise ("--mode must be one of:   parallel / preprocess / predict")




if __name__ == "__main__":
    main()

