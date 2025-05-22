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

MODEL_KMER = MODELS_PATH + "kmer_model/model_5-mer"
PREPROCESS_CPP = SCRIPTS_PATH + "preprocess/SWARM_preprocess"
PREPROCESS_py = SCRIPTS_PATH + "preprocess/SWARM_preprocess.py"
split_script = SCRIPTS_PATH + "preprocess/split_bams.py"
regPredict = SCRIPTS_PATH + "predict/predict_model1_parallel.py"
modsamPredict = SCRIPTS_PATH + "predict/predict_model1_parallel_modsam.py"
onlyPredict = SCRIPTS_PATH + "predict/predict_model1_from_pickle.py"
CHECK_KIT_SCRIPT= SCRIPTS_PATH + "preprocess/check_RNA_kit"

mod_dct = {
    "pU": "T",
    "m5C": "C",
    "ac4C": "C",
    "m6A": "A"
}

arch_dct = {
    "pU_RNA002": "Mini",
    "pU_RNA004": "Mini",
    "m5C_RNA002": "Mid",
    "m5C_RNA004": "Mini",
    "m6A_RNA002": "Large",
    "m6A_RNA004": "Large"
}

label_dct = {
    "pU_RNA002": "1",
    "pU_RNA004": "2",
    "m5C_RNA002": "3",
    "m5C_RNA004": "4",
    "m6A_RNA002": "5",
    "m6A_RNA004": "6"
}

def run_script(cmd,error_event):
    try:
        print("starting process with command",cmd)
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("PROCESS successfully started")
        # Print stdout and stderr from the script
        print(f"stdout from {cmd[0]}\n")
        print(process.stdout.decode("utf-8"))
        sys.stdout.flush()

        print(f"stderr from {cmd[0]}\n")
        print(process.stderr.decode("utf-8"))
        print(process.stderr.decode("utf-8"), file=sys.stderr)
        sys.stdout.flush()

    except subprocess.CalledProcessError as e:
        # The attribute to access stderr is e.stderr
        print(f"Error in {cmd[0]}: {e}\n\n")
        print(f"stderr:\n{e.stderr.decode('utf-8')}\n")
        print(f"stderr:\n{e.stderr.decode('utf-8')}\n", file=sys.stderr)
        sys.stdout.flush()

        # Raise the error again to propagate it
        error_event.set()
        raise e
    
def run_script2(cmd, error_event):
    try:
        # Start the subprocess with stdout and stderr pipes
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True  # Assuming you're using Python 3.7 or later for text output
        )

        # Read and print the output while the process is running
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(f"stdout from {cmd[0]}: {output.strip()}")
                # Flush stdout to ensure immediate visibility
                sys.stdout.flush()

            error = process.stderr.readline()
            if error:
                print(f"stderr from {cmd[0]}: {error.strip()}")
                sys.stdout.flush()

        # Wait for the process to complete and get the final output
        stdout, stderr = process.communicate()

        # Print the final output
        print(f"Final stdout from {cmd[0]}: {stdout.strip()}")
        print(f"Final stderr from {cmd[0]}: {stderr.strip()}")

        # Check if the process exited with an error
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)

    except subprocess.CalledProcessError as e:
        # Raise the error again to propagate it
        error_event.set()
        raise e

def run_check_kit(binary_path, blow5):
    """
    Runs a C++ binary and captures its stdout as a string.
    
    :param binary_path: Path to the C++ binary.
    :param args: List of arguments to pass to the binary.
    :return: stdout output from the binary as a string.
    """
    args = ["--raw", blow5]
    try:
        # Construct the command to run
        command = [binary_path] + args

        # Run the command and capture output
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,  # Ensures the output is captured as a string
            check=True  # Raises an exception if the binary exits with a non-zero code
        )

        # Return the captured stdout
        return result.stdout.strip()

    except subprocess.CalledProcessError as e:
        print(f"Error: The binary exited with code {e.returncode}.")
        print(f"stderr:\n{e.stderr.strip()}")
        return None
    except FileNotFoundError:
        print("Error: The specified binary could not be found.")
        return None



def main():
    parser = argparse.ArgumentParser(description="Detect RNA modifications at each parsed read/base\n")
    parser.add_argument("-m", "--RNAmod", required=False, help="\nTarget RNA modification\n")
    parser.add_argument("-o", "--out", required=True, help="\nPrefix to the output file. Outputs <prefix>.pred.tsv\n")

    parser.add_argument("--mode", required=False, default="parallel", help="\nparallel / preprocess / predict \n")
    parser.add_argument("-s", "--sam", required=False, help="\nPath to the input sam event align\n")
    parser.add_argument("-f", "--fasta", required=False, help="\nPath to the input fasta reference genome\n")
    parser.add_argument("-r", "--raw", required=False, help="\nPath to the input signals in blow5 format\n")
    
    parser.add_argument("--model1", required=False,default=None, help="\nPath to the trained model1\n")
    parser.add_argument("--kmer", required=False,default=None, help="\nPath to the kmer model\n")
    parser.add_argument("--cpp", required=False,default=None, help="\nPath to compiled c++ preprcessing\n")
    parser.add_argument("--kit", required=False,default="RNA002", help="\nRNA sequencing kit [RNA004/RNA002]\n")
    parser.add_argument("--temp", required=False, help="\nDirectory for temp files\n")
    parser.add_argument("--arch", required=False,default=None, help="\nModel1 network from Mini/Mid/Large. Mini is default\n")
    parser.add_argument("--modsam", help="\nMakes sam file (MM/ML tags) as <prefix>.mod.sam\n", action="store_true")
    parser.add_argument("--nworkers", help="\nNumber of C++ preprocessing workers\n", default=4,type=int)
    parser.add_argument("--limit", help="\nNumber of signals to preprocess\n", default=-1,type=int)

    parser.add_argument("-n", "--nanopolish", required=False, help="\nPath to the input nanopolish event align\n")
    parser.add_argument("-b", "--bam", required=False, help="\nPath to the input bam file\n")
    parser.add_argument("--base", required=False, help="\nTarget base to preprocess\n")
    parser.add_argument("-t", "--threads", required=False, default=1,type=int, help="\nNumber of threads for preprocessing\n")

    parser.add_argument("-p", "--pickle", required=False, help="\nPath to the input pickle file (model1 output)\n")

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
            

            CHECK_KIT_PATH = os.path.join(script_dir,CHECK_KIT_SCRIPT)
            # check for the sequencing kit if blow5 provided
            kit = run_check_kit(CHECK_KIT_PATH,args.raw)
            if not kit:
                print("Error detecting kit from the blow5. Using default kit.\n", file=sys.stderr)
                KIT="RNA002"
            else:
                print("detected kit is",kit, file=sys.stderr)
                if kit == "sqk-rna004":
                    KIT="RNA004"
                elif kit == "sqk-rna002" or kit == "sqk-rna001":
                    KIT="RNA002"
                else:
                    print("Unknown kit detected from the blow5. Using default kit.\n", file=sys.stderr)
                    KIT="RNA002"
            print("Using kit:",KIT, "\n")
            
            if not args.arch:
                ARCH = arch_dct[f"{RNAmod}_{KIT}"]
            else:
                ARCH= args.arch
            
            print("Using model arch:",ARCH, "\n")

            if args.model1:
                MODEL1_PATH = args.model1
            else:
                MODEL1_PATH = os.path.join(script_dir,MODELS_PATH + f"Model1/{KIT}/{RNAmod}/Model_100_epoch_relu.h5")
            print("Model1 =",MODEL1_PATH)
            if args.temp: # if temporary directory provided (ideally on SSD)
                TMP = args.temp
                if TMP[-1] == "/":
                    TMP = TMP[:-1]

                TEMP = f"{TMP}/t{hash(args.out)}{random.randint(0, 10000)}_tmp"

            else: # put temp files in provided output dir
                TMP=f"{os.path.dirname(args.out)}/.SWARM_tmp/"
                TEMP = f"{TMP}/t{hash(args.out)}{random.randint(0, 10000)}_tmp"
                if not os.path.exists(TMP):
                    os.mkdir(f"{os.path.dirname(args.out)}/.SWARM_tmp/")
            if args.modsam: # run prediction script which makes regular+modsam output
                SCRIPT_py = os.path.join(script_dir,modsamPredict)
                print("running modsam prediction")

                args_script_py = ["python3",os.path.join(script_dir,SCRIPT_py),
                              "-i", TEMP,
                              "-o", args.out,
                              "-m", MODEL1_PATH,
                              "--sam", args.sam,
                              "--arch", ARCH,
                              "-l", label_dct[f"{RNAmod}_{KIT}"]]

            else: # run prediction script which makes just regular tsv output
                SCRIPT_py = os.path.join(script_dir,regPredict)
                print("running regular prediction")
                args_script_py = ["python3",os.path.join(script_dir,SCRIPT_py),
                              "-i", TEMP,
                              "-o", args.out,
                              "-m", MODEL1_PATH,
                              "-l", label_dct[f"{RNAmod}_{KIT}"],
                              "--arch", ARCH,
                              "--nworkers", str(args.nworkers),
                              "--limit",str(args.limit)]
            


            if args.kmer:
                model_kmer_path = args.kmer
            else:
                model_kmer_path=os.path.join(script_dir,f"{MODEL_KMER}.{KIT}.csv")

            if args.cpp:
                cpp_bin = args.cpp
            else:
                cpp_bin = os.path.join(script_dir,PREPROCESS_CPP)


            ## split bam into nworkers chunks
            sam_filename= f"{os.path.dirname(args.out)}/{args.sam.split('/')[-1]}"
            args_split_script = ["python3",os.path.join(script_dir,split_script),
                              "-i", args.sam,
                              "-o", sam_filename,
                              "-n", str(args.nworkers)]

            error_event_split = multiprocessing.Event()
            process_split = multiprocessing.Process(target=run_script, args=[args_split_script, error_event_split])
            
            print("starting split")
            split_start_time = time.time()
            process_split.start()
            while process_split.is_alive():
                # Check if there was an error in either process
                if error_event_split.is_set():
                    print("\nTerminating due to an error in split processes.\n",file=sys.stderr)
                    # Terminate the other process
                    if process_split.is_alive():
                        process_split.terminate()
                        process_split.join()
                    exit(1)
                # Add a delay to avoid excessive checking
                time.sleep(2)
            process_split.join()
            print("Split done in ",time.time() - split_start_time)

            error_event = multiprocessing.Event()
            cpp_proc_list = []
            # Define arguments for each script
            for cpp_proc_i in range(args.nworkers):
                args_script_cpp = [cpp_bin,
                                   "--sam", f"{sam_filename}_{cpp_proc_i+1}.sam",
                                   "--raw", args.raw,
                                   "--fasta", args.fasta,
                                   "-m", model_kmer_path,
                                   "-o", f"{TEMP}_{cpp_proc_i}",
                                   "--base", BASE]



                process_script_cpp = multiprocessing.Process(target=run_script, args=[args_script_cpp,error_event])
                process_script_cpp.start()
                cpp_proc_list.append(process_script_cpp)

            process_script_py = multiprocessing.Process(target=run_script, args=[args_script_py,error_event])
            process_script_py.start()

            while process_script_py.is_alive() or any([process_script_cpp.is_alive() for process_script_cpp in cpp_proc_list]):
                # Check if there was an error in either process
                if error_event.is_set():
                    print("Terminating due to an error in one of the processes.")
                    # Terminate the other process
                    if process_script_py.is_alive():
                        process_script_py.terminate()
                        process_script_py.join()

                    for process_script_cpp in cpp_proc_list:
                        if process_script_cpp.is_alive():
                            process_script_cpp.terminate()
                            process_script_cpp.join()
                    exit(1)

                # Add a delay to avoid excessive checking
                #print("CPP",any([process_script_cpp.is_alive() for process_script_cpp in cpp_proc_list]))
                #print("py",process_script_py.is_alive())
                if not process_script_py.is_alive():
                    time.sleep(1)
                    for process_script_cpp in cpp_proc_list:
                        if process_script_cpp.is_alive():
                            print("Killing CPP worker as python prediction is no longer alive")
                            process_script_cpp.terminate()
                time.sleep(20)
            process_script_py.join()
            for process_script_cpp in cpp_proc_list:
                process_script_cpp.join()
            
            time.sleep(10)
            TMP_prefix = TEMP.split("/")[-1]
            for f in os.listdir(TMP):          # cleanup the temp files if any remained
                if TMP_prefix in str(f):
                    os.remove(TMP +"/" + f)

            for cpp_proc_i in range(args.nworkers):
                os.remove(f"{sam_filename}_{cpp_proc_i+1}.sam")

            print("All processes have completed in",time.time() - split_start_time , "seconds")



    elif args.mode == "preprocess":
        if not(args.nanopolish and args.bam):
            raise ("Missing args for preprocess mode. --nanopolish and --bam    must both be provided.")
        else:
            args_script_py = ["python3", os.path.join(script_dir,PREPROCESS_py),
                              "-i", args.nanopolish,
                              "-o", args.out,
                              "-m", os.path.join(script_dir,f"{MODEL_KMER}.{KIT}.csv"),
                              "-b", args.bam,
                              "-n", args.threads]
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
                              "-l", "0"]

            subprocess.run(args_script_py)
    else:
        raise ("--mode must be one of:   parallel / preprocess / predict")




if __name__ == "__main__":
    main()

