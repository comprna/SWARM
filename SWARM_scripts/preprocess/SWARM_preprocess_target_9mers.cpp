#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring>
#include <dirent.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include "argagg.hpp"
#include "htslib/sam.h"
#include <htslib/faidx.h>
#include "slow5/slow5.h"


#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))
using namespace std;


/// To fix: dependency on absolute path for using htslib


struct kmerStruct {
    array<int, 2> eventIndex;
    float* smoothCurrents;
    int dwell;
};


vector<string> split_string(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;
    while (getline(ss, tok, delimiter))
    {
        internal.push_back(tok);
    }
    return internal;
}


std::random_device rd;
std::mt19937 gen(42);

float* interpolateSignals(float* signals, int original_len, int target_scale=36){
    float* out = new float[target_scale];
    int interp_op = original_len-1; int interp_target = target_scale-1;
    int shuffle_steps[interp_op];
    int min_step_size= interp_target / interp_op;
    int max_step_size = min_step_size+1;
    int remainder = interp_target % interp_op;
    for (int i=0; i<remainder; i++){shuffle_steps[i]=max_step_size;}
    for (int i=remainder; i<interp_op; i++){shuffle_steps[i]=min_step_size;}
    shuffle(&shuffle_steps[0],&shuffle_steps[interp_op],gen);
    int out_index=0;
    for (int index=0; index<interp_op;index++){
        int window_size = shuffle_steps[index];
        float low_value = signals[index];
        float up_value = signals[index+1];
        float increment = (up_value-low_value) / window_size;
        float add=0;
        for (int j=0; j<window_size; j++){
            out[out_index++]= low_value+add;
            add+=increment;
        }
    }
    out[interp_target]=signals[interp_op];
    return out;
}


float* compresSignal(float* signals, int original_len, int target_len = 36){
    int min_step_size = original_len / target_len;
    int max_step_size = min_step_size+1;
    int remainder = original_len % target_len;
    int window_values[target_len];
    for (int i=0; i< remainder; i++){window_values[i]=max_step_size;}
    for (int i=remainder; i< target_len; i++){window_values[i]=min_step_size;}
    shuffle(&window_values[0], &window_values[target_len],gen);
    int index = 0;
    int out_index = 0;
    float* out = new float[target_len];
    for (int window_size : window_values){
        float suma = 0.0;
        for (int j = 0; j< window_size;j++){
            suma+= signals[index++];
        }
        out[out_index++] = suma / window_size;
    }
    return out;
}

unordered_map<string, int> get_target_umap(const string& filename) {
    unordered_map<string, int> umap;

    // Open the file
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return umap;
    }

    // Read file line by line
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string col1;
        int col2;

        // Extract col1 and col2 from each line
        if (iss >> col1 >> col2) {
            // Insert col2 into the unordered_set corresponding to col1
            umap[col1] = col2;
        } else {
            cerr << "Error: Invalid line format in file " << filename << endl;
        }
    }

    // Close the file
    file.close();

    return umap;
}



unordered_map<string, float> readModelKmer(const char* filename) {
    unordered_map<string, float> csvMap;

    ifstream inputFile(filename);
    if (!inputFile) {
        cerr << "Error opening file: " << filename << endl;
        return csvMap; // Return an empty map on error
    }
    string line;
    while (getline(inputFile, line)) {
        istringstream lineStream(line);
        string key, value;
        if (getline(lineStream, key, ',') && getline(lineStream, value, ',')) {
            // Use the first column as the key and the second column as the value
            try {
                float num_value = stod(value);
            csvMap[key] = num_value;
            // cout << key << "  " << value << endl;
            } catch (const std::invalid_argument& e){
                // skip first line 
            }
            
        } else {
            cerr << "Skipping invalid line: " << line << endl;
        }
    }

    inputFile.close();
    return csvMap;
}



std::string getDirectoryFromPath(const std::string& filePath) {
    size_t found = filePath.find_last_of("/"); // Platform-independent path separator
    if (found != std::string::npos) {
        return filePath.substr(0, found);
    }
    std:: cout << "FUCK";
    return ""; // Return an empty string if the path doesn't contain a directory part
}


std::string getFileNameFromPath(const std::string& filePath) {
    size_t lastSlash = filePath.find_last_of("/"); // Platform-independent path separator
    if (lastSlash != std::string::npos) {
        return filePath.substr(lastSlash + 1);
    }
    return filePath; // If there's no directory separator, the input is already a file name
}


int main(int argc, char *argv[]) {

    if (argc < 2)
    {
        std::cout << "Run with: ./SWARM_preprocess  --sam <events.sam> --fasta <ref.fa> --raw <signals.blow5> -o <outpath> --base <A/C/G/T> -m <model_kmer_path>" << std::endl;
        return EXIT_FAILURE;
    }

    argagg::parser argparser{{
        {"help", {"-h", "--help"}, "shows this help message\n", 0},
        {"samFileName", {"--sam"}, "SAM file from f5c. <sample.f5c.sam> (required)\n", 1},
        {"fastaFileName", {"--fasta"}, "Fasta file. <ref.fa> (required)\n", 1},
        {"slow5FileName", {"--raw"}, "Slow5/blow5 file. <sample.blow5> (required)\n", 1},
        {"model_kmer_path", {"-m", "--kmer-model"}, "file containing all the expected signal k-mer means (required)\n", 1},
        {"name_out", {"-o", "--out"}, "output file prefix (required)\n", 1},
        {"INPUT_BASE", {"--base"}, "Base to preprocess (required)\n", 1},
        {"threads", {"-n", "--CPU"}, "Number of CPUs (threads) to use (default: 1)\n", 1},
        {"TARGETS", {"--targets"}, "Target poisitions to preprocess where 1st column is contig and 2nd is 0-based position", 1},
        
    }};

    argagg::parser_results args;

    try
    {
        args = argparser.parse(argc, argv);
    }
    catch (const exception &e)
    {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }

    if (args["help"])
    {
        cerr << argparser;
        return EXIT_SUCCESS;
    }

    if (!args["samFileName"] | !args["fastaFileName"] | !args["slow5FileName"] | !args["model_kmer_path"] | !args["name_out"] |!args["INPUT_BASE"] |!args["TARGETS"])
    {
        cerr << "ERROR: Please provide all of required args: --sam --fasta --raw -m -o -i --targets" << endl;
        cerr << argparser;
        return EXIT_FAILURE;
    }

    string samFileName_s = args["samFileName"].as<string>(""); const char* samFileName = samFileName_s.c_str();
    string slow5FileName_s = args["slow5FileName"].as<string>(""); const char* slow5FileName = slow5FileName_s.c_str();
    string fastaFileName_s = args["fastaFileName"].as<string>(""); const char* fastaFileName = fastaFileName_s.c_str();
    string model_kmer_path_s = args["model_kmer_path"].as<string>(""); const char* model_kmer_path = model_kmer_path_s.c_str();
    string name_out = args["name_out"].as<string>("");
    string targets_path = args["TARGETS"].as<string>("");
    const char INPUT_BASE = args["INPUT_BASE"].as< string>("")[0];
    int n_threads = args["threads"].as<int>(1);
    
    cout << name_out << endl ;
    
    array<char,4> BASES = {'A','C','G','T'};

    int BASE_INDEX;
    if (INPUT_BASE == 'A'){
        BASE_INDEX = 0;
    } else if (INPUT_BASE == 'C'){
        BASE_INDEX = 1;
    } else if (INPUT_BASE == 'G'){
        BASE_INDEX = 2;
    } else if (INPUT_BASE == 'T'){
        BASE_INDEX = 3;
    } else {
        cerr << "input_base must be one of A/C/G/T" << endl;
        return 1;
    }

    int FILE_NUM = 0;
    int BATCH_SIZE = 32768;
    int PREV_PARSES = 0;
    

    string directoryPathStr = getDirectoryFromPath(name_out);
    const char* directoryPath = directoryPathStr.c_str();
    DIR* dir;
    struct dirent* entry;

    // string PREP_FILE = name_out + "_" + to_string(FILE_NUM);
    string PREP_FILE = name_out;
    std::ofstream outFile(PREP_FILE, std::ios::binary | std::ios::out);
    if (!outFile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
    }

    faidx_t* fai = fai_load(fastaFileName);

    if (!fai){
        cerr << "Error opening fasta." << endl;
        return 1;
    }

    samFile* samFile = hts_open(samFileName, "r"); // Open SAM file

    if (!samFile) {
        cerr << "Error opening SAM file." << endl;
        return 1;
    }


    // load blow5/slow5 file
    slow5_file_t* slow5File = slow5_open(slow5FileName, "r");

    if (!slow5File) {
        cerr << "Error opening slow5/blow5 file" <<endl;
        return 1;
    }

    // load blow5/slow5 index
    int ret = 0;
    ret = slow5_idx_load(slow5File);
    if (ret<0){
        cerr << "Error loading slow5/blow5 index" <<endl;
        return 1;
    }


    bam_hdr_t* header = sam_hdr_read(samFile); // Read SAM header

    if (!header) {
        cerr << "Error reading SAM header." << endl;
        sam_close(samFile);
        return 1;
    }

    bam1_t* samRecord = bam_init1(); // Initialize a SAM record

    int c_reads = 0;
    int c_parses = 0;

    // read kmer model map
    unordered_map<string, float> kmer_map = readModelKmer(model_kmer_path);
    if (kmer_map.empty()){
        cerr << "Error loading kmer file." << endl;
        return 1;
    }

    unordered_map<string, int> targets_map = get_target_umap(targets_path);
    unordered_map<string, int> preprocessed_targets_count;


    if (targets_map.empty()){
        cerr << "Error loading targets file." << endl;
        return 1;
    }


    while (sam_read1(samFile, header, samRecord) >= 0) {
        // Access the read name
        string readName = bam_get_qname(samRecord);

        // Access the contig (reference sequence) name
        const char* contigName = header->target_name[samRecord->core.tid];
        string contigNameString = contigName;

        // Check if contig is in target contigs, skip to next if not
        // if (targets_map.count(contigNameString)){
        //     cout << "processing " << contigNameString << endl;
        // } else { continue;}

        // unordered_set<int> target_pos_set = targets_map[contigNameString];
        
        
        int contigNameLength = strlen(contigName);



        // cout << "HERE   " << contigName[0] << endl;

        // Check if the read is aligned
        if (!(samRecord->core.flag & BAM_FUNMAP)) {

            c_reads++;
            unordered_map<int, int> qscore_umap;
            unordered_map<int, char> base_umap;
            unordered_map<int, int> readpos_umap;
            
            // Access the CIGAR string
            uint32_t* cigar = bam_get_cigar(samRecord);

            // Access the aligned bases and quality scores
            uint8_t* seq = bam_get_seq(samRecord);
            uint8_t* qual = bam_get_qual(samRecord);

            int readPosition = 0;
            int refPosition = samRecord->core.pos;

            for (int i = 0; i < samRecord->core.n_cigar; i++) {
                int cigarOp = bam_cigar_op(cigar[i]);
                int cigarLen = bam_cigar_oplen(cigar[i]);
                
                switch (cigarOp) {
                    case BAM_CMATCH:
                        for (int j = 0; j < cigarLen; j++) {
                            char base = seq_nt16_str[bam_seqi(seq, readPosition)]; // Get aligned base at this position
                            char quality = qual[readPosition]; // Get quality score at this position
                            int qscore = static_cast<int>(quality);

                            base_umap[refPosition] = base;
                            qscore_umap[refPosition] = qscore;
                            readpos_umap[refPosition] = readPosition;

                            readPosition++;
                            refPosition++;
                        }
                        break;

                    case BAM_CINS:
                        // Handle insertions
                        readPosition += cigarLen;
                        break;

                    case BAM_CDEL:
                        // Handle deletions
                        for (int j=0; j< cigarLen; j++){
                            qscore_umap[refPosition]=0;
                            refPosition++;

                        }
                        // refPosition += cigarLen;
                        break;

                    case BAM_CSOFT_CLIP:
                        // Handle soft-clipped bases
                        readPosition += cigarLen;
                        break;

                    case BAM_CREF_SKIP:
                        // Handle skip on the reference (e.g. spliced alignment)
                        for (int j=0; j< cigarLen; j++){
                            qscore_umap[refPosition]=0;
                            refPosition++;

                        }
                        // refPosition += cigarLen;
                        break;

                    case BAM_CHARD_CLIP:
                        // Handle hard-clipped bases
                        break;

                    default:
                        // Print if the operation is not handled
                        cout << "Operation not handled " << cigarOp << "  " << bam_cigar_opchr(cigarOp) << endl;
                }
            }
            // return 0;

            // Access the si tag
            uint8_t* siTag = bam_aux_get(samRecord, "si");

            // Access the ss tag
            uint8_t* ssTag = bam_aux_get(samRecord, "ss");

            // Access the sc tag
            uint8_t* scTag = bam_aux_get(samRecord, "sc");

            // Access the sh tag
            uint8_t* shTag = bam_aux_get(samRecord, "sh");


            if (siTag && ssTag && scTag && shTag) {
                //get text for si tag
                string siText = bam_aux2Z(siTag);

                //get text for ss tag
                string ssText = bam_aux2Z(ssTag);

                //get value for sc tag
                float scValue = bam_aux2f(scTag);

                //get value for sh tag
                float shValue = bam_aux2f(shTag);

                // Get  start_raw, end_raw, start_kmer, end_kmer from siText
                vector<string> siVector = split_string(siText, ',');
                int start_raw = stoi(siVector[0]);
                int end_raw = stoi(siVector[1]);
                int end_kmer = stoi(siVector[2]);   // end > start kmer as RNA is sequenced 3' -> 5'
                int start_kmer = stoi(siVector[3]);

                int cur_kmer_i = start_kmer;
                int cur_raw_i = start_raw;

                char buff_array[11];
                int buff_i = 0;
                int st_raw_arr[end_kmer+1];
                int end_raw_arr[end_kmer+1];

                for (int i=0; i< end_kmer+1; i++) {st_raw_arr[i] = end_raw_arr[i] = -1;}

                for (char ssChar : ssText){
                    if (ssChar == ',' || ssChar == 'I' || ssChar == 'D'){
                        if (buff_i < 0){
                            cerr << "Bad ss in sam: Preceding digit missing" << endl;
                            return 1;
                            } 

                        buff_array[buff_i]=0;      // null terminate buff
                        int num = atoi(buff_array); // get number preceeding operation
                        
                        if (num < 0) {
                            cerr << "Bad ss in sam: Cannot have a negative number" << endl;
                            return 1;  
                        }

                        buff_i = 0; buff_array[0] = 0;   // reset buff

                        if (ssChar == ','){

                            st_raw_arr[cur_kmer_i] = cur_raw_i;    // store end of old kmer event
                            cur_raw_i+=num;                         // increment start of new kmer event
                            end_raw_arr[cur_kmer_i] = cur_raw_i;     // store start of new kmer event  
                            cur_kmer_i++;                           // increment kmer 

                        } else if (ssChar == 'I'){
                            cur_raw_i+=num;
                        } else if (ssChar == 'D'){
                            cur_kmer_i+=num;
                        } 
                    } else {

                        if (!isdigit(ssChar)){
                            cout << "Bad ss in sam:" << ssChar << "not a digit or valid ss operation. " << endl;
                            cerr << "Bad ss in sam." << endl;
                            return 1;
                        }
                        buff_array[buff_i++] = ssChar;         
                    }
                }

                /// load the raw signals from fast5
                slow5_rec_t* rec=nullptr;
                int ret = 0;

                ret = slow5_get(readName.c_str(), &rec, slow5File);
                int readID_len; char* readID;
                if (ret < 0){
                    cerr << "Error fetching the read from slow5 in read " << readName << endl;
                    continue;
                } else {
                    readID = rec->read_id;
                    readID_len = rec->read_id_len;

                    // cout << rec->read_id << " " << c_reads << endl;
                    uint64_t len_raw_signal = rec->len_raw_signal;
                }
                int c_test = 0;
                int consecutive = 0;
                //start empty 9mer array-queue
                kmerStruct ninemer[5];
                //set the 5mer struct
                kmerStruct kmer;
                //start empty qscore array-queue
                float ninemer_qscores[9];
                float q_sum = 0;
                int q_count = 0;

                int refPos = start_kmer;

                for (int kmer_index = end_kmer-1; kmer_index >= start_kmer; kmer_index--){  // loop through 5mer events in the read
                    
                    int start_idx = st_raw_arr[kmer_index];
                    int end_idx = end_raw_arr[kmer_index];

                    if (start_idx!=-1){    // If event aligned to kmer
                        consecutive++;

                        if (consecutive==5){ // If first 9-mer built

                            // Fill 9mer = push event start and end
                            kmer.eventIndex = {start_idx,end_idx};
                            kmer.dwell = end_idx - start_idx;
                            ninemer[consecutive-1] = kmer;

                            // Build 9mer
                            for (int kmer_i=0; kmer_i < 5; kmer_i++){  //loop through 5 5mers
                                int event_start = ninemer[kmer_i].eventIndex[0];
                                int event_end = ninemer[kmer_i].eventIndex[1];
                                int dwell_time = event_end-event_start;
                                float event_arr[dwell_time];
                                int arr_index = 0;   
                   
                                // get squiggles from slow5
                                for (int raw_slow5_i = event_end -1; raw_slow5_i >= event_start; raw_slow5_i--){

                                    float pA = (TO_PICOAMPS(rec->raw_signal[raw_slow5_i],rec->digitisation,rec->offset,rec->range) -shValue) / scValue;
                                    event_arr[arr_index++]=pA;

                                // interpolate squiggles into 36 vector len

                                }
                                float* smooth_signals=nullptr;
                                if (dwell_time < 36){
                                    smooth_signals = interpolateSignals(event_arr, dwell_time);
                                
                                } else if (dwell_time >= 36){
                                    smooth_signals = compresSignal(event_arr, dwell_time);
                                }

                                ninemer[kmer_i].smoothCurrents = smooth_signals;

                            }
                            
                            // get 4DDQ
                            // get 9-mer sequence from fasta
                            int refBaseLen =9;
                            int ninemer_pos = refPos-4;
                            char* refBases = faidx_fetch_seq(fai, contigName, ninemer_pos, refPos+4, &refBaseLen);

                            q_sum = 0;
                            q_count = 0;
                            for (int i=0; i<9; i++){
        
                                char refBase = refBases[i];
                                char readBase = base_umap[ninemer_pos];
                                float qscore = qscore_umap[ninemer_pos];
                                if (qscore == 0){
                                    ninemer_qscores[i] = qscore;

                                } else if (refBase == readBase){
                                    ninemer_qscores[i] = qscore;
                                    q_sum += qscore;
                                    q_count++;
                                } else {
                                    float p = pow(10, (-1 * qscore)/10);
                                    float qadj = -10 * log10(1-p);
                                    ninemer_qscores[i] = qadj;
                                    q_sum += qadj;
                                    q_count++;
                                }
                                ninemer_pos++;
                            }

                            if (q_count > 0){
                                if (refBases[4] == INPUT_BASE){
                                    string ninemer_str = refBases;
                                    if (targets_map.count(ninemer_str)){
                                        if (preprocessed_targets_count.find(ninemer_str) == preprocessed_targets_count.end() ){
                                            preprocessed_targets_count[ninemer_str]=1;
                                        } else {
                                            preprocessed_targets_count[ninemer_str]++;
                                        }
                                        if (preprocessed_targets_count[ninemer_str]<=targets_map[ninemer_str]){
                                            float out_q[9];
                                            float avg_q = q_sum / static_cast<float> (q_count);
                                            for (int i=0; i<9; i++){
                                                float qscore = ninemer_qscores[i];
                                                if (qscore == 0){
                                                    out_q[i]=avg_q;
                                                } else {
                                                    out_q[i] = qscore;
                                                }
                                            }
                                            
                                            array<string,4> comp_ninemers;
                                            for (int i=0; i<4; i++){
                                                comp_ninemers[i] = ninemer_str.substr(0,4) + BASES[i] + ninemer_str.substr(5,10);
                                            }
                                            // Output 9mer
                                            
                                            // Write the integer len(string) followed by string in little-endian byte order
                                            
                                            // Write contig
                                            outFile.write(reinterpret_cast<const char*>(&contigNameLength), sizeof(int));
                                            outFile.write(contigName,contigNameLength);
                                            
                                            // Write position
                                            int ninemer_position = refPos-4;
                                            outFile.write(reinterpret_cast<const char*>(&ninemer_position), sizeof(int));

                                            // Write 9mer 
                                            outFile.write(refBases,9);

                                            //Write readID
                                            outFile.write(reinterpret_cast<const char*>(&readID_len), sizeof(int));
                                            outFile.write(readID, readID_len);

                                            //Write read position
                                            int readPosition = readpos_umap[refPos];
                                            outFile.write(reinterpret_cast<const char*>(&readPosition), sizeof(int));

                                            //Write qscore called in the read at target position. 0 if not called
                                            int qscoreTarget = qscore_umap[refPos];
                                            outFile.write(reinterpret_cast<const char*>(&qscoreTarget), sizeof(int));
                                            //Write Base called in the read at target position
                                            if (qscoreTarget != 0){
                                                    char targetBase = base_umap[refPos];
                                                    outFile.write(reinterpret_cast<const char*>(&targetBase), sizeof(char));
                                            } else {
                                                    char targetBase = 'D';
                                                    outFile.write(reinterpret_cast<const char*>(&targetBase), sizeof(char));
                                            }



                                            
                                            for (int kmer_i=0; kmer_i < 5; kmer_i++){
                                                string fivemer_seq = ninemer_str.substr(kmer_i, 5);
                                                array<float,4> comp_fivemers_curr;


                                                for (int i=0; i<4; i++){
                                                    comp_fivemers_curr[i]= kmer_map[comp_ninemers[i].substr(kmer_i,5)];
                                                //    if (comp_fivemers_curr[i] < 1){ cout << comp_ninemers[i] << " here " << ninemer_str << " " <<  ninemer_str.substr(5,10)  << endl; }
                                                //    else {cout << "GOOD" << endl;}
                                                }

                                                int current_out_start = 36 * kmer_i;
                                                for (int j =0; j<36; j++){

                                                    int arr_row = current_out_start + j;
                                                    int qscore_index = arr_row / 20;
                                                    float current_val = ninemer[kmer_i].smoothCurrents[j];
                                                    outFile.write(reinterpret_cast<const char*>(&current_val), sizeof(float));

                                                    float value = comp_fivemers_curr[0] - current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[1]- current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[2]- current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[3]- current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = ninemer[kmer_i].dwell;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = out_q[qscore_index] ;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));


                                                }
                                            }
                                            c_parses += 1;

                                        }
                                    }
                                }
                            }
                            free(refBases);
                            
                             
                        } else if (consecutive>5){
                            // Slide 9mer = Get new event from slow5, interpolate new into 36

                            int dwell_time = end_idx-start_idx;
                            float event_arr[dwell_time];
                            int arr_index = 0;                        
                            // get squiggles from slow5
                            for (int raw_slow5_i = end_idx -1; raw_slow5_i >= start_idx; raw_slow5_i--){
                                float pA = (TO_PICOAMPS(rec->raw_signal[raw_slow5_i],rec->digitisation,rec->offset,rec->range) -shValue) / scValue;
                                event_arr[arr_index++]=pA;
                            }
                            // smoothen the raw signals
                            float* smooth_signals=nullptr;
                            if (dwell_time < 36){
                                smooth_signals = interpolateSignals(event_arr, dwell_time);
                            
                            } else if (dwell_time >= 36){
                                smooth_signals = compresSignal(event_arr, dwell_time);
                            } 
                            // get updated array after sliding 9mer start
                            int kmer_sliding_start = consecutive %5;
                            int kmer_sliding_end = (kmer_sliding_start+4)%5;

                            kmer.eventIndex = {start_idx,end_idx};
                            kmer.smoothCurrents = smooth_signals;
                            kmer.dwell = dwell_time;
                            delete [] ninemer[kmer_sliding_end].smoothCurrents;
                            ninemer[kmer_sliding_end] = kmer;

                            // get 4DDQ
                            // get 9-mer sequence from fasta
                            int refBaseLen =9;
                            int ninemer_pos = refPos-4;
                            char* refBases = faidx_fetch_seq(fai, contigName, ninemer_pos, refPos+4, &refBaseLen);
                                                       
                            
                            // get new qscore 
                            char refBase = refBases[8];
                            char readBase = base_umap[refPos+4];
                            float qscore = qscore_umap[refPos+4];

                            // get new qscore index
                            int q_start_index = (consecutive +4)%9;
                            int q_end_index = (consecutive+3)%9;

                             
                            // slide qscores
                            float old_q = ninemer_qscores[q_end_index];

                            

                            if (old_q != 0){
                                q_count--;
                                q_sum -= old_q;
                                
                            }

                            if (qscore == 0){
                                ninemer_qscores[q_end_index] = qscore;

                            } else if (refBase == readBase){
                                ninemer_qscores[q_end_index] = qscore;
                                q_sum += qscore;
                                q_count++;
                            } else {
                                
                                float p = pow(10, (-1 * qscore)/10);
                                float qadj = -10 * log10(1-p);
                                ninemer_qscores[q_end_index] = qadj;
                                q_sum += qadj;
                                q_count++;
                                
                            }


                            if (q_count > 0){
                                if (refBases[4] == INPUT_BASE){
                                    string ninemer_str = refBases;
                                    if (targets_map.count(ninemer_str)){
                                        if (preprocessed_targets_count.find(ninemer_str) == preprocessed_targets_count.end() ){
                                            preprocessed_targets_count[ninemer_str]=1;
                                        } else {
                                            preprocessed_targets_count[ninemer_str]++;
                                        }
                                        if (preprocessed_targets_count[ninemer_str]<=targets_map[ninemer_str]){
                                        
                                            float out_q[9];
                                            float avg_q = q_sum / static_cast<float> (q_count);
                                            for (int i=0; i<9; i++){
                                                int q_index = (q_start_index +i) %9;
                                                float qscore = ninemer_qscores[q_index];
                                                if (qscore == 0){
                                                    out_q[i]=avg_q;
                                                } else {
                                                    out_q[i] = qscore;
                                                }
                                            }

                                            string ninemer_str = refBases;
                                            array<string,4> comp_ninemers;
                                            for (int i=0; i<4; i++){
                                                comp_ninemers[i] = ninemer_str.substr(0,4) + BASES[i] + ninemer_str.substr(5,4);
                                            }
                                            // Output 9mer
                                            
                                            // Write the integer len(string) followed by string in little-endian byte order
                                            
                                            // Write contig
                    
                                        
                                            outFile.write(reinterpret_cast<const char*>(&contigNameLength), sizeof(int));
                                            outFile.write(contigName,contigNameLength);
                                            
                                            // Write position
                                            int ninemer_position = refPos-4;
                                            outFile.write(reinterpret_cast<const char*>(&ninemer_position), sizeof(int));

                                            // Write 9mer 
                                            outFile.write(refBases,9);

                                            //Write readID
                                            outFile.write(reinterpret_cast<const char*>(&readID_len), sizeof(int));
                                            outFile.write(readID, readID_len);

                                            //Write read position
                                            int readPosition = readpos_umap[refPos];
                                            outFile.write(reinterpret_cast<const char*>(&readPosition), sizeof(int));


                                            //Write qscore called in the read at target position. 0 if not called
                                            int qscoreTarget = qscore_umap[refPos];
                                            outFile.write(reinterpret_cast<const char*>(&qscoreTarget), sizeof(int));
                                            //Write Base called in the read at target position
                                            if (qscoreTarget != 0){
                                                    char targetBase = base_umap[refPos];
                                                    outFile.write(reinterpret_cast<const char*>(&targetBase), sizeof(char));
                                            } else {
                                                    char targetBase = 'D';
                                                    outFile.write(reinterpret_cast<const char*>(&targetBase), sizeof(char));
                                            }



                                            for (int kmer_i=0; kmer_i < 5; kmer_i++){
                                                int kmer_arr_i = (kmer_i + kmer_sliding_start) % 5;
                                                string fivemer_seq = ninemer_str.substr(kmer_i, 5);
                                                array<float,4> comp_fivemers_curr;


                                                for (int i=0; i<4; i++){

                                                    comp_fivemers_curr[i]= kmer_map[comp_ninemers[i].substr(kmer_i,5)];

                                                }
                                                int current_out_start = 36 * kmer_i;
                                                for (int j =0; j<36; j++){

                                                    int arr_row = current_out_start + j;
                                                    int qscore_index = arr_row / 20;
                                                    float current_val = ninemer[kmer_arr_i].smoothCurrents[j];
                                                    outFile.write(reinterpret_cast<const char*>(&current_val), sizeof(float));

                                                    float value = comp_fivemers_curr[0] - current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[1] - current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[2] - current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = comp_fivemers_curr[3] - current_val;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = ninemer[kmer_arr_i].dwell;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                    value = out_q[qscore_index] ;
                                                    outFile.write(reinterpret_cast<const char*>(&value), sizeof(float));

                                                }
                                            }
                                            c_parses += 1;
                                        }
                                    }
                                }
                            }    
                            // } else { cout << "SKIP 0 Q" << endl;}
                            free(refBases);
                            
                            // if (c_parses % BATCH_SIZE == 0 && c_parses > PREV_PARSES){

                            //     FILE_NUM++;
                            //     PREV_PARSES = c_parses;
                            //     outFile.close();


                            //     // Check every 20 dumps if preprocess is too much ahead
                            //     if (FILE_NUM % 20 == 0 && FILE_NUM > 0){
                            //         // Look if 10 files pilled up, wait for predict to catch up
                            //         string PREP_FILE_OLD = getFileNameFromPath(name_out) + "_" + to_string(FILE_NUM-10);
                            //         int proceed = 0;
                            //         int waited = 0;
                            //         while (proceed == 0){
                            //             proceed = 1;
                            //             dir = opendir(directoryPath);
                            //             if (dir != NULL) {
                            //                 while ((entry = readdir(dir))) {
                            //                     if (entry->d_type == DT_REG) { // Regular file
                            //                         if (entry->d_name ==  PREP_FILE_OLD){
                            //                             proceed = 0;
                            //                         }
                            //                     }
                            //                 }

                            //                 if (proceed == 0){
                            //                     cout << "WAITING FOR " << PREP_FILE_OLD << " TO BE REMOVED, WAITED " << waited << endl;
                            //                     sleep(1);
                            //                     waited++;
                            //                 }

                            //                 closedir(dir);
                            //             } else {
                            //                 std::cerr << "Failed to open the directory." << std::endl;
                            //             }
                            //         }                                   
                            //     }
                            //     PREP_FILE = name_out + "_" + to_string(FILE_NUM);
                                

                            //     //open new file for new batch
                            //     outFile.open(PREP_FILE, std::ios::binary | std::ios::out);
                            //     if (!outFile) {
                            //         std::cerr << "Error opening file for writing." << std::endl;
                            //         return 1;
                            //     }
                            // }
                        } else {
                            // Fill 9mer = push event start and end
                            kmer.eventIndex = {start_idx,end_idx};
                            kmer.dwell = end_idx - start_idx;
                            ninemer[consecutive-1] = kmer;
                        }

                    } else {
                        // reset 9mer
                        //Free up memory for the currents pointers from the last parsed 9mer
                        if (consecutive > 4){
                            if (consecutive > 5){consecutive=5;}
                            for (int i=0; i<consecutive; i++){
                                delete [] ninemer[i].smoothCurrents;
                            }
                        }
                        
                        consecutive=0;
                    }
                    c_test++;
                    refPos++;                    
                }
            slow5_rec_free(rec);
            //Free up memory for the currents pointers from the last parsed 9mer
            if (consecutive > 4){
                if (consecutive > 5){consecutive=5;}
                for (int i=0; i<consecutive; i++){
                    delete [] ninemer[i].smoothCurrents;
                }
            }
            }
            else {
                cout << "One of -si -ss -sc -sh tags not in the sam entry " << readName << endl;
                cout << "Skipping" << endl;
            }

            // Parse only 10k reads during testing
        //    if (c_reads > 9) {
          //      cout << c_reads << " reads and " << c_parses << " bases parsed" << endl;
            //    break;
             //}
        
        }
    }

    // Clean up
    outFile.close();

    // PREP_FILE = name_out + "_DONE";

    // outFile.open(PREP_FILE, std::ios::binary | std::ios::out);
    //                             if (!outFile) {
    //                                 std::cerr << "Error opening file for writing." << std::endl;
    //                                 return 1;
    //                             }
    // outFile.close();
    bam_destroy1(samRecord);
    sam_hdr_destroy(header);
    sam_close(samFile);
    fai_destroy(fai);
    slow5_idx_unload(slow5File);
    slow5_close(slow5File);

    cout << c_reads << " reads and " << c_parses << " bases parsed" << endl;

    return 0;
}

