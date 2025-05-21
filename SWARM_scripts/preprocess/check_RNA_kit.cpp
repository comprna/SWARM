#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <cstring>
#include <dirent.h>
#include <unistd.h>
#include "argagg.hpp"
#include "slow5/slow5.h"

using namespace std;



int main(int argc, char *argv[]) {

    if (argc < 2)
    {
        std::cout << "Run with: ./SWARM_preprocess  --sam <events.sam> --fasta <ref.fa> --raw <signals.blow5> -o <outpath> --base <A/C/G/T> -m <model_kmer_path>" << std::endl;
        return EXIT_FAILURE;
    }

    argagg::parser argparser{{
        {"help", {"-h", "--help"}, "shows this help message\n", 0},
        {"slow5FileName", {"--raw"}, "Slow5/blow5 file. <sample.blow5> (required)\n", 1},
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

    if ( !args["slow5FileName"] )
    {
        cerr << "ERROR: Please provide all of required args: --sam --fasta --raw -m -o -i" << endl;
        cerr << argparser;
        return EXIT_FAILURE;
    }

    string slow5FileName_s = args["slow5FileName"].as<string>(""); const char* slow5FileName = slow5FileName_s.c_str();

    //cout << slow5FileName_s << endl;

    slow5_file_t* slow5File = slow5_open(slow5FileName, "r");


    const slow5_hdr_t* hdr = slow5File->header;

    // Check the first read group's sequencing kit
    char* kit = slow5_hdr_get("sequencing_kit", 0, hdr);
    if (kit == nullptr) {
        cerr << "ERROR: 'sequencing_kit' not found in SLOW5 header." << endl;
        slow5_close(slow5File);
        return EXIT_FAILURE;
    }

    cout << kit << endl;

    //bool is_rna004 = false;
    //if (strstr(kit, "rna004") != nullptr) {
    //    is_rna004 = true;
    //    cout << "Detected rna004 sequencing kit." << endl;
    //} else {
    //    cout << "Sequencing kit is not rna004." << endl;
    //}




    slow5_close(slow5File);

}
     
