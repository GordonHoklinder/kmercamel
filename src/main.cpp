#include <iostream>
#include <string>

#include "unistd.h"
#include "version.h"
#include "ac/global_ac.h"
#include "global.h"
#include "local.h"
#include "ac/local_ac.h"
#include "parser.h"
#include "ac/parser_ac.h"
#include "ac/streaming.h"
#include "khash_utils.h"
#include "conversions.h"

#include <iostream>
#include <string>
#include "unistd.h"
#include "version.h"
#include "masks.h"
#include "lower_bound.h"

int usage() {
    std::cerr << std::endl;
    std::cerr << "Program: KmerCamel - a tool for computation of masked superstrings." << std::endl;
    std::cerr << "Version: " << VERSION << std::endl;
    std::cerr << "Contact: Ondrej Sladky (ondra.sladky@gmail.com)" << std::endl << std::endl;
    std::cerr << "Usage:   kmercamel <command> [options]" << std::endl << std::endl;
    std::cerr << "Command:" << std::endl;
    std::cerr << "    ms         - Compute an approximately shortest masked superstring from input FASTA." << std::endl;
    std::cerr << "    optimize   - Optimize a given masked superstring." << std::endl;
    std::cerr << "    ms2msfa    - Split a masked superstring into a superstring and a mask." << std::endl;
    std::cerr << "    msfa2ms    - Join a masked superstring from a superstring and a mask." << std::endl;
    std::cerr << "    msfa2spss  - Compute rSPSS from a masked superstring." << std::endl;
    std::cerr << "    spss2msfa  - Compute masked superstring corresponding to (r)SPSS." << std::endl;
    std::cerr << "    lowerbound - Compute the lower bound on masked superstring size of a k-mer set." << std::endl;
    std::cerr << std::endl;
    return 1;
}

int usage_subcommand(std::string subcommand) {
    std::cerr << std::endl;
    std::cerr << "Usage:   kmercamel " << subcommand << " [options]";
    if (subcommand != "ms2msfa")
    std::cerr << " <fasta>";

    std::cerr  << std::endl << std::endl;
    std::cerr << "Options:" << std::endl;
    if (subcommand != "ms2msfa" && subcommand != "msfa2ms")
    std::cerr << "  -k INT   - k-mer size [required; up to 127]" << std::endl;

    if (subcommand == "ms")
    std::cerr << "  -a STR   - the algorithm to be run [global (default), globalAC, local, localAC, streaming]" << std::endl;

    else if (subcommand == "optimize")
    std::cerr << "  -a STR   - the algorithm to be run [maxone (default), minone, minruns, approxminruns]" << std::endl;

    if (subcommand != "lowerbound" && subcommand != "msfa2ms")
    std::cerr << "  -o FILE  - output, if not specified, the output is printed to stdout" << std::endl;
    
    if (subcommand == "ms")
    std::cerr << "  -M FILE  - if specified, print also mask maximizing ones to a separate file (possible only with global)" << std::endl;

    if (subcommand == "ms")
    std::cerr << "  -d INT   - d_max for local algorithm; default 5" << std::endl;
    if (subcommand == "ms" || subcommand == "optimize" || subcommand == "lowerbound")
    std::cerr << "  -u       - treat k-mer and its reverse complement as distinct" << std::endl;

    if (subcommand == "ms" || subcommand == "lowerbound")
    std::cerr << "  -x       - turn off the memory optimizations for global algorithm" << std::endl;

    if (subcommand == "ms2msfa" || subcommand == "msfa2ms") {
    std::cerr << "  -m FILE  - file with mask" << std::endl;
    std::cerr << "  -s FILE  - file with superstring" << std::endl;
    }
    std::cerr << "  -h       - print help" << std::endl;
    std::cerr << std::endl;
    return 1;
}

constexpr int MAX_K = 127;

void Version() {
    std::cerr << VERSION << std::endl;
}

/// Run KmerCamel with the given parameters.
template <typename kmer_t, typename kh_wrapper_t>
int kmercamel(kh_wrapper_t wrapper, kmer_t kmer_type, std::string path, int k, int d_max, std::ostream *of, std::ostream *maskf, bool complements, bool masks,
                    std::string algorithm, bool optimize_memory, bool lower_bound) {
    if (masks) {
        int ret = Optimize(wrapper, kmer_type, algorithm, path, *of, k, complements);
        if (ret) usage_subcommand("optimize");
        return ret;
    }

    /* Handle streaming algorithm separately. */
    if (algorithm == "streaming") {
        WriteName(k, *of);
        Streaming(path, *of,  k , complements);
    }
    /* Handle hash table based algorithms separately so that they consume less memory. */
    else if (algorithm == "global" || algorithm == "local") {
        auto *kMers = wrapper.kh_init_set();
        ReadKMers(kMers, wrapper, kmer_type, path, k, complements);
        if (!kh_size(kMers)) {
            std::cerr << "Path '" << path << "' contains no k-mers." << std::endl;
            return usage_subcommand("ms");
        }
        d_max = std::min(k - 1, d_max);
        if (!lower_bound) WriteName(k, *of);
        if (algorithm == "global") {
            auto kMerVec = kMersToVec(kMers, kmer_type);
            wrapper.kh_destroy_set(kMers);
            /* Turn off the memory optimizations if optimize_memory is set to false. */
            if(optimize_memory) PartialPreSort(kMerVec, k);
            else MEMORY_REDUCTION_FACTOR = 1;
            if (lower_bound) std::cout << LowerBoundLength(wrapper, kMerVec, k, complements);
            else Global(wrapper, kMerVec, *of, maskf, k, complements);
        }
        else Local(kMers, wrapper, kmer_type, *of, k, d_max, complements);
    } else {
        auto data = ReadFasta(path);
        if (data.empty()) {
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
            return usage_subcommand("ms");
        }
        d_max = std::min(k - 1, d_max);

        auto kMers = ConstructKMers(data, k, complements);
        WriteName(k, *of);
        if (algorithm == "globalAC") {
            GlobalAC(kMers, *of, complements);
        }
        else if (algorithm == "localAC") {
            LocalAC(kMers, *of, k, d_max, complements);
        }
        else {
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl;
            return usage_subcommand("ms");
        }
    }
    *of << std::endl;
    return 0;
}

int camel_compute(int argc, char **argv) {
    std::string subcommand = "ms";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    int d_max = 5;
    std::ofstream output;
    std::ostream *of = &std::cout;
    std::ofstream maskOutput;
    std::ostream *maskf = nullptr;
    std::string algorithm = "global";
    bool complements = true;
    bool optimize_memory = true;
    bool d_set = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:d:a:o:huxM:"))  != -1) {
            switch(opt) {
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  'd':
                    d_set = true;
                    d_max = std::stoi(optarg);
                    break;
                case  'a':
                    algorithm = optarg;
                    break;
                case  'u':
                    complements = false;
                    break;
                case 'x':
                    optimize_memory = false;
                    break;
                case 'M':
                    maskOutput.open(optarg);
                    maskf = &maskOutput;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    } else if (d_max < 0) {
        std::cerr << "d must be non-negative." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k > MAX_K && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > " << MAX_K << " not supported for the algorithm '" + algorithm + "'. Use the  AC version of the algorithm instead." << std::endl;
        return usage_subcommand(subcommand);
    } else if (d_set && (algorithm == "globalAC" || algorithm == "global" || algorithm == "streaming")) {
        std::cerr << "Unsupported argument d for algorithm '" + algorithm + "'." << std::endl;
        return usage_subcommand(subcommand);
    } else if (!optimize_memory && algorithm != "global") {
        std::cerr << "Memory optimization turn-off only supported for hash table global." << std::endl;
        return usage_subcommand(subcommand);
    } else if (maskf != nullptr && algorithm != "global") {
        std::cerr << "Outputting mask maximizing number of ones is possible only with global. For other algorithms, use separately kmercamel optimize." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, d_max, of, maskf, complements, false, algorithm, optimize_memory, false);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, d_max, of, maskf, complements, false, algorithm, optimize_memory, false);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, d_max, of, maskf, complements, false, algorithm, optimize_memory, false);
    }
}

int camel_optimize(int argc, char **argv) {
    std::string subcommand = "optimize";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    std::ofstream output;
    std::ostream *of = &std::cout;
    std::string algorithm = "maxone";
    bool complements = true;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:a:o:hu"))  != -1) {
            switch(opt) {
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  'a':
                    algorithm = optarg;
                    break;
                case  'u':
                    complements = false;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k > MAX_K && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > " << MAX_K << " not supported for the algorithm '" + algorithm + "'. Use the  AC version of the algorithm instead." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    }
}

int camel_lowerbound(int argc, char **argv) {
    std::string subcommand = "lowerbound";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    std::ostream *of = &std::cout;
    bool complements = true;
    bool optimize_memory = true;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:hu"))  != -1) {
            switch(opt) {
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case 'u':
                    complements = false;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                case 'm':
                    optimize_memory = false;
                    break;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, 0, of, nullptr, complements, false, "global", optimize_memory, true);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, 0, of, nullptr, complements, false, "global", optimize_memory, true);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, 0, of, nullptr, complements, false, "global", optimize_memory, true);
    }
}

int camel_split_ms(int argc, char **argv) {
    std::string subcommand = "msfa2ms";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    std::ofstream mask;
    std::ostream *maskf = &std::cout;

    std::ofstream superstring;
    std::ostream *superstringf = &std::cout;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "m:s:h"))  != -1) {
            switch(opt) {
                case  'm':
                    mask.open(optarg);
                    maskf = &mask;
                    break;
                case 's':
                    superstring.open(optarg);
                    superstringf = &superstring;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (!mask && !superstring) {
        std::cerr << "Cannot have both superstring and mask redirected to stdout." << std::endl;
        return usage_subcommand(subcommand);
    }
    split_ms(*superstringf, *maskf, path);
    return 0;    
}

int camel_join_ms(int argc, char **argv) {
    std::string subcommand = "ms2msfa";

    std::ifstream mask;
    std::istream *maskf = &std::cin;

    std::ifstream superstring;
    std::istream *superstringf = &std::cin;

    std::ofstream output;
    std::ostream *of = &std::cout;
    
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "m:s:o:h"))  != -1) {
            switch(opt) {
                case  'm':
                    mask.open(optarg);
                    maskf = &mask;
                    break;
                case 's':
                    superstring.open(optarg);
                    superstringf = &superstring;
                    break;
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (!mask && !superstring) {
        std::cerr << "Cannot have both superstring and mask redirected from stdin." << std::endl;
        return usage_subcommand(subcommand);
    }
    join_ms(*superstringf, *maskf, *of);
    return 0;    
}

int camel_ms_to_spss(int argc, char **argv) {
    std::string subcommand = "msfa2spss";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    std::ofstream output;
    std::ostream *of = &std::cout;

    int k = 0;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "o:k:h"))  != -1) {
            switch(opt) {
                case  'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'k':
                    k = std::stoi(optarg);
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    ms_to_spss(path, *of, k);
    return 0;
}

int camel_spss_to_ms(int argc, char **argv) {
    std::string subcommand = "spss2msfa";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    std::ofstream output;
    std::ostream *of = &std::cout;

    int k = 0;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "o:k:h"))  != -1) {
            switch(opt) {
                case  'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'k':
                    k = std::stoi(optarg);
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    spss_to_ms(path, *of, k);
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 1 || std::string(argv[1]) == "-h") {
        usage();
        return 0;
    } else if (std::string(argv[1]) == "ms") {
        return camel_compute(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "optimize") {
        return camel_optimize(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "lowerbound") {
        return camel_lowerbound(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "ms2msfa") {
        return camel_join_ms(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "msfa2ms") {
        return camel_split_ms(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "msfa2spss") {
        return camel_ms_to_spss(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "spss2msfa") {
        return camel_spss_to_ms(argc - 1, argv + 1);
    } else {
        std::cerr << "Command '" << argv[1] << "' not recognized." << std::endl;
        return usage(); 
    }
}
