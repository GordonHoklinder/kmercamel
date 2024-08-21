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
#include "masks.h"

int Help() {
    std::cerr << "KmerCamel version " << VERSION << std::endl;
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algorithm     - the algorithm to be run [global (default), globalAC, local, localAC, streaming]" << std::endl;
    std::cerr << "  -o output_path   - if not specified, the output is printed to stdout" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max; default 5" << std::endl;
    std::cerr << "  -c               - treat k-mer and its reverse complement as equal" << std::endl;
    std::cerr << "  -m               - turn off the memory optimizations for global" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "  -v               - print version" << std::endl;
    std::cerr << "Example usage:       ./kmercamel -p path_to_fasta -k 13 -d 5 -a local" << std::endl;
    std::cerr << "Possible algorithms: global globalAC local localAC streaming" << std::endl;
    std::cerr << std::endl;
    std::cerr << "For optimization of masks use `kmercamel optimize`."  << std::endl;
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algorithm     - the algorithm to be run [ones (default), runs, runsapprox, zeros]" << std::endl;
    std::cerr << "  -o output_path   - if not specified, the output is printed to stdout" << std::endl;
    std::cerr << "  -c               - treat k-mer and its reverse complement as equal" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "  -v               - print version" << std::endl;
    return 1;
}

constexpr int MAX_K = 63;

void Version() {
    std::cerr << VERSION << std::endl;
}

// This macro cannot be avoided without bloating the code too much.
#define INIT_KMERCAMEL(type) \
int kmercamel##type(std::string path, int k, int d_max, std::ostream *of, bool complements, bool masks, std::string algorithm, bool optimize_memory) { \
    kmer_dict##type##_t wrapper;                                                                                                                       \
    kmer##type##_t kmer_type; \
    if (masks) { \
        int ret = Optimize(wrapper, kmer_type, algorithm, path, *of, k, complements); \
        if (ret) Help(); \
        return ret; \
    } \
    \
    /* Handle streaming algorithm separately. */ \
    if (algorithm == "streaming") { \
        WriteName(k, *of); \
        Streaming(path, *of,  k , complements); \
    } \
    /* Handle hash table based algorithms separately so that they consume less memory. */ \
    else if (algorithm == "global" || algorithm == "local") { \
        auto *kMers = kh_init_S##type(); \
        ReadKMers(kMers, wrapper, kmer_type, path, k, complements); \
        if (!kh_size(kMers)) { \
            std::cerr << "Path '" << path << "' contains no k-mers." << std::endl; \
            return Help(); \
        } \
        d_max = std::min(k - 1, d_max); \
        WriteName(k, *of); \
        if (algorithm == "global") { \
            auto kMerVec = kMersToVec(kMers, kmer_type); \
            kh_destroy_S##type(kMers); \
            /* Turn off the memory optimizations if optimize_memory is set to false. */ \
            if(optimize_memory) PartialPreSort(kMerVec, k); \
            else MEMORY_REDUCTION_FACTOR = 1; \
            Global(wrapper, kMerVec, *of, k, complements); \
        } \
        else Local(kMers, wrapper, kmer_type, *of, k, d_max, complements); \
    } else { \
        auto data = ReadFasta(path); \
        if (data.empty()) { \
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl; \
            return Help(); \
        } \
        d_max = std::min(k - 1, d_max); \
        \
        auto kMers = ConstructKMers(data, k, complements); \
        WriteName(k, *of); \
        if (algorithm == "globalAC") { \
            GlobalAC(kMers, *of, complements); \
        } \
        else if (algorithm == "localAC") { \
            LocalAC(kMers, *of, k, d_max, complements); \
        } \
        else { \
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl; \
            return Help(); \
        } \
    } \
    *of << std::endl; \
    return 0; \
}

INIT_KMERCAMEL(64)
INIT_KMERCAMEL(128)

int main(int argc, char **argv) {
    std::string path;
    int k = 0;
    int d_max = 5;
    std::ofstream output;
    std::ostream *of = &std::cout;
    bool masks = false;
    std::string algorithm = "global";
    if (argc > 1 && std::string(argv[1]) == "optimize") {
        masks = true;
        argv++;
        argc--;
        algorithm = "ones";
    }
    bool complements = false;
    bool optimize_memory = true;
    bool d_set = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:o:hcvm"))  != -1) {
            switch(opt) {
                case  'p':
                    path = optarg;
                    break;
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
                    // Backwards compatability.
                    if (algorithm == "greedy") algorithm = "global";
                    if (algorithm == "greedyAC") algorithm = "globalAC";
                    if (algorithm == "pseudosimplitigs") algorithm = "local";
                    if (algorithm == "pseudosimplitigsAC") algorithm = "localAC";
                    break;
                case  'c':
                    complements = true;
                    break;
                case 'm':
                    optimize_memory = false;
                    break;
                case 'v':
                    Version();
                    return 0;
                case 'h':
                default:
                    Help();
                    return 0;
            }
        }
    } catch (std::invalid_argument&) {
        return Help();
    }
    if (path.empty()) {
        std::cerr << "Required parameter p not set." << std::endl;
        return Help();
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return Help();
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return Help();
    } else if (d_max < 0) {
        std::cerr << "d must be non-negative." << std::endl;
        return Help();
    } else if (k > MAX_K && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > " << MAX_K << " not supported for the algorithm '" + algorithm + "'. Use the 128bit version of KmerCamel or the AC version of the algorithm instead." << std::endl;
        return Help();
    } else if (d_set && (algorithm == "globalAC" || algorithm == "global" || algorithm == "streaming")) {
        std::cerr << "Unsupported argument d for algorithm '" + algorithm + "'." << std::endl;
        return Help();
    } else if (!optimize_memory && algorithm != "global") {
        std::cerr << "Memory optimization turn-off only supported for hash table global." << std::endl;
        return Help();
    } else if (masks && (d_set || !optimize_memory)) {
        std::cerr << "Not supported flags for optimize." << std::endl;
        return Help();
    }
    if (k < 32) {
        return kmercamel64(path, k, d_max, of, complements, masks, algorithm, optimize_memory);
    } else {
        return kmercamel128(path, k, d_max, of, complements, masks, algorithm, optimize_memory);
    }
}
