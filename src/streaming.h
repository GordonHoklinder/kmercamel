#pragma once
#include <unordered_set>
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>

#include "parser.h"
#include "kmers.h"
#include "khash_utils.h"

template <typename kmer_t, typename kh_wrapper_t>
void Streaming(kh_wrapper_t wrapper, kmer_t kmer_type, std::string &path, std::ostream &of, int k, bool complements) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    auto kMers = wrapper.kh_init_set();
    of << ">superstring" << std::endl;

    kmer_t mask = (kmer_t(1)) << (2 * k - 1);
    mask |= mask - 1;


    while (kseq_read(seq) >= 0) {
        kmer_t kMer = 0;
        size_t firstIndex = k - 1;
        int64_t lastOne = -k;
        for (size_t i = 0; i < seq->seq.l + k - 1; ++i) {
            kmer_t c = 4;
            if (i < seq->seq.l) c = nucleotideToInt[(uint8_t)seq->seq.s[i]];
            if (c >= 4) {
                kMer = 0;
                firstIndex = i + k;
            }
            kMer <<= 2;
            kMer |= c;
            kMer &= mask;
            if (i >= firstIndex && !containsKMer(kMers, wrapper, kMer, k, complements)) {
                int ret;
                wrapper.kh_put_to_set(kMers, kMer, &ret);
                of << Masked(seq->seq.s[i - k + 1], true);
                lastOne = i;
            } else if ((int64_t (i)) <= lastOne + k - 1) {
                of << Masked(seq->seq.s[i - k + 1], false);
            }
        }
    }

    of << std::endl;

    kseq_destroy(seq);
    gzclose(fp);
}