#pragma once

#include "parser.h"
#include <iostream>

inline bool is_upper(char c) {
    return c <= 'Z';
}

inline char to_upper(char c) {
    constexpr int alphabet_shift = 'a' - 'A';
    bool upper = is_upper(c);
    return upper ? c : c - alphabet_shift;
}

void split_ms(std::ostream &superstring, std::ostream &mask, std::string path) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    int64_t l = kseq_read(seq);

    for (int64_t i = 0; i < l; ++i) {
        bool bitmask = is_upper(seq->seq.s[i]);
        mask << (bitmask ? '1' : '0');
        superstring << to_upper(seq->seq.s[i]);
    }
    mask << std::endl;
    superstring << std::endl;

    kseq_destroy(seq);
    gzclose(fp);
}

void join_ms(std::istream &superstringf, std::istream &maskf, std::ostream &of) {
    std::string superstring, mask;
    superstringf >> superstring;
    maskf >> mask;
    of << ">superstring" << std::endl;
    for (size_t i = 0; i < superstring.length(); ++i) {
        of << Masked(superstring[i], mask[i] == '1');
    }
    of << std::endl;
}

void ms_to_spss(std::string path, std::ostream &of, int k) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    int64_t l = kseq_read(seq);

    bool masked = false;
    int counter = 0;
    for (int64_t i = 0; i < l; ++i) {
        if (is_upper(seq->seq.s[i])) {
            if (!masked) {
                of << ">" << counter++ << std::endl;
            }
            masked = true;
            of << seq->seq.s[i];
        } else {
            if (!masked) continue;
            masked = false;
            for (int64_t j = 0; j < k - 1; ++j) {
                if (i + j < l) of << to_upper(seq->seq.s[i + j]);
            }
            of << std::endl;
        }
    }

    kseq_destroy(seq);
    gzclose(fp); 
}

void spss_to_ms(std::string path, std::ostream &of, int k) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    int64_t l = 0;

    of << ">superstring " << path << std::endl;

    while ((l = kseq_read(seq)) >= 0) {
        if (l < k) continue;
        for (int64_t i = 0; i < l; ++i) of << Masked(seq->seq.s[i], i <= l - k); 
    }

    of << std::endl;

    kseq_destroy(seq);
    gzclose(fp); 
}

