#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "models.h"
#include "kmers.h"
#include "khash_utils.h"

/// Record of one fasta sequence.
struct FastaRecord {
    // The name of the record. Starts with '>'.
    std::string name;
    // The genetic sequence.
    std::string sequence;
};

/// Read fasta file with given path.
std::vector<FastaRecord> ReadFasta(std::string &path) {

    gzFile fp;
    kseq_t *seq;
    std::vector<FastaRecord> records;
    int l;

    fp = gzopen(path.c_str(), "r");
    if (fp == 0) {
        throw std::invalid_argument("couldn't open file " + path);
    }
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        records.push_back(FastaRecord{
            seq->name.s,
            seq->seq.s,
        });
    }
    kseq_destroy(seq);
    gzclose(fp);
    return records;
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
void AddKMersFromSequence(std::unordered_set<std::string> &kMers, std::string &data, int k) {
    // Convert the sequence to uppercase letters.
    std::transform(data.begin(), data.end(), data.begin(), toupper);
    size_t possibleKMerEnd = k;
    for (size_t i = 1; i <= data.size(); ++i) {
        if (data[i-1] != 'A' && data[i-1] != 'C' && data[i-1] != 'G' && data[i-1] != 'T') {
            // Skip this and the next k-1 k-mers.
            possibleKMerEnd = i + k;
        }
        if (i >= possibleKMerEnd) {
            kMers.insert(data.substr(i - k, k));
        }
    }
}

/// Return only the subset of k-mers where no two k-mers are complements of one another.
/// The k-mers are chosen with no particular property.
std::unordered_set<std::string> FilterKMersWithComplement(std::unordered_set<std::string> &kMers) {
    std::unordered_set<std::string> ret;
    for (auto &&x : kMers) {
        auto kMer = KMer{x};
        if (ret.count(ReverseComplement(kMer).value) == 0) ret.insert(x);
    }
    return ret;
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
std::vector<KMer> ConstructKMers(std::vector<FastaRecord> &data, int k, bool complements) {
    std::unordered_set<std::string> uniqueKMers;
    for (auto &&record : data) {
        AddKMersFromSequence(uniqueKMers, record.sequence, k);
    }
    if (complements) uniqueKMers = FilterKMersWithComplement(uniqueKMers);
    std::vector<KMer> result;
    for (auto &kMer : uniqueKMers) {
        result.push_back(KMer{kMer});
    }
    return result;
}

/// Read encoded k-mers from the given fasta file.
/// Return unique k-mers in no particular order.
/// If complements is set to true, the result contains only one of the complementary k-mers - it is not guaranteed which one.
/// This runs in O(sequence length) expected time.
void ReadKMers(kh_S64_t *kMers, std::string &path, int k, bool complements, bool case_sensitive = false) {
    std::ifstream fasta(path);
    if (fasta.is_open()) {
        char c;
        int beforeKMerEnd = k;
        kmer_t currentKMer = 0;
        kmer_t cases = 0;
        kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
        bool readingHeader = false;
        while (fasta >> std::noskipws >> c) {
            if (c == '>') {
                readingHeader = true;
                currentKMer = 0;
                beforeKMerEnd = k;
            }
            else if (c == '\n') readingHeader = false;
            if (readingHeader) continue;
            auto data = NucleotideToInt(c);
            // Disregard white space.
            if (c == '\n' || c == '\r' || c == ' ') continue;
            if (data == -1) {
                currentKMer = 0;
                beforeKMerEnd = k;
                continue;
            }
            currentKMer <<= 2;
            // K-mer is present if it is upper case or case-insensitive.
            cases |= !case_sensitive || c <= 'Z';
            cases <<= 1;
            currentKMer &= mask;
            currentKMer |= data;
            if(beforeKMerEnd > 0) --beforeKMerEnd;
            if (beforeKMerEnd == 0 && (!complements || kh_get_S64(kMers, ReverseComplement(currentKMer, k)) == kh_end(kMers))) {
                int ret;
                // If the k-mer was masked as present.
                if (cases & (1 << k)) kh_put_S64(kMers, currentKMer, &ret);
            }
        }
        fasta.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
}
