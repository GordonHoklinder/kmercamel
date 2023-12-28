#pragma once
#include "models.h"
#include "kmers.h"
#include "khash_utils.h"

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>
#include <cstdint>
#include <algorithm>
#include <fstream>


/// Find the right extension to the provided last k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<int64_t, int64_t> RightExtension(int64_t last, kh_S64_t *kMers, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (int64_t ext = 0; ext < (1 << (d << 1)); ++ext) {
        int64_t next = BitSuffix(last, k - d) << (d << 1) | ext;
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the left extension to the provided first k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<int64_t, int64_t> LeftExtension(int64_t first, kh_S64_t *kMers, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (int64_t ext = 0; ext < (1 << (d << 1)); ++ext) {
        int64_t next = ext << ((k - d) << 1) | BitPrefix(first, k, k - d);
        if (containsKMer(kMers, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the next generalized simplitig.
/// Update the provided superstring and the mask.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
void NextGeneralizedSimplitig(kh_S64_t *kMers, int64_t begin, std::ostream& of,  int k, int d_max, bool complements) {
     // Maintain the first and last k-mer in the simplitig.
    int64_t last = begin, first = begin;
    std::list<char> simplitig {NucleotideAtIndex(first, k, 0)};
    eraseKMer(kMers, last, k, complements);
    int d_l = 1, d_r = 1;
    while (d_l <= d_max || d_r <= d_max) {
        if (d_r <= d_l) {
            auto extension = RightExtension(last, kMers, k, d_r, complements);
            int64_t ext = extension.first;
            if (ext == -1) {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the generalized simplitig to the right.
                eraseKMer(kMers, extension.second, k, complements);
                for (int i = 1; i < d_r; ++i) simplitig.emplace_back((char)std::tolower(NucleotideAtIndex(last, k, i)));
                simplitig.emplace_back(NucleotideAtIndex(last, k, d_r));
                last = extension.second;
                d_r = 1;
            }
        } else {
            auto extension = LeftExtension(first, kMers, k, d_l, complements);
            int64_t ext = extension.first;
            if (ext == -1) {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                eraseKMer(kMers, extension.second, k, complements);
                for (int i = d_l - 1; i > 0; --i) simplitig.emplace_front((char)std::tolower(NucleotideAtIndex(extension.second, k, i)));
                simplitig.emplace_front(NucleotideAtIndex(extension.second, k, 0));
                first = extension.second;
                d_l = 1;
            }
        }
    }
    for (int i = 1; i < k; ++i) simplitig.emplace_back((char)std::tolower(NucleotideAtIndex(last, k, i)));
    of << std::string(simplitig.begin(), simplitig.end());
    of.flush();
}

/// Compute the generalized simplitigs greedily.
/// This runs in O(n d_max ^ k), where n is the number of k-mers, but for practical uses it is fast.
void GreedyGeneralizedSimplitigs(kh_S64_t *kMers, std::ostream& of, int k, int d_max, bool complements) {
    size_t lastIndex = 0;
    while(true) {
        int64_t begin = nextKMer(kMers, lastIndex);
        // No more k-mers.
        if (begin == -1) return;
        NextGeneralizedSimplitig(kMers, begin, of,  k, d_max, complements);
    }
}