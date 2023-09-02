#pragma once
#include "models.h"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <cstdint>

#include "kmers.h"
#include "khash.h"


/// Provide possibility to access reverse complements as if they were in the field.
#define access(field, index) (((field).size() > (index)) ? (field)[(index)] : \
        ReverseComplement((field)[(index) - (field).size()], k))

KHASH_MAP_INIT_INT64(P64, size_t)

constexpr int MEMORY_REDUCTION_FACTOR = 16;

/// Represents oriented edge in the overlap graph.
struct OverlapEdge {
    // Index of the first k-mer.
    size_t firstIndex;
    // Index of the second k-mer.
    size_t secondIndex;
    // Length of the overlap of the two k-mers.
    int overlapLength;
};

typedef std::pair<std::vector<size_t>, std::vector<unsigned char>> overlapPath;

/// Greedily find the approximate Hamiltonian path with longest overlaps.
/// k is the size of one k-mer and n is the number of distinct k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected to contain only one k-mer from a complement pair.
/// Moreover, if so, the resulting Hamiltonian path contains two superstrings which are reverse complements of one another.
overlapPath OverlapHamiltonianPath (std::vector<int64_t> &kMers, int k, bool complements) {

    size_t n = kMers.size();
    size_t kMersCount = n * (1 + complements);
    size_t batch_size = kMersCount / MEMORY_REDUCTION_FACTOR + 1;
    std::vector<size_t> edgeFrom(kMersCount, -1);
    std::vector<unsigned char> overlaps(kMersCount, -1);
    std::vector<bool> suffixForbidden(kMersCount, false);
    std::vector<bool> prefixForbidden(kMersCount, false);
    std::vector<size_t> first(kMersCount);
    std::vector<size_t> last(kMersCount);
    std::vector<size_t> next(kMersCount, -1);
    for (size_t i = 0; i < kMersCount; ++i) {
        first[i] = last[i] = i;
    }
    khash_t(P64)  *prefixes = kh_init(P64);
    kh_resize(P64, prefixes, (kMersCount / MEMORY_REDUCTION_FACTOR + 1 ) * 100 / 77 );
    for (int d = k - 1; d >= 0; --d) {
        // In order to reduce memory requirements, the prefixes are not processed at once, but in batches.
        // As a cost, this slows down the algorithm.
        for (int part = 0; part < MEMORY_REDUCTION_FACTOR; part++) {
            kh_clear(P64, prefixes);
            size_t to = std::min(kMersCount, (part + 1) * batch_size);
            for (size_t i = part * batch_size; i < to; ++i)
                if (!prefixForbidden[i]) {
                    next[i] = -1;
                    int64_t prefix = BitPrefix(access(kMers,i), k, d);
                    auto prefix_key = kh_get(P64, prefixes, prefix);
                    if (prefix_key != kh_end(prefixes)) {
                        next[i] = kh_val(prefixes, prefix_key);
                    } else {
                        int ret;
                        prefix_key = kh_put(P64, prefixes, prefix, &ret);

                    }
                    kh_value(prefixes, prefix_key) = i;
                }
            for (size_t i = 0; i < kMersCount; ++i)
                if (!suffixForbidden[i]) {
                    int64_t suffix = BitSuffix(access(kMers, i), d);
                    auto suffix_key = kh_get(P64, prefixes, suffix);
                    if (suffix_key == kh_end(prefixes)) continue;
                    size_t previous, j;
                    previous = j = kh_val(prefixes, suffix_key);
                    // If the path forms a cycle, or is between k-mer and its reverse complement, or the k-mers complement was already selected skip this path.
                    while (j != size_t(-1) &&
                           (first[i] % n == j % n || first[i] % n == last[j] % n || prefixForbidden[j])) {
                        size_t new_j = next[j];
                        // If the k-mer is forbidden, remove it to keep the complexity linear.
                        // This is not done with the first k-mer but that is not a problem.
                        if (prefixForbidden[j]) next[previous] = new_j;
                        else previous = j;
                        j = new_j;
                    }
                    if (j == size_t(-1)) {
                        continue;
                    }
                    std::vector<std::pair<size_t, size_t>> new_edges({{i, j}});
                    // Add also the edge between complementary k-mers in the opposite direction.
                    if (complements) new_edges.emplace_back((j + n) % kMersCount, (i + n) % kMersCount);
                    for (auto [x, y]: new_edges) {
                        edgeFrom[x] = y;
                        overlaps[x] = d;
                        prefixForbidden[y] = true;
                        first[last[y]] = first[x];
                        last[first[x]] = last[y];
                        suffixForbidden[x] = true;
                    }
                    next[previous] = next[j];
                }
        }
    }

    kh_destroy(P64, prefixes);
    return {edgeFrom, overlaps};
}

/// Construct the superstring and its mask from the given overlapPath path in the overlap graph.
/// If reverse complements are considered and the overlapPath path contains two paths which are reverse complements of one another,
/// return only one of them.
void SuperstringFromPath(const overlapPath &hamiltonianPath, const std::vector<int64_t> &kMers, std::ostream& of, const int k, const bool complements) {
    size_t kMersCount = kMers.size() * (1 + complements);
    auto edgeFrom = hamiltonianPath.first;
    auto overlaps = hamiltonianPath.second;

    // Find the vertex in the overlap graph with in-degree 0.
    std::vector<bool> isStart(kMersCount, true);
    for (auto edge : edgeFrom) {
        if (edge != size_t(-1)) isStart[edge] = false;
    }
    size_t start = 0;
    for (; start < kMersCount && !isStart[start]; ++start);

    int64_t last = BitSuffix(access(kMers, start), k-1);
    of << letters[BitPrefix(access(kMers, start), k, 1)];

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start] != size_t(-1)) {
        int overlapLength = overlaps[start];
        if (overlapLength != k - 1) {
            std::string unmaskedNucleotides = NumberToKMer(BitPrefix(last, k-1, k-1-overlapLength), k-1-overlapLength);
            std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
            of << unmaskedNucleotides;
        }
        last = BitSuffix(access(kMers, edgeFrom[start]), k-1);
        of << letters[BitPrefix(access(kMers, edgeFrom[start]), k, 1)];
        start = edgeFrom[start];
    }

    // Print the trailing k-1 characters.
    std::string unmaskedNucleotides = NumberToKMer(last, k-1);
    std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
    of << unmaskedNucleotides;
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
void Greedy(std::vector<int64_t> &kMers, std::ostream& of, int k, bool complements) {
    if (kMers.empty()) {
        throw std::invalid_argument("input cannot be empty");
    }
    auto hamiltonianPath = OverlapHamiltonianPath(kMers, k, complements);
    SuperstringFromPath(hamiltonianPath, kMers, of, k, complements);
}

// Undefine the access macro, so it does not interfere with other files.
#undef access
