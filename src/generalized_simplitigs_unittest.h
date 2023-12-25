#pragma once
#include "generalized_simplitigs.h"

#include "gtest/gtest.h"

namespace {
    TEST(RightExtensionTest, RightExtension) {
        struct TestCase {
            int64_t last;
            std::unordered_set<int64_t> kMers;
            int k;
            int d;
            bool complements;
            int64_t wantExt;
            int64_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, CTA, ACT, CCT}; A; CTA
                {0b000111, {0b110101, 0b011100, 0b000111, 0b010111}, 3, 1, false, 0b00, 0b011100},
                // ACT; {TCC, ACT, CCT}; CC; TCC
                {0b111111, {0b110101, 0b000111, 0b010111}, 3, 2,   false, 0b0101, 0b110101},
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, 1,  false,  -1, -1},
        };

        for (auto t: tests) {
            auto got = RightExtension(t.last, t.kMers, t.k, t.d, t.complements);
            auto gotExt = got.first;
            auto gotNext = got.second;
            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }

    TEST(LeftExtensionTest, LeftExtension) {
        struct TestCase {
            int64_t first;
            std::unordered_set<int64_t> kMers;
            int k;
            int d;
            bool complements;
            int64_t wantExt;
            int64_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, 1, false,   -1, -1},
                // TAC; {TCC, CTA, ACT, CCT}; C; CTA
                {0b110001, {0b110101, 0b011100, 0b000111, 0b010111}, 3, 1, false, 0b01, 0b011100},
        };

        for (auto t: tests) {
            auto got = LeftExtension(t.first, t.kMers, t.k, t.d, t.complements);
            auto gotExt = got.first;
            auto gotNext = got.second;
            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }


    TEST(NextGeneralizedSimplitigTest, NextGeneralizedSimplitig) {
        struct TestCase {
            std::unordered_set<int64_t> kMers;
            int k;
            int d_max;
            bool complements;
            std::vector<std::string> wantSuperstring;
            std::unordered_set<int64_t> wantKMers;
        };
        std::vector<TestCase> tests = {
                // {ACAA, AACA}
                {{0b00010000,  0b00000100}, 4, 2, false,
                 {"AAcaa"}, {},},
                // {ACAA, ATTT, TGTT, AAAT, TTGT, AACA}
                {{0b00010000, 0b00111111, 0b11101111, 0b00000011, 0b11111011, 0b00000100}, 4, 2, true,
                 {"AAcAaat", "AtTTgtt" }, {},},
        };

        for (auto &&t: tests) {
            std::stringstream of;

            NextGeneralizedSimplitig(t.kMers, of, t.k, t.d_max, t.complements);
            auto gotSuperstring = of.str();

            EXPECT_EQ(t.wantKMers, t.kMers);
            // Check that at least one valid superstring was returned.
            bool valid = false;
            for (size_t i = 0; i < t.wantSuperstring.size(); ++i) {
                valid |= t.wantSuperstring[i] == gotSuperstring;
            }
            EXPECT_TRUE(valid);
        }
    }

    TEST(GreedyGeneralizedSimplitigsTest, GreedyGeneralizedSimplitigs) {
        struct TestCase {
            std::vector<int64_t> kMers;
            int k;
            int d_max;
            bool complements;
            std::string wantSuperstring;
        };
        std::vector<TestCase> tests = {
                {{KMerToNumber(KMer{"GCT"}), KMerToNumber(KMer{"TAA"}), KMerToNumber(KMer{"AAA"})}, 3, 2, false, "GcTAaa"},
                {{KMerToNumber(KMer{"TAA"}), KMerToNumber(KMer{"AAA"}), KMerToNumber(KMer{"GCT"})}, 3, 2, false, "GcTAaa"},
                {{KMerToNumber(KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber(KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, 31, 5, false,
                 "TTtcttttttttttttttttttttttttttga"},
        };

        for (auto t: tests) {
            std::stringstream of;

            GreedyGeneralizedSimplitigs(t.kMers, of, t.k, t.d_max, t.complements);

            EXPECT_EQ(t.wantSuperstring, of.str());
        }
    }
}
