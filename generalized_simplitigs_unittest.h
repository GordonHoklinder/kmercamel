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
            int64_t wantExt;
            int64_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, CTA, ACT, CCT}; A; CTA
                {0b000111, {0b110101, 0b011100, 0b000111, 0b010111}, 3, 1, 0b00, 0b011100},
                // ACT; {TCC, ACT, CCT}; CC; TCC
                {0b111111, {0b110101, 0b000111, 0b010111}, 3, 2,   0b0101, 0b110101},
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, 1,    -1, -1},
        };

        for (auto t: tests) {
            auto got = RightExtension(t.last, t.kMers, t.k, t.d);
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
            int64_t wantExt;
            int64_t  wantNext;
        };
        std::vector<TestCase> tests = {
                // ACT; {TCC, ACT, CCT}
                {0b000111, {0b110101, 0b000111, 0b010111}, 3, 1,    -1, -1},
                // TAC; {TCC, CTA, ACT, CCT}; C; CTA
                {0b110001, {0b110101, 0b011100, 0b000111, 0b010111}, 3, 1,    0b01, 0b011100},
        };

        for (auto t: tests) {
            auto got = LeftExtension(t.first, t.kMers, t.k, t.d);
            auto gotExt = got.first;
            auto gotNext = got.second;
            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }


    TEST(NextGeneralizedSimplitigTest, NextGeneralizedSimplitig) {
        struct TestCase {
            std::string superstring;
            std::unordered_set<int64_t> kMers;
            std::vector<bool> mask;
            int k;
            int d_max;
            bool complements;
            std::string wantSuperstring;
            std::unordered_set<int64_t> wantKMers;
            std::vector<bool> wantMask;
        };
        std::vector<TestCase> tests = {
                // Behavior of the unordered_set.begin() is not specified. It may happen that it takes "ATTT" first and then this test fails.
                // {ACAA, ATTT, AACA}
                {"GGGG", {0b00010000, 0b00111111, 0b00000100}, {1,0,0, 0}, 4, 2, false,
                 "GGGGAACAA", {0b00111111}, {1,0,0,0,1,1,0,0, 0}},
                // {ACAA, ATTT, TGTT, AAAT, TTGT, AACA}
                {"GGGG", {0b00010000, 0b00111111, 0b11101111, 0b00000011, 0b11111011, 0b00000100}, {1,0,0, 0}, 4, 2, true,
                        "GGGGAACAAAT", {}, {1,0,0,0, 1,1,0,1,0,0, 0}},
        };

        for (auto t: tests) {
            NextGeneralizedSimplitig(t.kMers, t.superstring, t.mask, t.k, t.d_max, t.complements);

            EXPECT_EQ(t.wantSuperstring, t.superstring);
            EXPECT_EQ(t.wantMask, t.mask);
            EXPECT_EQ(t.wantKMers, t.kMers);
        }
    }

    TEST(GreedyGeneralizedSimplitigsTest, GreedyGeneralizedSimplitigs) {
        struct TestCase {
            std::vector<KMer> kMers;
            int k;
            int d_max;
            bool complements;
            std::string wantSuperstring;
            std::vector<bool> wantMask;
        };
        std::vector<TestCase> tests = {
                // Behavior of the unordered_set.begin() is not specified. It may happen that it takes "CCCC" first and then this test fails.
                { {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, 4, 3, false,
                 "AACAATTTCCCC", {1,1,0,0,1, 0,0,0, 1, 0,0,0}},
                { {KMer{"GCT"}, KMer{"TAA"}, KMer{"AAA"}}, 3, 2, false,
                        "GCTAAA", {1,0, 1,1, 0,0}},
                { {KMer{"TAA"}, KMer{"AAA"}, KMer{"GCT"}}, 3, 2, false,
                        "GCTAAA", {1,0, 1,1, 0,0}},
                { {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, 4, 3, true,
                        "ATTTGTTGGGG", {1,0,1, 1,0,0,0, 1, 0,0,0}},
                {{KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}, KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"}}, 31, 5, false,
                 "TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA", std::vector<bool> {1,1,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0}},
        };

        for (auto t: tests) {
            auto gotResult = GreedyGeneralizedSimplitigs(t.kMers, t.k, t.d_max, t.complements);

            EXPECT_EQ(t.wantSuperstring, gotResult.superstring);
            EXPECT_EQ(t.wantMask, gotResult.mask);
            EXPECT_EQ(t.k, gotResult.k);
        }
    }


}