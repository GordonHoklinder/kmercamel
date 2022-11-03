#include "kmers.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(BitSuffixTest, BitSuffix) {
        EXPECT_EQ(0b011110, BitSuffix(0b1100011110, 3));
        EXPECT_EQ(0b1110, BitSuffix(0b1110, 2));
        EXPECT_EQ(0b0, BitSuffix(0b1110, 0));
        EXPECT_EQ(0b111111'11111111'11111111'11111111'11111111'11111110LL, BitSuffix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL,  23));
        EXPECT_EQ(0b111111'11111110LL, BitSuffix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 7));
    }

    TEST(BitPrefixTest, BitPrefix) {
        EXPECT_EQ(0b110001, BitPrefix(0b1100011110, 5, 3));
        EXPECT_EQ(0b1110, BitPrefix(0b1110, 2, 2));
        EXPECT_EQ(0b0, BitPrefix(0b1110, 2, 0));
        EXPECT_EQ(0b111111'01111111'11111111'11111111'11111111'11111111LL, BitPrefix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, 23));
        EXPECT_EQ(0b111111'01111111LL, BitPrefix(0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, 7));
    }

    TEST(NucleotideToIntTest, NucleotideToInt) {
        struct TestCase {
            char nucleotide;
            int wantResult;
            bool wantSuccess;
        };
        std::vector<TestCase> tests = {
                {'A', 0, true},
                {'C', 1, true},
                {'G', 2, true},
                {'T', 3, true},
                {'B', 0, false},
        };

        for (auto t : tests) {
            int gotResult = 0;
            bool gotSuccess = true;
            try {
                gotResult = NucleotideToInt(t.nucleotide);
            } catch(std::invalid_argument) {
                gotSuccess = false;
            }
            EXPECT_EQ(t.wantSuccess, gotSuccess);
            if (t.wantSuccess) EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(NumberToKMerTest, NumberToKMer) {
        struct TestCase {
            long long encoded;
            int d;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, "GC"},
                {0b1011LL, 3, "AGT"},
                {0b111LL, 1, "T"},
                {0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, "TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"},
        };

        for (auto t: tests) {
            std::string gotResult = NumberToKMer(t.encoded, t.d);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(KMerToNumberTest, KMerToNumber) {
        struct TestCase {
            long long wantResult;
            KMer kMer;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, "GC"},
                {0b1011LL, "AGT"},
                {0b11LL, "T"},
                {0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, "TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"},
        };

        for (auto t: tests) {
            long long gotResult = KMerToNumber(t.kMer);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(ReverseComplementTest, ReverseComplement) {
        struct TestCase {
            long long input;
            int k;
            long long wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, 0b1001LL},
                {0b101111LL, 3, 0b000001LL},
                {0b11LL, 1, 0b00LL},
                {0b111111'01111111'11111111'11111111'11111111'11111111'11111111'11111110LL, 31, 0b010000'00000000'00000000'00000000'00000000'00000000'00000000'10000000LL },
        };

        for (auto t: tests) {
            long long gotResult = ReverseComplement(t.input, t.k);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }
}
