
#pragma once

#include "../src/conversions.h"
#include "gtest/gtest.h"

namespace {

// Retrieving current path on Windows does not work as on linux
// therefore the following unittests are linux-specific.
#ifdef __unix__
    TEST(Conversions, SplitMS) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/runstest.fa";

        std::stringstream mask;
        std::stringstream superstring;

        split_ms(superstring, mask, path); 

        EXPECT_EQ("10011100010110\n", mask.str());
        EXPECT_EQ("ACGCGTTACGTATT\n", superstring.str());
    }

    TEST(Conversions, MSToSPSS) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/runstest.fa"; 

        struct TestCase {
            int k;
            std::string wantResult;
        };

        auto tests = std::vector<TestCase> {
            {2, ">0\nAC\n>1\nCGTT\n>2\nGT\n>3\nATT\n"},
            {3, ">0\nACG\n>1\nCGTTA\n>2\nGTA\n>3\nATT\n"}
        };

        for (auto t : tests) {
            std::stringstream of;

            ms_to_spss(path, of, t.k);

            EXPECT_EQ(t.wantResult, of.str());       
        }
    }

    TEST(Conversions, SPSSToMS) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/spss.fa"; 

        struct TestCase {
            int k;
            std::string wantResult;
        };

        auto tests = std::vector<TestCase> {
            {1, ">superstring " + path + "\nACCCGAACCGTAATGCACCCGTTTAACGA\n"},
            {2, ">superstring " + path + "\nACCCGAAcCGTaATGcACCCGTTTAACg\n"},
            {3, ">superstring " + path + "\nACCCGAacCGtaATgcACCCGTTTAAcg\n"},
            {4, ">superstring " + path + "\nACCCGaacCgtaAtgcACCCGTTTAacg\n"},
            {5, ">superstring " + path + "\nACCCgaacACCCGTTTaacg\n"}
        };

        for (auto t : tests) {
            std::stringstream of;

            spss_to_ms(path, of, t.k);

            EXPECT_EQ(t.wantResult, of.str());       
        }
    }
#endif

    TEST(Conversions, JoinMS) {
        struct TestCase {
            std::string mask;
            std::string superstring;
            std::string wantResult;
        };

        auto tests = std::vector<TestCase> {
            {"1001101\n", "ACGTTTA\n", ">superstring\nAcgTTtA\n"},
            {"0111110000\n", "CCCGGTATTA\n", ">superstring\ncCCGGTatta\n"},
        };

        for (auto t : tests) {
            std::istringstream im (t.mask);
            std::istringstream is (t.superstring);
            std::stringstream of;

            join_ms(is, im, of);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
