#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <vector>
#include <deque>
#include "cmdline.h"
#include "gzstream.h"
#include "ThreadPool/ThreadPool.h"

u_int64_t revComp(u_int64_t& x, int& sizeKmer)
{
    u_int64_t res = x;
    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;
    return (res >> (2 * (32 - sizeKmer)));
}

std::vector<std::string> generateAllKmers(int n)
{
    std::vector<std::string> ret;
    if (n == 1)
    {
        ret.emplace_back("A");
        ret.emplace_back("T");
        ret.emplace_back("C");
        ret.emplace_back("G");
    }
    else
    {
        std::vector<std::string> ret_last = std::move(generateAllKmers(n - 1));
        for (auto&& kmer : ret_last)
        {
            ret.emplace_back("A" + kmer);
            ret.emplace_back("T" + kmer);
            ret.emplace_back("C" + kmer);
            ret.emplace_back("G" + kmer);
        }
    }
    return ret;
}

void calculateKmerMatrix(const std::string& reads_seq, int k, int w, std::vector<std::vector<double>>& freqMatrix, std::vector<std::vector<double>>& distMatrix, unsigned long bases)
{
    std::deque<std::pair<u_int64_t, size_t>> window; // Stores k-mers and their positions in the current window
    u_int64_t val = 0;
    std::size_t len = 0;
    
    for (std::size_t i = 0; i < reads_seq.length(); ++i)
    {
        if (!(reads_seq[i] == 'A' || reads_seq[i] == 'C' || reads_seq[i] == 'G' || reads_seq[i] == 'T'))
        {
            val = 0;
            len = 0;
            window.clear();
            continue;
        }
        val = (val << 2);
        val = val & bases;
        val += (reads_seq[i] >> 1 & 3);
        ++len;

        if (len == k)
        {
            // u_int64_t kmerIndex = std::min(val, revComp(val, k));
            u_int64_t kmerIndex = val;

            // Update the matrix for the current k-mer against the previous k-mers in the window
            for (const auto& prevKmer : window)
            {
                u_int64_t prevKmerIndex = prevKmer.first;
                size_t prevPos = prevKmer.second;
                size_t distance = i - prevPos;

                freqMatrix[prevKmerIndex][kmerIndex]++;
                distMatrix[prevKmerIndex][kmerIndex] += distance;
            }
            
            // cal revComp kmer distance and freq
            u_int64_t revkmerIndex = revComp(val, k);

            for (const auto& prevKmer : window)
            {
                u_int64_t    prevKmerIndex = prevKmer.first;
                u_int64_t revprevKmerIndex = revComp(prevKmerIndex,k);
                size_t prevPos = prevKmer.second;
                size_t distance = i - prevPos;

                freqMatrix[revprevKmerIndex][revkmerIndex]++;
                distMatrix[revprevKmerIndex][revkmerIndex] += distance;
            }           

            // Add the current k-mer to the window
            window.emplace_back(kmerIndex, i);

            // Remove k-mers that are out of the 5000 bp window
            if (window.size() > w)
            {
                window.pop_front();
            }

            --len; // Move to the next k-mer
        }
    }

    // Calculate average distances
    for (int i = 0; i < freqMatrix.size(); ++i)
    {
        for (int j = 0; j < freqMatrix.size(); ++j)
        {
            if (freqMatrix[i][j] > 0)
            {
                distMatrix[i][j] /= freqMatrix[i][j];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads", 'r', "Path to reads (gzipped FASTA).", true);
    argParser.add<std::string>("output", 'o', "Output (gzipped).", true);
    argParser.add<int>("kmer", 'k', "Length of k.", false, 4);
    argParser.add<int>("winsize", 'w', "Sliding window size.", false, 2048);
    argParser.add<int>("thread", 't', "Number of thread to use.", false, 16);
    argParser.parse_check(argc, argv);

    std::string reads = argParser.get<std::string>("reads");
    std::string output = argParser.get<std::string>("output");
    int kmer = argParser.get<int>("kmer");
    int winsize = argParser.get<int>("winsize");
    int thread = argParser.get<int>("thread");

    ThreadPool thread_pool(thread);
    std::vector<std::future<void>> results;

    std::vector<std::string> all_kmers = std::move(generateAllKmers(kmer));
    unsigned long bases = (unsigned long)std::pow(2, kmer * 2) - 1;

    // Initialize frequency and distance matrices
    int numKmers = all_kmers.size();
    std::vector<std::vector<double>> freqMatrix(numKmers, std::vector<double>(numKmers, 0));
    std::vector<std::vector<double>> distMatrix(numKmers, std::vector<double>(numKmers, 0));

    igzstream readsf(reads.c_str());
    unsigned long long cnt_line = 0;
    std::string line, reads_seq;

    while (getline(readsf, line))
    {
        if (line[0] == '>')
        {
            if (!reads_seq.empty())
            {
                results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, reads_seq, kmer, winsize, std::ref(freqMatrix), std::ref(distMatrix), bases));
                reads_seq.clear();
            }
        }
        else
        {
            reads_seq += line;
        }
    }
    if (!reads_seq.empty())
    {
        results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, reads_seq, kmer, winsize, std::ref(freqMatrix), std::ref(distMatrix), bases));
    }
    readsf.close();

    for (auto&& result : results)
    {
        result.get();
    }

    ogzstream outputf(output.c_str());

    // Write the header line
    outputf << "KMER";
    for (const auto& kmer : all_kmers)
    {
        outputf << "," << kmer;
    }
    outputf << std::endl;

    // Write the frequency and distance matrices
    for (int i = 0; i < numKmers; ++i)
    {
        outputf << all_kmers[i];
        for (int j = 0; j < numKmers; ++j)
        {
            if (i <= j)
            {
                outputf << "," << freqMatrix[i][j];
            }
            else
            {
                outputf << "," << distMatrix[i][j];
            }
        }
        outputf << std::endl;
    }

    outputf.close();
    return 0;
}
