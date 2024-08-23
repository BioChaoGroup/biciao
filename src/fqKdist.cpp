#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <vector>
#include <deque>
#include <tuple>
#include "cmdline.h"
#include "gzstream.h"
#include "ThreadPool/ThreadPool.h"

// Function to reverse complement a k-mer
u_int64_t revComp(u_int64_t x, int sizeKmer)
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

// Function to extract barcode and read name
std::pair<std::string, std::string> getBarcode(std::string& line, std::string& read_type)
{
    if (read_type.empty())
    {
        if (line.find("BX:Z") != std::string::npos)
            read_type = "10x";
        else if (line.find_first_of('#') != std::string::npos)
            read_type = "stLFR";
        else
            read_type = "other";
    }

    std::string read_name, barcode;
    if (read_type == "stLFR")
    {
        std::size_t pos1 = line.find_first_of('#');
        std::size_t pos2 = line.find_first_of('/', pos1 + 1);
        read_name = line.substr(0, pos1);
        barcode = line.substr(pos1 + 1, pos2 - pos1 - 1);
        if (barcode == "0_0_0")
            barcode.clear();
    }
    else if (read_type == "10x")
    {
        read_name = line.substr(0, line.find_first_of(" \r\t\n"));
        std::size_t pos1 = line.find("BX:Z");
        if (pos1 != std::string::npos)
        {
            std::size_t pos2 = line.find_first_of('-', pos1 + 5);
            barcode = line.substr(pos1 + 5, pos2 - pos1 - 5);
        }
    }
    else
    {
        std:size_t read_length = line.find_first_of("/ \r\t\n");
        read_name = line.substr(0, read_length);
        barcode.clear();
    }
    return std::make_pair(read_name, barcode);
}

// Function to generate all possible k-mers of length n
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
        std::vector<std::string> ret_last = generateAllKmers(n - 1);
        for (const auto& kmer : ret_last)
        {
            ret.emplace_back("A" + kmer);
            ret.emplace_back("T" + kmer);
            ret.emplace_back("C" + kmer);
            ret.emplace_back("G" + kmer);
        }
    }
    return ret;
}

// Function to calculate k-mer matrices
std::tuple<std::string, std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
calculateKmerMatrix(const std::string& reads_bc, const std::string& reads_seq, int k, int w,
                    std::vector<std::vector<int>> freqMatrix, std::vector<std::vector<int>> distMatrix, std::vector<std::vector<int>> distMitrix,
                    unsigned long bases)
{
    std::deque<std::pair<u_int64_t, size_t>> window;
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
            u_int64_t kmerIndex = val;

            for (const auto& prevKmer : window)
            {
                u_int64_t prevKmerIndex = prevKmer.first;

                size_t prevPos = prevKmer.second;
                size_t distance = i - prevPos;
                size_t dist_max = distMatrix[prevKmerIndex][kmerIndex];
                size_t dist_min = distMitrix[prevKmerIndex][kmerIndex];

                freqMatrix[prevKmerIndex][kmerIndex]++;

                // got max distance
                distMatrix[prevKmerIndex][kmerIndex] = std::max(distance, dist_max);
                // got min distance
                if (dist_min == 0 || dist_min > distance)
                    distMitrix[prevKmerIndex][kmerIndex] = distance;


                // for reverse complement kmers

                u_int64_t rev_kmerIndex = revComp(val, k);

                u_int64_t rev_prevKmerIndex = revComp(prevKmerIndex, k);

                dist_max = distMatrix[rev_kmerIndex][rev_prevKmerIndex];
                dist_min = distMitrix[rev_kmerIndex][rev_prevKmerIndex];

                freqMatrix[rev_kmerIndex][rev_prevKmerIndex]++;
                // got max distance
                distMatrix[rev_kmerIndex][rev_prevKmerIndex] = std::max(distance, dist_max);
                // got min distance                
                if (dist_min == 0 || dist_min > distance )
                    distMitrix[rev_kmerIndex][rev_prevKmerIndex] = distance;

                // stop when encounter the same kmer
                if (prevKmerIndex == kmerIndex)
                    break;
            }

            window.emplace_front(kmerIndex, i);
            if (window.size() > w)
            {
                window.pop_back();
            }

            --len;

        }
    }
    return std::make_tuple(reads_bc, freqMatrix, distMatrix, distMitrix);
}

int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads1", '1', "Path to reads1 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("reads2", '2', "Path to reads2 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("output", 'o', "Output (gzipped).", true);
    argParser.add<int>("kmer", 'k', "Length of k.", false, 4);
    argParser.add<int>("winsize", 'w', "Sliding window size.", false, 2048);
    argParser.add<int>("thread", 't', "Number of thread to use.", false, 8);
    argParser.add<std::string>("interleaved", 'i', "Path to interleaved reads. If specified, --reads1 and --reads2 are not needed.", false);
    argParser.parse_check(argc, argv);

    std::string reads1 = argParser.get<std::string>("reads1");
    std::string reads2 = argParser.get<std::string>("reads2");
    std::string output = argParser.get<std::string>("output");
    int kmer = argParser.get<int>("kmer");
    int winsize = argParser.get<int>("winsize");
    int thread = argParser.get<int>("thread");
    std::string interleaved = argParser.get<std::string>("interleaved");

    ThreadPool thread_pool(thread);
    std::vector<std::future<std::tuple<std::string, std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>>> results;

    std::vector<std::string> all_kmers = generateAllKmers(kmer);
    unsigned long bases = (unsigned long)std::pow(2, kmer * 2) - 1;

    int numKmers = all_kmers.size();
    std::vector<std::vector<int>> freqMatrix(numKmers, std::vector<int>(numKmers, 0));
    std::vector<std::vector<int>> distMatrix(numKmers, std::vector<int>(numKmers, 0));
    std::vector<std::vector<int>> distMitrix(numKmers, std::vector<int>(numKmers, 0));

    std::string header = "Barcode";
    for (const auto& kmer0 : all_kmers)
    {
        for (const auto& kmer1 : all_kmers)
        {
            header += "," + kmer0 + "-" + kmer1 + "-FREQ";
            header += "," + kmer0 + "-" + kmer1 + "-MAXD";
            header += "," + kmer0 + "-" + kmer1 + "-MIND";
        }
    }
    header += "\n";

    ogzstream outputf(output.c_str());
    outputf << header;

    if (interleaved.empty())
    {
        igzstream read1_file(reads1.c_str());
        igzstream read2_file(reads2.c_str());

        if (!read1_file.good() || !read2_file.good())
        {
            std::cerr << "Error opening reads files." << std::endl;
            return 1;
        }

        std::string head1, head2, seq1, seq2, seqs, line1, line2, read_name1, read_name2, read_bc1, read_bc2;
        std::string read_type;
        std::string prev_name, this_name;

        while (std::getline(read1_file, head1))
        {
            std::getline(read1_file, seq1);
            std::getline(read1_file, line1);
            std::getline(read1_file, line1);

            std::getline(read2_file, head2);
            std::getline(read2_file, seq2);
            std::getline(read2_file, line2);
            std::getline(read2_file, line2);

            std::tie(read_name1, read_bc1) = getBarcode(head1, read_type);
            std::tie(read_name2, read_bc2) = getBarcode(head2, read_type);

            if (read_bc1.empty() && read_bc2.empty())
            {
                this_name = read_name1;
            }
            else
            {
                this_name = read_bc1;
            }

            if (!prev_name.empty() && prev_name.compare(this_name))
            {
                results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, prev_name, seqs, kmer, winsize, freqMatrix, distMatrix, distMitrix, bases));
                seqs.clear();
            }

            seqs += "N" + seq1 + "N" + seq2;

            prev_name = this_name;
        }
        results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, prev_name, seqs, kmer, winsize, freqMatrix, distMatrix, distMitrix, bases));
        read1_file.close();
        read1_file.close();
    }
    else
    {
        igzstream read_file(interleaved.c_str());

        if (!read_file.good())
        {
            std::cerr << "Error opening interleaved reads file." << std::endl;
            return 1;
        }

        std::string line, read_name1, read_name2, read_bc1, read_bc2, seq1, seq2, seqs;
        std::string read_type;
        std::string prevBarcode = "";
        int cnt = 0;
        while (std::getline(read_file, line))
        {
            cnt++;
            std::getline(read_file, line);

            if (cnt % 8 == 1)
            {
                std::tie(read_name1, read_bc1) = getBarcode(line, read_type);
            }
            else if (cnt % 8 == 2)
            {
                seq1 = line;
            }
            else if (cnt % 8 == 6)
            {
                seq2 = line;
            }
            else if (cnt % 8 == 5)
            {
                std::tie(read_name2, read_bc2) = getBarcode(line, read_type);
            }
            else if (cnt % 8 == 0)
            {
                if (read_bc1.empty() && read_bc2.empty())
                {
                    results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, read_name1, seqs, kmer, winsize, freqMatrix, distMatrix, distMitrix, bases));
                    prevBarcode = read_bc1;
                    seqs.clear();
                }

                if (!prevBarcode.empty() && prevBarcode.compare(read_bc1))
                {
                    results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, prevBarcode, seqs, kmer, winsize, freqMatrix, distMatrix, distMitrix, bases));
                    prevBarcode = read_bc1;
                    seqs.clear();
                }
                else
                {
                    seqs += "N" + seq1 + "N" + seq2;
                }
                prevBarcode = read_bc1;
            }
        }
        results.emplace_back(thread_pool.enqueue(calculateKmerMatrix, prevBarcode, seqs, kmer, winsize, freqMatrix, distMatrix, distMitrix, bases));
        read_file.close();
    }

    for (auto&& result : results)
    {
        auto [barcode, freqMatrix, distMatrix, distMitrix] = result.get();

        if (barcode.empty())
        {
            continue;
        }

        std::string outputLine = barcode;

        for (int i = 0; i < numKmers; ++i)
        {
            for (int j = 0; j < numKmers; ++j)
            {
                outputLine += "," + std::to_string(freqMatrix[i][j]);
                outputLine += "," + std::to_string(distMatrix[i][j]);
                outputLine += "," + std::to_string(distMitrix[i][j]);
            }
        }

        outputLine += "\n";
        outputf << outputLine;
    }

    return 0;
}
