#ifndef BENCHMARK_HEADER_PAIRWISE_ALIGNMENT_BASE_H
#define BENCHMARK_HEADER_PAIRWISE_ALIGNMENT_BASE_H

#include <iostream>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

using namespace seqan;

struct Options {
    unsigned threads;
    CharString input;
    CharString output;
};

seqan::ArgumentParser::ParseResult parseCommandLine(std::string infoText, Options & options, int argc, char** const argv)
{
    seqan::ArgumentParser parser(infoText);
    setCategory(parser, "Demo");
    setVersion(parser, "0.1");
    setDate(parser, "Nov 2015");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN\\fP \\fIOUT\\fP ");

    addDescription(parser,
    "\\fIIN\\fP is a .fa / fasta input file.  \\fIOUT\\fP is a txt output file.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "IN"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, 0, "fasta, fa");
    setValidValues(parser, 1, "txt");

    addSection(parser, "Settings");
    addOption(parser, seqan::ArgParseOption("tc", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));

    #ifdef _OPENMP
    setDefaultValue(parser, "threads", omp_get_max_threads());
    #else
    setDefaultValue(parser, "threads", 1);
    #endif

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if(res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.threads, parser, "threads");
        getArgumentValue(options.input, parser, 0);
        getArgumentValue(options.output, parser, 1);
    }
    return res;
}

template <typename TSequence, typename TLambda>
inline int benchmark_pairwise_alignment_main(std::string infoText, TLambda computeAlignment, int argc, char **argv)
{
    Options options;
    ArgumentParser::ParseResult res = parseCommandLine(infoText, options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    #ifdef _OPENMP
    omp_set_num_threads(options.threads);
    #endif

    StringSet<CharString> ids;
    StringSet<TSequence> sequences;
    SeqFileIn seqFileIn(toCString(options.input));
    readRecords(ids, sequences, seqFileIn);
    clear(ids);

    unsigned const threads = options.threads;
    unsigned const len = length(sequences);
    unsigned const maxAlignments = len * (len - 1) / 2;
    unsigned const alignmentsPerThread = maxAlignments / threads + 1;
    unsigned const lastBucketSize = maxAlignments - alignmentsPerThread * (threads - 1);

    // allocate data structures
    std::vector<StringSet<Gaps<TSequence>>> gapsSetVectorH;
    std::vector<StringSet<Gaps<TSequence>>> gapsSetVectorV;
    String<String<int>> scores;

    resize(gapsSetVectorH, threads);
    resize(gapsSetVectorV, threads);
    resize(scores, threads);

    for (int jobId = 0; jobId < threads; jobId++)
    {
        auto bucketSize = jobId < threads - 1 ? alignmentsPerThread : lastBucketSize;
        resize(gapsSetVectorH[jobId], bucketSize);
        resize(gapsSetVectorV[jobId], bucketSize);
        // don't allocate scores[jobId], because this will be done in
        // computeAlignment by every thread
    }

    // initialize data structures
    for (int m = 0, k = 0; m < len; m++)
    for (int n = m + 1; n < len; n++, k++)
    {
        auto jobId = k / alignmentsPerThread;
        auto currentJob = k % alignmentsPerThread;

        auto & gapsH = gapsSetVectorH[jobId][currentJob];
        auto & gapsV = gapsSetVectorV[jobId][currentJob];
        assignSource(gapsH, sequences[m]);
        assignSource(gapsV, sequences[n]);
    }

    // compute data structures
    SEQAN_OMP_PRAGMA(parallel for)
    for (int jobId = 0; jobId < threads; ++jobId)
    {
        scores[jobId] = computeAlignment(gapsSetVectorH[jobId], gapsSetVectorV[jobId]);
    }

    //serial output scores
    std::ofstream ofs(toCString(options.output), std::ofstream::out);
    for (int m = 0, k = 0; m < len; m++)
    for (int n = m + 1; n < len; n++, k++)
    {
        auto jobId = k / alignmentsPerThread;
        auto currentJob = k % alignmentsPerThread;
        ofs << m << ", " << n << ": " << scores[jobId][currentJob] << std::endl;
    }

    return 0;
}

#endif
