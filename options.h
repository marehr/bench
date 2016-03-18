// ==========================================================================
//                             SeqAn Benchmark
// ==========================================================================
// Copyright (c) 2012-2014, Enrico Siragusa, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_BENCH_OPTIONS_H_
#define APP_BENCH_OPTIONS_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class BaseOptions
// ----------------------------------------------------------------------------

struct BaseOptions
{
    typedef std::string             TString;
    typedef std::vector<TString>    TList;

    enum AlphabetType
    {
        ALPHABET_DNA, ALPHABET_PROTEIN, ALPHABET_CHAR
    };

    CharString      textFile;
    unsigned        textCount;
    unsigned        textSum;
    unsigned        textLen;

    AlphabetType    alphabetType;
    TList           alphabetTypeList;

    bool            tsv;

    BaseOptions() :
        textCount(32),
        textSum(32),
        textLen(32),
        alphabetType(ALPHABET_DNA),
        tsv(false)
    {
        alphabetTypeList.push_back("dna");
        alphabetTypeList.push_back("protein");
        alphabetTypeList.push_back("char");
    }
};

// ----------------------------------------------------------------------------
// Class StreeOptions
// ----------------------------------------------------------------------------

struct StreeOptions : BaseOptions
{
    enum IndexType
    {
        INDEX_SA, INDEX_ESA, INDEX_LST, INDEX_QGRAM, INDEX_FMTL, INDEX_FMWT
    };

    CharString      textIndexFile;
    IndexType       textIndexType;
    TList           textIndexTypeList;

    StreeOptions() :
        BaseOptions(),
        textIndexType(INDEX_SA)
    {
        textIndexTypeList.push_back("sa");
        textIndexTypeList.push_back("esa");
        textIndexTypeList.push_back("lst");
        textIndexTypeList.push_back("qgram");
        textIndexTypeList.push_back("fm-tl");
        textIndexTypeList.push_back("fm-wt");
    }
};

// ----------------------------------------------------------------------------
// Class QGramOptions
// ----------------------------------------------------------------------------

struct QGramOptions : BaseOptions
{
    enum IndexType
    {
        INDEX_QGRAM_DIRECT, INDEX_QGRAM_OPEN
    };

    CharString      textIndexFile;
    IndexType       textIndexType;
    TList           textIndexTypeList;

    unsigned        q;

    QGramOptions() :
        BaseOptions(),
        textIndexType(INDEX_QGRAM_DIRECT),
        q(5)
    {
        textIndexTypeList.push_back("direct");
        textIndexTypeList.push_back("open");
    }
};

// ============================================================================
// Getters & Setters
// ============================================================================

template <typename TOption, typename TString, typename TOptionsList>
void getOptionEnum(TOption & option,
                   TString const & optionStr,
                   TOptionsList & optionsList)
{
    typedef typename Iterator<TOptionsList, Standard>::Type TOptionsIterator;

    TOptionsIterator optionsBegin = begin(optionsList, Standard());
    TOptionsIterator optionsEnd = end(optionsList, Standard());

    TOptionsIterator optionPos = std::find(optionsBegin, optionsEnd, optionStr);

    option = (optionPos != optionsEnd) ? TOption(optionPos - optionsBegin) : TOption();
}

template <typename TOption, typename TString, typename TOptionsList>
void getOptionValue(TOption & option,
                    ArgumentParser const & parser,
                    TString const & optionName,
                    TOptionsList & optionsList)
{
    typedef typename Value<TOptionsList>::Type              TOptionString;

    TOptionString optionStr;
    getOptionValue(optionStr, parser, optionName);

    return getOptionEnum(option, optionStr, optionsList);
}

template <typename TOptions>
void setTextLimits(ArgumentParser & parser, TOptions const & options)
{
    addSection(parser, "Limit Options");
    addOption(parser, ArgParseOption("tc", "text-count", "Limit the number of texts in the collection.", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("ts", "text-sum", "Limit the total length of the text collection.", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("tl", "text-length", "Limit the length of any text in the collection.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "text-count", options.textCount);
    setDefaultValue(parser, "text-sum", options.textSum);
    setDefaultValue(parser, "text-length", options.textLen);
}

template <typename TOptions>
void getTextLimits(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.textCount, parser, "text-count");
    getOptionValue(options.textSum, parser, "text-sum");
    getOptionValue(options.textLen, parser, "text-length");
}

template <typename TOptions>
void setIndexType(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("i", "index", "Select the full-text index type.", ArgParseOption::STRING));
    setValidValues(parser, "index", options.textIndexTypeList);
    setDefaultValue(parser, "index", options.textIndexTypeList[options.textIndexType]);
}

template <typename TOptions>
void getIndexType(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.textIndexType, parser, "index", options.textIndexTypeList);
}

template <typename TOptions>
void setAlphabetType(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("a", "alphabet", "Select the alphabet type.", ArgParseOption::STRING));
    setValidValues(parser, "alphabet", options.alphabetTypeList);
    setDefaultValue(parser, "alphabet", options.alphabetTypeList[options.alphabetType]);
}

template <typename TOptions>
void getAlphabetType(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.alphabetType, parser, "alphabet", options.alphabetTypeList);
}

template <typename TOptions>
void setAlgorithmType(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("g", "algorithm", "Select the search algorithm.", ArgParseOption::STRING));
    setValidValues(parser, "algorithm", options.algorithmTypeList);
    setDefaultValue(parser, "algorithm", options.algorithmTypeList[options.algorithmType]);
}

template <typename TOptions>
void getAlgorithmType(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.algorithmType, parser, "algorithm", options.algorithmTypeList);
}

#endif  // #ifndef APP_BENCH_OPTIONS_H_
