#include "pairwise_alignment_base.h"

using namespace seqan;

int main(int argc, char **argv)
{
    typedef Dna5String TSequence;
    typedef StringSet<Gaps<TSequence>> TGapSet;

    auto compute = [] (TGapSet & gapsSetH, TGapSet & gapsSetV){
        // match: 2, mismatch: -3, gap: -2
        return globalAlignment(gapsSetH, gapsSetV, Score<int>(2, -3, -2), AlignConfig<>(), NeedlemanWunsch());
    };
    return benchmark_pairwise_alignment_main<TSequence>("Benchmark: Global Pairwise Alignment with Linear Gap Model for DNA.", compute, argc, argv);
}
