#include "pairwise_alignment_base.h"

using namespace seqan;

int main(int argc, char **argv)
{
    typedef Dna5String TSequence;
    typedef StringSet<Gaps<TSequence>> TGapSet;

    auto compute = [] (TGapSet & gapsSetH, TGapSet & gapsSetV){
        // match: 2, mismatch: -3, gap: -2
        return localAlignment(gapsSetH, gapsSetV, Score<int>(2, -3, -2), NeedlemanWunsch());
    };
    return benchmark_pairwise_alignment_main<TSequence>("Benchmark: Local Pairwise Alignment with Linear Gap Model for DNA.", compute, argc, argv);
}
