#include "pairwise_alignment_base.h"

using namespace seqan;

int main(int argc, char **argv)
{
    typedef Dna5String TSequence;
    typedef StringSet<Gaps<TSequence>> TGapSet;

    auto compute = [] (TGapSet & gapsSetH, TGapSet & gapsSetV){
        // match: 2, mismatch: -3, gap_extend: -3, gap_open: -1
        return localAlignment(gapsSetH, gapsSetV, Score<int>(2, -3, -3, -1), Gotoh());
    };
    return benchmark_pairwise_alignment_main<TSequence>("Benchmark: Local Pairwise Alignment with Affine Gap Model for DNA.", compute, argc, argv);
}
