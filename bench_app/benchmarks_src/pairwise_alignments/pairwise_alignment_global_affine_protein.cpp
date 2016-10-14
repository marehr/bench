#include "pairwise_alignment_base.h"

using namespace seqan;

int main(int argc, char **argv)
{
    typedef String<AminoAcid> TSequence;
    typedef StringSet<Gaps<TSequence>> TGapSet;

    auto compute = [] (TGapSet & gapsSetH, TGapSet & gapsSetV){
        // Blosum62, gap_extend: -3, gap_open: -1
        return globalAlignment(gapsSetH, gapsSetV, Blosum62(-3, -1), AlignConfig<>(), Gotoh());
    };
    return benchmark_pairwise_alignment_main<TSequence>("Benchmark: Global Pairwise Alignment with Affine Gap Model for Protein.", compute, argc, argv);
}
