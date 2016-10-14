#include "pairwise_alignment_base.h"

using namespace seqan;

int main(int argc, char **argv)
{
    typedef String<AminoAcid> TSequence;
    typedef StringSet<Gaps<TSequence>> TGapSet;

    auto compute = [] (TGapSet & gapsSetH, TGapSet & gapsSetV){
        // Blosum62, gap: -1
        return localAlignment(gapsSetH, gapsSetV, Blosum62(), NeedlemanWunsch());
    };
    return benchmark_pairwise_alignment_main<TSequence>("Benchmark: Local Pairwise Alignment with Linear Gap Model for Protein.", compute, argc, argv);
}
