#!/usr/bin/env python
"""Convert from the miRDeep2 BED output file to SAF for quantifying expression
with featureCounts. Takes one argument:
    1) miRDeep2 BED file"""

import sys

try:
    bed_in = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def main(bed):
    """Main function."""
    with open(bed, 'rt') as f:
        for line in f:
            tmp = line.strip().split('\t')
            # Unpack the list for the fields we want:
            #   0: chromosome
            #   1: Start (0-based)
            #   2: End
            #   3: feature ID
            #   5: strand
            # SAF is printed in this column order:
            #   0: feature ID
            #   1: chromosome
            #   2: Start (1-based)
            #   3: End
            #   4: Strand
            chrom = tmp[0]
            start = int(tmp[1])
            end = tmp[2]
            feat_id = tmp[3]
            strand = tmp[5]
            toprint = [
                feat_id,
                chrom,
                str(start+1),
                end,
                strand]
            print('\t'.join(toprint))
    return


main(bed_in)
