"""
some helper functions to split simulated data
"""

import numpy as np

def first_half_of_chromosomes(pos, debug=False):
    """
    given position structure, return
    the indices of first half of chromosomes
    """

    num_snps = len(pos)
    all_snp_idx = np.array(range(num_snps))
    unique_chrom = list(set(pos[:,0]))

    total_snps = 0

    keeper_idx = []

    for chrom in unique_chrom:
        chrom_mask = (pos[:,0] == chrom)
        num_snps_chrom = sum(chrom_mask)
        total_snps += num_snps_chrom
        assert num_snps_chrom > 0

        # keep first half on each chromosome
        keeper_idx_chr = all_snp_idx[chrom_mask][0:(num_snps_chrom/2)]
        keeper_idx.extend(keeper_idx_chr)

    assert total_snps == num_snps
    assert 0.4 * num_snps < len(keeper_idx) < 0.6 * num_snps

    if debug:
        import pylab
        pylab.figure()
        #pylab.plot(pos[keeper_idx][:,2], pos[keeper_idx][:,0], "x")
        pylab.plot(keeper_idx, pos[keeper_idx][:,0], "x")
        pylab.show()

    return keeper_idx


def all_except_first_chr(pos, debug=False):
    """
    given position structure, return
    the indices of first half of chromosomes
    """

    num_snps = len(pos)
    all_snp_idx = np.array(range(num_snps))
    unique_chrom = list(set(pos[:,0]))

    except_chr = 1
    assert except_chr in unique_chrom
    keeper_mask = (pos[:,0] != except_chr)
    except_mask = (pos[:,0] == except_chr)

    keeper_idx = all_snp_idx[keeper_mask]
    except_idx = all_snp_idx[except_mask]

    assert len(keeper_idx) + len(except_idx) == len(all_snp_idx)

    if debug:
        import pylab
        pylab.figure()
        #pylab.plot(pos[keeper_idx][:,2], pos[keeper_idx][:,0], "x")
        pylab.plot(keeper_idx, pos[keeper_idx][:,0], "x")
        pylab.show()

    return keeper_idx, except_idx
    

def split_chr1_chr2_rest(pos, debug=False):
    """
    split position into three sets:
    chr1_idx --> used for t1 evaluation
    chr2_idx --> used for power evaluation
    rest_idx --> all other 
    """

    pos.flags.writeable = False


    num_snps = len(pos)
    all_snp_idx = np.array(range(num_snps))
    unique_chrom = list(set(pos[:,0]))

    assert 1 in unique_chrom
    assert 2 in unique_chrom
    chr1_mask = (pos[:,0] == 1)
    chr2_mask = (pos[:,0] == 2)
    rest_mask = np.bitwise_and(~chr1_mask, ~chr2_mask)
    

    chr1_idx = all_snp_idx[chr1_mask]
    chr2_idx = all_snp_idx[chr2_mask]
    rest_idx = all_snp_idx[rest_mask]
    
    chr1_idx.flags.writeable = False
    chr2_idx.flags.writeable = False
    rest_idx.flags.writeable = False

    assert len(chr1_idx) + len(chr2_idx) + len(rest_idx) == len(all_snp_idx)
    assert len(chr1_idx) > 0
    assert len(chr2_idx) > 0
    assert len(rest_idx) > 0

    if debug:
        import pylab
        pylab.figure()
        #pylab.plot(pos[keeper_idx][:,2], pos[keeper_idx][:,0], "x")
        pylab.plot(chr1_idx, pos[chr1_idx][:,0], "x")
        pylab.plot(chr2_idx, pos[chr2_idx][:,0], "x")
        pylab.plot(rest_idx, pos[rest_idx][:,0], "x")
        pylab.show()

    assert max(chr1_idx) < min(rest_idx)
    assert max(chr2_idx) < min(rest_idx)
    assert max(chr1_idx) == min(chr2_idx) - 1
    assert max(chr2_idx) == min(rest_idx) - 1

    return chr1_idx, chr2_idx, rest_idx
