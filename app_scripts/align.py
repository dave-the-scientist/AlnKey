from app_scripts.BLOSUM import get_matrix
import itertools
from math import sqrt
from collections import Counter

# Implement a sliding window for the variation calculation.
# Implement a threaded version of the variation function. Not sure if I want this module to implement a threaded version, or just provide a version of the function(s) that can be used by a thread pool managed in the user's other code.
# Implement the quality function, even if i'm not using it for this
# I think I might be able to vectorize the identity functions. Especially for all-against-all calculations.
#  - ALREADY DONE! I basically did this in .../work/heme_slp_diversity/tom_abaumanii/calculate_alignment_identity.py and it is crazy fast (300 sequences, ~100k comparisons in ~1 sec)
#    - Incorporate that into this file. Make sure it's working as intended. The method is slightly different from what i ahve below, likely better, but test. That one will be more memory-intensive; check if it gets too big with alignments >1000 or >4000 sequences. Maybe the per-column approach below will be fast enough and more efficient.
#    - I like the various public methods defined there, use them. Keep the below functions as fallbacks in case numpy can't be imported
#  - maintain 2 nxn matrices, 1 for non-gap matches, 1 for total comparisons. take sequences, go through column at a time. numpy to turn the column vector x its transpose into matrices, add to maintained matrices.
#    - or maybe 1 matrix for matches, 1 for gap-gap matches; ident is then (matches-gaps) / (total_length-gaps). don't forget to remove the impact of self-self comparisons
# - Could likely make use of triangular objects to store data. they'll get big for large alignments.
# - Likely will need 1 function for 1x1 identity, 1 for 1xall, and 1 for allxall

# variation() takes ~30 sec for lenient_old.aln; vectorize it.

# For convenience, import some of the major functions from sequ into this namespace. Like the load functions; don't want user to have to import sequ just for that.


def consensus(seqs, allow_gaps=True):
    """Expects seqs is a list of aligned sequences as strings/Sequences or a SeqList. At each position returns the most frequent character, the first encountered if there is a tie. If allow_gaps is False, returns the most frequent non-gap character, with arbitrary tie-breaking. Will only return a gap if the entire column is gaps."""
    consensus = []
    for column in zip(*seqs):
        cntr = Counter(column)
        if allow_gaps:
            c = max(column, key=cntr.get)
        else:
            cs = cntr.most_common(2)
            c = cs[1][0] if len(cs)>1 and cs[0][0]=='-' else cs[0][0]
        consensus.append(c)
    return ''.join(consensus)

def variants(seqs, number=None):
    """Returns a list, where each column in the alignment generates a sublist of the form: [(c1, count1), (c2, count2), ...]. Here, c1 is the most common character at with count1 occurances. Ties are broken arbitrarily. `number` indicates the maximum number of top hits to return; returns all if None."""
    vnts = []
    for column in zip(*seqs):
        cntr = Counter(column)
        vnts.append(cntr.most_common(number))
    return vnts

# # #  Alignment identity
def identity(seq1, seq2):
    """Does not compute an alignment, and expects that the sequences are already aligned. `seq1` and `seq2` can be Sequence objects or strings. Computes the number of any non-gap characters in agreement divided by the total number of columns that are not only gaps; that is, when both sequences do have a gap at the same position, it does not affect the numerator or denominator of the identity calculation. Returned values are percentages from 0 - 100%."""
    if len(seq1) != len(seq2):
        raise MolecbioAlignmentLengthError("cannot compute the identity between two sequences of different lengths. They should be aligned before calling the identity() function.")
    matches, total = 0, 0
    for c1, c2 in zip(seq1, seq2):
        if c1 == '-' == c2:
            continue
        if c1 == c2:
            matches += 1
        total += 1
    return matches / total * 100.0
def identities(seqs1, seqs2=[], average=True):
    """Normally both `seqs1` and `seqs2` should be containers (SeqList or list) of sequences (Sequence objects or strings); however each argument can also be a single sequence outside of a container. The sequences need to already be aligned. If `seqs2` is empty, will return all pairwise identities within `seqs1` (skipping self & redundant comparisons). Otherwise will calculate identities between each sequence in `seqs1` and each sequence in `seqs2`, without checking for self or redundant comparisons. If `average` is False will return a list of all raw identities, otherwise will return the average."""
    idents = []
    if seqs2:
        if len(seqs1[0]) == 1: # if seqs1 is a single sequence instead of a container of sequences
            seqs1 = [seqs1]
        if len(seqs2[0]) == 1: # if seqs2 is a single sequence instead of a container of sequences
            seqs2 = [seqs2]
        for seq1, seq2 in itertools.product(seqs1, seqs2):
            idents.append(identity(seq1, seq2))
    elif len(seqs1[0]) != 1 and len(seqs1) > 1:
        for seq1, seq2 in itertools.combinations(seqs1, 2):
            idents.append(identity(seq1, seq2))
    else:
        raise ValueError("if only one group of sequences is given, it must contain more than 1 sequence to calculate percent identity.")
    if average:
        return sum(idents) / len(idents)
    else:
        return idents


# # #  Alignment variation
def variation(seqs, matrix='BLOSUM62', stdev=False):
    """
    An adaptation of the quality calculation originally described in ClustalX (http://www.clustal.org/download/clustalx_help.html), implemented as x(). The default uses absolute mean deviation instead of standard deviation, which is more robust to outliers plus the interpretation does not depend on the scaling factors used in the particular BLOSUM matrix used (doubling the values will change stdev, but not absmeandev; TEST THAT).
    The `matrix` argument allows other scoring matrices to be used. Must be one of: 'BLOSUM30', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', '30', '45', '50', '62', or '80'."""
    blosum = get_matrix(matrix)

    dev_fxn = column_std_deviation if stdev else column_absmean_deviation
    num_seqs = len(seqs)
    devs = []
    for col in zip(*seqs):
        # try some normalization, don't think i really want it though
        #dev = (dev_fxn(col, blosum) - bl_min) / bl_range
        dev = dev_fxn(col, blosum)
        devs.append(dev)
    return devs
def column_absmean_deviation(column, blosum):
    if len(set(column)) == 1:
        return 0.0
    dists, _ = column_dists_and_mean(column, blosum)
    return sum(dists) / len(dists)
def column_std_deviation(column, blosum):
    if len(set(column)) == 1:
        return 0.0
    dists, _ = column_dists_and_mean(column, blosum)
    return sqrt(sum(d*d for d in dists) / len(dists))
def column_dists_and_mean(column, blosum):
    """Returns the distances from the mean vector for each residue in `column`, and the mean vector itself."""
    col_vecs = [blosum[c] for c in column]
    mean_vec = [sum(dim_vec)/len(col_vecs) for dim_vec in zip(*col_vecs)]
    dists = [sqrt(sum((dim-mean)**2 for dim in dim_vec)) for dim_vec, mean in zip(zip(*col_vecs), mean_vec)]
    return dists, mean_vec


# # #  Errors
class MolecbioAlignmentLengthError(ValueError):
    """Raised when sequences in an alignment are not the same length."""


# # #  TESTING
if __name__ == '__main__':
    from app_scripts import sequ
    seqs = sequ.load('test_dir\confirmed.clustal')
    vars = variation(seqs)

    per_line = 15
    for b_ind in range(0, len(seqs[0]), per_line):
        c_buff, v_buff = [], []
        for i in range(per_line):
            ind = b_ind + i
            if ind==len(seqs[0]): break
            c_buff.append('{:5}'.format(seqs[0][ind]))
            v_buff.append('{:5}'.format(str(round(vars[ind], 2))))
        print('{} - {}'.format(b_ind, b_ind+per_line))
        print(''.join(c_buff))
        print(''.join(v_buff))
