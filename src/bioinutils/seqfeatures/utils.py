import re

import edlib
from Bio import Seq
from primer3.bindings import (calc_hairpin, calc_heterodimer, calc_homodimer,
                              calcEndStability, calcTm)


def edlib_pw_align(seq1: str, seq2: str, mode: str = "HW"):

    # compute the reverse complement of the sequence and align
    rc_seq1 = Seq(seq1).__str__()
    out = edlib.align(rc_seq1, seq2, mode=mode, task="locations")
    # parse
    ed = out["editDistance"]
    sl = len(rc_seq1)
    aln_score = float((sl - ed) / sl)
    aln_start = out["locations"][0][0]
    aln_end = out["locations"][0][1]

    return aln_score, aln_start, aln_end


def homopolymer_counter(seq: str, limit: int = 4) -> int:
    pat = f"A{{{limit},}}|T{{{limit},}}|G{{{limit},}}|C{{{limit},}}|U{{{limit},}}"
    r = re.findall(pat, seq.upper())
    hit_num = len(r)
    correction = hit_num * (limit - 1)
    score = sum([len(h) for h in r]) - correction

    return -1 * score


def count_ambiguous_letters(seq: str) -> int:
    ambig = ["N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"]
    count = sum([seq.count(x) for x in ambig])

    return count


def calculate_ambiguous_percentage(seq: str) -> float:
    seq_len = len(seq)
    seq = seq.upper()
    ambig_count = count_ambiguous_letters(seq)
    perc = ambig_count / (seq_len + 1) * 100

    return perc


def calculate_lc(seq: str, k: int = 4) -> float:
    """Function to calculate the linguistic complexity of a sequence.

    :param seq: nucleotide sequence
    :type seq: str
    :param k: ksize, defaults to 4
    :type k: int, optional
    :return: linguistic complexity on a scale from 0 to 1
    :rtype: float
    """

    seq = seq.upper().replace("U", "T")
    max_denom = 4**k
    n = len(seq)
    denom = n - (k - 1)
    if denom > max_denom:
        denom = max_denom
    unique_kmers = len(set([seq[ind: ind + k] for ind in range(n - (k - 1))]))

    score = unique_kmers / denom

    return score


def calc_tm(seq: str):
    """calculates the sequence TM according to an average human cell

    :param seq:   seq
    :type seq: str
    :return: melting temp
    :rtype: int
    """
    tm = calcTm(
        seq, mv_conc=200, dv_conc=50, dntp_conc=50, dna_conc=50, max_nn_length=60
    )
    return tm


def calc_hetero(a: str, b: str):
    """calculate the Delta G for two sequences heterodimerization.
    The output returns in kcal/mol

    :param a: sequence A
    :type a: str
    :param b: sequence B
    :type b: str
    :return: change in gibbs free energy kcal/mol
    :rtype: int
    """
    dg = (
        calc_heterodimer(
            a, b, mv_conc=200, dv_conc=50, dntp_conc=50, dna_conc=50, temp_c=37
        ).dg
        / 1000
    )
    return dg


def calc_homo(seq: str):
    """calculate the homodimerization delta G

    :param seq: primer
    :type seq: str
    :return: change in gibbs free energy kcal/mol
    :rtype: int
    """
    dg = (
        calc_homodimer(
            seq, mv_conc=200, dv_conc=50, dntp_conc=50, dna_conc=50, temp_c=37
        ).dg
        / 1000
    )
    return dg


def calc_3_prime_stability(seq: str, signature: str):
    """Calculate the 3' stability of a primer compared to the surrounding seq

    :param seq: primer seq
    :type seq: str
    :param signature: surrounding region seq
    :type signature: str
    :return: change in gibbs free energy kcal/mol
    :rtype: int
    """
    dg = (
        calcEndStability(
            seq,
            signature,
            mv_conc=200,
            dv_conc=100,
            dntp_conc=50,
            dna_conc=50,
            temp_c=37,
        ).dg
        / 1000
    )

    return dg

def calc_hairpin_dg(seq: str):
    """calculate the homodimerization delta G

    :param seq: primer
    :type seq: str
    :return: change in gibbs free energy kcal/mol
    :rtype: int
    """
    dg = (
        calc_hairpin(
            seq, mv_conc=200, dv_conc=50, dntp_conc=50, dna_conc=50, temp_c=37
        ).dg
        / 1000
    )
    return dg
