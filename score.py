#!/usr/bin/env python3
"""
A module to calculate the deleterious score of an SNP
"""

__author__ = "Jorick Baron"


def score(alignments, pos, key):
    """
    :param alignments: The alignment object of all the sequences
    :param pos: The amino acid position of the mutation
    :param key:The key of the original sequence
    :return: the score
    """
    count = {}  # stores the amount of times a certain amino acid is found
    for alignment in alignments:
        if not alignment.name == key:
            if alignment.seq[pos] in count:  # if the aa already exists in the dict
                count[alignment.seq[pos]] += 1
            else:
                count[alignment.seq[pos]] = 1  # first encounter
        else:
            mutant = alignment.seq[pos]
    if mutant not in count.keys():
        count[mutant] = 1
    unique = len(count.keys()) - 1
    return 1 + 9*(unique/count[mutant])
