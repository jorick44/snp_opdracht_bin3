#!/usr/bin/env python3
"""
Module documentation.
"""

__author__ = "Jorick Baron"

# Imports
import sys
import argparse

import MSA_retriever
import score
from genes import Genes


# Global variables

# Class declarations

# Function declarations
def main():
    """
    Main function documentation
    """
    genes = Genes()
    key = genes.read_fasta(args.sequence)
    blasted = genes.retrieve_genes(key)
    genes.parse_blast(blasted)
    genes.mutate(args.position, args.snp, key)
    MSA_retriever.mk_file(genes.genes)
    alignment = MSA_retriever.align()
    out = score.score(alignment, args.position//3, key)
    return f"effect of the SNP is {out}"


# Main body
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given an SNP and a genetic sequence predict how "
                                                 "deleterious it is where 1 is harmless and 10 is very bad")
    parser.add_argument("-s", "--snp", help="The SNP you wish to calculate the effect of")
    parser.add_argument("-f", "--sequence",
                        help="A FASTA file of the sequence you wish to calculate the SNP's effect on")
    parser.add_argument("-p", "--position", help="The position of the SNP", type=int)
    args = parser.parse_args()

    sys.exit(main())
