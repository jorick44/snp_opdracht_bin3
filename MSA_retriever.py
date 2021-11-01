#!/usr/bin/env python3
"""
makes unaligned fastas and aligns them using Clustalo
"""

__author__ = "Jorick Baron"

import sys
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import subprocess
import os

# globals
in_file = "data/unalinged.fasta"
out_file = "data/alinged.fasta"


def align():
    """
    Takes the unaligned file and aligns them using clustalo and stores it in an alignment object
    :return: an alignment object
    """
    with open(out_file, "w") as FILE:
        print("Beginning alignment...")
        # makes a commandline argument to use Clustalo
        c_line = ClustalOmegaCommandline(infile=in_file, outfile=out_file, force=True)
        os.system(str(c_line))  # executes Clustalo
        print("Finished alignment...")
    with open(out_file) as FILE:
        alignment = AlignIO.read(FILE, "fasta")
    return alignment


def mk_file(genes):
    """
    makes the required files and fills the unaligned file with the found genes
    """
    if not os.path.isfile(out_file):
        subprocess.run(["touch", out_file])  # make sure the alignment file exists
    print("making unaligned fasta...")
    with open(in_file, "w") as FILE:
        for gene in genes:
            print(f">{gene}\n{genes[gene].replace('-', '')}", file=FILE)  # create the unaligned file
