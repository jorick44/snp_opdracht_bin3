#!/usr/bin/env python3
"""
reads a FASTA file for the SNP project
"""

__author__ = "Jorick Baron"

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq


def main():
    """
    tests the ability to read the fasta
    """
    genes = Genes()
    print(genes.read_fasta(test_file))
    return 0


class Genes:
    def __init__(self):
        self.genes = dict()

    def read_fasta(self, file):
        """
        reads a FASTA file and puts it in a dictionary
        :param file: a FASTA file
        :return: the key to access this dictonary entry
        """
        records = SeqIO.parse(open(file), format="fasta")
        for record in records:
            self.genes[record.name] = record.seq
            key = record.name
        return key

    def retrieve_genes(self, key):
        """
        :param key: the key of the gene you wish to blast for
        takes the given sequence and finds related genes trough a blast
        """
        print("Begin blasting...")
        blast_handle = NCBIWWW.qblast("blastx", "swissprot", self.genes.get(key), ncbi_gi=False, hitlist_size=15,
                                          format_type="XML")
        print("Done blasting...")
        result_handle = NCBIXML.parse(blast_handle)
        return result_handle

    def parse_blast(self, results):
        """
        parses the results of the blast and adds them to the gene dictionary
        """
        for result in results:
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    self.genes[alignment.accession] = hsp.query

    def mutate(self, pos, mutant, key):
        """

        """
        original = Seq(str(self.genes.get(key)))
        proper = original[original.find("ATG"):]
        proper = MutableSeq(str(proper))
        proper[pos] = mutant
        proper = Seq(str(proper))
        self.genes[key] = str(proper.translate(to_stop=True))


if __name__ == '__main__':
    test_file = "test_files/CFRT_DNA_sequence.fasta"
    sys.exit(main())
