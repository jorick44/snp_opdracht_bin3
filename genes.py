#!/usr/bin/env python3
"""
Processes the input from fasta files into a gene object and blasts for similar genes
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
        # parses the fasta file for the sequence and saves it in a record object
        records = SeqIO.parse(open(file), format="fasta")
        for record in records:
            self.genes[record.name] = record.seq
            key = record.name  # this accession number can be useful later to recognise the original sequence
            return key

    def retrieve_genes(self, key):
        """
        takes the given sequence and finds related genes trough a blast
        :param key: the key of the gene you wish to blast for
        """
        print("Begin blasting...")
        # does an online blastx search to find closely related sequences
        blast_handle = NCBIWWW.qblast("blastx", "swissprot", self.genes.get(key), ncbi_gi=False, hitlist_size=15,
                                      format_type="XML")
        print("Done blasting...")
        result_handle = NCBIXML.parse(blast_handle)
        return result_handle

    def parse_blast(self, results):
        """
        parses the results of the blast and adds them to the gene dictionary
        :param results: the results of the blast
        """
        for result in results:
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    self.genes[alignment.accession] = hsp.query  # uses the accession as a key and the query as value

    def mutate(self, pos, mutant, key):
        """
        Alters the gene according to the SNP
        :param pos: position of mutation
        :param mutant: the different nucleotide
        :param key: accession number of the gene to be mutated
        """
        original = Seq(str(self.genes.get(key)))
        proper = original[original.find("ATG"):]  # find the start codon
        proper = MutableSeq(str(proper))  # allows mutation
        proper[pos] = mutant  # mutate
        proper = Seq(str(proper))  # MutableSeq does not support translation
        self.genes[key] = str(proper.translate(to_stop=True))  # translate until the first stop codon


if __name__ == '__main__':
    test_file = "test_files/CFRT_DNA_sequence.fasta"
    sys.exit(main())
