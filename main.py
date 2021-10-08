#!/usr/bin/env python3
"""
Module documentation.
"""

__author__ = "Jorick Baron"

# Imports
import sys
import argparse
import Bio


# Global variables

# Class declarations

# Function declarations

def main():
    """
    Main function documentation
    """
    pass


# Main body
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given an SNP and a genetic sequence\n"
                                                 "predict how deleterious it is where 1 is harmless and 10 is very bad")
    parser.add_argument(metavar="S", help="the SNP you wish to calculate the effect of")
    parser.add_argument(metavar="P", help="The position of the SNP", type=int)
    parser.add_argument(metavar="M", help="a MSA file containing the sequences you compare the SNP to")

    sys.exit(main())
