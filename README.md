# SNP_anotation
Calculates the effect of a SNP and rates it from 1(harmless) to 10(harmful)

## requirements
Requires connection to the internet and the clustalo alignment tool

## usage 
Some examples:  
for a low score: ``python3 main.py -f testfiles/CFRT_DNA_sequence -p 108 -s A``  
for a high score: ``python3 main.py -f test_files/CFRT_DNA_sequence.fasta -p 1 -s G``  
-f is for a FASTA file of the sequence you wish to calculate the SNP's effect on.  
-p is for the position of the SNP.  
-s is for the SNP you wish to calculate the effect of.

Note: blasting might take a while.

## contact
For any questions mail me at:  
j.baron@st.hanze.nl
