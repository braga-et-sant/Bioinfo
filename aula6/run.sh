conda create -n bioinfo_lab6 -c conda-forge -c bioconda \
samtools minimap2 bowtie2 blast emboss seqkit python=3.11 -y

needle -asequence "asis:HEAGAWGHEE" -bsequence "asis:PAWHEAE" -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -outfile needle_out.txt -auto



needle -asequence "asis:HEAGAWGHEE" -bsequence "asis:PAWHEAE" -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -outfile needle_out.txt -auto


#4.3 Biological significance
 #
 #The global alignment score is much lower because Needleman–Wunsch forces the entire lengths of both sequences to align, including poorly matching ends and gaps.
 #
 #The local alignment score is higher because Smith–Waterman finds only the best matching subsequence and ignores unrelated flanking regions.
 #
 #Meaning biologically
 #A low global score suggests the two full proteins are not strongly similar over their entire length.
 #A higher local score suggests they may still share a conserved motif, domain, or short functional region.
 #
 #So the difference tells you these sequences likely have local similarity, but not strong overall/global similarity.

