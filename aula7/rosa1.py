from Bio.Seq import Seq

seq = Seq("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")

print(str(seq.count("A")) +" " + str(seq.count("C")) +" " + str(seq.count("G"))+" " + str(seq.count("T")))

