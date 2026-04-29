from Bio import SeqIO
from io import StringIO

fastq_data = """@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!*((((***+))%%%++)(%%%%).1***-+*****))**55CCF>>>>>>CCCCCCC65
"""

input_handle = StringIO(fastq_data)
output_handle = StringIO()

count = SeqIO.write(
    SeqIO.parse(input_handle, "fastq"),
    output_handle,
    "fasta"
)

print(output_handle.getvalue())
print(f"Converted {count} records.")