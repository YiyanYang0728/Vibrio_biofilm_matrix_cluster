from Bio import SeqIO
import sys

def filter_sequences(input_fasta, output_fasta):
    # Allowing standard amino acids plus ambiguous characters B, J, O, U, Z, X
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYBJOUZX")
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = record.seq.upper()
            if set(sequence).issubset(valid_amino_acids):
                SeqIO.write(record, output_handle, "fasta")

# Replace 'input.fasta' with your input fasta file path
# Replace 'output.fasta' with your desired output file path
input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

filter_sequences(input_fasta, output_fasta)