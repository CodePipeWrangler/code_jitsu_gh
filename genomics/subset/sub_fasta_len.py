from Bio import SeqIO

input_fasta = "input.fasta"  # Replace with your actual input file
output_fasta = "output.fasta"  # The file to save trimmed sequences

start, end = 60, 213  # Define the range (1-based indexing)

with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        trimmed_seq = record.seq[start-1:end]  # Extract the range
        record.seq = trimmed_seq
        SeqIO.write(record, output_handle, "fasta")

print(f"Trimmed sequences saved to {output_fasta}")

