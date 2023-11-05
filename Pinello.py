""" Briefly we aim to quantify the frequency of the top 10 deletions in sequencing reads from a PacBio run, using the human reference as a reference (hg38). 
 
This is the potential breakdown in steps:

1) Download this small FASTQ file from the NCBI SRA archive: https://www.ncbi.nlm.nih.gov/sra/SRX15881277[accn] 
This file is from this publication, if you want to learn more about the provenance.

2) Align the reads to the reference genome using minimap2.

3) Parse the CIGAR strings in the obtained BAM file to report the top 10 deletions, ranked by their frequency in descending order. """


import os
 
os.system("echo GeeksForGeeks")
#./minimap2 -c /Users/benvyshedskiy/Downloads/GCF_000001405.26/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna /Users/benvyshedskiy/Downloads/SRR19837561.fastq.g >mayvework.paf

import re

# Open the PAF file
paf_file = open("/Users/benvyshedskiy/Documents/GitHub/minimap2/mayvework.paf", "r")

# Initialize a dictionary to store the deletion lengths and their frequencies
deletions = {}

# Iterate over the PAF file and parse the CIGAR strings
for line in paf_file:
    # Split the PAF line into fields
    paf_fields = line.split("\t")

    if len(paf_fields)<23:
        continue
    # Get the CIGAR string from the PAF line
    cigar_string = paf_fields[22]
    querystart = paf_fields[2]
    currentbase = int(querystart)


    # Parse the CIGAR string
    cigar_operations = re.findall(r"(\d+)([MIDNSH])", cigar_string)

    # Calculate the total deletion length
    for operation in cigar_operations:
        if operation[1] != "D":
            currentbase += int(operation[0])
        if operation[1] == "D":
            for i in range(0,int(operation[0])):
                if currentbase in deletions:
                    deletions[currentbase] += 1
                else:
                    deletions[currentbase] = 1
                currentbase+=1

# Close the PAF file
paf_file.close()

# Sort the deletion lengths by frequency in descending order
sorted_deletion_lengths = sorted(deletions, key=deletions.get, reverse=True)

# Report the top 10 deletions
print("Top 10 deletions, ranked by frequency in descending order:")
for i in range(10):
    deletion_length = sorted_deletion_lengths[i]
    deletion_frequency = deletions[deletion_length]
    print(f"position {deletion_length} deleted: {deletion_frequency} times")


