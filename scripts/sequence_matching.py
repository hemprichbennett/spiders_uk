# Script matches the names of species-level assignments from RDP classifier to
# long sequences that were manually downloaded from BOLD

# Setup, open files
import sys
import re

del sys.argv[0]

if not len(sys.argv) == 3:
    sys.exit("Error, please provide three filepaths: first the input species\
     name, then the input fasta file,\
      and finally the path for the output file")

spFile = sys.argv[0]
# Make a list of all the species whilst reading in the file
spNames = open(spFile).readlines()

fastaFile = sys.argv[1]
fastaData = open(fastaFile, mode='r')

outFileName = sys.argv[2]
outFile = open(outFileName, mode='w')

# Replace all of the newline characters at the end of spNames
spNames = [re.sub("\n", "", x) for x in spNames]
print(spNames)

accessionPattern = re.compile(">.+")

# loop through the fasta file, returning the sequence if needed
for line in fastaData:
    if accessionPattern.match(line):
        ID = line
        species = line.split('|')[1]  # May not actually be a species
        if species in spNames:  # If its a species matching one of our inputs
            seq = next(fastaData)
            outFile.write('>' + species + '\n' + seq + '\n')
