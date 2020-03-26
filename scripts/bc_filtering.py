# Script formats the outputs from bcdatabaser, to get it into the correct
# style for RDP classifier.

# Setup, open files
import sys
import re

del sys.argv[0]

if not len(sys.argv) == 4:
    sys.exit("Error, please provide four filepaths: first the input fasta\
     name, then the output fasta file,\
      then the path for the output taxonomy file,\
      and finally the path for where bad sequences go to die")

inFile = sys.argv[0]
# Make a list of all the species whilst reading in the file
inFasta = open(inFile, mode='r')

outFastaPath = sys.argv[1]
outFasta = open(outFastaPath, mode='w')

outTaxaPath = sys.argv[2]
outTaxa = open(outTaxaPath, mode='w')

outTaxa.write('0*Root*-1*0*rootrank\n')

badTaxaPath = sys.argv[3]
badTaxa = open(badTaxaPath, mode='w')

# Pattern for finding accessions
accessionPattern = re.compile(">.+")

# patterns for checking sequence has full taxonomy
taxonomicPatterns = ['k:', ',p:', ',c:', ',o:', ',f:', ',g:', ',s:']

# things which show a non-accepted species
badSpecies = ['sp.', 'nr.', 'aff.', 'cf.', '/']


#  iterator for when making the taxonomy file
i = 1

# Dictionary which we populate below in the loop, to track taxonomy for each
# taxonomic group
taxonomy_dict = {}

# a second, simpler dictionary which is used to track the placement of each
# taxonomic level, to provide the final number needed in the taxonomy file
level_dict = {}


# loop through the fasta file, returning the sequence if needed
for line in inFasta:
    line = line.rstrip()
    if accessionPattern.match(line):
        ID_line = line
        ID = line.split(';')[0]  # May not actually be a species
        # check the sequence has full taxonomy
        if all(x in ID_line for x in taxonomicPatterns):

            if any(x in ID_line for x in badSpecies):

                badTaxa.write(ID_line + '\n')

                continue
            else:
                # write the sequences to filtered fasta
                # the ID line will need some formatting before writing, # to
                # allow it to work with RDP
                seq = next(inFasta)
                outFasta.write(ID_line + '\n' + seq)

                #  Now get the taxonomic information from ID_line

                # first, for kingdom. This ones simple to do

                # basic regex to isolate the kingdom from ID_line
                tax_kingdom = re.sub('.+k\\:', '', ID_line)
                tax_kingdom = re.sub(',.+', '', tax_kingdom)

                # add kingdom to level_dict, if its not there already
                if('kingdom' not in level_dict):
                    level_dict['kingdom'] = str(i)

                # if its not in the big dictionary of taxonomy, create a value
                # for that kingdom (it's linenumber) and a dictionary for all
                # taxa that are in that kingdom
                if tax_kingdom not in taxonomy_dict:
                    taxonomy_dict[tax_kingdom] = {}
                    taxonomy_dict[tax_kingdom]['number'] = i
                    # the writing of taxonomy for kingdom is simpler, as there
                    # will presumably only be one kingdom per reference database
                    outTaxa.write(str(i) + '*' + tax_kingdom + '*0*1*kingdom\n')

                    i += 1

                # now do the same thing for phylum
                tax_phylum = re.sub('.+p\\:', '', ID_line)
                tax_phylum = re.sub(',.+', '', tax_phylum)

                if('phylum' not in level_dict):
                    level_dict['phylum'] = str(i)
                if tax_phylum not in taxonomy_dict[tax_kingdom]:
                    taxonomy_dict[tax_kingdom][tax_phylum] = {}
                    taxonomy_dict[tax_kingdom][tax_phylum]['number'] = i

                    outTaxa.write(str(i) + '*' + tax_phylum + '*' +
                        str(taxonomy_dict[tax_kingdom]['number']) + '*' +
                        level_dict['phylum'] + '*phylum\n')

                    i += 1


                tax_class = re.sub('.+c\\:', '', ID_line)
                tax_class = re.sub(',.+', '', tax_class)
                #print('tax_class is ' + tax_class)

                if('class' not in level_dict):
                    level_dict['class'] = str(i)

                if tax_class not in taxonomy_dict[tax_kingdom][tax_phylum]:
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class] = {}
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class]['number'] = i

                    outTaxa.write(str(i) + '*' + tax_class + '*' +
                        str(taxonomy_dict[tax_kingdom][tax_phylum]['number'])
                        + '*' + level_dict['class'] + '*class\n')

                    i += 1


                tax_order = re.sub('.+o\\:', '', ID_line)
                tax_order = re.sub(',.+', '', tax_order)

                if('order' not in level_dict):
                    level_dict['order'] = str(i)

                if tax_order not in taxonomy_dict[tax_kingdom][tax_phylum][tax_class]:
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order] = {}
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order]['number'] = i

                    outTaxa.write(str(i) + '*' + tax_order + '*' +
                        str(taxonomy_dict[tax_kingdom][tax_phylum][tax_class]['number'])
                        + '*' + level_dict['order'] + '*order\n')

                    i += 1

                tax_family = re.sub('.+f\\:', '', ID_line)
                tax_family = re.sub(',.+', '', tax_family)

                if('family' not in level_dict):
                    level_dict['family'] = str(i)

                if tax_family not in taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order]:
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family] = {}
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family]['number'] = i

                    outTaxa.write(str(i) + '*' + tax_family + '*' +
                        str(taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order]['number'])
                        + '*' + level_dict['family'] + '*family\n')

                    i += 1

                if('genus' not in level_dict):
                    level_dict['genus'] = str(i)

                tax_genus = re.sub('.+g\\:', '', ID_line)
                tax_genus = re.sub(',.+', '', tax_genus)
                if tax_genus not in taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family]:
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family][tax_genus] = {}
                    taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family][tax_genus]['number'] = i

                    outTaxa.write(str(i) + '*' + tax_genus + '*' +
                        str(taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family]['number'])
                        + '*' + level_dict['genus'] + '*genus\n')

                    i += 1

                if('species' not in level_dict):
                    level_dict['species'] = str(i)
                tax_species = re.sub('.+s\\:', '', ID_line)
                tax_species = re.sub(',.+', '', tax_species)

                if tax_species not in taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family][tax_genus]:

                    outTaxa.write(str(i) + '*' + tax_species + '*'
                        + str(taxonomy_dict[tax_kingdom][tax_phylum][tax_class][tax_order][tax_family][tax_genus]['number'])
                        + '*' + level_dict['species'] + '*species\n')

                    i += 1

        else:
            badTaxa.write(ID_line + '\n')
            #  add writing out badseqs step here

#  print(taxonomy_dict)
