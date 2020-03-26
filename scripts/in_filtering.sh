# input directory
unfiltered=$1
echo $unfiltered

date_str="$(echo $unfiltered | sed 's/.*uncurated\///')"
echo $date_str

# the unfiltered files
files=$unfiltered/*/*sequences.tax.fa

# Combine all those files into a single one
cat $files > $unfiltered/concatenated.tax.fa

# now remove the newlines in the sequences, as otherwise the python script only
# gives one sequence file
# (code from https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $unfiltered/concatenated.tax.fa > $unfiltered/better_concatenated.tax.fa

outdir=data/reference_library/curated/$date_str
mkdir $outdir

python scripts/bc_filtering.py $unfiltered/better_concatenated.tax.fa $outdir/arthropods.fa $outdir/arthropods.tax $outdir/bad.fa
