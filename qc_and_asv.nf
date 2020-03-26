params.indir = ''

toSearch = params.indir + '**/*_R{1,2}*.fastq.gz'

input_files = Channel
  .fromFilePairs(toSearch)

/* the root directory of the project */
base = file('./')

/* because onedrive is awful, we need to escape any spaces in filepaths */
String base = base
basedir = base.replace(" ", "\\ ")
basedir = file(basedir)

process cutadapt_channel {
  input:
  set sampleID, file(infiles) from input_files

  output:
  set sampleID, file('f_clipped'), file('r_clipped') into clipped_channel

  shell:
  """
  # separate the infiles variable, into the forward and reverse
  fwd_file=\$(echo $infiles | sed 's/ .*//')
  rev_file=\$(echo $infiles | sed 's/.* //')
  echo \$fwd_file
  primerpair=\$(echo $sampleID | sed 's/-.*//')
  fwd_primer=\$(sed '1q;d' $basedir/data/primers/\$primerpair) #forward primer is the first line of the primer-file
  rev_primer=\$(sed '2q;d' $basedir/data/primers/\$primerpair) #reverse primer is the second line of the primer-file

  cutadapt -g \$fwd_primer -G \$rev_primer -o f_clipped -p r_clipped \$fwd_file \$rev_file --minimum-length 1

  """
}


process dada2_quality_filtering{
  input:
  set id, file(f_clipped), file(r_clipped) from clipped_channel

  output:
  file "${id}_f" into filtered_fwds
  file "${id}_r" into filtered_revs

  shell:
  """
  Rscript $basedir/scripts/dada2_filtering.R $f_clipped $r_clipped ${id}_f ${id}_r
  """
}

process input_dir_for_asvs{

  input:
  file f_file from filtered_fwds.collect()
  file r_file from filtered_revs.collect()

  output:
  file f_dir into f_dir_channel
  file r_dir into r_dir_channel

  shell:
  """
  mkdir -p f_dir
  mv $f_file f_dir
  mkdir -p r_dir
  mv $r_file r_dir
  """
}


process dada2_merge{
  input:
  file f_dir from f_dir_channel
  file r_dir from r_dir_channel

  shell:
  """
  Rscript $basedir/scripts/dada2_errors_and_merge.R $f_dir $r_dir $params.indir $basedir
  """
}
