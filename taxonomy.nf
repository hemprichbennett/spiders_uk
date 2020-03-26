/* the root directory of the project */
base = file('./')

/* because onedrive is awful, we need to escape any spaces in filepaths */
String base = base
basedir = base.replace(" ", "\\ ")
basedir = file(basedir)

SeqFile_channel  = Channel
  .fromPath('data/processed_data/interactions/**/all_asvs.fasta')
  .collectFile(name: 'bigfa', newLine: false)


  bigEdgeList_channel = Channel.fromPath('data/processed_data/interactions/**/simple_edgelist.csv')
    .collectFile(name: 'bigedge', keepHeader: true)


process classify_sequences{

  input:
  file bigfa from SeqFile_channel

  output:
  file classified_sequences into classified_to_plot, classified_to_query, classified_for_edgelist

  shell:
  """

  java -Xmx8g -jar ~/RDPTools/classifier.jar classify \
    -t $basedir/data/reference_library/rRNAClassifier.properties \
    -o classified_sequences $bigfa
  cp classified_sequences $basedir/data/processed_data/classifications.tsv
  """
}



process plot_assignments{

  input:
  file classified_sequences from classified_to_plot

  shell:
  """
  Rscript $basedir/scripts/asv_assignment_confidences.R $classified_sequences $basedir
  """
}

process query_seqlibrary{
  input:
  file classified_sequences from classified_to_query

  shell:
  """
  Rscript $basedir/scripts/sp_isolating.R $classified_sequences splist
  python $basedir/scripts/sequence_matching.py splist $basedir/data/all_arthropods.fas queried_longreads
  cp queried_longreads $basedir/data/processed_data/queried_longreads.fa
  """
}



process combine_edgelist_and_assignments{
  input:
  file asv_edgeList from bigEdgeList_channel
  file classified_sequences from classified_for_edgelist

  shell:
  """
  Rscript $basedir/scripts/edgelist_taxonomy.R $asv_edgeList $classified_sequences $basedir
  """
}
