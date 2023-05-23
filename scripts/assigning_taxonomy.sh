# script uses classifiers made in the CRAGinformatics pipeline (separate project)
# to assing taxonomy to our ASVs

qiime tools import --input-path data/processed_data/interactions/spiders_uk/all_asvs.fasta \
--output-path data/processed_data/interactions/spiders_uk/all_asvs.qza  \
--type 'FeatureData[Sequence]'


qiime feature-classifier classify-sklearn \
  --i-classifier data/reference_library/2023-04/anml_classifier.qza \
  --i-reads data/processed_data/interactions/spiders_uk/all_asvs.qza \
  --o-classification data/processed_data/classifications/anml_classifications.qza
  
qiime feature-classifier classify-sklearn \
  --i-classifier data/reference_library/2023-04/zbj_classifier.qza \
  --i-reads data/processed_data/interactions/spiders_uk/all_asvs.qza \
  --o-classification data/processed_data/classifications/zbj_classifications.qza
  
qiime rescript merge-taxa \
	--i-data data/processed_data/classifications/zbj_classifications.qza \
	data/processed_data/classifications/anml_classifications.qza \
	--o-merged-data data/processed_data/classifications/both_classifications.qza \
	--p-mode score
