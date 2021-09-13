# kernel bash
# conda info --envs
source activate qiime2-2020.6

cd ~/sync/lichens/bio_info/data_sequencages/2021_07

qiime tools import \
--input-path ./dada2/seqtabnochim_mosaique_16S_18S_ITS.fna \
--type 'FeatureData[Sequence]' \
--output-path ./qiime2/rep-seqs_mosaique_210710_16S_18S_ITS.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07

biom convert -i ./dada2/seqtabnochim_mosaique_16S_18S_ITS.txt -o ./qiime2/seqtabnochim_mosaique_16S_18S_ITS.biom --table-type="OTU table" --to-hdf5

qiime tools import \
--input-path ./qiime2/seqtabnochim_mosaique_16S_18S_ITS.biom \
--type 'FeatureTable[Frequency]' \
--output-path ./qiime2/table_mosaique_210710_16S_18S_ITS.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2
mkdir ./OTU-ITS2

qiime feature-classifier extract-reads \
  --i-sequences rep-seqs_mosaique_final_16S_18S_ITS.qza \
  --p-f-primer GTGAATCATCGAATCTTTGAA \
  --p-r-primer TCCTCCGCTTATTGATATGC \
  --p-identity 0.9 \
  --p-n-jobs 32 \
  --o-reads ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza
  
qiime feature-table filter-features \
  --i-table table_mosaique_210710_16S_18S_ITS.qza \
  --m-metadata-file ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza \
  --o-filtered-table ./OTU-ITS2/table_mosaique_210710_ITS.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2
mkdir ./OTU-16S

qiime feature-classifier extract-reads \
  --i-sequences rep-seqs_mosaique_210710_16S_18S_ITS.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer CCGYCAATTYMTTTRAGTTT \
  --p-identity 0.7 \
  --p-n-jobs 32 \
  --o-reads ./OTU-16S/rep-seqs_mosaique_210710_16S.qza

qiime feature-table filter-features \
  --i-table table_mosaique_210710_16S_18S_ITS.qza \
  --m-metadata-file ./OTU-16S/rep-seqs_mosaique_210710_16S.qza \
  --o-filtered-table ./OTU-16S/table_mosaique_210710_16S.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2
mkdir ./OTU-18S

qiime feature-classifier extract-reads \
  --i-sequences rep-seqs_mosaique_210710_16S_18S_ITS.qza \
  --p-f-primer CCCTGCCHTTTGTACACAC \
  --p-r-primer CCTTCYGCAGGTTCACCTAC \
  --p-identity 0.7 \
  --p-n-jobs 32 \
  --o-reads ./OTU-18S/rep-seqs_mosaique_210710_18S.qza

qiime feature-table filter-features \
  --i-table table_mosaique_210710_16S_18S_ITS.qza \
  --m-metadata-file ./OTU-18S/rep-seqs_mosaique_210710_18S.qza \
  --o-filtered-table ./OTU-18S/table_mosaique_210710_18S.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2

qiime vsearch cluster-features-de-novo \
  --i-table ./OTU-ITS2/table_mosaique_210710_ITS.qza \
  --i-sequences ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza \
  --p-perc-identity 0.97 \
  --p-threads 32 \
  --o-clustered-table ./OTU-ITS2/table-97_mosaique_210710_ITS.qza \
  --o-clustered-sequences ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza

qiime vsearch cluster-features-de-novo \
  --i-table ./OTU-16S/table_mosaique_210710_16S.qza \
  --i-sequences ./OTU-16S/rep-seqs_mosaique_210710_16S.qza \
  --p-perc-identity 0.97 \
  --p-threads 32 \
  --o-clustered-table ./OTU-16S/table-97_mosaique_210710_16S.qza \
  --o-clustered-sequences ./OTU-16S/rep-seqs-97_mosaique_210710_16S.qza
  
qiime vsearch cluster-features-de-novo \
  --i-table ./OTU-18S/table_mosaique_210710_18S.qza \
  --i-sequences ./OTU-18S/rep-seqs_mosaique_210710_18S.qza \
  --p-perc-identity 0.97 \
  --p-threads 32 \
  --o-clustered-table ./OTU-18S/table-97_mosaique_210710_18S.qza \
  --o-clustered-sequences ./OTU-18S/rep-seqs-97_mosaique_210710_18S.qza

cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2

# ITS
qiime tools export \
    --input-path ./OTU-ITS2/table-97_mosaique_210710_ITS.qza \
    --output-path ./OTU-ITS2/
    
biom convert \
    -i ./OTU-ITS2/feature-table.biom \
    -o ./OTU-ITS2/feature-table-97_mosaique_210710_ITS.tsv --to-tsv

# 16S
qiime tools export \
    --input-path ./OTU-16S/table-97_mosaique_210710_16S.qza \
    --output-path ./OTU-16S
    
biom convert \
    -i ./OTU-16S/feature-table.biom \
    -o ./OTU-16S/feature-table-97_mosaique_210710_16S.tsv --to-tsv

# 18S
qiime tools export \
    --input-path ./OTU-18S/table-97_mosaique_210710_18S.qza \
    --output-path ./OTU-18S
    
biom convert \
    -i ./OTU-18S/feature-table.biom \
    -o ./OTU-18S/feature-table-97_mosaique_210710_18S.tsv --to-tsv

########### ITS2
# BLAST+ Unite
cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2

qiime feature-classifier classify-consensus-blast \
    --i-query ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza \
    --i-reference-reads ~/sync/references_metabarcoding/references_taxons/sh_refs_qiime_ver8_97_s_04.02.2020.fasta.qza \
    --i-reference-taxonomy  ~/sync/references_metabarcoding/references_taxons/sh_taxonomy_qiime_ver8_97_s_04.02.2020.qza \
    --o-classification ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza

# exporter en .biom
qiime tools export \
    --input-path ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza \
    --output-path ./OTU-ITS2/Unite

########### ITS2
# BLAST+ PLANiTS
cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2

qiime feature-classifier classify-sklearn \
   --i-reads ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza \
   --i-classifier ~/sync/references_metabarcoding/ITS2/PLANiTS_29-03-2020/ITS2_PLANiTS_20-03-2020_classifier_q2020-6.qza \
   --p-n-jobs -1 \
   --o-classification ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza
   
# exporter en .biom
qiime tools export \
    --input-path ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza \
    --output-path ./OTU-ITS2/PLANiTS

########### PROKARYOTES
# source activate qiime2-2020.6
# procédure de curation des références 16S
# import and training 16S classifiers (SILVA_138)
# https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

# source activate qiime2-2020.6
cd /home/tony/sync/mangroves/references_taxons/SILVA_138_QIIME_release

qiime rescript get-silva-data \
    --p-version '138' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-ssu-nr99-seqs.qza \
    --o-silva-taxonomy silva-138-ssu-nr99-tax.qza

qiime rescript cull-seqs \
    --i-sequences silva-138-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138-ssu-nr99-seqs-cleaned.qza
    
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138-ssu-nr99-seqs-discard.qza
    
qiime rescript dereplicate \
    --i-sequences silva-138-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138-ssu-nr99-tax.qza \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138-ssu-nr99-tax-derep-uniq.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138-ssu-nr99-seqs-derep-uniq.qza \
  --i-reference-taxonomy silva-138-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138-ssu-nr99-classifier.qza

qiime feature-classifier extract-reads \
    --i-sequences silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACNVGGGTWTCTAAT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138-ssu-nr99-seqs-515f-806r.qza

qiime rescript dereplicate \
    --i-sequences silva-138-ssu-nr99-seqs-515f-806r.qza \
    --i-taxa silva-138-ssu-nr99-tax-derep-uniq.qza \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138-ssu-nr99-seqs-515f-806r-uniq.qza \
    --o-dereplicated-taxa  silva-138-ssu-nr99-tax-515f-806r-derep-uniq.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138-ssu-nr99-seqs-515f-806r-uniq.qza \
    --i-reference-taxonomy silva-138-ssu-nr99-tax-515f-806r-derep-uniq.qza \
    --o-classifier silva-138-ssu-nr99-515f-806r-classifier.qza

# train the classifier
source activate qiime2-2020.6
cd /home/tony/sync/mangroves/
qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads ./references_taxons/16S_silva138_NR99_otus.qza \
   --i-reference-taxonomy /media/tony/DATA2/tax/SILVA_138/silva-138-ssu-nr99-tax.qza \
   --o-classifier ./16S_classifier_qiime_NR99_silva138.qza

# 16S BLAST
# le pipeline recript qui a permis d'entraîner le classificateur a été lancé avec sci-sklearn 0.21.2
# en upgradant qiime2, on upgrade la version de sci-sklearn et le classify-sklearn n'est plus compatible
# il faut alors réinstaller la version 0.21.2 de sci-sklearn sous environnement qiime2 2020.6 :
# python -m pip install scikit-learn==0.21.2
# il faut activer la mémoire swap (128Go)
# ou alors descendre --p-n-jobs à 8

cd ~/sync/lichens/bio_info/data_sequencages/2021_06/210608/qiime2

qiime feature-classifier classify-sklearn \
   --i-reads ./OTU-16S/rep-seqs-97_mosaique_210608_16S.qza \
   --i-classifier /media/tony/DATA2/tax/SILVA_138/silva-138-ssu-nr99-515f-926r-classifier.qza \
   --p-n-jobs 32 \
   --o-classification ./OTU-16S/classification-97_Silva138_mosaique_210608_16S.qza
   
# exporter en .biom
qiime tools export \
    --input-path ./OTU-16S/classification-97_Silva138_mosaique_210608_16S.qza \
    --output-path ./OTU-16S/Silva_138

# 18S BLAST

cd ~/sync/lichens/bio_info/data_sequencages/2021_06/210608/qiime2

qiime feature-classifier classify-sklearn \
   --i-reads ./OTU-18S/rep-seqs-97_mosaique_210608_18S.qza \
   --i-classifier /media/tony/DATA2/tax/SILVA_138/silva-138-ssu-nr99-18SV9-classifier.qza \
   --p-n-jobs 32 \
   --o-classification ./OTU-18S/classification-97_Silva138_mosaique_210608_18S.qza
   
# exporter en .biom
qiime tools export \
    --input-path ./OTU-18S/classification-97_Silva138_mosaique_210608_18S.qza \
    --output-path ./OTU-18S/Silva_138
