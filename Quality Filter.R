# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("dada2")
library(dada2)

# On traite les deux marqueurs (16SV4 et ITS2) ensemble, on sort des seqtab_nochim communes
# puis on sépare les marqueurs dans QIIME2 avec $ qiime feature-classifier extract-reads


setwd(path)
# path <- "/Volumes/BOUGREDANE/191224_mang1"
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1_final.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_final.fastq", full.names = TRUE))

fnFs <- sort(fnFs)
fnRs <- sort(fnRs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
sample.names
# plotQualityProfile(fnFs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=100, matchIDs=TRUE, maxN=0, maxEE=c(3,3), rm.phix=TRUE, compress=TRUE, multithread=TRUE)  #trimLeft=15, trimRight=15

out

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# plot(table(nchar(getSequences(seqtab))))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)
# plot(table(nchar(getSequences(seqtab.nochim))))
write.csv(t(seqtab.nochim),"~/sync/lichens/bio_info/data_sequencages/2021_07/dada2/seqtabnochim_mosaique_16S_18S_ITS.csv")
uniquesToFasta(seqtab.nochim,"~/sync/lichens/bio_info/data_sequencages/2021_07/dada2/seqtabnochim_mosaique_16S_18S_ITS.fasta")

write.table(t(seqtab.nochim), "~/sync/lichens/bio_info/data_sequencages/2021_07/dada2/seqtabnochim_mosaique_16S_18S_ITS.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, "~/sync/lichens/bio_info/data_sequencages/2021_07/dada2/seqtabnochim_mosaique_16S_18S_ITS.fna", ids=colnames(seqtab.nochim))


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,"~/sync/lichens/bio_info/data_sequencages/2021_07/dada2/seqtabnochim_mosaique_16S_18S_ITS_track.txt")




##QIIME2 
#Bash

source activate qiime2-2020.6 cd ~/
  
  {
    "cells": [
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "## Pipeline qiime2 BLAST des OTUs sur amplicons Illumina\n",
          "\n",
          "* importation des données DADA2\n",
          "Séquences illumina avec plusieurs marqueurs par index, ici 16SV4 et ITS2.\n",
          "On débruite et déréplique tout en bloc avec DADA2 sous R car c'est beaucoup plus rapide que DADA2 sous environnement qiime2. Puis on importe la table des séquences uniques de DADA2 en artefacts qiime2.\n",
          "\n",
          "**sous R (DADA2):**\n",
          "write.table(t(seqtab.nochim), \"~/sync/mangroves/data_sequencages/guyane/dada2-analysis/seqtab-nochim.txt\", sep=\"\\t\", row.names=TRUE, col.names=NA, quote=FALSE)\n",
          "\n",
          "uniquesToFasta(seqtab.nochim, \"~/sync/mangroves/data_sequencages/guyane/dada2-analysis/rep-seqs.fna\", ids=colnames(seqtab.nochim))\n",
          "\n",
          "* sous bash, activation environnement conda qiime2"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 3,
        "metadata": {
          "scrolled": true
        },
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "# kernel bash\n",
          "# conda info --envs\n",
          "source activate qiime2-2020.6"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* création de l'artefact rep-seqs.qza issu de seqtab.nochim (renommé en rep-seqs.fna sous DADA2/R)"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 2,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mImported ./dada2/seqtabnochim_mosaique_16S_18S_ITS.fna as DNASequencesDirectoryFormat to ./qiime2/rep-seqs_mosaique_210710_16S_18S_ITS.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07\n",
          "\n",
          "qiime tools import \\\n",
          "--input-path ./dada2/seqtabnochim_mosaique_16S_18S_ITS.fna \\\n",
          "--type 'FeatureData[Sequence]' \\\n",
          "--output-path ./qiime2/rep-seqs_mosaique_210710_16S_18S_ITS.qza"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* création de l'artefact table.qza par biom"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 3,
        "metadata": {
          "scrolled": true
        },
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mImported ./qiime2/seqtabnochim_mosaique_16S_18S_ITS.biom as BIOMV210DirFmt to ./qiime2/table_mosaique_210710_16S_18S_ITS.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07\n",
          "\n",
          "biom convert -i ./dada2/seqtabnochim_mosaique_16S_18S_ITS.txt -o ./qiime2/seqtabnochim_mosaique_16S_18S_ITS.biom --table-type=\"OTU table\" --to-hdf5\n",
          "\n",
          "qiime tools import \\\n",
          "--input-path ./qiime2/seqtabnochim_mosaique_16S_18S_ITS.biom \\\n",
          "--type 'FeatureTable[Frequency]' \\\n",
          "--output-path ./qiime2/table_mosaique_210710_16S_18S_ITS.qza"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* séparation des reads mergés par les amorces ITS2 (ITS86f-next et ITS4r-next)\n",
          "** quid des NNN de décalage pré-amorces ?  -> pas de prob avec --p-identity 0.9"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 4,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Sequence] to: ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-ITS2/table_mosaique_210710_ITS.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "# mkdir ./OTU-ITS2\n",
          "\n",
          "qiime feature-classifier extract-reads \\\n",
          "  --i-sequences rep-seqs_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --p-f-primer GTGAATCATCGAATCTTTGAA \\\n",
          "  --p-r-primer TCCTCCGCTTATTGATATGC \\\n",
          "  --p-identity 0.9 \\\n",
          "  --p-n-jobs 32 \\\n",
          "  --o-reads ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza\n",
          "  \n",
          "qiime feature-table filter-features \\\n",
          "  --i-table table_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --m-metadata-file ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza \\\n",
          "  --o-filtered-table ./OTU-ITS2/table_mosaique_210710_ITS.qza"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* séparation des reads mergés par les amorces 16SV4 (926r-next et 515f-next), avec 30% de flexibilité pour les bases dégénérées"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 5,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) mkdir: cannot create directory ‘./OTU-16S’: File exists\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Sequence] to: ./OTU-16S/rep-seqs_mosaique_210710_16S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-16S/table_mosaique_210710_16S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "mkdir ./OTU-16S\n",
          "\n",
          "qiime feature-classifier extract-reads \\\n",
          "  --i-sequences rep-seqs_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --p-f-primer GTGYCAGCMGCCGCGGTAA \\\n",
          "  --p-r-primer CCGYCAATTYMTTTRAGTTT \\\n",
          "  --p-identity 0.7 \\\n",
          "  --p-n-jobs 32 \\\n",
          "  --o-reads ./OTU-16S/rep-seqs_mosaique_210710_16S.qza\n",
          "\n",
          "qiime feature-table filter-features \\\n",
          "  --i-table table_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --m-metadata-file ./OTU-16S/rep-seqs_mosaique_210710_16S.qza \\\n",
          "  --o-filtered-table ./OTU-16S/table_mosaique_210710_16S.qza"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 6,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Sequence] to: ./OTU-18S/rep-seqs_mosaique_210710_18S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-18S/table_mosaique_210710_18S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "mkdir ./OTU-18S\n",
          "\n",
          "qiime feature-classifier extract-reads \\\n",
          "  --i-sequences rep-seqs_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --p-f-primer CCCTGCCHTTTGTACACAC \\\n",
          "  --p-r-primer CCTTCYGCAGGTTCACCTAC \\\n",
          "  --p-identity 0.7 \\\n",
          "  --p-n-jobs 32 \\\n",
          "  --o-reads ./OTU-18S/rep-seqs_mosaique_210710_18S.qza\n",
          "\n",
          "qiime feature-table filter-features \\\n",
          "  --i-table table_mosaique_210710_16S_18S_ITS.qza \\\n",
          "  --m-metadata-file ./OTU-18S/rep-seqs_mosaique_210710_18S.qza \\\n",
          "  --o-filtered-table ./OTU-18S/table_mosaique_210710_18S.qza"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* VSEARCH produit la table des OTUs (<3%), pour chaque marqueur"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 7,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-ITS2/table-97_mosaique_210710_ITS.qza\u001b[0m\r\n",
              "\u001b[32mSaved FeatureData[Sequence] to: ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-16S/table-97_mosaique_210710_16S.qza\u001b[0m\r\n",
              "\u001b[32mSaved FeatureData[Sequence] to: ./OTU-16S/rep-seqs-97_mosaique_210710_16S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureTable[Frequency] to: ./OTU-18S/table-97_mosaique_210710_18S.qza\u001b[0m\r\n",
              "\u001b[32mSaved FeatureData[Sequence] to: ./OTU-18S/rep-seqs-97_mosaique_210710_18S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "qiime vsearch cluster-features-de-novo \\\n",
          "  --i-table ./OTU-ITS2/table_mosaique_210710_ITS.qza \\\n",
          "  --i-sequences ./OTU-ITS2/rep-seqs_mosaique_210710_ITS.qza \\\n",
          "  --p-perc-identity 0.97 \\\n",
          "  --p-threads 32 \\\n",
          "  --o-clustered-table ./OTU-ITS2/table-97_mosaique_210710_ITS.qza \\\n",
          "  --o-clustered-sequences ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza\n",
          "\n",
          "qiime vsearch cluster-features-de-novo \\\n",
          "  --i-table ./OTU-16S/table_mosaique_210710_16S.qza \\\n",
          "  --i-sequences ./OTU-16S/rep-seqs_mosaique_210710_16S.qza \\\n",
          "  --p-perc-identity 0.97 \\\n",
          "  --p-threads 32 \\\n",
          "  --o-clustered-table ./OTU-16S/table-97_mosaique_210710_16S.qza \\\n",
          "  --o-clustered-sequences ./OTU-16S/rep-seqs-97_mosaique_210710_16S.qza\n",
          "  \n",
          "qiime vsearch cluster-features-de-novo \\\n",
          "  --i-table ./OTU-18S/table_mosaique_210710_18S.qza \\\n",
          "  --i-sequences ./OTU-18S/rep-seqs_mosaique_210710_18S.qza \\\n",
          "  --p-perc-identity 0.97 \\\n",
          "  --p-threads 32 \\\n",
          "  --o-clustered-table ./OTU-18S/table-97_mosaique_210710_18S.qza \\\n",
          "  --o-clustered-sequences ./OTU-18S/rep-seqs-97_mosaique_210710_18S.qza"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* exporter les tables d'OTUs"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 8,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-ITS2/table-97_mosaique_210710_ITS.qza as BIOMV210DirFmt to directory ./OTU-ITS2/\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-16S/table-97_mosaique_210710_16S.qza as BIOMV210DirFmt to directory ./OTU-16S\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-18S/table-97_mosaique_210710_18S.qza as BIOMV210DirFmt to directory ./OTU-18S\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "# ITS\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-ITS2/table-97_mosaique_210710_ITS.qza \\\n",
          "    --output-path ./OTU-ITS2/\n",
          "    \n",
          "biom convert \\\n",
          "    -i ./OTU-ITS2/feature-table.biom \\\n",
          "    -o ./OTU-ITS2/feature-table-97_mosaique_210710_ITS.tsv --to-tsv\n",
          "\n",
          "# 16S\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-16S/table-97_mosaique_210710_16S.qza \\\n",
          "    --output-path ./OTU-16S\n",
          "    \n",
          "biom convert \\\n",
          "    -i ./OTU-16S/feature-table.biom \\\n",
          "    -o ./OTU-16S/feature-table-97_mosaique_210710_16S.tsv --to-tsv\n",
          "\n",
          "# 18S\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-18S/table-97_mosaique_210710_18S.qza \\\n",
          "    --output-path ./OTU-18S\n",
          "    \n",
          "biom convert \\\n",
          "    -i ./OTU-18S/feature-table.biom \\\n",
          "    -o ./OTU-18S/feature-table-97_mosaique_210710_18S.tsv --to-tsv"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "* feature-classifier BLAST\n",
          "* de préférence sur un cluster / un calculateur"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 9,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Taxonomy] to: ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza as TSVTaxonomyDirectoryFormat to directory ./OTU-ITS2/Unite\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "########### ITS2\n",
          "# BLAST+ Unite\n",
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "qiime feature-classifier classify-consensus-blast \\\n",
          "    --i-query ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza \\\n",
          "    --i-reference-reads ~/sync/references_metabarcoding/references_taxons/sh_refs_qiime_ver8_97_s_04.02.2020.fasta.qza \\\n",
          "    --i-reference-taxonomy  ~/sync/references_metabarcoding/references_taxons/sh_taxonomy_qiime_ver8_97_s_04.02.2020.qza \\\n",
          "    --o-classification ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza\n",
          "\n",
          "# exporter en .biom\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-ITS2/classification_mosaique_210710_ITS_Unite_qiime_ver8_97_s_04.02.2020.qza \\\n",
          "    --output-path ./OTU-ITS2/Unite"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 10,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Taxonomy] to: ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza as TSVTaxonomyDirectoryFormat to directory ./OTU-ITS2/PLANiTS\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "########### ITS2\n",
          "# BLAST+ PLANiTS\n",
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "qiime feature-classifier classify-sklearn \\\n",
          "   --i-reads ./OTU-ITS2/rep-seqs-97_mosaique_210710_ITS.qza \\\n",
          "   --i-classifier ~/sync/references_metabarcoding/ITS2/PLANiTS_29-03-2020/ITS2_PLANiTS_20-03-2020_classifier_q2020-6.qza \\\n",
          "   --p-n-jobs -1 \\\n",
          "   --o-classification ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza\n",
          "   \n",
          "# exporter en .biom\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-ITS2/classification_mosaique_210710_ITS_PLANiTS_qiime_20-03-2020.qza \\\n",
          "    --output-path ./OTU-ITS2/PLANiTS"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 11,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Taxonomy] to: ./OTU-16S/classification-97_Silva138_mosaique_210710_16S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-16S/classification-97_Silva138_mosaique_210710_16S.qza as TSVTaxonomyDirectoryFormat to directory ./OTU-16S/Silva_138\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "# 16S BLAST\n",
          "# le pipeline recript qui a permis d'entraîner le classificateur a été lancé avec sci-sklearn 0.21.2\n",
          "# en upgradant qiime2, on upgrade la version de sci-sklearn et le classify-sklearn n'est plus compatible\n",
          "# il faut alors réinstaller la version 0.21.2 de sci-sklearn sous environnement qiime2 2020.6 :\n",
          "# python -m pip install scikit-learn==0.21.2\n",
          "# il faut activer la mémoire swap (128Go)\n",
          "# ou alors descendre --p-n-jobs à 8\n",
          "\n",
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "qiime feature-classifier classify-sklearn \\\n",
          "   --i-reads ./OTU-16S/rep-seqs-97_mosaique_210710_16S.qza \\\n",
          "   --i-classifier /media/tony/DATA2/tax/SILVA_138/silva-138-ssu-nr99-515f-926r-classifier.qza \\\n",
          "   --p-n-jobs 32 \\\n",
          "   --o-classification ./OTU-16S/classification-97_Silva138_mosaique_210710_16S.qza\n",
          "   \n",
          "# exporter en .biom\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-16S/classification-97_Silva138_mosaique_210710_16S.qza \\\n",
          "    --output-path ./OTU-16S/Silva_138"
        ]
      },
      {
        "cell_type": "code",
        "execution_count": 12,
        "metadata": {},
        "outputs": [
          {
            "name": "stdout",
            "output_type": "stream",
            "text": [
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mSaved FeatureData[Taxonomy] to: ./OTU-18S/classification-97_Silva138_mosaique_210710_18S.qza\u001b[0m\r\n",
              "(qiime2-2020.6) (qiime2-2020.6) (qiime2-2020.6) \u001b[32mExported ./OTU-18S/classification-97_Silva138_mosaique_210710_18S.qza as TSVTaxonomyDirectoryFormat to directory ./OTU-18S/Silva_138\u001b[0m\r\n",
              "(qiime2-2020.6) \n"
            ]
          }
        ],
        "source": [
          "# 18S BLAST\n",
          "\n",
          "cd ~/sync/lichens/bio_info/data_sequencages/2021_07/qiime2\n",
          "\n",
          "qiime feature-classifier classify-sklearn \\\n",
          "   --i-reads ./OTU-18S/rep-seqs-97_mosaique_210710_18S.qza \\\n",
          "   --i-classifier /media/tony/DATA2/tax/SILVA_138/silva-138-ssu-nr99-18SV9-classifier.qza \\\n",
          "   --p-n-jobs 32 \\\n",
          "   --o-classification ./OTU-18S/classification-97_Silva138_mosaique_210710_18S.qza\n",
          "   \n",
          "# exporter en .biom\n",
          "qiime tools export \\\n",
          "    --input-path ./OTU-18S/classification-97_Silva138_mosaique_210710_18S.qza \\\n",
          "    --output-path ./OTU-18S/Silva_138"
        ]
      },
      {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
          "** On a ainsi créé, par marqueur :\n",
          "- un tableau des OTUs à 3% (séquences et abondances des reads après filtres qualités) [feature-table_guyane_ITS2.tsv]\n",
          "- un tableau des assignations taxonomiques des OTUs [taxonomy.tsv]"
        ]
      }
    ],
    "metadata": {
      "kernelspec": {
        "display_name": "Calysto Bash",
        "language": "bash",
        "name": "calysto_bash"
      },
      "language_info": {
        "file_extension": ".sh",
        "help_links": [
          {
            "text": "MetaKernel Magics",
            "url": "https://metakernel.readthedocs.io/en/latest/source/README.html"
          }
        ],
        "mimetype": "text/x-sh",
        "name": "bash",
        "version": "0.2.2"
      }
    },
    "nbformat": 4,
    "nbformat_minor": 4
  }

