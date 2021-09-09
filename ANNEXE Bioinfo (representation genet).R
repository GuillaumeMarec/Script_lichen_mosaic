library(ggplot2)
library(wesanderson)
library(grDevices)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(ggrepel)
library(ade4)
library(gridExtra)
library(ampvis2)
library(vegan)
library(tidyr)
library(cowplot)


#################################
###   Microbial communities   ###
#################################

#Choisir Dossier esp�ce 
setwd("/Users/guill/OneDrive/Documents/Stage M2/Sequencage/Core microbiome/")

### Bio Mol ###

library(phyloseq)
library(dplyr)
library(microbiome)
#citation("microbiome")

samples_mosaique <- read.table("./sample_data_mosaique_210710.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique) <- samples_mosaique[,1]
samples_mosaique %>% as_tibble()


####16S
abond16S <- read.csv2("./table-97_OTUs_mosaique_210710_16S_abseules.csv", sep=",", header= T,dec = ".")
rownames(abond16S) <- paste("p",(1:nrow(abond16S)))
as_tibble(abond16S)
abond16S=as.matrix(abond16S)
dim(abond16S)
OTU_16S = otu_table(abond16S, taxa_are_rows = TRUE)

tax16S <- read.csv2("./table-97_OTUs_mosaique_210710_16S_taxa.csv", sep=";", header= F)
dim(tax16S)
rownames(tax16S) <- rownames(abond16S)
colnames(tax16S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax16S %>% as_tibble()
tax16S <- as.matrix(tax16S)

TAX_16S = tax_table(as.matrix(tax16S))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_16S = phyloseq(OTU_16S, TAX_16S, sampledata)
mosaique_16S_rel <- microbiome::transform(mosaique_16S, "compositional") # passe en abondances relatives
mosaique_16S_rel <- core(mosaique_16S_rel, detection = 0.01, prevalence = 0) # seuil ab rel = 1%
mosaique_16S
rel_16S<-taxa_names(mosaique_16S_rel)

mosaique_16S_final<-prune_taxa(rel_16S,mosaique_16S)




##### ITS
abondITS <- read.csv2("./table-97_OTUs_mosaique_210710_ITS_abseules.csv", sep=",", header= T)
rownames(abondITS) <- paste("i",(1:nrow(abondITS)))
as_tibble(abondITS)
dim(abondITS)
abondITS=as.matrix(abondITS)
OTU_ITS = otu_table(abondITS, taxa_are_rows = TRUE)

# je ne sais pas pourquoi mais il a fallu réouvrir le ITS_taxa.csv sur libreoffice et le ré-enregistrer avec sep=";"
taxITS <- read.csv2("./table-97_OTUs_mosaique_210710_ITS_taxa.csv", sep=";", header= F)
dim(taxITS)
as_tibble(taxITS)
tail(taxITS)
rownames(taxITS) <- rownames(abondITS)
colnames(taxITS) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
taxITS %>% as_tibble()
taxITS <- as.matrix(taxITS)

TAX_ITS = tax_table(as.matrix(taxITS))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_ITS = phyloseq(OTU_ITS, TAX_ITS, sampledata)
mosaique_ITS_rel <- microbiome::transform(mosaique_ITS, "compositional") # passe en abondances relatives
mosaique_ITS_rel <- core(mosaique_ITS_rel, detection = 0.01, prevalence = 0) # seuil ab rel = 5%
mosaique_ITS

rel_ITS<-taxa_names(mosaique_ITS_rel)
mosaique_ITS_final<-prune_taxa(rel_ITS,mosaique_ITS)


##### 18S
abond18S <- read.csv2("./table-97_OTUs_mosaique_210710_18S_abseules.csv", sep=",", header= T)
rownames(abond18S) <- paste("e",(1:nrow(abond18S)))
as_tibble(abond18S)
dim(abond18S)
abond18S=as.matrix(abond18S)
OTU_18S = otu_table(abond18S, taxa_are_rows = TRUE)

# je ne sais pas pourquoi mais il a fallu réouvrir le ITS_taxa.csv sur libreoffice et le ré-enregistrer avec sep=";"
tax18S <- read.csv2("./table-97_OTUs_mosaique_210710_18S_taxa.csv", sep=";", header= F)
dim(tax18S)
as_tibble(tax18S)
tail(tax18S)
rownames(tax18S) <- rownames(abond18S)
colnames(tax18S) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
tax18S %>% as_tibble()
tax18S <- as.matrix(tax18S)

TAX_18S = tax_table(as.matrix(tax18S))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_18S = phyloseq(OTU_18S, TAX_18S, sampledata)
mosaique_18S_rel <- microbiome::transform(mosaique_18S, "compositional") # passe en abondances relatives
mosaique_18S_rel <- core(mosaique_18S_rel, detection = 0.01, prevalence = 0) # seuil ab rel = 5%
mosaique_18S

rel_18S<-taxa_names(mosaique_18S_rel)
mosaique_18S_final<-prune_taxa(rel_18S,mosaique_18S)


mosaique_16S_ITS_18S <- merge_phyloseq(mosaique_16S_final, mosaique_ITS_final, mosaique_18S_final)
mosaique_16S_ITS_18S_rel<-microbiome::transform(mosaique_16S_ITS_18S, "compositional")




##########################
##                      ##
##    Ordinations       ##
##                      ##
##########################
# via phyloseq
nMDS_mosaique_all <- ordinate(mosaique_16S_ITS_18S_rel, "PCoA", "bray")
(ordplot <- plot_ordination(mosaique_16S_ITS_18S_rel, nMDS_mosaique_all, type = "samples", color="espece_code"))

#orddata <- ordplot$data
#?stat_ellipse
#png(file="./nMDS_otus_mosaique1.png", width=800, height=800)
ordplot + 
  stat_ellipse(type = "t") +
  theme_bw() + ggtitle("nMDS sur abondances relatives des OTUs, tous marqueurs")
#dev.off()

# via ampvis2
# d'abord enlever les échantillons qui n'ont pas d'otus
mosaique_16S_seul <- prune_samples(sample_sums(mosaique_16S_final)>0, mosaique_16S_final)
d_16S <- amp_load(otutable = otu_table(mosaique_16S_seul), taxonomy = tax_table(mosaique_16S_seul), 
                  metadata = samples_mosaique)

mosaique_ITS_seul <- prune_samples(sample_sums(mosaique_ITS_final)>0, mosaique_ITS_final)
d_ITS <- amp_load(otutable = otu_table(mosaique_ITS_seul), taxonomy = tax_table(mosaique_ITS_seul), 
                  metadata = samples_mosaique)

mosaique_18S_seul <- prune_samples(sample_sums(mosaique_18S_final)>0, mosaique_18S_final)
d_18S <- amp_load(otutable = otu_table(mosaique_18S_seul), taxonomy = tax_table(mosaique_18S_seul), 
                  metadata = samples_mosaique)

d_all <- amp_load(otutable = otu_table(mosaique_16S_ITS_18S_rel), taxonomy = tax_table(mosaique_16S_ITS_18S_rel), 
                  metadata = samples_mosaique)
d_16S
d_ITS
d_18S
d_all

pcoa_16S <- amp_ordinate(d_16S, 
                         type = "pcoa",
                         distmeasure = "bray",
                         sample_color_by = "espece_code",
                         sample_colorframe = TRUE,
                         sample_colorframe_label = "espece_code") + theme(legend.position = "blank") + ggtitle("Procaryotes 16SV4")+theme(plot.title = element_text(size = 20,face = "bold"))

pcoa_ITS <- amp_ordinate(d_ITS, 
                         type = "pcoa",
                         distmeasure = "bray",
                         sample_color_by = "espece_code",
                         sample_colorframe = TRUE,
                         sample_colorframe_label = "espece_code") + theme(legend.position = "blank") + ggtitle("ITS2")+theme(plot.title = element_text(size = 20,face = "bold"))

pcoa_18S <- amp_ordinate(d_18S, 
                         type = "pcoa",
                         distmeasure = "bray",
                         sample_color_by = "espece_code",
                         sample_colorframe = TRUE,
                         sample_colorframe_label = "espece_code") + theme(legend.position = "blank") + ggtitle("Eukaryotes 18SV9")+theme(plot.title = element_text(size = 20,face = "bold"))

pcoa_all <- amp_ordinate(d_all, 
                         type = "pcoa",
                         distmeasure = "bray",
                         sample_color_by = "espece_code",
                         sample_colorframe = TRUE,
                         sample_colorframe_label = "espece_code") + theme(legend.position = "blank") + ggtitle("All markers")+theme(plot.title = element_text(size = 20,face = "bold"))

#png(file="./pcoa_otus_mosaique1.png", width=800, height=800)
grid.arrange(pcoa_16S,pcoa_18S,pcoa_ITS, pcoa_all, ncol=2)
#dev.off()



############################
###                      ###
###   Core micrombiome   ###
###                      ###
############################


samples_mosaique <- read.table("./sample_data_mosaique.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique) <- samples_mosaique[,1]
samples_mosaique %>% as_tibble()
samples_mosaique_thalle <- read.table("./sample_data_mosaique_thalle.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique_thalle) <- samples_mosaique_thalle[,1]
samples_mosaique_thalle %>% as_tibble()

####16S
abond16S <- read.csv2("./table-97_OTUs_mosaique_210710_16S_abseules.csv", sep=",", header= T,dec = ".")
rownames(abond16S) <- paste("p",(1:nrow(abond16S)))
as_tibble(abond16S)
abond16S=as.matrix(abond16S)
dim(abond16S)
OTU_16S = otu_table(abond16S, taxa_are_rows = TRUE)

tax16S <- read.csv2("./table-97_OTUs_mosaique_210710_16S_taxa_unk.csv", sep=";", header= F)
dim(tax16S)
rownames(tax16S) <- rownames(abond16S)
colnames(tax16S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax16S %>% as_tibble()
tax16S <- as.matrix(tax16S)

TAX_16S = tax_table(as.matrix(tax16S))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_16S = phyloseq(OTU_16S, TAX_16S, sampledata)
#mosaique_16S <- rarefy_even_depth(mosaique_16S, sample.size = 1500,
#                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

mosaique_16S_rel <- microbiome::transform(mosaique_16S, "compositional") # passe en abondances relatives
# filtre qualité p/p 16S mock community 1%
mosaique_16S_rel <- core(mosaique_16S_rel, detection = 0.01, prevalence = 0) # seuil ab rel = 1%
mosaique_16S_rel

rel_16S<-taxa_names(mosaique_16S_rel)
mosaique_16S_final<-prune_taxa(rel_16S,mosaique_16S)


##### ITS
abondITS <- read.csv2("./table-97_OTUs_mosaique_210710_ITS_abseules.csv", sep=",", header= T)
rownames(abondITS) <- paste("i",(1:nrow(abondITS)))
as_tibble(abondITS)
dim(abondITS)
abondITS=as.matrix(abondITS)
OTU_ITS = otu_table(abondITS, taxa_are_rows = TRUE)

# je ne sais pas pourquoi mais il a fallu réouvrir le ITS_taxa.csv sur libreoffice et le ré-enregistrer avec sep=";"
taxITS <- read.csv2("./table-97_OTUs_mosaique_210710_ITS_taxa.csv", sep=";", header= F)
dim(taxITS)
as_tibble(taxITS)
tail(taxITS)
rownames(taxITS) <- rownames(abondITS)
colnames(taxITS) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
taxITS %>% as_tibble()
taxITS <- as.matrix(taxITS)

TAX_ITS = tax_table(as.matrix(taxITS))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_ITS = phyloseq(OTU_ITS, TAX_ITS, sampledata)
#mosaique_ITS <- rarefy_even_depth(mosaique_ITS, sample.size = 1500,
#                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
mosaique_ITS_rel <- microbiome::transform(mosaique_ITS, "compositional") # passe en abondances relatives
mosaique_ITS_rel <- core(mosaique_ITS_rel, detection = 0.05, prevalence = 0) # seuil ab rel = 5%
mosaique_ITS_rel

mosaique_16S_ITS <- merge_phyloseq(mosaique_16S,mosaique_ITS)


rel_ITS<-taxa_names(mosaique_ITS_rel)

mosaique_ITS_final<-prune_taxa(rel_ITS,mosaique_ITS)

mosaique_16S_ITS <- merge_phyloseq(mosaique_16S_final,mosaique_ITS_final)
otu_table(mosaique_16S_ITS) -> otu_table
tax_table(mosaique_16S_ITS) -> tax
paste(tax[,1], tax[,2],tax[,3],tax[,4],tax[,5],tax[,6],tax[,7]) -> tax
tax -> rownames(otu_table)
otu_table[1:2,1:2]



# Noyau microbiotique par espece de lichen

#Tephro
mosaique_16S_ITS_Tephro <- subset_samples(mosaique_16S_ITS, espece_nom=="Tephro")
Tephro_core_list <- core_members(mosaique_16S_ITS_Tephro, detection = 0, prevalence = 80/100,
                                 include.lowest = FALSE)
# mosaique_16S_ITS_Tephro_core <- prune_taxa(Tephro_core_list, mosaique_16S_ITS_Tephro)
mosaique_16S_ITS_Tephro_core <- subset_taxa(mosaique_16S_ITS_Tephro, rownames(tax_table(mosaique_16S_ITS_Tephro)) %in% Tephro_core_list)
mosaique_16S_ITS_Tephro_rare <- subset_taxa(mosaique_16S_ITS_Tephro, !(rownames(tax_table(mosaique_16S_ITS_Tephro)) %in% Tephro_core_list))
mosaique_16S_ITS_Tephro_autres <- merge_taxa(mosaique_16S_ITS_Tephro_rare, taxa_names(mosaique_16S_ITS_Tephro_rare)[1:length(taxa_names(mosaique_16S_ITS_Tephro_rare))], 2)
taxa_names(mosaique_16S_ITS_Tephro_autres) <- "autres"

otu_table(mosaique_16S_ITS_Tephro_autres)
mosaique_16S_ITS_Tephro_core_autres <- merge_phyloseq(mosaique_16S_ITS_Tephro_core,mosaique_16S_ITS_Tephro_autres)
mosaique_16S_ITS_Tephro_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Tephro_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Tephro_core_autres_rel)

#Ochro
mosaique_16S_ITS_Ochro <- subset_samples(mosaique_16S_ITS, espece_nom=="Ochro")
Ochro_core_list <- core_members(mosaique_16S_ITS_Ochro, detection = 0, prevalence = 80/100,
                                include.lowest = FALSE)
mosaique_16S_ITS_Ochro_core <- subset_taxa(mosaique_16S_ITS_Ochro, rownames(tax_table(mosaique_16S_ITS_Ochro)) %in% Ochro_core_list)
mosaique_16S_ITS_Ochro_rare <- subset_taxa(mosaique_16S_ITS_Ochro, !(rownames(tax_table(mosaique_16S_ITS_Ochro)) %in% Ochro_core_list))
mosaique_16S_ITS_Ochro_autres <- merge_taxa(mosaique_16S_ITS_Ochro_rare, taxa_names(mosaique_16S_ITS_Ochro_rare)[1:length(taxa_names(mosaique_16S_ITS_Ochro_rare))], 2)
taxa_names(mosaique_16S_ITS_Ochro_autres) <- "autres"

otu_table(mosaique_16S_ITS_Ochro_autres)
mosaique_16S_ITS_Ochro_core_autres <- merge_phyloseq(mosaique_16S_ITS_Ochro_core,mosaique_16S_ITS_Ochro_autres)
mosaique_16S_ITS_Ochro_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Ochro_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Ochro_core_autres_rel)


#Anapt
mosaique_16S_ITS_Anapt <- subset_samples(mosaique_16S_ITS, espece_nom=="Anapt")
Anapt_core_list <- core_members(mosaique_16S_ITS_Anapt, detection = 0, prevalence = 80/100,
                                include.lowest = FALSE)
mosaique_16S_ITS_Anapt_core <- subset_taxa(mosaique_16S_ITS_Anapt, rownames(tax_table(mosaique_16S_ITS_Anapt)) %in% Anapt_core_list)
mosaique_16S_ITS_Anapt_rare <- subset_taxa(mosaique_16S_ITS_Anapt, !(rownames(tax_table(mosaique_16S_ITS_Anapt)) %in% Anapt_core_list))
mosaique_16S_ITS_Anapt_autres <- merge_taxa(mosaique_16S_ITS_Anapt_rare, taxa_names(mosaique_16S_ITS_Anapt_rare)[1:length(taxa_names(mosaique_16S_ITS_Anapt_rare))], 2)
taxa_names(mosaique_16S_ITS_Anapt_autres) <- "autres"

otu_table(mosaique_16S_ITS_Anapt_autres)
mosaique_16S_ITS_Anapt_core_autres <- merge_phyloseq(mosaique_16S_ITS_Anapt_core,mosaique_16S_ITS_Anapt_autres)
mosaique_16S_ITS_Anapt_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Anapt_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Anapt_core_autres_rel)


#Xanth
mosaique_16S_ITS_Xanth <- subset_samples(mosaique_16S_ITS, espece_nom=="Xanth")
Xanth_core_list <- core_members(mosaique_16S_ITS_Xanth, detection = 0, prevalence = 80/100,
                                include.lowest = FALSE)
mosaique_16S_ITS_Xanth_core <- subset_taxa(mosaique_16S_ITS_Xanth, rownames(tax_table(mosaique_16S_ITS_Xanth)) %in% Xanth_core_list)
mosaique_16S_ITS_Xanth_rare <- subset_taxa(mosaique_16S_ITS_Xanth, !(rownames(tax_table(mosaique_16S_ITS_Xanth)) %in% Xanth_core_list))
mosaique_16S_ITS_Xanth_autres <- merge_taxa(mosaique_16S_ITS_Xanth_rare, taxa_names(mosaique_16S_ITS_Xanth_rare)[1:length(taxa_names(mosaique_16S_ITS_Xanth_rare))], 2)
taxa_names(mosaique_16S_ITS_Xanth_autres) <- "autres"

otu_table(mosaique_16S_ITS_Xanth_autres)
mosaique_16S_ITS_Xanth_core_autres <- merge_phyloseq(mosaique_16S_ITS_Xanth_core,mosaique_16S_ITS_Xanth_autres)
mosaique_16S_ITS_Xanth_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Xanth_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Xanth_core_autres_rel)


#Diplo
mosaique_16S_ITS_Diplo <- subset_samples(mosaique_16S_ITS, espece_nom=="Diplo")
Diplo_core_list <- core_members(mosaique_16S_ITS_Diplo, detection = 0, prevalence = 80/100,
                                include.lowest = FALSE)
mosaique_16S_ITS_Diplo_core <- subset_taxa(mosaique_16S_ITS_Diplo, rownames(tax_table(mosaique_16S_ITS_Diplo)) %in% Diplo_core_list)
mosaique_16S_ITS_Diplo_rare <- subset_taxa(mosaique_16S_ITS_Diplo, !(rownames(tax_table(mosaique_16S_ITS_Diplo)) %in% Diplo_core_list))
mosaique_16S_ITS_Diplo_autres <- merge_taxa(mosaique_16S_ITS_Diplo_rare, taxa_names(mosaique_16S_ITS_Diplo_rare)[1:length(taxa_names(mosaique_16S_ITS_Diplo_rare))], 2)
taxa_names(mosaique_16S_ITS_Diplo_autres) <- "autres"

otu_table(mosaique_16S_ITS_Diplo_autres)
mosaique_16S_ITS_Diplo_core_autres <- merge_phyloseq(mosaique_16S_ITS_Diplo_core,mosaique_16S_ITS_Diplo_autres)
mosaique_16S_ITS_Diplo_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Diplo_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Diplo_core_autres_rel)


#Buel
mosaique_16S_ITS_Buel <- subset_samples(mosaique_16S_ITS, espece_nom=="Buel")
Buel_core_list <- core_members(mosaique_16S_ITS_Buel, detection = 0, prevalence = 80/100,
                               include.lowest = FALSE)
mosaique_16S_ITS_Buel_core <- subset_taxa(mosaique_16S_ITS_Buel, rownames(tax_table(mosaique_16S_ITS_Buel)) %in% Buel_core_list)
mosaique_16S_ITS_Buel_rare <- subset_taxa(mosaique_16S_ITS_Buel, !(rownames(tax_table(mosaique_16S_ITS_Buel)) %in% Buel_core_list))
mosaique_16S_ITS_Buel_autres <- merge_taxa(mosaique_16S_ITS_Buel_rare, taxa_names(mosaique_16S_ITS_Buel_rare)[1:length(taxa_names(mosaique_16S_ITS_Buel_rare))], 2)
taxa_names(mosaique_16S_ITS_Buel_autres) <- "autres"

otu_table(mosaique_16S_ITS_Buel_autres)
mosaique_16S_ITS_Buel_core_autres <- merge_phyloseq(mosaique_16S_ITS_Buel_core,mosaique_16S_ITS_Buel_autres)
mosaique_16S_ITS_Buel_core_autres_rel <- microbiome::transform(mosaique_16S_ITS_Buel_core_autres, "compositional")
otu_table(mosaique_16S_ITS_Buel_core_autres_rel)



plot_mosaique_rel<-merge_phyloseq(mosaique_16S_ITS_Tephro_core_autres_rel, mosaique_16S_ITS_Ochro_core_autres_rel, mosaique_16S_ITS_Anapt_core_autres_rel, mosaique_16S_ITS_Xanth_core_autres_rel, mosaique_16S_ITS_Diplo_core_autres_rel, mosaique_16S_ITS_Buel_core_autres_rel)

mosaique_16S_ITS_core_tout_absolue<-merge_phyloseq(mosaique_16S_ITS_Tephro_core_autres, mosaique_16S_ITS_Ochro_core_autres, mosaique_16S_ITS_Anapt_core_autres, mosaique_16S_ITS_Xanth_core_autres, mosaique_16S_ITS_Diplo_core_autres, mosaique_16S_ITS_Buel_core_autres)
mosaique_16S_ITS_core_tout_absolue
core_otu_finale<-as.data.frame(mosaique_16S_ITS_core_tout_absolue@otu_table)


mosaique_16S_ITS_core_tout_absolue_rel <- microbiome::transform(mosaique_16S_ITS_core_tout_absolue,"compositional")

SD <- merge_samples(mosaique_16S_ITS_core_tout_absolue_rel, "thalle")
SD_rel <- microbiome::transform(SD,"compositional")
sampledata_tha = sample_data(data.frame(samples_mosaique_thalle, row.names=rownames(samples_mosaique_thalle), stringsAsFactors=FALSE))
sample_variables(SD_rel)
sample_names(sampledata_tha)
sample_data(SD_rel) <- sampledata_tha

# par Ordre
speciesListOrdre = unique(tax_table(SD_rel)[,"Order"])
#speciesPaletteOrdre = getPalette(length(speciesListOrdre))
speciesPaletteOrdre = c("steelblue","steelblue1","steelblue2","steelblue3","seagreen","orange","orange1","black","orange2","orange3","steelblue4","royalblue","sienna1","royalblue1","royalblue2","sienna2","royalblue3")
names(speciesPaletteOrdre) = speciesListOrdre

ordre_legende_Order <- c(",o__Acetobacterales",",o__Rhizobiales",",o__Acidobacteriales",",o__Armatimonadales", ",o__Chitinophagales",",o__Bryobacterales","",",o__Cyanobacteriales",",o__Isosphaerales","Trebouxiaceae","o__Caliciales","o__Lecanorales","o__Ostropales","o__Teloschistales","o__Pertusariales","o__Saccharomycetales", NA)

#pdf("main_symbiontes_histogrammes_Ordres.pdf")
plot_bar(SD_rel, fill="Order") + facet_grid(~espece_nom, scale="free_x", drop=TRUE) +
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") +
  scale_color_manual(values=speciesPaletteOrdre, guide="none") +
  scale_fill_manual(values= speciesPaletteOrdre, breaks=ordre_legende_Order, name = "Main symbiontes (Order)", labels = c("Acetobacterales (B)","Rhizobiales (B)","Acidobacteriales (B)","Armatimonadales (B)","Chitinophagales (B)","Bryobacterales (B)","inconnu (B)","Cyanobacteriales (B)","Isosphaerales (B)","Trebouxiaceae (P)","Caliciales (M)","Lecanorales (M)","Ostropales (M)","Teloschistales (M)","Pertusariales (M)","Saccharomycetales (M)", "non-core taxa (B,P,M)"))
#dev.off()

# par Classe
speciesListClasse = unique(tax_table(SD_rel)[,"Class"])
speciesPaletteClasse = c("steelblue","steelblue1","steelblue2","limegreen","orange","wheat","orange3","steelblue3","royalblue2","steelblue4","royalblue")
names(speciesPaletteClasse) = speciesListClasse
ordre_legende_Class <-  c(",c__Alphaproteobacteria",",c__Acidobacteriae",",c__Armatimonadia",",c__Bacteroidia",",c__Cyanobacteriia",",c__Planctomycetes","","Trebouxiales","c__Lecanoromycetes","c__Saccharomycetes",NA)

#pdf("main_symbiontes_histogrammes_Classes.pdf")
plot_bar(SD_rel, fill="Class") + facet_grid(~espece_nom, scale="free_x", drop=TRUE) +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  scale_color_manual(values=speciesPaletteClasse, guide="none") +
  scale_fill_manual(values= speciesPaletteClasse, breaks=ordre_legende_Class, name = "Main symbiontes (Class)", labels = c("Alphaproteobacteria (B)","Acidobacteriae (B)","Armatimonadia (B)","Bacteroidia (B)","Cyanobacteriia (B)","Planctomycetes (B)","unknown (B)","Trebouxiales (P)","Lecanoromycetes (M)","Saccharomycetes (M)","non-core taxa (B,P,M)"))
#dev.off()


##  Tree  ##
library("ape")
Tax_core<-tax_table(mosaique_final_sans_autres_rel)


Tree_Core<-rtree(ntaxa(mosaique_final_sans_autres_rel),rooted=TRUE,tip.label = taxa_names(Tax_core))
plot(Tree_Core)
mosaique_Core_tree<-merge_phyloseq(mosaique_final_sans_autres_rel,sampledata,Tree_Core)
ntaxa(mosaique_Core_tree)

plot_tree(mosaique_Core_tree, nodelabf=nodeplotblank, label.tips="Order", ladderize="left",color = "espece_nom", 
          title = "Tree core OTU 16S_ITS", size="abundance", base.spacing=0.03)
#+ coord_polar(theta="y")



###############
### Boxplot ###
###############

#Abondance relative par �chantillons
otu_table_core<-otu_table(mosaique_16S_ITS_core_tout_absolue_rel)
otu_table_core<-as.data.frame(otu_table_core)
otu_table_core<-t(otu_table_core)
otu_table_core<-as.data.frame(otu_table_core)

sp_lichen<-substr(rownames(otu_table_core),1,3)
th_lichen<-substr(rownames(otu_table_core),1,4)
samp_lichen<-substr(rownames(otu_table_core),1,5)


ggplot(otu_table_core,aes(x=sp_lichen,y=otu_table_core$`i 15`,fill=sp_lichen))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Set1")+ theme_classic()


ggplot(otu_table_core,aes(x=th_lichen,y=`i 93`,fill=sp_lichen))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Set1")



