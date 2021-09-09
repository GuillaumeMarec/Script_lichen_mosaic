library(ggplot2)
library(wesanderson)
library(dplyr)
library(grDevices)
#install.packages("RColorBrewer")
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library(corrplot)
#install.packages("ggrepel")
library(ggrepel)
library(ade4)
library(wesanderson)
library(ExPosition)
library(stats)
library(phyloseq)
library(dplyr)
library(microbiome)

setwd("/Users/guill/OneDrive/Documents/Stage M2/GC-MS/Representation finale/")

samples_mosaique <- read.table("./sample_data_mosaique_representation_thalle.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique) <- samples_mosaique[,1]
samples_mosaique %>% as_tibble()

#################
###   CHCL3   ###
#################

abondCHCL3<-read.csv("./PCA CHCL3_thalle.csv",sep = ",",header = TRUE,dec = ".")
rownames(abondCHCL3) <- paste("C",(1:nrow(abondCHCL3)))
as_tibble(abondCHCL3)
abondCHCL3=as.matrix(abondCHCL3)
dim(abondCHCL3)
OTU_CHCL3 = otu_table(abondCHCL3, taxa_are_rows = TRUE)

taxCHCL3 <- read.csv2("./Tax CHCL3.csv", sep=",", header= F)
dim(taxCHCL3)
rownames(taxCHCL3) <- rownames(abondCHCL3)
colnames(taxCHCL3) <- c("Code", "NIST", "family")
taxCHCL3 %>% as_tibble()
taxCHCL3 <- as.matrix(taxCHCL3)

TAX_CHCL3 = tax_table(as.matrix(taxCHCL3))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_CHCL3 = phyloseq(OTU_CHCL3, TAX_CHCL3, sampledata)
mosaique_CHCL3_rel <- microbiome::transform(mosaique_CHCL3, "compositional") # passe en abondances relatives
mosaique_CHCL3_rel



plot_bar(mosaique_CHCL3_rel, x="thalle", fill = "NIST")+ facet_grid(~espece_nom,scales = "free_x") +
  theme_minimal()+
  theme(axis.text.x = element_text(size=16,angle=90),
        axis.text.y = element_text(size=14))+ theme(legend.position='none')

+theme(legend.text = element_text(size = 18) )+theme(legend.title = element_text(size = 18))











#################
###   MEOH   ###
#################

abondMEOH<-read.csv("./PCA MEOH_thalle.csv",sep = ",",header = TRUE,dec = ".")
rownames(abondMEOH) <- paste("M",(1:nrow(abondMEOH)))
as_tibble(abondMEOH)
abondMEOH=as.matrix(abondMEOH)
dim(abondMEOH)
OTU_MEOH = otu_table(abondMEOH, taxa_are_rows = TRUE)

taxMEOH <- read.csv2("./Tax MEOH.csv", sep=",", header= F)
dim(taxMEOH)
rownames(taxMEOH) <- rownames(abondMEOH)
colnames(taxMEOH) <- c("Code", "NIST", "family")
taxMEOH %>% as_tibble()
taxMEOH <- as.matrix(taxMEOH)

TAX_MEOH = tax_table(as.matrix(taxMEOH))
sampledata = sample_data(data.frame(samples_mosaique, row.names=rownames(samples_mosaique), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_MEOH = phyloseq(OTU_MEOH, TAX_MEOH, sampledata)
mosaique_MEOH_rel <- microbiome::transform(mosaique_MEOH, "compositional") # passe en abondances relatives
mosaique_MEOH_rel



plot_bar(mosaique_MEOH_rel, x="thalle", fill = "NIST") +
  theme_minimal()+
  theme(axis.text.x = element_text(size=14,angle=45),
        axis.text.y = element_text(size=14))

+ theme(legend.position='none')




#############################################################################################
############ TOUT LES POINTS D ECHANTILLONNAGE ##############################################
#############################################################################################

samples_mosaique_tot <- read.table("./sample_data_mosaique_210710.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique_tot) <- samples_mosaique_tot[,1]
samples_mosaique_tot %>% as_tibble()


abondMEOH_tot<-read.csv("./PCA MEOH.csv",sep = ",",header = TRUE,dec = ".")
rownames(abondMEOH_tot) <- paste("M",(1:nrow(abondMEOH_tot)))
as_tibble(abondMEOH_tot)
abondMEOH_tot=as.matrix(abondMEOH_tot)
dim(abondMEOH_tot)
OTU_MEOH_tot = otu_table(abondMEOH_tot, taxa_are_rows = TRUE)

taxMEOH <- read.csv2("./Tax MEOH.csv", sep=",", header= F)
dim(taxMEOH)
rownames(taxMEOH) <- rownames(abondMEOH)
colnames(taxMEOH) <- c("Code", "NIST", "family")
taxMEOH %>% as_tibble()
taxMEOH <- as.matrix(taxMEOH)

TAX_MEOH = tax_table(as.matrix(taxMEOH))
sampledata = sample_data(data.frame(samples_mosaique_tot, row.names=rownames(samples_mosaique_tot), stringsAsFactors=FALSE))
sample_names(sampledata)
mosaique_MEOH_tot = phyloseq(OTU_MEOH_tot, TAX_MEOH, sampledata)
mosaique_MEOH_tot_rel <- microbiome::transform(mosaique_MEOH_tot, "compositional") # passe en abondances relatives
mosaique_MEOH_tot_rel



plot_bar(mosaique_MEOH_tot_rel, x="sample", fill = "NIST") +
  theme_minimal()+
  theme(axis.text.x = element_text(size=14,angle=90),
        axis.text.y = element_text(size=14))+ theme(legend.position='none')









#############################################################################################


# via phyloseq
nMDS_mosaique_core <- ordinate(mosaique_CHCL3_rel, "PCoA", "bray")
(ordplot <- plot_ordination(mosaique_CHCL3_rel, nMDS_mosaique_core, type = "samples", color="espece_code"))

ordplot + 
  stat_ellipse(type = "t") +
  theme_bw() + ggtitle("nMDS core OTU relative abundance (16S_ITS)")
#dev.off()

# via ampvis2
# d'abord enlever les échantillons qui n'ont pas d'otus
mosaique_core_seul <- prune_samples(sample_sums(mosaique_CHCL3_rel)>0, mosaique_CHCL3_rel)
d_core <- amp_load(otutable = otu_table(mosaique_core_seul), taxonomy = tax_table(mosaique_core_seul), 
                   metadata = samples_mosaique)

d_core


pcoa_core_CHCL3 <- amp_ordinate(d_core, 
                          type = "pcoa",
                          distmeasure = "bray",
                          sample_color_by = "espece_code",
                          sample_colorframe = TRUE,
                          sample_colorframe_label = "espece_code",
                          transform = "none") + theme(legend.position = "blank") + ggtitle("CHCL3")








# via phyloseq
nMDS_mosaique_core <- ordinate(mosaique_MEOH_rel, "PCoA", "bray")
(ordplot <- plot_ordination(mosaique_MEOH_rel, nMDS_mosaique_core, type = "samples", color="espece_code"))

ordplot + 
  stat_ellipse(type = "t") +
  theme_bw() + ggtitle("nMDS core OTU relative abundance (16S_ITS)")
#dev.off()

# via ampvis2
# d'abord enlever les échantillons qui n'ont pas d'otus
mosaique_core_seul <- prune_samples(sample_sums(mosaique_MEOH_rel)>0, mosaique_MEOH_rel)
d_core <- amp_load(otutable = otu_table(mosaique_core_seul), taxonomy = tax_table(mosaique_core_seul), 
                   metadata = samples_mosaique)

d_core


pcoa_core_MEOH <- amp_ordinate(d_core, 
                          type = "pcoa",
                          distmeasure = "bray",
                          sample_color_by = "espece_code",
                          sample_colorframe = TRUE,
                          sample_colorframe_label = "espece_code",
                          transform = "none") + theme(legend.position = "blank") + ggtitle("MeOH")



grid.arrange(pcoa_core_CHCL3,pcoa_core_MEOH, ncol=2)






