########################### SCRIPT ################################
# by Cédric Hubas 
########################### SCRIPT ################################

#####################
# PACKAGES
#####################
library(ade4)
library(ggplot2)
library(scales)
library(cowplot)
#library(devtools)

#####################
# AESTHETICS
#####################
col.pal <- colorRampPalette(c("red3","orange","green3","royalblue"))
row.pal <- colorRampPalette(c("#DB2218","#77ABD6","#8FCF3C","#853894","#FC7F3C"))
# please modify col.pal and row.pal to change scores and loadings plot label colors

#####################
# Description
#####################

# The folowing function performs a specific supervised analysis.
# It performs a Multiple Factor Analysis as described by Escofier and Pages [1] by using the dudi.pca function of package ade4.
# The generated objects of class pca and dudi. are then used to perform a Between Class Analysis as described by Dolédec and Chessel [2].
# The final result is a supervised MFA called BC-MFA

#####################
# Usage
#####################

#bc.mfa(df,bloc,fac,spcos=0,X=1,Y=2)

#####################
# Arguments
#####################
# note: The used must identify several groups of variables within the data frame in odert to build argument bloc.

# df => a data frame with n rows (individuals) and p columns (numeric variables).
# bloc => a vector or factor object giving the groups for the corresponding groups of variable of df.
# fac => an external factor used for the supervised analysis (BCA). 
# spcos => a numerical value giving the cos2 by which variables text and symbols should be magnified in the variable plot.
# X and Y => the dimension of the principal components to be plotted.

#####################
# Note
#####################

# The function returns a plot of the Multiple factor analysis ordination (top-left pannel)
# The function returns also a plot of the supervised MFA (i.e. BC-MFA) in the bottom-left pannel
# Total Inertia Explained (T.I.E) by the chosen factor is given in the second plot
# The function returns also a plot of the loadings (i.e. variables) of the BC-MFA
# percentages of inertia explained by each axis is also given in axis titles
# The functions also returns the result of the generic function randtest() whch performs a Monte-Carlo test of the BC-MFA
# Please note that the BCA allows the extraction of a number k-1 of components which will be a function of the number of modalities (k) of the factor.
# If the fac factor has less than 3 modalities, the analysis will not be possible. 

#####################
# geom_convexhull function
# by : Charles Martin https://github.com/cmartin/ggConvexHull
#####################

# devtools::install_github("cmartin/ggConvexHull",force=T)

StatConvexHull <- ggplot2::ggproto(
  "StatConvexHull",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  compute_group = function(self, data, scales, params) {
    data[chull(data$x, data$y), ]
  }
)

stat_convexhull <- function(mapping = NULL, data = NULL, geom = "polygon",
                            position = "identity", show.legend = NA,
                            inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatConvexHull,
    data = data, mapping = mapping, geom = geom, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, params = list(...)
  )
}

geom_convexhull <- function (mapping = NULL, data = NULL, stat = "convex_hull", position = "identity",
                             ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = ggplot2::GeomPolygon,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#####################
# bc.mfa function
# by Cédric Hubas : https://github.com/Hubas-prog/BC-MFA
#####################

bc.mfa<-function(df,bloc,fac,spcos=0,X=1,Y=2,...){
  v<-NULL
  for(i in 1:length(bloc)){
    v[i]<-rep(paste("Group",i,sep=""))
  }
  splitfac<-factor(rep(v,bloc))
  
  LI<-NULL
  for(i in 1:length(bloc)){
    LI[i]<-list(df[,splitfac==paste("Group",i,sep="")])
  }
  names(LI)<-v
  
  # MFA
  eig<-unlist(lapply(LI,function(x){dudi.pca(x,scannf=F,nf=2)$eig[1]}))
  res.mfa<-dudi.pca(df,col.w=rep(1/eig,bloc),scannf=F,nf=ncol(df))
  varexp1<-res.mfa$eig*100/sum(res.mfa$eig)
  
  # BC-MFA
  k<-length(levels(factor(fac)))
  res.bcmfa<-bca(res.mfa,fac,scannf=F,nf=k-1)
  varexp2<-res.bcmfa$eig*100/sum(res.bcmfa$eig)
  
  #plots
  
  mfaX<-res.mfa$li[,X]
  mfaY<-res.mfa$li[,Y]
  
  MFA.ind<-ggplot(res.mfa$li,aes(x=res.mfa$li[,X],y=res.mfa$li[,Y],col=fac))+
    geom_point()+
    geom_convexhull(alpha = 0.3,aes(fill = fac))+
    scale_fill_manual(values=row.pal(length(levels(fac))))+
    scale_colour_manual(values=row.pal(length(levels(fac))))+
    ggtitle("MFA scores")+
    xlab(paste("Axis ",X," : ",round(varexp1[X],2),"%"))+
    ylab(paste("Axis ",Y," : ",round(varexp1[Y],2),"%"))+
    labs(fill="Grouping variable",col="Grouping variable")+
    theme_bw()
  
  BCMFA.ind<-ggplot(res.bcmfa$ls,aes(x=res.bcmfa$ls[,X],y=res.bcmfa$ls[,Y],col=fac))+
    geom_point()+
    scale_fill_manual(values=row.pal(length(levels(fac))))+
    scale_colour_manual(values=row.pal(length(levels(fac))))+
    geom_convexhull(alpha = 0.3,aes(fill = fac))+
    ggtitle(paste("BC-MFA scores"," | T.I.E=",
                  round(res.bcmfa$ratio,2)*100,
                  "%",
                  sep=""))+
    xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
    ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
    labs(fill="Explanatory variable",col="Explanatory variable")+
    theme_bw()
  
  cos2 <- as.matrix(res.bcmfa$co[,c(X,Y)])*as.matrix(res.bcmfa$co[,c(X,Y)])
  var.group<-factor(rep(v,bloc),levels=v)
  res.bcmfa$co$col<-col.pal(length(bloc))[var.group]
  col.filter<-cos2[,1]>spcos | cos2[,2]>spcos
  
  if(spcos==0) {
    BCMFA.var<-ggplot(res.bcmfa$co,aes(x=res.bcmfa$co[,X],
                                       y=res.bcmfa$co[,Y],col=var.group))+
      scale_colour_manual(values=col.pal(length(levels(var.group))))+
      geom_segment(aes(x = 0, y = 0,
                       xend = res.bcmfa$co[,X],
                       yend = res.bcmfa$co[,Y]),
                   arrow = arrow(length = unit(0.5, "cm")),
                   col="lightgrey")+
      annotate(geom="text",
               x=res.bcmfa$co[,X],
               y=res.bcmfa$co[,Y],
               label=rownames(res.mfa$co),
               size=4,
               col=res.bcmfa$co$col)+
      ggtitle("BC-MFA variables")+
      xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
      ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
      xlim(c(-1,1))+
      ylim(c(-1,1))+
      theme_bw()
  }else {
    BCMFA.var<-ggplot(res.bcmfa$co,aes(x=res.bcmfa$co[,X],
                                       y=res.bcmfa$co[,Y],col=var.group))+
      geom_segment(aes(x = 0, y = 0,
                       xend = res.bcmfa$co[,X],
                       yend = res.bcmfa$co[,Y]),
                   arrow = arrow(length = unit(0.5, "cm")),
                   col="lightgrey")+
      annotate(geom="text",
               x=res.bcmfa$co[,X],
               y=res.bcmfa$co[,Y],
               label=rownames(res.mfa$co),
               size=4,
               col=alpha(res.bcmfa$co$col,c(0.15,1)[factor(col.filter)]))+
      ggtitle("BC-MFA variables")+
      xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
      ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
      xlim(c(-1,1))+
      ylim(c(-1,1))+
      theme_bw()
  }
  
  left<-plot_grid(MFA.ind,BCMFA.ind,labels=c("a","b"),ncol=1)
  PLOT<-plot_grid(left,BCMFA.var,labels=c("","c"),ncol=2)
  return(list(PLOT,perm.test=randtest(res.bcmfa),
              BCMFAcos2=cos2,
              BCMFAco=res.bcmfa$co,
              BCMFAind=res.bcmfa$ls))
}

# col=alpha(res.bcmfa$co$col,c(0.15,1)[factor(col.filter)])

#############################################################################################################

#######################
###   Import Data   ###
#######################

###############
##  Bio mol  ## 
###############

library(phyloseq)
library(dplyr)
library(microbiome)

setwd("/Users/guill/OneDrive/Documents/Stage M2/Sequencage/Core microbiome/")

samples_mosaique <- read.table("./sample_data_mosaique_210710.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mosaique) <- samples_mosaique[,1]
samples_mosaique %>% as_tibble()


#### Core OTU 16S_ITS ### 

#Faire d'abord tourner script pour obtenir objet phyloseq contenant core OTU
  
otu_table_core<-otu_table(mosaique_final_sans_autres)
otu_table_core <- t(otu_table_core)
otu_table_core <- as.data.frame(otu_table_core)
rownames(otu_table_core) -> ID
ID <- gsub('_', ' ', ID)
ID -> otu_table_core$ID
ID
dim(otu_table_core) 
  



##################
##  Metabolite  ## 
##################

setwd("/Users/guill/OneDrive/Documents/Stage M2/BC-MFA/")

CHCL3<-read.csv("./MFA/t_CHCL3_code.csv",sep = ",", header = TRUE)
MeOH<-read.csv("./MFA/t_MeOH_code.csv",header=T,sep = ",")


Metabo<-merge(MeOH,CHCL3,by="ID") 



##################
###   BC-MFA   ###
##################

DATA<-merge(otu_table_core,Metabo,by="ID")
rownames(DATA)<-DATA[,1]

DATA <- DATA[, colSums(DATA != 0) > 0] # on enlève les colonnes à somme zero
length(DATA[, colSums(DATA != 0) > 0])
dim(DATA) # 58 121

#De 2 a 31=OTU, de 32 a 121=Metabo (90)

bloc<-c(dim(DATA[2:31])[2],
        dim(DATA[32:120])[2])

DATA$specie<-substr(DATA$ID,1,3)
DATA$thalle<-substr(DATA$ID,1,4)


res<-bc.mfa(df=DATA[,2:120],
       bloc=bloc,
       fac=factor(DATA$specie),
       spcos=0)

Cos2_metabo_coreOTU<-as.data.frame(res$BCMFAcos2)






