#libraries
library(clipr)
library(ggplot2)
library(readxl)
library(cowplot)
library(pheatmap)
library(nVennR)
library(tidyr)
library(dplyr)
library(viridis)
library(readxl)
library(clusterProfiler)

##########################################
########## DATA TO IMPORT ################
##########################################
#modify the working directory below accordingly
#setwd("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji")

#correspondance between probes and AGI for ATH1, modify the path accordingly
Probe.AGI.Gene.Atlas.correspondant.sand.may.2016<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/Probe.AGI.Gene.Atlas.correspondant.sand.may.2016.txt",header=T, row.names=1, sep="\t",fill=TRUE)

#list of treatments
colonnes<-c("hydro","NaCl25","NaCl50","NaCl75","NaCl100","Sorb50","Sorb100","Sorb150","PEG75",
            "PEG105","PEG125","PEG150","EG50","EG100","EG150","EG200")

#Import of transcriptome data, gcrma file from Gene Atlas system, modify the path accordingly. This file is available from GEO, GSE223207
a<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/Osmo.data.gcrma.txt",header=T, row.names=1, sep="\t",fill=TRUE)
a.omit<-na.omit(a)

#physiological data, modify the path accordingly
physio<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/New-13-11-2020-CPPdata .csv",sep = ",",header = T,)
physio<-physio[,-c(1,4)]
rownames(physio)<-colnames(a)

#Initial data import
#setwd("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji")
Probe.AGI.Gene.Atlas.correspondant.sand.may.2016<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/Probe.AGI.Gene.Atlas.correspondant.sand.may.2016.txt",header=T, row.names=1, sep="\t",fill=TRUE)

#list of treatments
colonnes<-c("hydro","NaCl25","NaCl50","NaCl75","NaCl100","Sorb50","Sorb100","Sorb150","PEG75",
            "PEG105","PEG125","PEG150","EG50","EG100","EG150","EG200")
#original gcrma file
a<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/Osmo.data.gcrma.txt",header=T, row.names=1, sep="\t",fill=TRUE)
a.omit<-na.omit(a)

#physiological data
physio<-read.table("C:/Yann B/AQUA2/marqueurs stress/Agropolis water potential markers/Watermarkers transcriptome/220121-files anova Yunji/New-13-11-2020-CPPdata .csv",sep = ",",header = T,)
physio<-physio[,-c(1,4)]
rownames(physio)<-colnames(a)

#CPP measurements for Figure 1
require(readxl)
Fig1_kinetic <- read_excel("Figure_1_data.xlsx")
#View(Fig1_kinetic)
Fig1_PiP <- read_excel("Figure_1_data.xlsx", 
                       sheet = "summary")

#Factors used in the ANOVA model
Na.Fact<-as.numeric(c(0,25,50,75,100,0,0,0,0,0,0,0,0,0,0,0,0,25,50,75,100,0,0,0,0,0,0,0,0,0,0,0))
Sorb.Fact<-as.numeric(c(0,0,0,0,0,50,100,150,0,0,0,0,0,0,0,0,0,0,0,0,0,50,100,150,0,0,0,0,0,0,0,0))
Eg.Fact<-as.numeric(c(0,0,0,0,0,0,0,0,0,0,0,0,50,100,150,200,0,0,0,0,0,0,0,0,0,0,0,0,50,100,150,200))
PEG.Fact<-as.numeric(c(0,0,0,0,0,0,0,0,75,100,125,150,0,0,0,0,0,0,0,0,0,0,0,0,75,100,125,150,0,0,0,0))
all.Fact<-matrix(c(Na.Fact,Sorb.Fact,PEG.Fact,Eg.Fact),nrow=32, ncol=4)

require(clipr)
#copying data from supplemental material of Sorenson et al. 2018, sov alpha values because Col is sov mutant
#columns A and C of sheet S2-Modeling Results
refdecay<-read_clip_tbl()

#Importing the RT-qPCR data
#start from gene list "220928-20_genes_frist_paper_QPCR.xls"
#copied only the AGI numbers
yjh.tested.agi<-read_clip_tbl(x=read_clip(), header= F)
#then the data
yjh_QPCR_data <- read_excel("C:/Yann B/AQUA2/marqueurs stress/2019-2023 Yunji Huang PhD/QPCRs/220928-20_genes_first_paper_QPCR.xlsx", 
                            sheet = "transposed for R first paper")
colnames(yjh_QPCR_data)[1]<-"AGI"

##########################################
########## FUNCTIONS #####################
##########################################

#ANOVA model with the 4 solutes factors altogether
Anova.over.data.v.treat.car.pckg<-function (data) 
{
  require(car)
  out <- 1:(length(data[, 1]) * 4)
  dim(out) <- c(length(data[, 1]),4)
  colnames(out) <- c("NaCl.pvalue","Sorb.pvalue","PEG.pvalue","Eg.pvalue")
  for (i in 1:length(data[, 1])) {
    Anova(lm(as.numeric(a[i ,]) ~ Na.Fact+Sorb.Fact+PEG.Fact+Eg.Fact))[1:4,4]->yo
    #print(yo)
    out[i, ] <- yo
  }
  row.names(out) <- row.names(data)
  return(out)
}

#ANOVA model using only 1 factor, used for Pi and P
Anova.over.1treat<-function (data, fact) 
{
  out <- 1:(length(data[, 1]))
  dim(out) <- c(length(data[, 1]),1)
  #colnames(out) <- c(".pvalue")
  for (i in 1:length(data[, 1])) {
    summary(aov(as.numeric(a[i ,]) ~ fact))[[1]][1,5]->yo
    #print(yo)
    out[i, ] <- yo
  }
  row.names(out) <- row.names(data)
  return(out)
}
#working function for plotting affy gene expression as a function of Pi or P
ggplot.lowess.physio.agi<-function(x,donneesphy, donnees){
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(cowplot)
  
  dfPi<-data.frame(matrix(ncol = dim(donneesphy)[1], nrow = 1))
  colnames(dfPi)<-sort(donneesphy$Psi)
  dfP<-data.frame(matrix(ncol = dim(donneesphy)[1], nrow = 1))
  colnames(dfP)<-sort(donneesphy$P)
  dfPi[1,]<-lowess(donneesphy$Psi, y= as.numeric(donnees[x,]),f=2/3)$y
  dfP[1,]<-lowess(donneesphy$P, y= as.numeric(donnees[x,]),f=2/3)$y
  
  dfPi<-cbind(x,dfPi)
  colnames(dfPi)[1]<-"AGI"
  dfP<-cbind(x,dfP)
  colnames(dfP)[1]<-"AGI"
  
  df<-gather(dfPi,key="osmotic_potential",value="centered_expression",2:dim(donneesphy)[1]+1)
  gPi<-ggplot(df,aes(x=as.numeric(osmotic_potential), y=centered_expression, group=AGI)) +
    geom_line(aes(colour=AGI), size= 0.2)+
    scale_colour_viridis_d() +
    theme_classic(base_size = 12) +
    ggtitle(x)+
    theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
          axis.text.x = element_text(hjust=1),
          axis.text=element_text(size=6),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size=8, hjust = 0.5),
    )+
    xlab("osmotic potential")+ ylab("centered expression level")
  #+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
  
  df<-gather(dfP,key="turgor_pressure",value="centered_expression",2:dim(donneesphy)[1]+1)
  gP<-ggplot(df,aes(x=as.numeric(turgor_pressure), y=centered_expression, group=AGI)) +
    geom_line(aes(colour=AGI), size= 0.2)+
    scale_colour_viridis_d() +
    theme_classic(base_size = 12) +
    #ggtitle("downregulated")+
    theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
          axis.text.x = element_text(hjust=1),
          axis.title.y = element_blank(),
          axis.text=element_text(size=6),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size=8, hjust = 0.5),
    )+
    xlab("turgor pressure (bars)")#+ ylab("centered expression level")
  #+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
  
  gall<-plot_grid(gPi,NULL,gP,ncol = 3,rel_widths=c(1,-0.05,0.95))
  rm(df)
  return(gall)
}

#calculate the SEM over QPCR data
SEMcalc<-function(x){
  sd(x, na.rm=T)/sqrt(sum(!is.na(x))-1)
}

#plot Affymetrix data (averaged and by agi) with RT-q-PCR data
ggplot.AffyQP<-function(agi){
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(ggpubr)
  longueurdf<-sum(!is.na(yjh.average.data[agi,]))
  df<-data.frame(matrix(ncol = 4, nrow = longueurdf))
  colnames(df)<-c("Affy","QPCR","SEMAffy","SEMQP")
  df[,1]<-unlist(average_agi[agi,colnames(yjh.average.data)[1:longueurdf]])
  df[,2]<-unlist(yjh.average.data[agi,1:longueurdf])
  df[,3]<-unlist(average_agi.SEM[agi,colnames(yjh.average.data)[1:longueurdf]])
  df[,4]<-unlist(yjh.average.data.SEM[agi,1:longueurdf])
  corcoef<-format(round(cor.test(df$Affy,df$QPCR)$p.value, 3), nsmall = 3)
  gAffyQp<-ggplot(df,aes(x=Affy, y=QPCR)) +
    theme_classic()+
    labs(title = paste(agi, ", cor. test. pvalue =",corcoef),x="Affymetrix signal",y="Q-PCR signal")+
    theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
          axis.text.x = element_text(hjust=1),
          axis.text=element_text(size=6),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size=6, hjust = 0.5),
    )+
    stat_regline_equation( aes(label = ..rr.label..), size=2)+
    geom_errorbarh(aes(xmin=Affy-SEMAffy, xmax=Affy+SEMAffy), height=.02, color="dark grey", size=0.1)+
    geom_errorbar(aes(ymin=QPCR-SEMQP, ymax=QPCR+SEMQP), width=.05, color="dark grey", size=0.1)+
    geom_point(size=0.5)+
    geom_smooth(method='lm', se=F, size=0.2, color="dark red")
  
  return(gAffyQp)
}

#creating a heatmap and clustering, saving the information
heatmapcluster<-function(datas,nbclusters){
  packages = c("pheatmap")
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  #print("rows containing NAs are removed from the data")
  datas<-na.omit(datas)
  
  result <- pheatmap(datas,cluster_row=T, scale="row", show_rownames = T,
                     show_colnames = T,#cluster_cols = col_cluster,
                     color=colorRampPalette(c("navy","white","firebrick3"))(100),
                     treeheight_row=80,treeheight_col=30,
                     border_color=NA,fontsize=12,fontsize_row=6,fontsize_col=10,
                     #annotation_col=annotation_col,
                     cutree_rows = nbclusters,angle_col = 45,cellwidth = 30,
                     silent = T)
  geneClust<-result$tree_row
  geneCluster<-cutree(geneClust,k=nbclusters)
  clusters<-list()
  for (i in 1:nbclusters){
    clusters[i]<-as.data.frame(names(geneCluster[geneCluster==i]))
  }
  return(clusters)
}

#plot average agi centered data with colors but with nolabels
plotcolor.osmo.nolabels<-function(datas,probelist,viridisoptionAtoH,lwdth){
  n<-length(probelist)
  require(viridis)
  colpal<-viridisMap(n, option= viridisoptionAtoH)
  for (i in 1:n){
    if (i<2){
      plot(as.numeric(datas[probelist[i],]),xlim = c(0.5,length(datas)+0.5), type = 'l', 
           ylim = c(-3,3), xlab="", ylab="",col=rgb(colpal[i,]),xaxt = "n",yaxt="n",
           lwd=0.1, frame=F)
      axis(1,at=1:16, labels = F, tck = 0.04, lwd = lwdth, lwd.ticks = lwdth)
      axis(2,at = c(-2,0,2),tck=0.04, labels = F, lwd = lwdth,lwd.ticks = lwdth)
      box(lwd=lwdth)
    } else {
      points(as.numeric(datas[probelist[i],]), type = 'l', col = rgb(colpal[i,]),xaxt = "n")
    }
  }
}

#last version of correspondence from probe to AGI in ATH2
Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif<-function(y)
{
  intersect(y,row.names(Probe.AGI.Gene.Atlas.correspondant.sand.may.2016))->x
  as.character(x)->x
  Probe.AGI.Gene.Atlas.correspondant.sand.may.2016[x,1]->out
  as.vector(out)->out
  as.character(out)->out
  out[grep(c('RESCUE'),out)]<-NA
  out[grep(c('UNIDENtIFIED'),out)]<-NA
  out[grep(c('AMBIgUOUS'),out)]<-NA
  out[grep(c('RRNA'),out)]<-NA
  sort(out)->out
  return(out)
}

#function to calculate the maximal p-value (not adjusted) from an anova that respects the FDR given by q_star
p.false.discovery <-function (p, q_star) 
{
  n <- length(p)
  p <- sort(p)
  J = 1:n
  q <- n * p/J
  {
    if (any(q <= q_star)) {
      I <- which(q <= q_star, TRUE)
      p.false.discovery <- p[max(I)]
    }
    else p.false.discovery <- 0
  }
  p.false.discovery
}

#miror function of the above,
#function to calculate the FDR of the anova genes selected with a given p-value (p.value.treshold)
FDR.calculation<-function(p,p.value.treshold) {
  n <- length(p)
  p <- sort(p)
  J = 1:n
  q <- n * p/J
  {
    if (any(p <= p.value.treshold)) {
      I <- which(p <= p.value.treshold, TRUE)
      FDR <- q[max(I)]
    }
    else FDR <- 0
  }
  FDR
}
####END OF FUNCTIONS####

##########################################
########## DATA MINING ###################
##########################################

#create a dataframe with pvalues for the transcriptome, based on probes ID
PvalueID<-Anova.over.data.v.treat.car.pckg(a)
out<-Anova.over.1treat(a,physio$Psi)
out2<-Anova.over.1treat(a,physio$P)
PvalueID<-cbind(PvalueID,out,out2)
colnames(PvalueID)<-c("NaCl.pvalue","Sorbitol.pvalue","PEG.pvalue","EG.pvalue","Pi.pvalue","P.pvalue")

#calculate the average value for each treatment and each probe
average_id<-data.frame(matrix(ncol = length(a)/2, nrow = length(a[,1])))
SD_id<-data.frame(matrix(ncol = length(a)/2, nrow = length(a[,1])))
names(SD_id)<-colonnes
for (j in 1:(length(a)/2)){
  average_id[,j]<-(a[,j]+a[,j+16])/2
  SD_id[,j]<-apply(cbind(a[,j],a[,j+16]),1,SEMcalc)
}
row.names(average_id)<-row.names(a)
row.names(SD_id)<-row.names(a)
colnames(average_id)<-colonnes
average_id<-na.omit(average_id)
average_id_centered<-as.data.frame(t(scale(t(average_id))))

#calculate the average value relative to hydropony
average_rel_id<-data.frame(matrix(ncol = length(a)/2, nrow = length(a[,1])))
for (j in 1:(length(a)/2)){
  average_rel_id[,j]<-(a[,j]/a[,1]+a[,j+16]/a[,17])/2
}
row.names(average_rel_id)<-row.names(a)
colnames(average_rel_id)<-colonnes
#then the log version of it, which should be different from the average of log (the latest is probably wrong)
average_rel.log_id<-log(average_rel_id)

# establish Probe/AGI correspondance and suppress redundancy
#récupération des listes d'AGI à partir des probes sur le fichier Osmo.data.gcrma.txt
agilist<-c()
for (i in 1:length(row.names(a))){
  print(i)
  temp<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(row.names(a)[i])
  if (length(temp)!=0)
  {agilist<-c(agilist,temp)}
  else {agilist<-c(agilist,0)}
}
agilist<-as.data.frame(agilist)
rm(temp)

#calculating the average of signal intensity for genes that are defined by multiple probes
n_occur<-data.frame(table(agilist))
#211 genes have from 2 to 16 probes, most of them show little discrepancy (at least in hydrop)
temp<-cbind(average_id,agilist)
average_agi<-aggregate.data.frame(temp[,1:16],agilist,mean)
average_agi<-subset(average_agi,average_agi$agilist!=0)
row.names(average_agi)<-average_agi$agilist
average_agi<-average_agi[,2:17]
rm(temp)
average_agi_centered<-as.data.frame(t(scale(t(average_agi))))

temp<-cbind(SD_id,agilist)
average_agi.SEM<-aggregate.data.frame(temp[,1:16],agilist,mean)
average_agi.SEM<-subset(average_agi.SEM,average_agi.SEM$agilist!=0)
row.names(average_agi.SEM)<-average_agi.SEM$agilist
average_agi.SEM<-average_agi.SEM[,2:17]
rm(temp)

#######################################
######### FIGURE 1  ###################
#######################################
require(clipr)
require(ggplot2)
require(cowplot)

#View(Fig1_PiP)

plot1a<-ggplot(data=Fig1_kinetic,aes(x=timeNaCl,y=PNaCl/10, group=treatmentNaCl, fill=treatmentNaCl))+
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5,color = "#FF33FF", face="bold"),
        )+
  ggtitle("NaCl")+
  labs(x= "time after treatment (min)", y="turgor pressure (MPa)")+
  geom_vline(xintercept=0, linetype=3, size=0.5, color="darkgrey")+
  geom_point(aes(shape= treatmentNaCl),size=1, stroke=0.2)+
  scale_shape_manual(values=c(21,23,24,25,22))+
  scale_fill_manual(values = c("#FFFFFF","#CC00CC","#FFCCFF","#FF66FF","#FF00FF"))+
  scale_color_manual(values = "#000000")+
  geom_smooth(data=Fig1_kinetic[Fig1_kinetic$timeNaCl>-2,],aes(color=treatmentNaCl),method="loess", se=F, size=0.5)+
  scale_color_manual(values = c("#000000CC","#CC00CCCC","#FFCCFFCC","#FF66FFCC","#FF00FFCC"))+
  scale_y_continuous(limits = c(0,0.58))

plot1b<-ggplot(data=Fig1_kinetic,aes(x=timesorbitol,y=Psorbitol/10, group=treatmentsorbitol, fill=treatmentsorbitol))+
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5,color = "#0000FF", face="bold"),
  )+
  ggtitle("Sorbitol")+
  labs(x= "time after treatment (min)", y="turgor pressure (MPa)")+
  geom_vline(xintercept=0, linetype=3, size=0.5, color="darkgrey")+
  geom_point(aes(shape= treatmentsorbitol),size=1, stroke=0.2)+
  scale_shape_manual(values=c(21,23,24,25,22))+
  scale_fill_manual(values = c("#FFFFFF","#3333FF","#0000FF","#6666FF"))+
  scale_color_manual(values = "#000000")+
  geom_smooth(data=Fig1_kinetic[Fig1_kinetic$timesorbitol>-2,],aes(color=treatmentsorbitol),method="loess", span=0.6,se=F, size=0.5)+
  scale_color_manual(values = c("#000000CC","#3333FFCC","#0000FFCC","#6666FFCC"))+
  scale_y_continuous(limits = c(0,0.58))

plot1c<-ggplot(data=Fig1_kinetic,aes(x=timePEG,y=PPEG/10, group=treatmentPEG, fill=treatmentPEG))+
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5,color = "#66FF66", face="bold"),
  )+
  ggtitle("PEG")+
  labs(x= "time after treatment (min)", y="turgor pressure (MPa)")+
  geom_vline(xintercept=0, linetype=3, size=0.5, color="darkgrey")+
  geom_point(aes(shape= treatmentPEG),size=1, stroke=0.2)+
  scale_shape_manual(values=c(21,23,24,25,22))+
  scale_fill_manual(values = c("#FFFFFF","#66CC66","#339933","#006600","#99FF99"))+
  scale_color_manual(values = "#000000")+
  geom_smooth(data=Fig1_kinetic[Fig1_kinetic$timePEG>-2,],aes(color=treatmentPEG),method="loess",span=0.8, se=F, size=0.5)+
  scale_color_manual(values = c("#000000CC","#66CC66CC","#339933CC","#006600CC","#99FF99CC"))+
  scale_y_continuous(limits = c(0,0.58))

plot1d<-ggplot(data=Fig1_kinetic,aes(x=timeEG,y=PEG/10, group=treatmentEG, fill=treatmentEG))+
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5,color = "#FF0000", face="bold"),
  )+
  ggtitle("EG")+
  labs(x= "time after treatment (min)", y="turgor pressure (MPa)")+
  geom_vline(xintercept=0, linetype=3, size=0.5, color="darkgrey")+
  geom_point(aes(shape= treatmentPEG),size=1, stroke=0.2)+
  scale_shape_manual(values=c(21,23,24,25,22))+
  scale_fill_manual(values = c("#FF6666","#FF3333","#FF0000","#FF9999","#FFFFFF"))+
  scale_color_manual(values = "#000000")+
  geom_smooth(data=Fig1_kinetic[Fig1_kinetic$timeEG>-2,],aes(color=treatmentEG),method="loess", span=0.6,se=F, size=0.5)+
  scale_color_manual(values = c("#FF6666CC","#FF3333CC","#FF0000CC","#FF9999CC","#000000CC"))+
  scale_y_continuous(limits = c(0,0.58))

plot1e<-ggplot(data = Fig1_PiP,aes(x=`Pi(Mpa)`,y=`P(Mpa)`, group=osmoticum, fill=osmoticum))+
  geom_point()+
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  geom_errorbar(aes(ymin= `P(Mpa)`-`SEMP(Mpa)`,ymax=`P(Mpa)`+`SEMP(Mpa)`,color=osmoticum),size=0.5, width=0.01)+
  geom_point(aes(shape=osmoticum),size=1.5, stroke=0.1)+
  scale_shape_manual(values=c(21,23,24,22))+
  scale_fill_manual(values = c("#FF0000","#99FF99","#FF00FF","#0000FF"))+
  #scale_color_manual(values = "#000000")+
  geom_line(aes(color=osmoticum))+
  scale_color_manual(values = c("#FF0000CC","#99FF99CC","#FF00FFCC","#0000FFCC"))
  
figure1<-plot_grid(plot1a,plot1b,plot1c,plot1d,plot1e, NULL, ncol=2)

save_plot("2211213-fig1.svg",figure1,base_height=9.5, base_asp=0.75)

#######################################
#######   FIGURE 2   ##################
#######################################
require(pheatmap)
require(nVennR)
require(cowplot)
require(ggplot2)
require(ggupset)
require(tidyr)
require(dplyr)
require(viridis)

#Figure 2a
d<-dist(t(average_agi_centered))
hc<-hclust(d)
plot2a<-function()
{ par(mar=c(2,3,2,1))
  plot(hc, xlab="",ylab="",main="")}

#p-value selection based on a FDR <0.2 on transcriptome data
#for the anova model on solutes:
thresholdsolutes<-p.false.discovery(unlist(PvalueID[,1:4]),0.2)#0.001
#for the one-way anova on Pi and P:
thresholdPi<-p.false.discovery(PvalueID[,5],0.2)#0.0004
thresholdP<-p.false.discovery(PvalueID[,6],0.2)#0.0012

#fig2data<-t(matlabPValue_ID[,1:4]<thresholdsolutes)
fig2data<-t(PvalueID[,1:4]<thresholdsolutes)
rownames(fig2data)<-c("NaCl","Sorbitol","PEG","EG")

tidy_fig2data<-fig2data %>%
  as_tibble(rownames = "Treatment") %>%
  gather(Gene, Member, -Treatment) %>%
  filter(Member) %>%
  select(- Member)
tidy_fig2data %>%
  group_by(Gene) %>%
  summarize(Treatments = list(Treatment))
#not a figure but allow to extract the number of probes per solute or shared by solutes
tfig2data<-PvalueID[,1:4]<thresholdsolutes
myvprobessolutes<-plotVenn(list(`NaCl (1)`=row.names(tfig2data[tfig2data[,1]==TRUE,]), `Sorbitol (2)`=row.names(tfig2data[tfig2data[,2]==TRUE,]), 
                                `PEG (3)`=row.names(tfig2data[tfig2data[,3]==TRUE,]), `EG (4)`=row.names(tfig2data[tfig2data[,4]==TRUE,])), nCycles = 15000, showPlot = F)
showSVG(myvprobessolutes, opacity=0.7, outFile = "230525-Venn_4_treatments_probes.svg")

#Figure 2b
plot2b<-tidy_fig2data %>%
  group_by(Gene) %>%
  summarize(Treatments = list(Treatment)) %>%
  ggplot(aes(x = Treatments)) +
  theme_classic()+
  theme(axis.ticks.length=unit(-0.051, "cm"))+
  labs(y="differentialy expressed probes")+
  geom_bar() +
  scale_x_upset()

#Figure 2C
# let's proceed to select the ID below the treshold, for NaCl, Sorbitol, PEG and EG
#and translate the information of Affy probes into AGI
solutes.probes<-c()
for (i in 1:4){
  solutes.probes[i]<-as.data.frame(subset(rownames(PvalueID), PvalueID[,i]<thresholdsolutes))
}
solutes.genes<-c()
for (i in 1:4){
  tempprobe<-as.data.frame(solutes.probes[[i]])
  templist<-c()
  for (j in 1:dim(tempprobe)[1]){
    temp<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(tempprobe[[1]][j])
    if (length(temp)!=0)
    {templist<-c(templist,temp)}
  }
  templist<-unique(templist)
  solutes.genes[[i]]<-templist
  rm(temp)
  rm(templist)
  rm(tempprobe)
}

#retrieve specific regions by venn diagram 
#cluster them in 2 for up- and down-regulation with pheatmap
myV2<-plotVenn(list(`NaCl (1)`=solutes.genes[[1]], `Sorbitol (2)`=solutes.genes[[2]], `PEG (3)`=solutes.genes[[3]], `EG (4)`=solutes.genes[[4]]), nCycles = 15000, showPlot = F)
showSVG(myV2, opacity=0.7, outFile = "230525-fig2c_Venn_4_treatments_new-tresholds.svg")#"221122-fig2c_Venn_4_treatments.svg"
liste<-data.frame(c(1,0,0,0), c(0,1,0,0),c(0,0,1,0), c(0,0,0,1))
plotlist<-{}
k<-1
for (i in 1:4){
  print(i)
  region<-getVennRegion(myV2,liste[[i]])
  if (length(region)>1){
    tempclusters<-heatmapcluster(average_agi_centered[region,],2)
    print(length(tempclusters[[1]]))
    print(length(tempclusters[[2]]))
    for (j in 1:2){
      n<-length(tempclusters[[j]])
      cluster<-toString(c("cluster ",j))
      ngenes<-toString((c(n," genes")))
      df1<-average_agi_centered[tempclusters[[j]],]
      df1<-cbind(rownames(df1),df1)
      colnames(df1)[1]<-"AGI"
      df<-gather(df1,key="treatment",value="centered_expression",2:17)
      df%>% mutate(treatment=factor(treatment,levels = colonnes))
      hydroNaCl<-df[1:max(grep("NaCl",df$treatment)),]
      sorbitol<-df[min(grep("Sorb",df$treatment)):max(grep("Sorb",df$treatment)),]
      PEG<-df[min(grep("PEG",df$treatment)):max(grep("PEG",df$treatment)),]
      EG<-df[(max(grep("PEG",df$treatment))+1):max(grep("EG",df$treatment)),]
      plotlist[[k]]<-ggplot(df,aes(x=factor(treatment,level = colonnes), y=centered_expression, group=AGI)) +
        geom_point(aes(colour=AGI),alpha=0)+
        geom_line(data= hydroNaCl,aes(colour=AGI), size= 0.2, alpha=1)+#x=factor(treatment,level = c("hydro","NaCl25","NaCl50","NaCl75","NaCl100","NaCl125")),
        geom_line(data= sorbitol,aes(colour=AGI), size= 0.2)+
        geom_line(data= PEG,aes(colour=AGI), size= 0.2)+
        geom_line(data= EG,aes(colour=AGI), size= 0.2)+
        scale_colour_viridis_d() +
        theme_classic(base_size = 12) +
        theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
          axis.text.x = element_blank(),#element_text(angle = 45, hjust=1),
          axis.text=element_text(size=5),
          axis.title = element_text(size = 7),
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          plot.title = element_text(size=8, hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
        )+
        xlab("") + ylab("")
      k<-k+1
    rm(df1)
    rm(hydroNaCl)
    rm(sorbitol)
    rm(PEG)
    rm(EG)
    }
    rm(tempclusters)
  }
  rm(region)
  }
plot2c<-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],
                  plotlist[[6]],plotlist[[7]],plotlist[[8]],ncol=2)
#plot2c
fig2ac<-plot_grid(plot2a,plot2c,ncol=1,rel_heights = c(0.7,1))
fig2bnull<-plot_grid(plot2b,NULL,ncol = 1,rel_heights = c(1,0.5))
figure2<-plot_grid(fig2ac,fig2bnull, ncol=2)
figure2
save_plot("230525-Fig2-new-tresholds.svg",figure2,base_height = 7.5, base_width = 7.5)#"230428-Fig2.svg"
rm(plotlist)


#########################################
#######   FIGURE 3     ##################
#########################################
require(tidyr)
require(dplyr)
require(ggplot2)
require(viridis)
require(cowplot)
require(nVennR)

#retrieving DEPs for Pi from the anova, and splitting into 2 for up- and down-regulation
temp<-rownames(PvalueID[PvalueID[,5]<thresholdPi,])
Pi.only.probe<-heatmapcluster(a[temp,],2)
rm(temp)
#retrieving DEPs for P from the anova, and splitting into 2 for up- and down-regulation
temp<-rownames(PvalueID[PvalueID[,6]<thresholdP,])
P.only.probe<-heatmapcluster(a[temp,],2)
rm(temp)
#transforming list of DEPs into list of DEGs
Pi.only.agi<-c()
Pi.only.agi[[1]]<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(Pi.only.probe[[1]])
Pi.only.agi[[2]]<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(Pi.only.probe[[2]])
P.only.agi<-c()
P.only.agi[[1]]<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(P.only.probe[[1]])
P.only.agi[[2]]<-Probe.to.AGI.Gene.Atlas.may.2016.cleaner.gb.no.notif(P.only.probe[[2]])

#lists edited to remove genes that have multiple probes which are giving different p-values
problematigenesPi<-c("At1g72850", "At1g78270", "At2g24540", "At2g33810", "At4g13920", "At4g24410", "At4g25880", "At4g26490" )
length(unlist(Pi.only.agi))
length(unlist(Pi.only.agi)[unlist(Pi.only.agi) %in% problematigenesPi==F])
Pi.only.agi<-unlist(Pi.only.agi)[unlist(Pi.only.agi) %in% problematigenesPi==F]

problematigenesP<-c("At1g07130","At1g07725","At1g08590","At1g51640","At1g56240","At1g72850","At2g11851",
                    "At2g22960", "At3g22070", "At3g54630", "At3g56770" ,"At4g24410", "At4g28650", "At4g36030",
                    "At4g38210","At4g38550", "At5g59730")
length(unlist(P.only.agi))
length(unlist(P.only.agi)[unlist(P.only.agi) %in% problematigenesP==F])
P.only.agi<-unlist(P.only.agi)[unlist(P.only.agi) %in% problematigenesP==F]

####FIGURE 3A
#venn diagram for these lists, as plot 3a
plot3a<-plotVenn(list("Pi (1)"=unlist(Pi.only.agi),"P (2)"=unlist(P.only.agi)), #bquote(Pi ~ " (1)")
               nCycles = 5500, showPlot = F)
showSVG(plot3a, opacity=0.8, outFile = "230525_Fig_3A-Venn_PPi_only_new-tresholds.svg")

# Retrieving specific genes from the Venn diagram above
region<-getVennRegion(plot3a,c(1,0))
tempclusters<-heatmapcluster(average_agi_centered[region,],2)
Pi.specific.agi<-c()
Pi.specific.agi[[1]]<-tempclusters[[1]]
Pi.specific.agi[[2]]<-tempclusters[[2]]
rm(region)
rm(tempclusters)
region<-getVennRegion(plot3a,c(0,1))
tempclusters<-heatmapcluster(average_agi_centered[region,],2)
P.specific.agi<-c()
P.specific.agi[[1]]<-tempclusters[[1]]
P.specific.agi[[2]]<-tempclusters[[2]]

# Figure 3B, gene expression in the different treatments
liste<-data.frame(c(1,0), c(0,1))
plotlist<-{}
k<-1
for (i in 1:2){
  print(i)
  region<-getVennRegion(plot3a,liste[[i]])
  if (length(region)>1){
    tempclusters<-heatmapcluster(average_agi_centered[region,],2)
    for (j in 1:2){
      n<-length(tempclusters[[j]])
      cluster<-toString(c("cluster ",j))
      ngenes<-toString((c(n," genes")))
      df1<-average_agi_centered[tempclusters[[j]],]
      df1<-cbind(rownames(df1),df1)
      colnames(df1)[1]<-"AGI"
      df<-gather(df1,key="treatment",value="centered_expression",2:17)
      df%>% mutate(treatment=factor(treatment,levels = colonnes))
      hydroNaCl<-df[1:max(grep("NaCl",df$treatment)),]
      sorbitol<-df[min(grep("Sorb",df$treatment)):max(grep("Sorb",df$treatment)),]
      PEG<-df[min(grep("PEG",df$treatment)):max(grep("PEG",df$treatment)),]
      EG<-df[(max(grep("PEG",df$treatment))+1):max(grep("EG",df$treatment)),]
      plotlist[[k]]<-ggplot(df,aes(x=factor(treatment,level = colonnes), y=centered_expression, group=AGI)) +
        geom_point(aes(colour=AGI),alpha=0)+
        geom_line(data= hydroNaCl,aes(colour=AGI), size= 0.2, alpha=1)+
        geom_line(data= sorbitol,aes(colour=AGI), size= 0.2)+
        geom_line(data= PEG,aes(colour=AGI), size= 0.2)+
        geom_line(data= EG,aes(colour=AGI), size= 0.2)+
        scale_colour_viridis_d() +
        theme_classic(base_size = 12) +
        theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
              axis.text.x = element_blank(),#element_text(angle = 45, hjust=1),
              axis.text=element_text(size=5),
              axis.title = element_text(size = 7),
              axis.line = element_line(size = 0.2),
              axis.ticks = element_line(size = 0.2),
              plot.title = element_text(size=8, hjust = 0.5),
              plot.margin = unit(c(0, 0, 0, 0), "cm"),
        )+
        xlab("") + ylab("")
      k<-k+1
      rm(df1)
      rm(hydroNaCl)
      rm(sorbitol)
      rm(PEG)
      rm(EG)
    }
    rm(tempclusters)
  }
  rm(region)
}
xaxe=c("0","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200")
plotlist[[3]]<-plotlist[[3]]+
  scale_x_discrete(labels = xaxe)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
plotlist[[4]]<-plotlist[[4]]+
  scale_x_discrete(labels = xaxe)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plot3B<-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],ncol=2 )
plot3B
save_plot("230525-Fig_3B_new-tresholds.svg",plot3B,base_height = 3.54,base_asp = 1)
rm(plotlist)

# Figure 3C, gene expression f(Pi or P) for cluster 1
liste1<-Pi.specific.agi[[1]]
df1Pi<-data.frame(matrix(ncol = 16, nrow = length(liste1)))
colnames(df1Pi)<-sort(physio$Psi[1:16])
df1P<-data.frame(matrix(ncol = 16, nrow = length(liste1)))
colnames(df1P)<-sort(physio$P[1:16])
for (i in 1:length(liste1)){
  df1Pi[i,]<-lowess(physio[1:16,1], y= as.numeric(average_agi_centered[liste1[i],]),f=2/3)$y
  df1P[i,]<-lowess(physio[1:16,2], y= as.numeric(average_agi_centered[liste1[i],]),f=2/3)$y
}

df1Pi<-cbind(liste1,df1Pi)
colnames(df1Pi)[1]<-"AGI"
df1P<-cbind(liste1,df1P)
colnames(df1P)[1]<-"AGI"

df<-gather(df1Pi,key="osmotic_potential",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g5<-ggplot(df,aes(x=as.numeric(osmotic_potential), y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "magma") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("osmotic potential (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g5

df<-gather(df1P,key="turgor_pressure",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g6<-ggplot(df,aes(x=as.numeric(turgor_pressure)/10, y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "magma") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.title.y = element_blank(),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("turgor pressure (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g6

liste2<-Pi.specific.agi[[2]]
df2Pi<-data.frame(matrix(ncol = 16, nrow = length(liste2)))
colnames(df2Pi)<-sort(physio$Psi[1:16])
df2P<-data.frame(matrix(ncol = 16, nrow = length(liste2)))
colnames(df2P)<-sort(physio$P[1:16])
for (i in 1:length(liste2)){
  df2Pi[i,]<-lowess(physio[1:16,1], y= as.numeric(average_agi_centered[liste2[i],]),f=2/3)$y
  df2P[i,]<-lowess(physio[1:16,2], y= as.numeric(average_agi_centered[liste2[i],]),f=2/3)$y
}

df2Pi<-cbind(liste2,df2Pi)
colnames(df2Pi)[1]<-"AGI"
df2P<-cbind(liste2,df2P)
colnames(df2P)[1]<-"AGI"

df<-gather(df2Pi,key="osmotic_potential",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g7<-ggplot(df,aes(x=as.numeric(osmotic_potential), y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "inferno") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("osmotic potential (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g7

df<-gather(df2P,key="turgor_pressure",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g8<-ggplot(df,aes(x=as.numeric(turgor_pressure)/10, y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "inferno") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("turgor pressure (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g8

plot3cup<-plot_grid(g5,g6)
title3cup<-ggdraw()+
  draw_label(bquote(Pi ~ "-specific cluster 1," ~ .(dim(df1Pi)[1]) ~ "genes"), size = 8)
plot3cup<-plot_grid(title3cup, plot3cup, ncol=1,rel_heights = c(0.06,1))
plot3cdown<-plot_grid(g7,g8)
title3cdown<-ggdraw()+
  draw_label(bquote(Pi ~ "-specific cluster 2," ~ .(dim(df2Pi)[1])~ "genes"), size = 8)
plot3cdown<-plot_grid(title3cdown, plot3cdown, ncol=1,rel_heights = c(0.06,1))
plot3c<-plot_grid(plot3cup,plot3cdown,ncol=1)
save_plot("230525-Fig_3C_new-tresholds.svg",plot3c,base_height = 3.54, base_asp = 1)
rm(df)


#Figure 3D
liste3<-P.specific.agi[[1]]
df3Pi<-data.frame(matrix(ncol = 16, nrow = length(liste3)))
colnames(df3Pi)<-sort(physio$Psi[1:16])
df3P<-data.frame(matrix(ncol = 16, nrow = length(liste3)))
colnames(df3P)<-sort(physio$P[1:16])
for (i in 1:length(liste3)){
  df3Pi[i,]<-lowess(physio[1:16,1], y= as.numeric(average_agi_centered[liste3[i],]),f=2/3)$y
  df3P[i,]<-lowess(physio[1:16,2], y= as.numeric(average_agi_centered[liste3[i],]),f=2/3)$y
}

df3Pi<-cbind(liste3,df3Pi)
colnames(df3Pi)[1]<-"AGI"
df3P<-cbind(liste3,df3P)
colnames(df3P)[1]<-"AGI"

df<-gather(df3Pi,key="osmotic_potential",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g9<-ggplot(df,aes(x=as.numeric(osmotic_potential), y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "mako") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("osmotic potential (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g9

df<-gather(df3P,key="turgor_pressure",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g10<-ggplot(df,aes(x=as.numeric(turgor_pressure)/10, y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "mako") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.title.y = element_blank(),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("turgor pressure (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g10

liste4<-P.specific.agi[[2]]
df4Pi<-data.frame(matrix(ncol = 16, nrow = length(liste4)))
colnames(df4Pi)<-sort(physio$Psi[1:16])
df4P<-data.frame(matrix(ncol = 16, nrow = length(liste4)))
colnames(df4P)<-sort(physio$P[1:16])
for (i in 1:length(liste4)){
  df4Pi[i,]<-lowess(physio[1:16,1], y= as.numeric(average_agi_centered[liste4[i],]),f=2/3)$y
  df4P[i,]<-lowess(physio[1:16,2], y= as.numeric(average_agi_centered[liste4[i],]),f=2/3)$y
}

df4Pi<-cbind(liste4,df4Pi)
colnames(df4Pi)[1]<-"AGI"
df4P<-cbind(liste4,df4P)
colnames(df4P)[1]<-"AGI"

df<-gather(df4Pi,key="osmotic_potential",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g11<-ggplot(df,aes(x=as.numeric(osmotic_potential), y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "rocket") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("osmotic potential (MPa)")+ ylab("centered expression level")
#+ scale_x_discrete(labels= c("hydro","25","50","75","100","50","100","150","75","100","125","150","50","100","150","200"))
g11

df<-gather(df4P,key="turgor_pressure",value="centered_expression",2:17)
#df%>% mutate(treatment=factor(treatment,levels = colonnes))
g12<-ggplot(df,aes(x=as.numeric(turgor_pressure)/10, y=centered_expression, group=AGI)) +
  geom_line(aes(colour=AGI), size= 0.2)+
  scale_colour_viridis_d(option = "rocket") +
  theme_classic(base_size = 12) +
  #ggtitle("downregulated")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=1),
        axis.text=element_text(size=6),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=8, hjust = 0.5),
  )+
  xlab("turgor pressure (MPa)")+ ylab("centered expression level")

g12

plot3dup<-plot_grid(g9,g10)
title3dup<-ggdraw()+
  draw_label(bquote("P -specific cluster 1," ~ .(dim(df3Pi)[1]) ~ "genes"), size = 8)
plot3dup<-plot_grid(title3dup, plot3dup, ncol=1,rel_heights = c(0.06,1))
plot3ddown<-plot_grid(g11,g12)
title3ddown<-ggdraw()+
  draw_label(bquote("P -specific cluster 2," ~ .(dim(df4Pi)[1])~ "genes"), size = 8)
plot3ddown<-plot_grid(title3ddown, plot3ddown, ncol=1,rel_heights = c(0.06,1))
plot3d<-plot_grid(plot3dup,plot3ddown,ncol=1)
save_plot("230525-Fig_3D_new-treshold.svg",plot3d,base_height = 3.54, base_asp = 1)


######################################################
#########Figure 4: mRNA decay halftime ###############
######################################################
###analysis of decay rate of the candidates based on Sorenson et al. 2018

#calculating the T1/2 according to Sorenson et al. 2018
refdecay$halftime<-log(2)/refdecay$alpha_sov
rownames(refdecay)<-refdecay$AGI
median(refdecay$halftime)

bootstrapvalues<-vector(length=4)
#statistics for T1/2 for each cluster
t12Picl1<-as.matrix(sapply(Pi.specific.agi[[1]], toupper))
median((unlist(refdecay[t12Picl1,3])), na.rm = TRUE)
bootstrapvalues[1]<-bootstrap.decay(t12Picl1)
t12Picl2<-as.matrix(sapply(Pi.specific.agi[[2]], toupper))
median((unlist(refdecay[t12Picl2,3])), na.rm = TRUE)
bootstrapvalues[2]<-bootstrap.decay(t12Picl2)
t12Pcl1<-as.matrix(sapply(P.specific.agi[[1]], toupper))
median((unlist(refdecay[t12Pcl1,3])), na.rm = TRUE)
bootstrapvalues[3]<-bootstrap.decay(t12Pcl1)
t12Pcl2<-as.matrix(sapply(P.specific.agi[[2]], toupper))
median((unlist(refdecay[t12Pcl2,3])), na.rm = TRUE)
bootstrapvalues[4]<-bootstrap.decay(t12Pcl2)

#comparing the T1/2 of Pi and P clusters to a set of random genes of the same study
#making 10^4 comparisons
#use the code above to provide to the variable "testupper" the proper list
bootstrap.decay<-function(atester){
  df<-refdecay[unlist(atester),3]
  df<-df[is.finite(df)]
    #median(df)
    #median(refdecay$halftime)
  compteur=0
  for (i in 1:10000){
    agitest<-sample(refdecay$AGI,length(df))
    dfsampledecay<-refdecay[agitest,3]
    dfsampledecay<-dfsampledecay[is.finite(dfsampledecay)]
    if (median(dfsampledecay)<median(df)){compteur<-compteur+1}
  }
  return(compteur/10000)
  }

testupper<-as.matrix(sapply(P.specific.agi[[1]], toupper))
df<-as.matrix(sapply(testupper, function(i) refdecay[(i),3]))
df<-as.data.frame(df[is.finite(df)])
df[,2]<-rep("P_cluster_1",dim(df)[1])
colnames(df)<-c("half_time","cluster")

testupper<-as.matrix(sapply(P.specific.agi[[2]], toupper))
df2<-as.matrix(sapply(testupper, function(i) refdecay[(i),3]))
df2<-as.data.frame(df2[is.finite(df2)])
df2[,2]<-rep("P_cluster_2",dim(df2)[1])
colnames(df2)<-c("half_time","cluster")

testupper<-as.matrix(sapply(Pi.specific.agi[[1]], toupper))
df3<-as.matrix(sapply(testupper, function(i) refdecay[(i),3]))
df3<-as.data.frame(df3[is.finite(df3)])
df3[,2]<-rep("Pi_cluster_1",dim(df3)[1])
colnames(df3)<-c("half_time","cluster")

testupper<-as.matrix(sapply(Pi.specific.agi[[2]], toupper))
df4<-as.matrix(sapply(testupper, function(i) refdecay[(i),3]))
df4<-as.data.frame(df4[is.finite(df4)])
df4[,2]<-rep("Pi_cluster_2",dim(df4)[1])
colnames(df4)<-c("half_time","cluster")

dftot<-rbind(df3,df4,df,df2)
dftot$half_time<-as.numeric(dftot$half_time)
dftot$cluster<-factor(dftot$cluster, c("Pi_cluster_1","Pi_cluster_2","P_cluster_1","P_cluster_2"))
dfbootstrap<-as.data.frame(cbind(c("Pi_cluster_1","Pi_cluster_2","P_cluster_1","P_cluster_2"),bootstrapvalues))
colnames(dfbootstrap)[1]<-"cluster"

plotfig4<-ggplot(dftot,aes(x=cluster,y=half_time))+
  geom_boxplot(lwd=0.3,outlier.shape=NA)+
  #geom_dotplot(binaxis='x',stackdir='center', dotsize=1)+
  geom_jitter(size=0.2, alpha=0.3)+
  coord_cartesian(ylim=c(0,600))+
  theme_classic()+
  labs(x="",y="mRNA decay half-time (min)")+
  theme(axis.ticks.length=unit(-0.051, "cm"),legend.position = "none",
        axis.text.x = element_text(hjust=0.5),
        axis.text=element_text(size=6),
        axis.title = element_text(size = 5),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(size=6, hjust = 0.5),
  )+
  scale_x_discrete(labels = c(expression(Pi ~" cluster 1"),expression(Pi~" cluster 2"),"P cluster 1", "P cluster 2"))+
  geom_label(data=dfbootstrap,aes(x=cluster,y=580,label = bootstrapvalues), size=2, label.size = NA, alpha=0.5)
plotfig4
save_plot("230526-fig4_mRNA_half-time_new-tresholds.svg",plotfig4,base_width = 3, base_height = 2)

#########################################
########## FIGURE 5 #####################
#########################################

#double check if data are loaded, if not go to the import section at the begining of the document
length(yjh.tested.agi[,1])
#quick graph for visualization of the Affy data for the genes tested by RT-QPCR
myplots<-lapply(yjh.tested.agi[,1],ggplot.lowess.physio.agi,physio[1:16,],average_agi_centered)
gall<-plot_grid(plotlist = myplots,ncol = 2)
gall
rm(gall)

#creating the dataframes with average and SEM for the RT-qPCR data
yjh.average.data<-data.frame(matrix(nrow = length(yjh_QPCR_data$AGI),ncol = 12))
colnames(yjh.average.data)<-c("hydro","NaCl25","NaCl50","NaCl75","NaCl100","Sorb50","Sorb100","Sorb150","EG50","EG100","EG150","EG200")
yjh.average.data$hydro<-rowSums(yjh_QPCR_data[,grep("control",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("control",colnames(yjh_QPCR_data))]))
yjh.average.data$NaCl25<-rowSums(yjh_QPCR_data[,grep("NaCl25",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("NaCl25",colnames(yjh_QPCR_data))]))
yjh.average.data$NaCl50<-rowSums(yjh_QPCR_data[,grep("NaCl50",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("NaCl50",colnames(yjh_QPCR_data))]))
yjh.average.data$NaCl75<-rowSums(yjh_QPCR_data[,grep("NaCl75",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("NaCl75",colnames(yjh_QPCR_data))]))
yjh.average.data$NaCl100<-rowSums(yjh_QPCR_data[,grep("NaCl100",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("NaCl100",colnames(yjh_QPCR_data))]))
yjh.average.data$Sorb50<-rowSums(yjh_QPCR_data[,grep("Sorb50",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("Sorb50",colnames(yjh_QPCR_data))]))
yjh.average.data$Sorb100<-rowSums(yjh_QPCR_data[,grep("Sorb100",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("Sorb100",colnames(yjh_QPCR_data))]))
yjh.average.data$Sorb150<-rowSums(yjh_QPCR_data[,grep("Sorb150",colnames(yjh_QPCR_data))], na.rm = T)/rowSums(!is.na(yjh_QPCR_data[,grep("Sorb150",colnames(yjh_QPCR_data))]))
yjh.average.data$EG50<-rowSums(yjh_QPCR_data[,grep("EG50",colnames(yjh_QPCR_data))], na.rm = T)/3
yjh.average.data$EG100<-rowSums(yjh_QPCR_data[,grep("EG100",colnames(yjh_QPCR_data))], na.rm = T)/3
yjh.average.data$EG150<-rowSums(yjh_QPCR_data[,grep("EG150",colnames(yjh_QPCR_data))], na.rm = T)/3
yjh.average.data$EG200<-rowSums(yjh_QPCR_data[,grep("EG200",colnames(yjh_QPCR_data))], na.rm = T)/3
rownames(yjh.average.data)<-yjh_QPCR_data$AGI

yjh.average.data.SEM<-data.frame(matrix(nrow = length(yjh_QPCR_data$AGI),ncol = 12))
colnames(yjh.average.data.SEM)<-c("hydro","NaCl25","NaCl50","NaCl75","NaCl100","Sorb50","Sorb100","Sorb150","EG50","EG100","EG150","EG200")
yjh.average.data.SEM$hydro<-apply(yjh_QPCR_data[,grep("control",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$NaCl25<-apply(yjh_QPCR_data[,grep("NaCl25",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$NaCl50<-apply(yjh_QPCR_data[,grep("NaCl50",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$NaCl75<-apply(yjh_QPCR_data[,grep("NaCl75",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$NaCl100<-apply(yjh_QPCR_data[,grep("NaCl100",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$Sorb50<-apply(yjh_QPCR_data[,grep("Sorb50",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$Sorb100<-apply(yjh_QPCR_data[,grep("Sorb100",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$Sorb150<-apply(yjh_QPCR_data[,grep("Sorb150",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$EG50<-apply(yjh_QPCR_data[,grep("EG50",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$EG100<-apply(yjh_QPCR_data[,grep("EG100",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$EG150<-apply(yjh_QPCR_data[,grep("EG150",colnames(yjh_QPCR_data))],1,SEMcalc)
yjh.average.data.SEM$EG200<-apply(yjh_QPCR_data[,grep("EG200",colnames(yjh_QPCR_data))],1,SEMcalc)
rownames(yjh.average.data.SEM)<-yjh_QPCR_data$AGI

#plots Figure 4 are saved individually because it is easier to assemble them in Inkscape
myplots<-lapply(c("At3g14440","At4g36550","At5g44060"),ggplot.AffyQP)
gall<-plot_grid(plotlist=myplots, ncol=3)
save_plot("230525-plotsAffy-QPCR_Picluster1_new-tresholds.svg",gall,base_width = 7, base_height = 2)

myplots<-lapply(c("At3g45210","At4g37580","At3g26510"),ggplot.AffyQP)
gall<-plot_grid(plotlist=myplots, ncol=3)
save_plot("220930-plotsAffy-QPCR_Picluster2.svg",gall,base_width = 7, base_height = 2)

myplots<-lapply(c("At1g64640","At1g78260","At1g55700"),ggplot.AffyQP)
gall<-plot_grid(plotlist=myplots, ncol=3)
save_plot("220930-plotsAffy-QPCR_Pcluster1.svg",gall,base_width = 7, base_height = 2)

myplots<-lapply(c("At3g50260","At5g46710","At4g17490"),ggplot.AffyQP)
gall<-plot_grid(plotlist=myplots, ncol=3)
save_plot("220930-plotsAffy-QPCR_Pcluster2.svg",gall,base_width = 7, base_height = 2)
rm(gall)
rm(myplots)

#################################################
######  FIGURE 6b, GO enrichment ################
#################################################
#note: figure 6a is the direct output of the genesecloud website
#note 2: figure 6c is the direct output of the genesect function from Virtual Plant

require(clusterProfiler)

#Pi correlated specific analysis, cluster 1
gene<-Pi.specific.agi[[1]]
print(length(gene))
ego <- enrichGO(gene = toupper(gene), OrgDb="org.At.tair.db", ont = "ALL", keyType = "TAIR",
                pAdjustMethod = "none", pvalueCutoff=0.05, qvalueCutoff = 0.05)
print(dim(ego)[1])
if (dim(ego)[1]!=0){
  plot6b1<-dotplot(ego, showCategory=10, font.size = 8, title = bquote(Pi ~ " cluster 1"))
}
rm(gene)

#Pi correlated specific analysis, cluster 2
gene<-Pi.specific.agi[[2]]
print(length(gene))
ego <- enrichGO(gene = toupper(gene), OrgDb="org.At.tair.db", ont = "ALL", keyType = "TAIR",
                pAdjustMethod = "none", pvalueCutoff=0.05, qvalueCutoff = 0.05)
print(dim(ego)[1])
if (dim(ego)[1]!=0){
  plot6b2<-dotplot(ego, showCategory=10, font.size = 8, title = bquote(Pi ~ " cluster 2"))
}
rm(gene)

#P correlated specific analysis, cluster 1
gene<-P.specific.agi[[1]]
print(length(gene))
ego <- enrichGO(gene = toupper(gene), OrgDb="org.At.tair.db", ont = "ALL", keyType = "TAIR",
                pAdjustMethod = "none", pvalueCutoff=0.05, qvalueCutoff = 0.05)
print(dim(ego)[1])
if (dim(ego)[1]!=0){
  plot6b3<-dotplot(ego, showCategory=10, font.size = 8, title = "P cluster 1")
}
rm(gene)

#P correlated specific analysis, cluster 2
gene<-P.specific.agi[[2]]
print(length(gene))
ego <- enrichGO(gene = toupper(gene), OrgDb="org.At.tair.db", ont = "ALL", keyType = "TAIR",
                pAdjustMethod = "none", pvalueCutoff=0.05, qvalueCutoff = 0.05)
print(dim(ego)[1])
if (dim(ego)[1]!=0){
  plot6b4<-dotplot(ego, showCategory=10, font.size = 8, title = "P cluster 2")
}
rm(gene)

plot6b<-plot_grid(plot6b1,plot6b4,ncol = 2)
save_plot("230526-fig6b-GO enrichment.svg",plot6b,base_width = 7, base_height = 2)



##################################################################
######supplemental file S4: Pi and P specific genes list ranking##
##################################################################
# rank by the highest absolute slope of the Pi or P to gene expression, the R square and the average expression
Pi.specific.ranking<-c()
for (i in 1:2){
  Pi.specific.ranking[[i]]<-data.frame(matrix(ncol = 3, nrow = length(Pi.specific.agi[[i]])))
  rownames(Pi.specific.ranking[[i]])<-Pi.specific.agi[[i]]
  colnames(Pi.specific.ranking[[i]])<-c("coeff","adj.r.sqr","avg.lvl")
  for (j in 1:length(Pi.specific.agi[[i]])){
    data<-average_agi_centered[Pi.specific.agi[[i]][j],]
    templm<-lm(unlist(data)~physio$Psi[1:16])
    Pi.specific.ranking[[i]][j,1]<-summary(templm)$coefficients[[2]]
    Pi.specific.ranking[[i]][j,2]<-summary(templm)$adj.r.squared
    Pi.specific.ranking[[i]][j,3]<-average_agi[Pi.specific.agi[[i]][j],1]
    rm(templm)
    rm(data)
  }
  Pi.specific.ranking[[i]]$score<-rank(abs(Pi.specific.ranking[[i]]$coeff))+rank(Pi.specific.ranking[[i]]$adj.r.sqr)+rank(Pi.specific.ranking[[i]]$avg.lvl) 
}

P.specific.ranking<-c()
for (i in 1:2){
  P.specific.ranking[[i]]<-data.frame(matrix(ncol = 3, nrow = length(P.specific.agi[[i]])))
  rownames(P.specific.ranking[[i]])<-P.specific.agi[[i]]
  colnames(P.specific.ranking[[i]])<-c("coeff","adj.r.sqr","avg.lvl")
  for (j in 1:length(matlab.P.specific.agi[[i]])){
    data<-average_agi_centered[P.specific.agi[[i]][j],]
    templm<-lm(unlist(data)~physio$P[1:16])
    P.specific.ranking[[i]][j,1]<-summary(templm)$coefficients[[2]]
    P.specific.ranking[[i]][j,2]<-summary(templm)$adj.r.squared
    P.specific.ranking[[i]][j,3]<-average_agi[P.specific.agi[[i]][j],1]
    rm(templm)
    rm(data)
  }
  P.specific.ranking[[i]]$score<-rank(abs(P.specific.ranking[[i]]$coeff))+rank(P.specific.ranking[[i]]$adj.r.sqr)+rank(P.specific.ranking[[i]]$avg.lvl) 
}
