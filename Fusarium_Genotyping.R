# fusarium clustering analysis
library(tidyverse)

path="inputs"

strain<- read.csv(paste0(path, "/strain_info2.csv"), na.strings = "")
genotypes<- read.csv(paste0(path, "/clusters_GACA4_M13.csv"))
suscep<- read.csv(paste0(path, "/antimicotic_suscep.csv"))
genotypes2<- merge(genotypes[genotypes$marker=="GACA4", ], genotypes[genotypes$marker=="M13", ], by="strain", all=T)
colnames(genotypes2)<- c("strain", "marker_GACA4", "genotype_GACA4","cluster_GACA4", "cluster_LR_GACA4",
                         "marker_M13", "genotype_M13", "cluster_M13", "cluster_LR_M13" )

genotypes2<- merge(genotypes2, strain, by="strain", all=T)
genotypes2<- merge(genotypes2, suscep, by="strain", all=T)

genotypes2$genotype_M13<- factor(genotypes2$genotype_M13, levels = str_sort(unique(as.character(genotypes2$genotype_M13)), numeric = T))
genotypes2$Voriconazol_phenotype<- ifelse(genotypes2$Voriconazol<=0.1, "S", ifelse(genotypes2$Voriconazol<=2, "SDD", "R"))
genotypes2$Amp_phenotype<- ifelse(genotypes2$Amfotericin_B<=0.1, "S", ifelse(genotypes2$Amfotericin_B<=4, "SDD", "R"))

pdf(paste0(path, "/figures_paper.pdf"), width = 18, height = 15, fonts = "ArialMT", pointsize = 14)

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_M13) & !is.na(genotypes2$TEF1_a),], aes(fill=TEF1_a, x=as.factor(cluster_LR_M13))) + 
  geom_bar(position="stack", stat="count") + ggtitle("TEF1-a distribution") +
  xlab("M13 cluster")+ theme_grey(base_size = 50)+ scale_fill_brewer(palette="Dark2")

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_M13) & !is.na(genotypes2$rDNA_28S),], aes(fill=rDNA_28S, x=as.factor(cluster_LR_M13))) + 
  geom_bar(position="stack", stat="count") + ggtitle("28S ID distribution") +
  xlab("M13 cluster")+ scale_fill_brewer(palette="Dark2")+ theme_grey(base_size = 50)

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_M13) & !is.na(genotypes2$TEF1_a),], aes(fill=Isolation_tissue,  x=as.factor(cluster_LR_M13))) + 
  geom_bar(position="stack", stat="count") + ggtitle("Isolation tissue TEF1-a ID") +
  xlab("M13 cluster")+ theme_grey(base_size = 50)+ scale_fill_brewer(palette="Dark2")

#GACA4
ggplot(genotypes2[!is.na(genotypes2$cluster_LR_GACA4) & !is.na(genotypes2$TEF1_a),], aes(fill=TEF1_a, x=as.factor(cluster_LR_GACA4))) + 
  geom_bar(position="stack", stat="count") + ggtitle("TEF1-a distribution") +
  xlab("GACA4 cluster")+ scale_fill_brewer(palette="Dark2")+ theme_grey(base_size = 50)

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_GACA4) & !is.na(genotypes2$rDNA_28S),], aes(fill=rDNA_28S, x=as.factor(cluster_LR_GACA4))) + 
  geom_bar(position="stack", stat="count") + ggtitle("28S ID distribution") +
  xlab("GACA4 cluster")+ scale_fill_brewer(palette="Dark2")+ theme_grey(base_size = 50)

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_GACA4) & !is.na(genotypes2$TEF1_a),], aes(fill=Isolation_tissue,  x=as.factor(cluster_LR_GACA4))) + 
  geom_bar(position="stack", stat="count") + ggtitle("Isolation tissue TEF1-a ID") +
  xlab("GACA4 cluster")+ theme_grey(base_size = 50)+ scale_fill_brewer(palette="Dark2")

#antimicotic resistance
colfunc<-colorRampPalette(c("blue","red"))
#plot(rep(1,50),col=(colfunc(50)), pch=19,cex=2)

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_M13) & !is.na(genotypes2$Voriconazol),], aes(fill=as.factor(Voriconazol), x=as.factor(cluster_LR_M13))) + 
  geom_bar(position="stack", stat="count") + ggtitle("Voriconazol response") +
  xlab("M13 cluster")+ theme_grey(base_size = 50)+scale_fill_manual(values=colfunc(12))

ggplot(genotypes2[!is.na(genotypes2$cluster_LR_M13) & !is.na(genotypes2$Amphotericin_B),], aes(fill=as.factor(Amphotericin_B), x=as.factor(cluster_LR_M13))) + 
  geom_bar(position="stack", stat="count") + ggtitle("Amphotericin B response") +
  xlab("M13 cluster")+ scale_fill_manual(values=colfunc(6))+ theme_grey(base_size = 50)

dev.off()

#DI calculations

calculate_HGDI<- function(genotypes){
  HGDI<- 1-(sum(mapply(function(x) sum(genotypes==x)*(sum(genotypes==x)-1), unique(genotypes))))/(length(genotypes)*(length(genotypes)-1))
  return(HGDI)
}

HGDI_M13_cluster<- calculate_HGDI(genotypes2$cluster_M13[!is.na(genotypes2$cluster_M13)])
HGDI_M13_cluster
HGDI_M13_cluster_LR<- calculate_HGDI(genotypes2$cluster_LR_M13[!is.na(genotypes2$cluster_LR_M13)])
HGDI_M13_cluster_LR
HGDI_M13_genotype<- calculate_HGDI(genotypes2$genotype_M13[!is.na(genotypes2$genotype_M13)])
HGDI_M13_genotype

HGDI_GACA4_cluster<- calculate_HGDI(genotypes2$cluster_GACA4[!is.na(genotypes2$cluster_GACA4)])
HGDI_GACA4_cluster
HGDI_GACA4_cluster_LR<- calculate_HGDI(genotypes2$cluster_LR_GACA4[!is.na(genotypes2$cluster_LR_GACA4)])
HGDI_GACA4_cluster_LR
HGDI_GACA4_genotype<- calculate_HGDI(genotypes2$genotype_GACA4[!is.na(genotypes2$genotype_GACA4)])
HGDI_GACA4_genotype


HGDI_PCR_molec<- calculate_HGDI(genotypes2$spp_molec[!is.na(genotypes2$spp_molec)])
HGDI_PCR_molec

HGDI_H_TFE1_a<- calculate_HGDI(genotypes2$haplotype_TEF1_a[!is.na(genotypes2$haplotype_TEF1_a)])
HGDI_H_TFE1_a
HGDI_H_28S<- calculate_HGDI(genotypes2$haplotype_28S_rDNA[!is.na(genotypes2$haplotype_28S_rDNA)])
HGDI_H_28S

HGDI_table<- data.frame(Method=c("M13_cluster_LR", "M13_cluster_HR", "M13_genotype", "(GACA)4_cluster_LR", "(GACA)4_cluster_HR","(GACA)4_genotype", "PCR_ID", "TEF1_a_Haplotypes", "rDNA_28S_Haplotypes"),
                        HGDI=c(HGDI_M13_cluster_LR, HGDI_M13_cluster, HGDI_M13_genotype, HGDI_GACA4_cluster_LR, HGDI_GACA4_cluster, HGDI_GACA4_genotype, HGDI_PCR_molec, HGDI_H_TFE1_a, HGDI_H_28S))

write.csv(HGDI_table, paste0(path, "/HGDI_table.csv"), row.names = F)

#HG test
HG_test<- function(trts_name, cluster_name, cluster_table){
  library(effsize)
  #table trts by clusters
  cluster_table_filtered<- cluster_table[!is.na(cluster_table[,trts_name]) & !is.na(cluster_table[,cluster_name]), c(cluster_name,trts_name)]
  trts<- unique(cluster_table_filtered[,trts_name])
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = length(unique(cluster_table_filtered[,cluster_name])))) 
  #hypergeom_test[,1]<- unique(cluster_table_filtered[,cluster_name])
  colnames(hypergeom_test)<- c(cluster_name, unlist(sapply(trts, function(x) c(paste0(x,"_HG_pvalue"), paste0(x,"_HG_padj")))))
  IDs_cluster<- data.frame(table(cluster_table_filtered[,cluster_name]))
  hypergeom_test[,1]<-IDs_cluster$Var1
  
  for(m in trts){
    IDs_trt<-data.frame(table(cluster_table_filtered[cluster_table_filtered[,trts_name]==m, cluster_name]))
    IDs_trt<- merge(IDs_trt, IDs_cluster, by="Var1", all.x=T)
    hypergeom_test[match(IDs_trt$Var1, hypergeom_test[,cluster_name]), 
                   paste0(m,"_HG_pvalue")]<- phyper(IDs_trt$Freq.x, sum(cluster_table_filtered[,trts_name]==m), 
                                                    sum(IDs_trt$Freq.y - IDs_trt$Freq.x), 
                                                    IDs_trt$Freq.y, lower.tail = FALSE)
    
    hypergeom_test[hypergeom_test[,cluster_name] %in% IDs_trt$Var1, paste0(m,"_HG_padj")]<- p.adjust(hypergeom_test[hypergeom_test[,cluster_name] %in% IDs_trt$Var1, paste0(m,"_HG_pvalue")], method = "BH")
    
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  return(hypergeom_test)
}

HG_Voriconazol<- HG_test(trts_name="Voriconazol", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_Voriconazol_phen<- HG_test(trts_name="Voriconazol_phenotype", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_Amfotericin_B<- HG_test(trts_name="Amphotericin_B", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_Amfotericin_B_phen<- HG_test(trts_name="Amp_phenotype", cluster_name="cluster_LR_M13", cluster_table=genotypes2)

HG_TEF1_a<- HG_test(trts_name="TEF1_a", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_28S<- HG_test(trts_name="rDNA_28S", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_tissue<- HG_test(trts_name="Isolation_tissue", cluster_name="cluster_LR_M13", cluster_table=genotypes2)
HG_gender<- HG_test(trts_name="Gender", cluster_name="cluster_LR_M13", cluster_table=genotypes2)

HG_TEF1_a_GACA4<- HG_test(trts_name="TEF1_a", cluster_name="cluster_GACA4", cluster_table=genotypes2)
HG_28S_GACA4<- HG_test(trts_name="rDNA_28S", cluster_name="cluster_GACA4", cluster_table=genotypes2)
HG_tissue_GACA4<- HG_test(trts_name="Isolation_tissue", cluster_name="cluster_GACA4", cluster_table=genotypes2)

write.csv(HG_TEF1_a_GACA4, paste0(path, "/HG_tests/HG_TEF1_a_GACA4.csv"))
write.csv(HG_28S_GACA4, paste0(path, "/HG_tests/HG_28S_GACA4.csv"))

write.csv(HG_Voriconazol, paste0(path, "/HG_tests/HG_Voriconazol.csv"))
write.csv(HG_Voriconazol_phen, paste0(path, "/HG_tests/HG_Voriconazol_phen.csv"))
write.csv(HG_Amfotericin_B, paste0(path, "/HG_tests/HG_Amfotericin_B.csv"))
write.csv(HG_TEF1_a, paste0(path, "/HG_tests/HG_TEF1_a.csv"))
write.csv(HG_28S, paste0(path, "/HG_tests/HG_28S.csv"))
write.csv(HG_tissue, paste0(path, "/HG_tests/HG_tissue.csv"))
write.csv(HG_gender, paste0(path, "/HG_tests/HG_gender.csv"))


#comparing matrices

#https://popgen.nescent.org/2015-12-15-microsatellite-differentiation.html
#http://userwww.sfsu.edu/efc/classes/biol710/amova/amova.htm
#https://stat.ethz.ch/pipermail/r-help/2005-October/080261.html


#AMOVA
library("adegenet")
library("pegas")
library("mmod")
library("reshape2")
library("ggplot2")
library("proxy")

genotypes2$haplotypes<- paste0(genotypes2$haplotype_28S_rDNA, "_",genotypes2$haplotype_TEF1_a)
M13_dist  <- read.csv(paste0(path, "/M13_matrix.csv"), row.names = 1)
M13_dist$X[duplicated(rownames(M13_dist))]
M13_stra  <- genotypes2[genotypes2$strain %in% rownames(M13_dist), ]
rownames(M13_dist)[!rownames(M13_dist)%in% genotypes2$strain]
colnames(M13_dist)<- rownames(M13_dist)
M13_dist<- as.matrix(M13_dist[M13_stra$strain, M13_stra$strain])
M13_stra$Isolation_tissue<- as.factor(M13_stra$Isolation_tissue)
M13_stra$spp_molec<- as.factor(M13_stra$spp_molec)
M13_stra$haplotypes<- as.factor(M13_stra$haplotypes)

M13_dist5<-M13_dist[M13_stra$strain[M13_stra$haplotypes!="NA_NA"],M13_stra$strain[M13_stra$haplotypes!="NA_NA"]]
#M13_dist4 = 1-M13_dist4
M13_dist5 = as.dist(M13_dist5)
M13_amova5 <- pegas::amova(M13_dist5 ~ haplotypes, data = M13_stra[M13_stra$haplotypes!="NA_NA", ], nperm = 100)
M13_amova5

# Analysis of Molecular Variance
# 
# Call: pegas::amova(formula = M13_dist5 ~ haplotypes, data = M13_stra[M13_stra$haplotypes != 
#                                                                        "NA_NA", ], nperm = 100)
# 
# SSD      MSD df
# haplotypes  7768914 323704.7 24
# Error       3396575 212285.9 16
# Total      11165489 279137.2 40
# 
# Variance components:
#   sigma2 P.value
# haplotypes  71100  0.0099
# Error      212286        
# 
# Phi-statistics:
#   haplotypes.in.GLOBAL 
# 0.2508945 

# Variance coefficients:
#   a 
# 1.567073 


#Mantel test
library("cultevo")
library("seqinr")
library("LncFinder")
library("Biostrings")
library("DECIPHER")
library("ade4")
path="inputs"

TEF1_a_stringset<-readDNAStringSet(paste0(path, "/fac_elon.fasta"))
TEF1_a_dist2<- DistanceMatrix(TEF1_a_stringset, type="matrix")
TEF1_a_dist3<- DistanceMatrix(TEF1_a_stringset, type="dist")              
is.euclid(TEF1_a_dist3)
is.euclid(small)
# cool the distances matrices are euclidian :)

# to compare the 2 matrices
M13_dist5<-as.dist(M13_dist[rownames(TEF1_a_dist2), rownames(TEF1_a_dist2)])
mt1 <- mantel.randtest(M13_dist5,TEF1_a_dist3,nrepet=10000)
mt1
plot(mt1)

