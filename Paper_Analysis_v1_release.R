library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(pheatmap)
library(reshape2)
set.seed(1)

saveext = "~/Amnioids_Figures/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#cType <- c("Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am")
#BaseCol <- c("lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF")
#cType <- c("Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","TB_CS5")
#BaseCol <- c("lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#921FE6")

#Now load in other datasets
#Key20307 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CRUK_SLX-20307/allQC20307.txt",sep="\t",header = T, row.names=1)
#raw_counts20307 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing data/CRUK_SLX-20307/featurecountsAll_extended_SLX-20307.csv",sep=",",header = T, row.names=1)
#Key20308 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CRUK_SLX-20308/allQC20308.txt",sep="\t",header = T, row.names=1)
#raw_counts20308 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing data/CRUK_SLX-20308/featurecountsAll_extended_SLX-20308.csv",sep=",",header = T, row.names=1)
#KeyCS5 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS5/CS5Key_210122.csv",sep=",",header = T, row.names=1)
#raw_countsCS5 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS5/featurecountsCS5.csv",sep=",",header = T, row.names=1)
#KeyCS6 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS6/CS6Key_210402.csv",sep=",",header = T, row.names=1)
#raw_countsCS6 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS6/featurecountsCS6.csv",sep=",",header = T, row.names=1)
#KeyCS7 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS7/CS7Key_210217.csv",sep=",",header = T, row.names=1) 
#raw_countsCS7 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS7/featurecountsCS7.csv",sep=",",header = T, row.names=1)
#KeyCRUK <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Annotations_Revision\ seq/files/CRUKKey.csv",sep=",",header = T, row.names=1) 
#raw_countsCRUK <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Annotations_Revision\ seq/files/featurecounts-CRUK.csv",sep=",",header = T, row.names=1)
#KeyCRUK2 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CRUK_210301_A00489_0804_AH3FMGDRXY/allQC.csv",sep=",",header = T, row.names=1) 
#raw_countsCRUK2 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing data/CRUK_210301_A00489_0804_AH3FMGDRXY/featurecounts.csv",sep=",",header = T, row.names=1)
#KeyCRUK3 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CRUK_210303_A00489_0808_BH3NYWDRXY/allQC.csv",sep=",",header = T, row.names=1) 
#raw_countsCRUK3 <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing data/CRUK_210303_A00489_0808_BH3NYWDRXY/featurecounts.csv",sep=",",header = T, row.names=1)
#KeyPre <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS1-3/CS1-3Key.csv",sep=",",header = T, row.names=1) 
#raw_countsPre <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/CS1-3/featurecountsCS1-3.csv",sep=",",header = T, row.names=1)

#marmoset_data_20307 <- CreateSeuratObject(counts = raw_counts20307[,which(Key20307$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_20307) <- Key20307$Primary.Annotation.Old[which(Key20307$QC>0)]
#marmoset_data_20307$Stage <- "CS6"
#marmoset_data_20307 <- subset(marmoset_data_20307, idents = c("Tb_CS6","EmDisc_Gast_CS6","ExMes_CS6","Am_CS6","SYS_CS6","PGC_CS6","VE_CS6","EmDisc_CS6","Am_CS6_EmDisc","EmDisc_Stalk_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_20307 <- NormalizeData(marmoset_data_20307, verbose = FALSE)
#marmoset_data_20307$Dataset <- "InVivo"
#marmoset_data_20308 <- CreateSeuratObject(counts = raw_counts20308[,which(Key20308$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_20308) <- Key20308$Primary.Annotation[which(Key20308$QC>0)]
#marmoset_data_20308$Stage <- "CS7"
#marmoset_data_20308 <- subset(marmoset_data_20308, idents = c("Tb_CS7","EmDisc_Gast_CS7","VE_CS7","EmDisc_CS7","Am_CS7","ExMes_stalk_CS7","PGC_CS7","SYS_CS7","ExMes_CS7","Am_CS7_EmDisc")) #,"ReStroma_CS7")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_20308 <- NormalizeData(marmoset_data_20308, verbose = FALSE)
#marmoset_data_20308$Dataset <- "InVivo"

#marmoset_data_Pre <- CreateSeuratObject(counts = raw_countsPre[,which(KeyPre$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_Pre) <- KeyPre$Primary.lineage[which(KeyPre$QC>0)]
#marmoset_data_Pre$LOC <- KeyPre$Location[which(KeyPre$QC>0)]
#marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Tb_CS3","Epi_CS3","ICM_CS3","Hyp_CS3","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3"))
#marmoset_data_Pre <- NormalizeData(marmoset_data_Pre, verbose = FALSE)
#marmoset_data_Pre$Stage <- "Pre"
#marmoset_data_Pre$Dataset <- "InVivo"

#marmoset_data_CRUK <- CreateSeuratObject(counts = raw_countsCRUK[,which(KeyCRUK$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CRUK) <- KeyCRUK$Primary.Annotation[which(KeyCRUK$QC>0)]
#marmoset_data_CRUK$LOC <- KeyCRUK$Loc[which(KeyCRUK$QC>0)]
#marmoset_data_CRUK$Stage <- KeyCRUK$Stage[which(KeyCRUK$QC>0)]
#marmoset_data_CRUK <- subset(marmoset_data_CRUK, idents = c("Tb_CS5","EmDisc_CS6","Am_CS7","VE_CS5","EmDisc_CS7","ExMes_CS6","ExMes_CS7","Am_CS6","SYS_CS7","EmDisc_CS5","EmDisc_CS7_Am","VE_CS6","ExMes_CS5","Am_CS5_PGC","Am_CS5","SYS_CS5")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_CRUK <- NormalizeData(marmoset_data_CRUK, verbose = FALSE)
#marmoset_data_CRUK$Dataset <- "InVivo"

#marmoset_data_CRUK2 <- CreateSeuratObject(counts = raw_countsCRUK2[,which(KeyCRUK2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CRUK2) <- KeyCRUK2$Primary.Annotation[which(KeyCRUK2$QC>0)]
#marmoset_data_CRUK2$LOC <- KeyCRUK2$Loc[which(KeyCRUK2$QC>0)]
#marmoset_data_CRUK2$Stage <- KeyCRUK2$Stage[which(KeyCRUK2$QC>0)]
#marmoset_data_CRUK2 <- subset(marmoset_data_CRUK2, idents = c("EmDisc_CS6","Am_CS6","EmDisc_CS5","Am_CS6_EmDisc","EmDisc_CS6_Gast","PGC_CS6","EmDisc_CS6_Am","ExMes_CS6","SYS_CS6","VE_CS6","Tb_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_CRUK2 <- NormalizeData(marmoset_data_CRUK2, verbose = FALSE)
#marmoset_data_CRUK2Dataset <- "InVivo"

#marmoset_data_CRUK3 <- CreateSeuratObject(counts = raw_countsCRUK3[,which(KeyCRUK3$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CRUK3) <- KeyCRUK3$Primary.Annotation[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3$LOC <- KeyCRUK3$Loc[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3$Stage <- KeyCRUK3$Stage[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3 <- subset(marmoset_data_CRUK3, idents = c("ExMes_CS6","Am_CS6","Am_CS7","EmDisc_CS7","EmDisc_CS6","EmDisc_CS7_Gast","PGC_CS7","VE_CS7","VE_CS6","SYS_CS7","SYS_CS6","Stroma_CS6")) #,"Stroma_CS6")) 
#marmoset_data_CRUK3 <- NormalizeData(marmoset_data_CRUK3, verbose = FALSE)
#marmoset_data_CRUK3$Dataset <- "InVivo"

#marmoset_data_CRUK3 <- CreateSeuratObject(counts = raw_countsCRUK3[,which(KeyCRUK3$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CRUK3) <- KeyCRUK3$Primary.Annotation[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3$LOC <- KeyCRUK3$Loc[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3$Stage <- KeyCRUK3$Stage[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3 <- subset(marmoset_data_CRUK3, idents = c("Am_CS7","EmDisc_CS7","EmDisc_CS7_Gast","PGC_CS7","VE_CS7","SYS_CS7")) #,"Stroma_CS6")) 
#marmoset_data_CRUK3 <- NormalizeData(marmoset_data_CRUK3, verbose = FALSE)
#marmoset_data_CRUK3$Dataset <- "InVivo"

#marmoset_data_CS5 <- CreateSeuratObject(counts = raw_countsCS5[,which(KeyCS5$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CS5) <- KeyCS5$Primary.Annotation[which(KeyCS5$QC>0)]
#marmoset_data_CS5$LOC <- KeyCS5$Location[which(KeyCS5$QC>0)]
#marmoset_data_CS5$Stage <- "CS5"
#marmoset_data_CS5 <- subset(marmoset_data_CS5, idents = c("Tb_CS5","ExMes_CS5","EmDisc_CS5","SYS_CS5","Am_CS5","VE_CS5","Am_CS5_PGC","EmDisc_CS5_Am")) #","Gland_CS5","ReGland_CS5","ReStroma_CS5")) 
#marmoset_data_CS5 <- NormalizeData(marmoset_data_CS5, verbose = FALSE)
#marmoset_data_CS5$Dataset <- "InVivo"

#marmoset_data_CS6 <- CreateSeuratObject(counts = raw_countsCS6[,which(KeyCS6$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CS6) <- KeyCS6$Primary.Anotation[which(KeyCS6$QC>0)]
#marmoset_data_CS6$LOC <- KeyCS6$Location[which(KeyCS6$QC>0)]
#marmoset_data_CS6$Stage <- "CS6"
#marmoset_data_CS6 <- subset(marmoset_data_CS6, idents = c("Tb_CS6","ExMes_CS6","EmDisc_CS6","SYS_CS6","Am_CS6","VE_CS6","PGC_CS6","EmDisc_gast_CS6","EmDisc_stalk_CS6")) #,"Gland_CS6","Stroma_CS6","ReStroma_CS6"))
#marmoset_data_CS6 <- NormalizeData(marmoset_data_CS6, verbose = FALSE)
#marmoset_data_CS6$Dataset <- "InVivo"

#marmoset_data_CS7 <- CreateSeuratObject(counts = raw_countsCS7[,which(KeyCS7$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_CS7) <- KeyCS7$Primary.Annotation[which(KeyCS7$QC>0)]
#marmoset_data_CS7$Stage <- "CS7"
#marmoset_data_CS7$LOC <- KeyCS7$Location[which(KeyCS7$QC>0)]
#marmoset_data_CS7 <- subset(marmoset_data_CS7, idents = c("Tb_CS7","ExMes_CS7","EmDisc_CS7","SYS_CS7","Am_CS7","Stalk_CS7_PGC","EmDisc_CS7_PGC","ExMes_stalk_CS7","EmDisc_stalk_CS7","Stalk_CS7")) #"Gland_CS7","ReStroma_CS7","ReGland_CS7","Myo_CS7")) 
#marmoset_data_CS7 <- NormalizeData(marmoset_data_CS7, verbose = FALSE)
#marmoset_data_CS7$Dataset <- "InVivo"

#marmoset_dataInVivo <- merge(marmoset_data_CS5, y = c(marmoset_data_Pre,marmoset_data_CS6,marmoset_data_CS7,marmoset_data_CRUK,marmoset_data_CRUK2,marmoset_data_CRUK3,marmoset_data_20307,marmoset_data_20308), project = "merged")
#marmoset_dataInVivo$Dataset <- "2) Marmoset in vivo"

#marmoset_dataInVivo2 <- marmoset_dataInVivo
#marmoset_dataInVivo2 <- FindVariableFeatures(marmoset_dataInVivo2, selection.method = "vst", nfeatures = 20000)
#marmoset_dataInVivo2 <- ScaleData(marmoset_dataInVivo2, verbose = FALSE)
#marmoset_dataInVivo2 <- RunPCA(marmoset_dataInVivo2, npcs = 20, verbose = FALSE)
#marmoset_dataInVivo2 <- RunUMAP(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
#marmoset_dataInVivo2 <- RunTSNE(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
#marmoset_dataInVivo2 <- FindNeighbors(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)

cType <- c("2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Stalk_CS6")
BaseCol <- c("lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#c49a00")

#colind <- integer( length( levels(Idents(marmoset_dataInVivo2)) )  )
#for (i in 1:length( levels(Idents(marmoset_dataInVivo2)) ) ) {
#  colind[i] <- which(cType==levels(Idents(marmoset_dataInVivo2))[i])
#}
#coluse <- BaseCol[colind]

#DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_All.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

#DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_All.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

#DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/TSNE_All.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)


#KeyInVitro <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/InVitro/InVitroKey_merged.csv",sep=",",header = T, row.names=1) 
#raw_countsInVitro <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/Chris_sequencing\ data/InVitro/featurecountsInVitro_merged.csv",sep=",",header = T, row.names=1)

#marmoset_data_InVitro <- CreateSeuratObject(counts = raw_countsInVitro[,which(KeyInVitro$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data_InVitro) <- KeyInVitro$Primary.lineage[which(KeyInVitro$QC>0)]
#marmoset_data_InVitro$LOC <- "Invitro"
#marmoset_data_InVitro <- subset(marmoset_data_InVitro, idents = c("EmD","EmDisc","Am","Amnoid_bead","CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"))
#marmoset_data_InVitro <- NormalizeData(marmoset_data_InVitro, verbose = FALSE)
#marmoset_data_InVitro$Stage <- "Invitro"
#marmoset_data_InVitro$Dataset <- "InVitro"

#marmoset_dataInVivo <- merge(marmoset_data_CS5, y = c(marmoset_data_Pre,marmoset_data_CS6,marmoset_data_CS7,marmoset_data_CRUK,marmoset_data_CRUK2,marmoset_data_CRUK3,marmoset_data_20307,marmoset_data_20308,marmoset_data_InVitro), project = "merged")
#marmoset_dataInVivo$Dataset <- "2) Marmoset in vivo"
#marmoset_dataInVivo2 <- marmoset_dataInVivo
#marmoset_dataInVivo2 <- FindVariableFeatures(marmoset_dataInVivo2, selection.method = "vst", nfeatures = 3000)
#marmoset_dataInVivo2 <- ScaleData(marmoset_dataInVivo2, verbose = FALSE)
#marmoset_dataInVivo2 <- RunPCA(marmoset_dataInVivo2, npcs = 20, verbose = FALSE)
#marmoset_dataInVivo2 <- RunUMAP(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
#marmoset_dataInVivo2 <- RunTSNE(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
#marmoset_dataInVivo2 <- FindNeighbors(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)

#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo, marmoset_data_InVitro), dims = 1:20, anchor.features = 5000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- merge(marmoset_data, y = c(marmoset_data2), project = "merged")
#mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 20000)

#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#Read in aligned datasets
mammal.combined <- readRDS("Data/Amnioids_aligned.rds")

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="Amnioid_#2")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Amnioid_#1")] <- 1
cluster_letters[which(Idents(mammal.combined)=="EmDisc_#2")] <- 2
cluster_letters[which(Idents(mammal.combined)=="EmDisc_#1")] <- 2
cluster_letters[which(Idents(mammal.combined)=="CHIR_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="SB43_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="ActA_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="BMP_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="FGF_noMEF")] <- 1
cluster_letters[which(Idents(mammal.combined)=="BMP_noMEF")] <- 1
cluster_letters[which(Idents(mammal.combined)=="ActA_noMEF")] <- 1

cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined, metadata = cluster_letters,col.name = 'cell.orig')

DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by = "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split4.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split4.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by = "cell.orig", reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge4.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge4.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

saveRDS(mammal.combined, file=paste(saveext,"/Amnioids_aligned.rds",sep=""))


mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split4b.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split4b.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge4b.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge4b.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

#Extract out EmDisc only.
library(e1071)
mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","Amnioid_#1","EmDisc_#2","EmDisc_#1"),invert=TRUE)
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="Am_CS5")) <- "Amnion"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="Am_CS6")) <- "Amnion"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="Am_CS7")) <- "Amnion"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="EmDisc_CS5")) <- "EmDisc"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="EmDisc_CS6")) <- "EmDisc"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="EmDisc_CS7")) <- "EmDisc"
Loadis <- as.data.frame(Loadings(mammal.combined, reduction = "pca")[, 1:2])
#Loadis_sub <- Loadis[ intersect( intersect(rownames(Loadis),TF), allMarkers$gene  ),]
#p <- ggplot(Loadis_sub, aes(x=PC_1, y=PC_2)) + geom_text(label=rownames(Loadis_sub),size=3)+ theme_classic() + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0))
#ggsave(filename=paste(saveext,"/DimRed/PCA_Loadings_1_2_TF",".pdf",sep=""),width = 40, height = 40,p)
dat <- as.data.frame(mammal.combined2 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined2)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "polynomial",degree=3, cost = 10, scale = FALSE)

make.grid = function(x, n = 75) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  expand.grid(X1 = x1, X2 = x2)
}
xgrid = make.grid(as.matrix(dat[,1:2]))
ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
plot(xgrid, col = as.numeric(ygrid), pch = 20, cex = .2)
points(as.matrix(dat[,1:2]), col = as.factor(dat$ID), pch = 19)
contour(seq(min(dat[,1]),max(dat[,1]),length.out=100), seq(min(dat[,2]),max(dat[,2]),length.out=100), matrix(func, 100, 100), level = 0.5, add = TRUE, col = "blue", lwd = 2)
contour(seq(0,1,length.out=100), seq(0,1,length.out=100), matrix(func, 100, 100), level = 0.5, add = TRUE, col = "blue", lwd = 2)

#Plot example expression
DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2,  reduction = "pca", features = "ISL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_ISL1.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_POU5F1.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_SOX2.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_TFAP2C.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_TFAP2A.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_VTCN1.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_T.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "EOMES", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_EOMES.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "MIXL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_expression_MIXL1.pdf",sep=""),width =16, height = 8, useDingbats = FALSE)


#Marker genes
FortyEight <- c("OOEP","TCL1A","WEE2",
                "IL1RN","NOV","ZNF80",
                "SPIC","ESRRB","STAT3",
                "KLF17","SOX15",
                "POU5F1","NANOG","SOX2","SFRP2","DNMT3B",
                "T","MIXL1","EOMES","LHX1","SNAI2","FOXA2",
                "PAX6","SOX1",
                "GATA3",
                "TFAP2C","TFAP2A","VTCN1","WNT6","ISL1",
                "SOX17","PRDM1","NANOS3","PRDM14",
                "GATA6","GATA4","PDGFRA","APOA1","TTR","APOB","HAND2","TBX4","HGF",
                "JAM2","GATA2","CGB3","CGA") 
DefaultAssay(mammal.combined2) <- "RNA"

Idents(mammal.combined2) <- factor(Idents(mammal.combined2), levels= c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7", "EmDisc_#1","Amnioid_#1","EmDisc_#2","Amnioid_#2")) 
for (i in 1:length( FortyEight ) ) {
  
  FeaturePlot(mammal.combined2,  reduction = "pca", features = FortyEight[i], combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/DimRed/PCA_expression_", FortyEight[i], ".pdf",sep=""),width =16, height = 8, useDingbats = FALSE)

  VlnPlot(mammal.combined2,  features = FortyEight[i], pt.size = 4)
  ggsave(filename=paste(saveext,"/DimRed/VIOLIN_expression_", FortyEight[i], ".pdf",sep=""),width =16, height = 8, useDingbats = FALSE)
  
}

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split5.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split5.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge5.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge5.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_SB43MEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_SB43MEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_SB43MEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_SB43MEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_BMPMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_BMPMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_BMPMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_BMPMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_CHIRMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_CHIRMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_CHIRMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_CHIRMEF.pdf",sep=""),width = 12, height = , useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_ActAMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_ActAMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_ActAMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_ActAMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_ActANoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_ActANoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_ActANoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_ActANoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_BMPNoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_BMPNoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_BMPNoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_BMPNoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","BMP_noMEF","ActA_noMEF","Amnioid_#2","EmDisc_#2"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_FGFnoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_FGFnoMEF.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_FGFnoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_FGFnoMEF.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

#Now do without preimp
DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2,  reduction = "umap", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_TFAP2A.pdf",sep=""),width =16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_TFAP2A.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_VTCN1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_VTCN1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_SOX2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_SOX2.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_POU5F1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_POU5F1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_T.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_T.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "umap", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_TFAP2C.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_TFAP2C.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_PDGFRA.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_PDGFRA.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "umap", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_SOX17.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_SOX17.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "umap", features = "NANOS3", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_NANOS3.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "NANOS3", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_NANOS3.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "umap", features = "NANOG", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_NANOG.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "NANOG", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_NANOG.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "umap", features = "FOXA2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_FOXA2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined2,  reduction = "pca", features = "FOXA2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_FOXA2.pdf",sep=""),width = 16, height = 8)

#Now redo volcano plots
AvExp <- AverageExpression(mammal.combined)
Ae <- AvExp$RNA
Ae$gene <- rownames(Ae)
SIGNAL<-read.table("Data/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("Data/TF.txt",header = F)
TF <- TF$V1

#Now do the volcano plots
DefaultAssay(mammal.combined) <- "RNA"
Idents(mammal.combined) <- mammal.combined$Cells

#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("EmDisc_CS5"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl1 <- FindMarkers(mammal.combined, ident.1 = c("EmDisc_CS5"), ident.2 = c("Epi_CS3"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
#write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDiscCS5_Epi_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Epi_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
ggsave(filename=paste(saveext,"EmDiscCS5_Epi_TF.eps",sep=""),width = 13, height = 13, plot = p1, useDingbats = FALSE)

dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Epi_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("Am_CS5"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("Am_CS5"), ident.2 = c("Epi_CS3"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/AmCS5_Epi_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS5_Epi_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS5_Epi_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5"), ident.1 = c("Am_CS5"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5"), ident.1 = c("Am_CS5"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
#Cl2 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5"), ident.1 = c("Am_CS5"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))


write.csv(as.data.frame(Cl1), file=paste(saveext,"/AmCS5_EmDiscCS5_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = rownames(Ae5)[pospos1]
genes.to.label1 = intersect(rownames(Ae5)[pospos1],TF)
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS5_EmDiscCS5_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))  + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS5_EmDiscCS5_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("Am_CS6"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("Am_CS6"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#1"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
Cl3 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#2"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))

write.csv(as.data.frame(Cl1), file=paste(saveext,"/AmCS6_EmDiscCS6_DEseq2.csv",sep=""))
write.csv(as.data.frame(Cl2), file=paste(saveext,"/Am_EmD_DEseq2.csv",sep=""))
write.csv(as.data.frame(Cl3), file=paste(saveext,"/Amnoid_EmDisc_DEseq2.csv",sep=""))


overlap1 <- intersect(rownames(Cl1)[which(Cl1$p_val_adj<0.05)],rownames(Cl2)[which(Cl2$p_val_adj<0.05)])

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.5) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label0 <- overlap1
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label1,genes.to.label2,genes.to.label3,genes.to.label4,as.character(TF) ))
genes.to.label5 <- intersect(genes.to.label5,overlap1)
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_EmDiscCS6_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_EmDiscCS6_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = overlap1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_EmDiscCS6_overlap.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl2),"Pval1"] <- -log(Cl2$p_val_adj)
Ae5[rownames(Cl2),"FC1"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.5) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label0 <- overlap1
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label1,genes.to.label2,genes.to.label3,genes.to.label4,as.character(TF) ))
genes.to.label5 <- intersect(genes.to.label5,overlap1)
Ae5 <- Ae5[rownames(Cl2),]
Ae5 <- Ae5[order(Ae5$Pval1),]

Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EmD_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EmD_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = overlap1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EmD_overlap.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#2"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#2"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))

write.csv(as.data.frame(Cl1), file=paste(saveext,"/Amnoidbead_EmDisc_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnoidbead_EmDisc_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnoidbead_EmDisc_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#1"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_#1"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Am_EmD_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EmD_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EmD_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Amnoidbead_EpiCS3_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioidbead_EpiCS3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioidbead_EpiCS3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Am_EpiCS3_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EpiCS3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_EpiCS3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("EmDisc_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("EmDisc_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDisc_EpiCS3_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_EpiCS3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_EpiCS3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("EmDisc_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Epi_CS3"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmD_EpiCS3_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmD_EpiCS3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmD_EpiCS3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("EmDisc_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("EmDisc_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDisc_EmDiscCS6_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_EmDiscCS6_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_EmDiscCS6_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_CS6"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_CS6"), ident.1 = c("Amnioid_#2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Amnoidbead_AmCS6_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnoidbead_AmCS6_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnoidbead_AmCS6_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmD_EmDiscCS6_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmD_EmDiscCS6_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmD_EmDiscCS6_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do the volcano plots
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_CS6"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_CS6"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Am_AmCS6_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_AmCS6_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_AmCS6_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Now cluster broadly into types
#mammal.combined$Cells <- Idents(mammal.combined)
AltID <- as.character(Idents(mammal.combined))
AltID[which(AltID=="EmDisc_CS5")] <- "EmDisc"
AltID[which(AltID=="EmDisc_CS6")] <- "EmDisc"
AltID[which(AltID=="EmDisc_CS7")] <- "EmDisc"
AltID[which(AltID=="Am_CS5")] <- "Am"
AltID[which(AltID=="Am_CS6")] <- "Am"
AltID[which(AltID=="Am_CS7")] <- "Am"
mammal.combined$AltID <- AltID

tbl <- table(Idents(mammal.combined))
tbl[order(-as.numeric(names(tbl)))]

#Now subcluster ...
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined$IDs <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.3)
mammal.combined$Cl3 <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.9)
mammal.combined$Cl9 <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$Cl5
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_Invit_Cl5.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 8, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_All_Invit_Cl5.pdf",sep=""),width = 20, height = 16,limitsize = FALSE, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_All_Invit_Cl5.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

Idents(mammal.combined) <- mammal.combined$Cl3
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_Invit_Cl3.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 8, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_All_Invit_Cl3.pdf",sep=""),width = 20, height = 16,limitsize = FALSE, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_All_Invit_Cl3.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

longID <- paste(mammal.combined$AltID,mammal.combined$Cl3,sep="_")
Idents(mammal.combined) <- longID
DimPlot(mammal.combined, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_Invit_Cl3_ID.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)
DimPlot(mammal.combined,  pt.size = 8, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_All_Invit_Cl3_ID.pdf",sep=""),width = 20, height = 16,limitsize = FALSE, useDingbats = FALSE)
DimPlot(mammal.combined,  pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_All_Invit_Cl3_ID.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined,idents = c("EmDisc_1","EmDisc_3","Am_0","Am_2"))
VlnPlot(mammal.combined2,features=c("SOX2","POU5F1","T","MIXL1","TFAP2A","VTCN1"))
ggsave(filename=paste(saveext,"/DimRed/Violin_Cl3_ID.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

DefaultAssay(mammal.combined) <- "RNA"
AvExp <- AverageExpression(mammal.combined)

list1 <- c("EmDisc_1", 
           "EmDisc_0",
           "EmDisc_3",  
           "Am_0",  
           "Am_2")

list2 <- c("EmDisc_#1_1",   
           "EmDisc_#2_1", 
           "Amnioid_#1_0",
           "Amnioid_#2_0")

C1 <- cor( log(as.matrix(AvExp$RNA[,list1]) +1), log(as.matrix(AvExp$RNA[,list2]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), gaps_row = c(3), gaps_col = c(2),display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterRNACross",".pdf",sep=""),width=3,height=3.5)

C1 <- cor( log(as.matrix(AvExp$integrated[,list1]) +1), log(as.matrix(AvExp$integrated[,list2]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), gaps_row = c(3), gaps_col = c(2), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterIntCross",".pdf",sep=""),width=3,height=3.5)

list3 <- c("EmDisc_#1_1",   
           "EmDisc_#2_1", 
           "ActA_MEF_1",
           "BMP_MEF_0",
           "SB43_MEF_0",
           "SB43_MEF_2",
           "CHIR_MEF_2",
           "Amnioid_#1_0",
           "Amnioid_#2_0",
           "BMP_noMEF_0",
           "ActA_noMEF_1",
           "ActA_noMEF_3",
           "FGF_noMEF_3")


C1 <- cor( log(as.matrix(AvExp$RNA[,list1]) +1), log(as.matrix(AvExp$RNA[,list3]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), gaps_row = c(3), gaps_col = c(3,5,7,10),display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterRNACross_screen",".pdf",sep=""),width=5,height=3.5)

C1 <- cor( log(as.matrix(AvExp$integrated[,list1]) +1), log(as.matrix(AvExp$integrated[,list3]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), gaps_row = c(3), gaps_col = c(3,5,7,10), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterIntCross_screen",".pdf",sep=""),width=5,height=3.5)


Idents(mammal.combined) <- mammal.combined$Cells
mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","BMP_noMEF","ActA_MEF","FGF_noMEF","ActA_noMEF"),invert=TRUE)
Idents(mammal.combined2) <- mammal.combined2$Cl3
DimPlot(mammal.combined2,  pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_Cl.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2,  pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_Cl.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2,  pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_Cl.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2,  pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_Cl.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


#Now do the volcano plots
mammal.combined$longID <- longID
DefaultAssay(mammal.combined) <- "RNA"
Idents(mammal.combined) <- mammal.combined$longID
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_0"), ident.1 = c("EmDisc_#1_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_1"), ident.1 = c("EmDisc_#1_1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDisc_1_Em_ht1_1_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_0_EmDisc_ht1_0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_0_EmDisc_ht1_0_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

DefaultAssay(mammal.combined) <- "RNA"
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_0"), ident.1 = c("EmDisc_#2_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_1"), ident.1 = c("EmDisc_#2_!"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDisc_1_Em_ht2_1_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_0_EmDisc_ht2_0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_0_EmDisc_ht2_0_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

DefaultAssay(mammal.combined) <- "RNA"
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_1"), ident.1 = c("Amnioid_#1_1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_0"), ident.1 = c("Amnioid_#1_0"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Am_0_Am_ht1_0_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_1_Am_ht1_1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_1_Am_ht1_1_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

DefaultAssay(mammal.combined) <- "RNA"
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_5"), ident.1 = c("Amnioid_#2_5"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Am_0"), ident.1 = c("Amnioid_#2_0"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Am_0_Am_ht2_0_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_1_Am_ht1_1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Am_1_Am_ht1_1_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
#
#
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("BMP_MEF"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("BMP_MEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/BMP4_EmDisc_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
#p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_BMP_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_BMP_MEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("CHIR_MEF"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("CHIR_MEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/CHIR_EmDisc_DEseq2.csv",sep=""))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.5) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
Ae5 <- Ae5[which(is.na(Ae5$Pval1)==FALSE),]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_CHIR_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_CHIR_MEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("SB43_MEF"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("SB43_MEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/SB43_EmDisc_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_SB43_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_SB43_MEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("ActA_MEF"), ident.1 = c("EmDisc_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("ActA_MEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/ActAMEF_EmDisc_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_ActAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_ActAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("ActA_noMEF"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("ActA_noMEF"), ident.2 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/ActA_Am_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_ActnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_ActnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("BMP_noMEF"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("BMP_noMEF"), ident.2 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/BMP_Am_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_BMPnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_BMPnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now do comparisons with various perturbations?
Idents(mammal.combined) <- mammal.combined$Cells
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("FGF_noMEF"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(mammal.combined, ident.1 = c("FGF_noMEF"), ident.2 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"FGF_Am_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_FGFnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_FGFnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.1 = c("FGF_noMEF"), ident.2 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"FGF_Am_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_FGFnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Amnioid_1_FGFnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()l1 <- FindMarkers(mammal.combined, ident.1 = c("FGF_noMEF"), ident.2 = c("Amnioid_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"FGF_Am_DEseq2.csv",sep=""))

Cl1 <- FindMarkers(mammal.combined, ident.1 = c("FGF_noMEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"FGF_EmDisc_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_FGFnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_FGFnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(mammal.combined, ident.1 = c("ActA_noMEF"), ident.2 = c("EmDisc_#1"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/ActA_EmDisc_DEseq2.csv",sep=""))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.1)) ,]
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_ActnoAMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.01),linetype="dashed") + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDisc_1_ActnoAMEF_lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

library("ggtern")
Idents(mammal.combined) <- mammal.combined$Cells
#Now whole compare labelled datasets


#Cl0 <- FindMarkers(mammal.combined, ident.2 = c("Amnioid_#2"), ident.1 = c("Amnioid_#1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc"), ident.1 = c("EmDisc_CS6"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#Cl2 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5"), ident.1 = c("Am_CS5"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#Cl3 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("Am_CS6"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#Cl4 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS6"), ident.1 = c("Epi_CS3"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#Cl5 <- FindMarkers(mammal.combined, ident.2 = c("Am_CS6"), ident.1 = c("Epi_CS3"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)

#genes.to.label <- intersect(c(rownames(Cl0),rownames(Cl1),rownames(Cl2),rownames(Cl3),rownames(Cl4),rownames(Cl5)),TF)

#TernX$colsi <- 0
#TernX[genes.to.label,c("colsi")] <- 1
#TernX$colsi <- as.factor(TernX$colsi)

#TernX3 <- TernX[order(TernX$colsi),]
#TernX2 <- TernX3[which(TernX3$colsi==1),]

#p1 <- ggtern(data=TernX3, aes_string(x="EmDisc_CS6",z="Am_CS6", y="Amnoid_bead")) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("lightgrey", "black")) #+ scale_color_manual(values=c("grey", "black")) 
#p1 <- p1 + annotate(geom  = 'text', x = TernX2$EmDisc_CS6/(TernX2$EmDisc_CS6+TernX2$Amnoid_bead+TernX2$Am_CS6),z= TernX2$Am_CS6/(TernX2$EmDisc_CS6+TernX2$Amnoid_bead+TernX2$Am_CS6),y= TernX2$Amnoid_bead/(TernX2$EmDisc_CS6+TernX2$Amnoid_bead+TernX2$Am_CS6),label = rownames(TernX2),color = c("red")) 
#ggsave(filename=paste(saveext,"GGTern1.pdf",sep=""),width = 40, height = 40, plot = p1)

#p1 <- ggtern(data=TernX3, aes_string(x="EmDisc_CS6",z="Am_CS6", y="EmDisc")) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("lightgrey", "black")) #+ scale_color_manual(values=c("grey", "black")) 
#p1 <- p1 + annotate(geom  = 'text', x = TernX2$EmDisc_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),z= TernX2$Am_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),y= TernX2$EmDisc/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),label = rownames(TernX2),color = c("red")) 
#ggsave(filename=paste(saveext,"GGTern2.pdf",sep=""),width = 40, height = 40, plot = p1)

#p1 <- ggtern(data=TernX3, aes_string(x="EmDisc_CS6",z="Am_CS6", y="EmDisc")) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("lightgrey", "black")) #+ scale_color_manual(values=c("grey", "black")) 
#p1 <- p1 + annotate(geom  = 'text', x = TernX2$EmDisc_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),z= TernX2$Am_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),y= TernX2$EmDisc/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),label = rownames(TernX2),color = c("red")) 
#ggsave(filename=paste(saveext,"GGTern2.pdf",sep=""),width = 40, height = 40, plot = p1)

#write.csv(as.data.frame(Cl1), file=paste(saveext,"/EmDisc_AmnoidBead.csv",sep=""))
#write.csv(as.data.frame(Cl2), file=paste(saveext,"/EmDiscCS5_AmCS5.csv",sep=""))
#write.csv(as.data.frame(Cl3), file=paste(saveext,"/EmDiscCS6_AmCS6.csv",sep=""))
#write.csv(as.data.frame(Cl4), file=paste(saveext,"/EmDisc_EmDiscCS6.csv",sep=""))
#write.csv(as.data.frame(Cl5), file=paste(saveext,"/Am_AmCS6.csv",sep=""))

uID <- as.character(Idents(mammal.combined))
uID[which(uID=="EmDisc_#1")] <- "EmDisc"
uID[which(uID=="Amnioid_#1")] <- "Amnoid_bead"
Idents(mammal.combined) <- uID

library("ggtern")
AvExp <- AverageExpression(mammal.combined)
TernX <- AvExp$RNA
#TernX2$IV <- max()
genes.to.label <- c("VTCN1","POU5F1","SOX2","WNT6A","TFAP2A","TFAP2C","NANOG","PDGFRA","MESP2","T","EOMES","SOX17","MIXL1","GATA2")
TernXl <- TernX[genes.to.label,]
p1 <- ggtern(data=TernXl,
               aes(x = EmDisc_CS6,
                   z= Am_CS6,
                   y= EmDisc,
                   xend = EmDisc_CS6,
                   zend= Am_CS6,
                   yend= BMP_MEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),y= TernXl$EmDisc/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_BMP_MEF.pdf",sep=""),width = 20, height = 20, plot = p1)


p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= EmDisc,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= ActA_MEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),y= TernXl$EmDisc/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_ActA_MEF.pdf",sep=""),width = 20, height = 20, plot = p1)


p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= EmDisc,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= SB43_MEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),y= TernXl$EmDisc/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_SB43_MEF.pdf",sep=""),width = 20, height = 20, plot = p1)


p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= Amnoid_bead,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= ActA_noMEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),y= TernXl$Amnoid_bead/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_ActA_noMEF.pdf",sep=""),width = 20, height = 20, plot = p1)


p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= Amnoid_bead,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= BMP_noMEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),y= TernXl$Amnoid_bead/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_BMP_noMEF.pdf",sep=""),width = 20, height = 20, plot = p1)


p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= Amnoid_bead,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= FGF_noMEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),y= TernXl$Amnoid_bead/(TernXl$EmDisc_CS6+TernXl$Amnoid_bead+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_FGF_noMEF.pdf",sep=""),width = 20, height = 20, plot = p1)

#genes.to.label <- c("VTCN1","POU5F1","SOX2","WNT6A","TFAP2A","TFAP2C","NANOG","PDGFRA","MESP2","T","EOMES")
TernXl <- TernX[genes.to.label,]
p1 <- ggtern(data=TernXl,
             aes(x = EmDisc_CS6,
                 z= Am_CS6,
                 y= EmDisc,
                 xend = EmDisc_CS6,
                 zend= Am_CS6,
                 yend= CHIR_MEF)) + geom_point() + geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1 + annotate(geom  = 'text', x = TernXl$EmDisc_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),z= TernXl$Am_CS6/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),y= TernXl$EmDisc/(TernXl$EmDisc_CS6+TernXl$EmDisc+TernXl$Am_CS6),label = rownames(TernXl),color = c("red")) 
ggsave(filename=paste(saveext,"GGTern_CHIR_MEF.pdf",sep=""),width = 20, height = 20, plot = p1)

Idents(mammal.combined) <- mammal.combined$Cells

exps <- c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7")
expe2 <- c("EmDisc_#2","EmDisc_#1","ActA_MEF","ActA_noMEF","FGF_noMEF","Amnioid_#1","Amnioid_#2","BMP_noMEF","BMP_MEF","SB43_MEF","CHIR_MEF")

#Finally generate correlations for identity mapping
subs1 <- which(Idents(mammal.combined) %in% exps)
subs2 <- which(Idents(mammal.combined) %in% expe2)
Mat1 <- GetAssayData(mammal.combined,assay = "RNA")
C1 <- cor(as.matrix(Mat1[,subs1]),as.matrix(Mat1[,subs2]), method = "pearson")
write.csv(as.data.frame(C1), file=paste(saveext,"/Correlation_alltissues_pearson_RNA.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[subs1]), file=paste(saveext,"/Idents_alltissues1_RNA.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[subs2]), file=paste(saveext,"/Idents_alltissues2_RNA.csv",sep=""))
Mat1 <- GetAssayData(mammal.combined,assay = "integrated")
C1 <- cor(as.matrix(Mat1[,subs1]),as.matrix(Mat1[,subs2]), method = "pearson")
write.csv(as.data.frame(C1), file=paste(saveext,"/Correlation_alltissues_pearson_Int.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[subs1]), file=paste(saveext,"/Idents_alltissues1_Int.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[subs2]), file=paste(saveext,"/Idents_alltissues2_Int.csv",sep=""))

#end game

raw_countsMEF <- read.table("/Volumes/GoogleDrive/Shared\ drives/Munger\ et\ al./Chris\'\ Sequencing\ Dungeon/MouseMEFs.txt",sep="\t",header = T, row.names=1)
KeyMEF <- read.table("/Volumes/GoogleDrive/Shared\ drives/Munger\ et\ al./Chris\'\ Sequencing\ Dungeon/MEFIDs.csv",sep=",",header = T, row.names=1)
rownames(raw_countsMEF) <- toupper(rownames(raw_countsMEF))

marmoset_data_MEF <- CreateSeuratObject(counts = raw_countsMEF[,which(KeyMEF$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_MEF) <- KeyMEF$Annotation[which(KeyMEF$QC>0)]
marmoset_data_MEF <- subset(marmoset_data_MEF, idents = c("ESC_1"),invert = TRUE) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_MEF$Stage <- "MEF"
marmoset_data_MEF <- NormalizeData(marmoset_data_MEF, verbose = FALSE)
marmoset_data_MEF$Dataset <- "MEF"

marmoset_data_MEF <- FindVariableFeatures(marmoset_data_MEF, selection.method = "vst", nfeatures = 1000)
marmoset_data_MEF <- ScaleData(marmoset_data_MEF, verbose = FALSE)
marmoset_data_MEF <- RunPCA(marmoset_data_MEF, npcs = 20, verbose = FALSE)
marmoset_data_MEF <- RunUMAP(marmoset_data_MEF, reduction = "pca", dims = 1:20)
marmoset_data_MEF <- RunTSNE(marmoset_data_MEF, reduction = "pca", dims = 1:20)
marmoset_data_MEF <- FindNeighbors(marmoset_data_MEF, reduction = "pca", dims = 1:20)

AVE <- AverageExpression(marmoset_data_MEF)

write.csv(as.data.frame(AVE$RNA), file=paste(saveext,"/MEF_average_couns_per_10k.csv",sep=""))


C1 <- cor( as.matrix( log(AVE$RNA[,]+1) ), method = "pearson")
mat_breaks <- seq(0.65, 0.85, length.out = 20)
pheatmap(C1,color =  redblue1(60), breaks=mat_breaks,  display_numbers = round(C1,2), fontsize_number = 5, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/MEFCorr",".pdf",sep="") ,width=3.4,height=3.2)


DimPlot(marmoset_data_MEF, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_MEF.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)
DimPlot(marmoset_data_MEF, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_MEF.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)
DimPlot(marmoset_data_MEF,  pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_MEF.pdf",sep=""),width = 20, height = 16, useDingbats = FALSE)

DEMEF <- FindMarkers(marmoset_data_MEF,ident.1 = "iMEF", ident.2 = "ESC_2",verbose = FALSE, test.use = "DESeq2", only.pos =TRUE)
write.csv(as.data.frame(DEMEF), file=paste(saveext,"/MEF_up_reg.csv",sep=""))

#FindClusters()
marmoset_data_InVitro <- CreateSeuratObject(counts = raw_countsInVitro[,which(KeyInVitro$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_InVitro) <- KeyInVitro$Primary.lineage[which(KeyInVitro$QC>0)]
marmoset_data_InVitro$LOC <- "Invitro"
marmoset_data_InVitro <- subset(marmoset_data_InVitro, idents = c("ESC_conv2"))
marmoset_data_InVitro <- NormalizeData(marmoset_data_InVitro, verbose = FALSE)
marmoset_data_InVitro$Stage <- "Invitro"
marmoset_data_InVitro$Dataset <- "InVitro"

mammal.combined <- merge(marmoset_data_MEF, y = c(marmoset_data_InVitro), project = "merged")
mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 3000)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

DimPlot(mammal.combined, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/MEF_UMAP_marm_split.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)

DimPlot(mammal.combined, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/MEF_PCA_marm_split.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)

DefaultAssay(mammal.combined) <- "RNA"
AVE_2 <- AverageExpression(mammal.combined)
write.csv(as.data.frame(AVE2$RNA), file=paste(saveext,"/MEF_with_marm_average_couns_per_10k.csv",sep=""))


#marmoset_data_MEF_sub <- subset(marmoset_data_MEF,idents = c("iMEF","ESC_2"))
#AVE <- AverageExpression(marmoset_data_MEF_sub)

DEMEF2 <- FindMarkers(mammal.combined,ident.1 = "iMEF", ident.2 = "ESC_conv2",verbose = FALSE, test.use = "DESeq2", only.pos =TRUE)
write.csv(as.data.frame(DEMEF2), file=paste(saveext,"/MEF_up_reg_vsMarm.csv",sep=""))

TF<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Leaving_package/Dimensionality\ reduction\ techniques\ smart-seq2\ -\ Boroviaklab\ data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

LIST1 <- intersect(toupper(rownames(DEMEF)),SIGNAL1)
LIST2 <- intersect(toupper(rownames(DEMEF)),SIGNAL2)
LIST3 <- intersect(toupper(rownames(DEMEF)),SIGNAL3)

AVE1 <- AVE$RNA[LIST1,]
AVE1$X <- as.factor(rownames(AVE1))
AVE2 <- AVE$RNA[LIST2,]
AVE2$X <- as.factor(rownames(AVE2))
AVE3 <- AVE$RNA[LIST3,]
AVE3$X <- as.factor(rownames(AVE3))

AVE4 <- AVE$RNA[sort(unique(c(LIST2,LIST3)) ),]
AVE4$X <- as.factor(rownames(AVE4))

D1 <- melt(AVE1)
D2 <- melt(AVE2)
D3 <- melt(AVE3)
D4 <- melt(AVE4)

p1 <- ggplot(D2, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/Ligand.pdf",sep=""),width = 30, height = 16, p1)

p1 <- ggplot(D3, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/ECM.pdf",sep=""),width = 30, height = 16, p1)

p1 <- ggplot(D4, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/LIGAND_ECM.pdf",sep=""),width = 30, height = 16, p1)


LIST1 <- intersect(toupper(rownames(DEMEF2)),SIGNAL1)
LIST2 <- intersect(toupper(rownames(DEMEF2)),SIGNAL2)
LIST3 <- intersect(toupper(rownames(DEMEF2)),SIGNAL3)

AVE1 <- AVE_2$RNA[LIST1,]
AVE1$X <- as.factor(rownames(AVE1))
AVE2 <- AVE_2$RNA[LIST2,]
AVE2$X <- as.factor(rownames(AVE2))
AVE3 <- AVE_2$RNA[LIST3,]
AVE3$X <- as.factor(rownames(AVE3))

AVE4 <- AVE_2$RNA[sort(unique(c(LIST2,LIST3)) ),]
AVE4$X <- as.factor(rownames(AVE4))

D1 <- melt(AVE1)
D2 <- melt(AVE2)
D3 <- melt(AVE3)
D4 <- melt(AVE4)

p1 <- ggplot(D2, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/Ligand_marm.pdf",sep=""),width = 30, height = 16, p1)

p1 <- ggplot(D3, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/ECM_marm.pdf",sep=""),width = 30, height = 16, p1)

p1 <- ggplot(D4, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + theme_bw() + theme(axis.text=element_text()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(filename=paste(saveext,"/DimRed/LIGAND_ECM_marm.pdf",sep=""),width = 30, height = 16, p1)



sdsakdsalkdsalmdlsamd

#CHIRM <- FindMarkers(mammal.combined,ident.2 = c("EmDisc"),ident.1 = "CHIR_MEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#write.csv(as.data.frame(CHIRM), file=paste(saveext,"/Chiron_vs_EmDisc.csv",sep=""))

#BMPM <- FindMarkers(mammal.combined,ident.2 = c("EmDisc"),ident.1 = "BMP_MEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#write.csv(as.data.frame(BMPM), file=paste(saveext,"/BMPMEF_vs_EmDisc.csv",sep=""))

#SBM <- FindMarkers(mammal.combined,ident.2 = c("EmDisc"),ident.1 = "SB43_MEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#write.csv(as.data.frame(SBM), file=paste(saveext,"/SB43MEF_vs_EmDisc.csv",sep=""))

#ActM <- FindMarkers(mammal.combined,ident.2 = c("EmDisc"),ident.1 = "ActA_MEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE)
#write.csv(as.data.frame(ActM), file=paste(saveext,"/ActAMEF_vs_EmDisc.csv",sep=""))

##p1 <- ggtern(data=TernX2, aes_string(x="EmDisc_CS6",z="Am_CS6", y="EmDisc")) + geom_point(aes(color=colsi)) +  geom_segment(aes(xend = 0.1+TernX2$EmDisc_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),zend= TernX2$Am_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),yend= TernX2$EmDisc/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6)), arrow = arrow(length = unit(0.1, "cm"))) + scale_color_manual(values=c("lightgrey", "black")) #+ scale_color_manual(values=c("grey", "black")) 
#p1 <- p1 + annotate(geom  = 'text', x = TernX2$EmDisc_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),z= TernX2$Am_CS6/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),y= TernX2$EmDisc/(TernX2$EmDisc_CS6+TernX2$EmDisc+TernX2$Am_CS6),label = rownames(TernX2),color = c("red")) 

#geom_point(size = 1) + 
#  geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.1, "cm"))) + 
  

ggsave(filename=paste(saveext,"GGTern2.pdf",sep=""),width = 40, height = 40, plot = p1)



+geom_segment(arrow = arrow())


#Mat1 <- GetAssayData(mammal.combined,assay = "RNA")
#C1 <- cor(as.matrix(Mat1[,subs1]),as.matrix(Mat1[,subs2]), method = "pearson")
#mat_breaks <- seq(0.65, 0.85, length.out = 20)
#pheatmap(C1,color =  redblue1(20),gaps_col=c(5), gaps_row=c(3,6,9,12,15), breaks=mat_breaks,border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Alltissueheamtap",".pdf",sep="") ,width=4,height=4)
#mycolors <- c("#0c9cf5",
#              "#0767DA",
#              "#0233BF",
#              "#c2a5cf","#7b3294","#4d4d4d",
#              "#0571b0","#92c5de",
#              "#877bd6",
#              "#5F54C7",
#              "#1A0873",
#              "#008837","#a6dba0",
#              "#ca0020","#f4a582",
#              "#d8b365","#5ab4ac")
#avexp  <- AverageExpression(object = mammal.combined, slot = "data")
#a <- avexp$RNA
##annotation_col = data.frame(Stage = factor(colnames(X)))
#rownames(annotation_col) <- colnames(X)
#names(mycolors) <- colnames(X)
#anno_colors <- list(Stage = mycolors)
#redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#mat_breaks <- seq(0, 20, length.out = 20)
#pheatmap(X*100,color =  redblue1(20),gaps_col=c(3,8,11), gaps_row=c(5,9,11,13),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/heamtap",".pdf",sep="") ,width=5,height=4)
#mat_breaks <- seq(-2, 2, length.out = 20)
#pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(3,8,11),gaps_row=c(5,9,11,13),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/heamtapscale",".pdf",sep=""),width=5,height=4 )

#exps <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","ICM_CS3","Epi_CS3","Hypo_CS3","Tb_CS3",
#  "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
#  "PGC_CS6","PGC_CS7",
#          "Tb_CS5","Tb_CS6","Tb_CS7",
#          "ExMes_CS5","ExMes_CS6","ExMes_CS7",
#          "VE_CS5","VE_CS6","VE_CS7",
#          "SYS_CS5","SYS_CS6","SYS_CS7",
#  "EmDisc","EmD","ActA_MEF","ActA_noMEF","FGF_noMEF","Am","Amnoid_bead","BMP_noMEF","BMP_MEF","SB43_MEF","CHIR_MEF")
#
#
#
#expe2 <- c("EmDisc","EmD","ActA_MEF","ActA_noMEF","FGF_noMEF","Am","Amnoid_bead","BMP_noMEF","BMP_MEF","SB43_MEF","CHIR_MEF")
#
#subs1 <- which(Idents(mammal.combined) %in% exps)
#subs2 <- which(Idents(mammal.combined) %in% expe2)

#C1 <- cor( as.matrix( log(AvExp$RNA[,exps]+1) ), as.matrix( log(AvExp$RNA[,expe2] +1) ), method = "pearson")
#Mat1 <- GetAssayData(mammal.combined,assay = "RNA")
#C1 <- cor(as.matrix(Mat1[,subs1]),as.matrix(Mat1[,subs2]), method = "pearson")
#
#mat_breaks <- seq(0.65, 0.85, length.out = 20)
#pheatmap(C1,color =  redblue1(20),gaps_col=c(5), gaps_row=c(3,6,9,12,15), breaks=mat_breaks,border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Alltissueheamtap",".pdf",sep="") ,width=4,height=4)
#marmoset_data_InVitro <- subset(marmoset_data_InVitro, idents = c("EmD","EmDisc","Am","Amnoid_bead","CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"))
#marmoset_data_InVitro <- NormalizeData(marmoset_data_InVitro, verbose = FALSE)
#marmoset_data_InVitro$Dataset <- "1) Marmoset in vitro"
#marmoset_data_InVitro$Stage <- "In vitro"

#Load marmoset in vivo data
#BS <- read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/Keycorrect_updatedAno.csv",sep=",",header = T, row.names=1)
#raw_counts <- read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
#marmoset_data <- CreateSeuratObject(counts = raw_counts[,which(BS$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_data) <- BS$FinalAnnotation[which(BS$QC>0)]
#marmoset_data$Dataset <- "2) Marmoset in vivo"

#Extract out the primed and naive as these make good controls. Can do seperte analyses with and without them
#marmoset_dataInVivo <- subset(marmoset_data, idents = c("Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_dataInVivo <- NormalizeData(marmoset_dataInVivo, verbose = FALSE)

#marmoset_dataDylan <- subset(marmoset_data, idents = c("Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Tb_CS3"))
#marmoset_dataDylan <- NormalizeData(marmoset_dataDylan, verbose = FALSE)

#marmoset_dataDylan1 <- NormalizeData(marmoset_dataDylan, verbose = FALSE)
#marmoset_dataDylan1 <- FindVariableFeatures(marmoset_dataDylan1, selection.method = "vst", nfeatures = 2000)
#marmoset_dataDylan1 <- ScaleData(marmoset_dataDylan1, verbose = FALSE)
#marmoset_dataDylan1 <- RunPCA(marmoset_dataDylan1, npcs = 20, verbose = FALSE)
#marmoset_dataDylan1 <- RunUMAP(marmoset_dataDylan1, reduction = "pca", dims = 1:20)
#marmoset_dataDylan1 <- RunTSNE(marmoset_dataDylan1, reduction = "pca", dims = 1:20)
#marmoset_dataDylan1 <- FindNeighbors(marmoset_dataDylan1, reduction = "pca", dims = 1:20)
#DimPlot(marmoset_dataDylan1, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_Tb.pdf",sep=""),width = 8, height = 8)

#marmoset_dataDylan1 <- FindClusters(marmoset_dataDylan1, resolution = 1.0)
#marmoset_dataDylan1$Cl5 <- Idents(marmoset_dataDylan1)
#DimPlot(marmoset_dataDylan1, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_ClTb.pdf",sep=""),width = 8, height = 8)

#Cl0 <- WhichCells(marmoset_dataDylan1, idents = "0")
#Cl1 <- WhichCells(marmoset_dataDylan1, idents = "1")
#Cl2 <- WhichCells(marmoset_dataDylan1, idents = "2")
#Cl3 <- WhichCells(marmoset_dataDylan1, idents = "3")

#marmoset_dataAP <- subset(marmoset_data, idents = c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7"))
#marmoset_dataAP <- NormalizeData(marmoset_dataAP, verbose = FALSE)
#marmoset_dataAP <- FindVariableFeatures(marmoset_dataAP, selection.method = "vst", nfeatures = 2000)
#marmoset_dataAP <- ScaleData(marmoset_dataAP, verbose = FALSE)
#marmoset_dataAP <- RunPCA(marmoset_dataAP, npcs = 20, verbose = FALSE)
#marmoset_dataAP <- RunUMAP(marmoset_dataAP, reduction = "pca", dims = 1:20)
#marmoset_dataAP <- RunTSNE(marmoset_dataAP, reduction = "pca", dims = 1:20)
#marmoset_dataAP <- FindNeighbors(marmoset_dataAP, reduction = "pca", dims = 1:20)
#FeaturePlot(marmoset_dataAP,  reduction = "umap", features = "POU5F1", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP_POU5F1.pdf",sep=""),width = 8, height = 8)
#FeaturePlot(marmoset_dataAP,  reduction = "umap", features = "SOX2", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP_SOX2.pdf",sep=""),width = 8, height = 8)
#DimPlot(marmoset_dataAP, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCAmarmoset_dataAP.pdf",sep=""),width = 8, height = 8)

#FeaturePlot(marmoset_dataAP,  reduction = "umap", features = "TFAP2C", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP_TFAP2C.pdf",sep=""),width = 8, height = 8)
#FeaturePlot(marmoset_dataAP,  reduction = "umap", features = "TFAP2A", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP_TFAP2A.pdf",sep=""),width = 8, height = 8)

#FeaturePlot(marmoset_dataAP,  reduction = "umap", features = "T", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP_T.pdf",sep=""),width = 8, height = 8)


#DimPlot(marmoset_dataAP, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAP.pdf",sep=""),width = 8, height = 8)
#marmoset_dataAP$ID <- Idents(marmoset_dataAP)
#marmoset_dataAP <- FindClusters(marmoset_dataAP, resolution = 1)
#DimPlot(marmoset_dataAP, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmosetCl_dataAP.pdf",sep=""),width = 8, height = 8)

#marmoset_dataAPam <- subset(marmoset_data, idents = c("Am_CS5","Am_CS6","Am_CS7"))
#marmoset_dataAPam <- NormalizeData(marmoset_dataAPam, verbose = FALSE)
#marmoset_dataAPam <- FindVariableFeatures(marmoset_dataAPam, selection.method = "vst", nfeatures = 2000)
#marmoset_dataAPam <- ScaleData(marmoset_dataAPam, verbose = FALSE)
#marmoset_dataAPam <- RunPCA(marmoset_dataAPam, npcs = 20, verbose = FALSE)
#marmoset_dataAPam <- RunUMAP(marmoset_dataAPam, reduction = "pca", dims = 1:20)
#marmoset_dataAPam <- RunTSNE(marmoset_dataAPam, reduction = "pca", dims = 1:20)
#marmoset_dataAPam <- FindNeighbors(marmoset_dataAPam, reduction = "pca", dims = 1:20)

#FeaturePlot(marmoset_dataAPam,  reduction = "umap", features = "POU5F1", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAm_POU5F1.pdf",sep=""),width = 8, height = 8)
#FeaturePlot(marmoset_dataAPam,  reduction = "umap", features = "VTCN1", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAm_VTCN1.pdf",sep=""),width = 8, height = 8)
#FeaturePlot(marmoset_dataAPam,  reduction = "umap", features = "TFAP2A", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAm_TFAP2A.pdf",sep=""),width = 8, height = 8)

#FeaturePlot(marmoset_dataAPam,  reduction = "umap", features = "TFAP2C", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAm_TFAP2C.pdf",sep=""),width = 8, height = 8)

#DimPlot(marmoset_dataAPam, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCAmarmoset_dataAm.pdf",sep=""),width = 8, height = 8)
#DimPlot(marmoset_dataAPam, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmoset_dataAm.pdf",sep=""),width = 8, height = 8)
#marmoset_dataAPam$ID <- Idents(marmoset_dataAPam)
#marmoset_dataAPam <- FindClusters(marmoset_dataAPam, resolution = 1)
#DimPlot(marmoset_dataAPam, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAPmarmosetCl_dataAm.pdf",sep=""),width = 8, height = 8)


#Extract other in vitro datasets 
#marmoset_dataInVitro1<- subset(marmoset_data, idents = c("ESC_primed"))
#marmoset_dataInVitro1 <- NormalizeData(marmoset_dataInVitro1, verbose = FALSE)
#marmoset_dataInVitro1$Dataset <- "1) Marmoset in vitro"

#Load the in vitro datasets
#/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/InVitro/InVitroKey_merged.csv 
#/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/InVitro/featurecountsInVitro_merged.csv 


#BS2 <- read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/Annot.csv",sep=",",header = T, row.names=1)
#raw_counts2<-read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/featureCounts_Ensembl92.csv",sep=",",header = T, row.names=1)
#marmoset_dataInVitro2 <- CreateSeuratObject(counts = raw_counts2[,which(BS2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_dataInVitro2) <- BS2$ID[which(BS2$QC>0)]
#marmoset_dataInVitro2$Dataset <- "1) Marmoset in vitro"
#marmoset_dataInVitro2 <- subset(marmoset_dataInVitro2, idents = c("EmD","EmDisc","Am","Amnoid_bead","CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"))
#marmoset_dataInVitro2 <- NormalizeData(marmoset_dataInVitro2, verbose = FALSE)
#marmoset_dataInVitro <- merge(marmoset_dataInVitro1, y = c(marmoset_dataInVitro2), project = "merged")



#BS2 <- read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/Annot.csv",sep=",",header = T, row.names=1)
#raw_counts2<-read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/featureCounts_Ensembl92.csv",sep=",",header = T, row.names=1)
##marmoset_dataInVitro2 <- CreateSeuratObject(counts = raw_counts2[,which(BS2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_dataInVitro2) <- BS2$ID[which(BS2$QC>0)]
#marmoset_dataInVitro2$Dataset <- "1) Marmoset in vitro"
#marmoset_dataDylanIV1 <- subset(marmoset_dataInVitro2, idents = c("newTSP3","OKAEP5esc"))
#marmoset_dataDylanIV1 <- NormalizeData(marmoset_dataDylanIV1, verbose = FALSE)
#marmoset_dataDylanIV1 <- FindVariableFeatures(marmoset_dataDylanIV1, selection.method = "vst", nfeatures = 2000)
#marmoset_dataDylanIV1 <- ScaleData(marmoset_dataDylanIV1, verbose = FALSE)
#marmoset_dataDylanIV1 <- RunPCA(marmoset_dataDylanIV1, npcs = 20, verbose = FALSE)
#marmoset_dataDylanIV1 <- RunUMAP(marmoset_dataDylanIV1, reduction = "pca", dims = 1:20)
#marmoset_dataDylanIV1 <- RunTSNE(marmoset_dataDylanIV1, reduction = "pca", dims = 1:20)
#marmoset_dataDylanIV1 <- FindNeighbors(marmoset_dataDylanIV1, reduction = "pca", dims = 1:20)
#DimPlot(marmoset_dataDylanIV1, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_TbInVitro.pdf",sep=""),width = 8, height = 8)




#Load the in vitro datasets
#BS2 <- read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/Annot.csv",sep=",",header = T, row.names=1)
#raw_counts2<-read.table("/Volumes/GoogleDrive/Shared\ drives/FH-TEB\ Collaboration/CRUK\ sequencing/featureCounts_Ensembl92.csv",sep=",",header = T, row.names=1)
#marmoset_dataInVitro3 <- CreateSeuratObject(counts = raw_counts2[,which(BS2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
#Idents(marmoset_dataInVitro3) <- BS2$ID[which(BS2$QC>0)]
#marmoset_dataInVitro3$Dataset <- "1) Marmoset in vitro"
#marmoset_dataInVitro3 <- subset(marmoset_dataInVitro3, idents = c("HYPO","newTSP3","OKAEP5esc"))
#marmoset_dataInVitro3 <- NormalizeData(marmoset_dataInVitro3, verbose = FALSE)

#marmoset_dataDylanMax <- merge(marmoset_dataInVitro3, y = c(marmoset_dataDylan), project = "merged")
#Idents(marmoset_dataDylanMax, cells=Cl0) <- "Cl0"
#Idents(marmoset_dataDylanMax, cells=Cl1) <- "Cl1"
#Idents(marmoset_dataDylanMax, cells=Cl2) <- "Cl2"
#Idents(marmoset_dataDylanMax, cells=Cl3) <- "Cl3"


#which.nonnum <- function(x) {
#  which(is.na(suppressWarnings(as.numeric(as.character(x)))))
#}
#which.nonnum(as.data.frame(raw_counts2)[[1]])


#mammal.combined <- merge(marmoset_data, y = c(marmoset_data2), project = "merged")
#marmoset_datax <- FindVariableFeatures(marmoset_datax, selection.method = "vst", nfeatures = 20000)

#marmoset_datax <- ScaleData(marmoset_datax, verbose = FALSE)
#marmoset_datax <- RunPCA(marmoset_datax, npcs = 20, verbose = FALSE)
#marmoset_datax <- RunUMAP(marmoset_datax, reduction = "pca", dims = 1:20)
#marmoset_datax <- RunTSNE(marmoset_datax, reduction = "pca", dims = 1:20)
#marmoset_datax <- FindNeighbors(marmoset_datax, reduction = "pca", dims = 1:20)

#Local run for marker ID
#marmoset_dataDylanMax <- FindVariableFeatures(marmoset_dataDylanMax, selection.method = "vst", nfeatures = 20000)
#marmoset_dataDylanMax <- ScaleData(marmoset_dataDylanMax, verbose = FALSE)
#marmoset_dataDylanMax <- RunPCA(marmoset_dataDylanMax, npcs = 20, verbose = FALSE)
#marmoset_dataDylanMax <- RunUMAP(marmoset_dataDylanMax, reduction = "pca", dims = 1:20)
#marmoset_dataDylanMax <- RunTSNE(marmoset_dataDylanMax, reduction = "pca", dims = 1:20)
#marmoset_dataDylanMax <- FindNeighbors(marmoset_dataDylanMax, reduction = "pca", dims = 1:20)

mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
Idents(mammal.combined2) <- mammal.combined2$Cl3
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by = "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cluster","_split_wPre.pdf",sep=""),width = 16, height = 8)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Cluster","_split_wPre.pdf",sep=""),width = 16, height = 8)

uID <- paste(mammal.combined2$Dataset,mammal.combined2$Cl3,sep="_")
Idents(mammal.combined2) <- uID

#Now do the volcano plots
DefaultAssay(mammal.combined2) <- "RNA"
Cl1 <- FindMarkers(mammal.combined2, ident.2 = c("2) Marmoset in vivo_1"), ident.1 = c("InVitro_1"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined2, ident.2 = c("2) Marmoset in vivo_0"), ident.1 = c("InVitro_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(mammal.combined2, ident.2 = c("2) Marmoset in vivo_2"), ident.1 = c("InVitro_2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl2_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"Cl2_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo, marmoset_data_InVitro), dims = 1:20, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- merge(marmoset_data, y = c(marmoset_data2), project = "merged")
#mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 20000)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_20dim_5k_alltissues_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","20dim_5k_alltissues_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_20dim_5k_alltissues_merge.pdf",sep=""),width = 8, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_20dim_5k_alltissues_merge.pdf",sep=""),width = 8, height = 8)





DefaultAssay(mammal.combined) <- "RNA"
AvExp <- AverageExpression(mammal.combined)

exps <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
          "Tb_CS5","Tb_CS6","Tb_CS7",
          "ExMes_CS5","ExMes_CS6","ExMes_CS7",
          "VE_CS5","VE_CS6","VE_CS7",
          "SYS_CS5","SYS_CS6","SYS_CS7")
                    
expe2 <- c("EmDisc","EmD","ActA_MEF","ActA_noMEF","FGF_noMEF","Am","Amnoid_bead","BMP_noMEF","BMP_MEF","SB43_MEF","CHIR_MEF")

subs1 <- which(Idents(mammal.combined) %in% exps)
subs2 <- which(Idents(mammal.combined) %in% expe2)

C1 <- cor( as.matrix( log(AvExp$RNA[,exps]+1) ), as.matrix( log(AvExp$RNA[,expe2] +1) ), method = "pearson")

Mat1 <- GetAssayData(mammal.combined,assay = "RNA")
C1 <- cor(as.matrix(Mat1[,subs1]),as.matrix(Mat1[,subs2]), method = "pearson")
  
  
mat_breaks <- seq(0.65, 0.85, length.out = 20)
pheatmap(C1,color =  redblue1(20),gaps_col=c(5), gaps_row=c(3,6,9,12,15), breaks=mat_breaks,border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Alltissueheamtap",".pdf",sep="") ,width=4,height=4)



write.csv(as.data.frame(C1), file=paste(saveext,"/Correlation_alltissues_pearson_RNA.csv",sep=""))

write.csv(as.data.frame(Idents(mammal.combined)[subs1]), file=paste(saveext,"/Idents_alltissues1.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[subs2]), file=paste(saveext,"/Idents_alltissues2.csv",sep=""))



#DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Stage", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_20dim_5k_splitstage.pdf",sep=""),width = 40, height = 8)

#DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Stage", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","20dim_5k_splitstage.pdf",sep=""),width = 40, height = 8)


#DefaultAssay(mammal.combined) <- "RNA"
#AvExp <- AverageExpression(mammal.combined)

#Which genes/conditons
exps <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc","EmD","ActA_MEF","ActA_noMEF","FGF_noMEF","Am_CS5","Am_CS6","Am_CS7","Am","Amnoid_bead","BMP_noMEF","BMP_MEF","SB43_MEF","CHIR_MEF")
#FortyEight <- c("OOEP","TCL1A","WEE2",
#                "IL1RN","NOV","ZNF80",
#                "SPIC","ESRRB","STAT3",
#                "KLF17","SOX15","POU5F1","NANOG","SOX2","SFRP2","DNMT3B","T",
#                "TFAP2C","TFAP2A","VTCN1","PRDM1","PRDM14","NANOS3",
#                "GATA6","GATA4","PDGFRA","SOX17","OTX2","APOA1","TTR","APOB","HAND2","TBX4","HGF",
#                "JAM2","GATA3","GATA2","CGB3","CGA") 


Cl1 <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7"), ident.1 = c("Am_CS5","Am_CS6","Am_CS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Cluster_EmvAm_CS567.csv",sep=""))

Cl1e <- FindMarkers(mammal.combined, ident.2 = c("EmDisc_CS5","EmDisc_CS6"), ident.1 = c("Am_CS5","Am_CS6"), verbose = FALSE, test.use = "MAST", only.pos = FALSE)
write.csv(as.data.frame(Cl1), file=paste(saveext,"/Cluster_EmvAm_CS56.csv",sep=""))



FortyEight <- c("SOX2","POU5F1","NANOG","SFRP2","DNMT3B",
                "T","MIXL1","EOMES",
                "LHX1",
                "SNAI2",
                "PDGFRA",
                "SOX17",
                "FOXA2",
                "GATA3","TFAP2C","TFAP2A","VTCN1","WNT6","ISL1")



X <- (AvExp$RNA[FortyEight,exps])

mat_breaks <- seq(0, 2, length.out = 20)

annotationL <- c("EmDisc_CS5",
           "EmDisc_CS6",
           "EmDisc_CS7",
           "EmD","EmDisc",
           "ActA_MEF","ActA_noMEF","FGF_noMEF",
           "Am_CS5",
           "Am_CS6",
           "Am_CS7",
           "Am","Amnoid_bead",
           "BMP_noMEF","BMP_MEF",
           "SB43_MEF","CHIR_MEF")
mycolors <- c("#0c9cf5",
             "#0767DA",
             "#0233BF",
             "#c2a5cf","#7b3294","#4d4d4d",
             "#0571b0","#92c5de",
             "#877bd6",
             "#5F54C7",
             "#1A0873",
             "#008837","#a6dba0",
             "#ca0020","#f4a582",
             "#d8b365","#5ab4ac")
   
#avexp  <- AverageExpression(object = mammal.combined, slot = "data")
#a <- avexp$RNA

annotation_col = data.frame(Stage = factor(colnames(X)))
rownames(annotation_col) <- colnames(X)

names(mycolors) <- colnames(X)
anno_colors <- list(Stage = mycolors)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

mat_breaks <- seq(0, 20, length.out = 20)
pheatmap(X*100,color =  redblue1(20),gaps_col=c(3,8,11), gaps_row=c(5,9,11,13),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/heamtap",".pdf",sep="") ,width=5,height=4)
mat_breaks <- seq(-2, 2, length.out = 20)
pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(3,8,11),gaps_row=c(5,9,11,13),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/heamtapscale",".pdf",sep=""),width=5,height=4 )




mat_breaks <- seq(-2, 2, length.out = 20)
X <- (AvExp$RNA[c(rownames(Cl1e)[1:100],rownames(Cl1e)[2232:2332]),exps])
X <- (AvExp$RNA[c(rownames(Cl1)[1:100],rownames(Cl1)[2801:2901]),exps])
pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(3,8,11),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/heamtapscale100",".pdf",sep="") ,width=8,height=60)

AvExp[genes_,treats_]
C1 <- cor( log(as.matrix(AvExp$integrated[,TypeKeep]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=TRUE,cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterInt",".pdf",sep=""),width=8,height=8)



Expe1 <- GetAssayData(mammal.combined, assay = "RNA")
Expe2 <- GetAssayData(mammal.combined, assay = "integrated")
ind1 <- which(mammal.combined$Dataset=="2) Marmoset in vivo")
ind2 <- which(mammal.combined$Dataset=="1) Marmoset in vitro")

C1 <- cor( as.matrix(Expe1[,ind1]), as.matrix(Expe1[,ind2]), method = "pearson")
C2 <- cor( as.matrix(Expe2[,ind1]), as.matrix(Expe2[,ind2]), method = "pearson")


order1 <- c(which(grepl("EmDisc_CS5",mammal.combined$IDs) ),which(grepl("EmDisc_CS6",mammal.combined$IDs) ),which(grepl("EmDisc_CS7",mammal.combined$IDs) ),which(grepl("Am_CS5",mammal.combined$IDs) ),which(grepl("Am_CS6",mammal.combined$IDs) ),which(grepl("Am_CS7",mammal.combined$IDs)) )

order1 <- c(
  which(grepl("EmD",mammal.combined$IDs) ),
  which(grepl("EmDisc",mammal.combined$IDs) ),
  which(grepl("ActA_MEF",mammal.combined$IDs) ),
  which(grepl("ActA_noMEF",mammal.combined$IDs) ),
  which(grepl("FGF_noMEF",mammal.combined$IDs) ),
  which(grepl("CHIR_MEF",mammal.combined$IDs) ),
  which(grepl("Am",mammal.combined$IDs) ),
  which(grepl("Amnoid_bead",mammal.combined$IDs) ),
  which(grepl("BMP_noMEF",mammal.combined$IDs) ),
  
  
  
  
  
  
write.csv(as.data.frame(C1), file=paste(saveext,"/Correlation_pearson_RNA.csv",sep=""))
write.csv(as.data.frame(C2), file=paste(saveext,"/Correlation_pearson_Int.csv",sep=""))


ind1 <- which(mammal.combined$IDs=="Am_CS6")
ind2 <- which(mammal.combined$IDs=="Amnoid_bead")

C1 <- cor( as.matrix(Expe2[,ind1]), as.matrix(Expe2[,ind2]), method = "pearson")


C1 <- cor( as.matrix(Expe1[,ind1]), as.matrix(Expe1[,ind2]), method = "pearson")
C2 <- cor( as.matrix(Expe2[,ind1]), as.matrix(Expe2[,ind2]), method = "pearson")


write.csv(as.data.frame(C1), file=paste(saveext,"/Correlation_pearson_RNA.csv",sep=""))
write.csv(as.data.frame(C2), file=paste(saveext,"/Correlation_pearson_Int.csv",sep=""))


write.csv(as.data.frame(Idents(mammal.combined)[ind1]), file=paste(saveext,"/Idents1.csv",sep=""))
write.csv(as.data.frame(Idents(mammal.combined)[ind2]), file=paste(saveext,"/Idents2.csv",sep=""))










#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_data2, marmoset_datax), dims = 1:20, anchor.features = 4000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#DefaultAssay(mammal.combined) <- "integrated"

#Remove most extraembryonic tissues.
marmoset_dataInVivo2 <- subset(marmoset_dataInVivo,idents = c("Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS5","VE_CS5","VE_CS6","SYS_CS5","SYS_CS7","SYS_CS6","ExMes_CS5","ExMes_CS7","ExMes_CS6"),invert=TRUE)
#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo2, marmoset_data_InVitro), dims = 1:20, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])[1]
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Lab","_20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Lab","20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Lab","_20dim_5k_merge.pdf",sep=""),width = 8, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Lab","_20dim_5k_merge.pdf",sep=""),width = 8, height = 8)

#mammal.combined$ID2 <- Idents(mammal.combined)
#Idents(mammal.combined) <- mammal.combined$ID2



DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Stage", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Lab","_20dim_5k_splitstage.pdf",sep=""),width = 32, height = 8)


DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "ID2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Labsplit","_20dim_5k_split.pdf",sep=""),width = 86, height = 8,limitsize = FALSE)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "ID2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Labsplit","20dim_5k_split.pdf",sep=""),width = 86, height = 8,limitsize = FALSE)










DefaultAssay(mammal.combined) <- "integrated"
mammal.combined$IDs <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)

#DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Lab","_20dim_5k_merge.pdf",sep=""),width = 16, height = 16)


#cType <- c("Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF")
#BaseCol <- c("lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d")



#DefaultAssay(mammal.combined) <- "integrated"


#mammal.combined <- merge(marmoset_data, y = c(marmoset_data2), project = "merged")
#mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 20000)
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

DefaultAssay(mammal.combined) <- "RNA"
FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_TFAP2A.pdf",sep=""),width =16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_TFAP2A.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_VTCN1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_VTCN1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_SOX2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_SOX2.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_POU5F1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_POU5F1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_T.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_T.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_TFAP2C.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_TFAP2C.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_PDGFRA.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_PDGFRA.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_SOX17.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_SOX17.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "NANOS3", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_NANOS3.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "NANOS3", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_NANOS3.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "NANOG", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_NANOG.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "NANOG", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_NANOG.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "FOXA2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids_FOXA2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "FOXA2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids_FOXA2.pdf",sep=""),width = 16, height = 8)


DefaultAssay(mammal.combined) <- "RNA"
AvExp <- AverageExpression(mammal.combined)

#AvExp <- AverageExpression(marmoset_dataInVitro)

#marmoset_dataInVitro

order1 <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7")
order2 <- c("EmD","EmDisc","Am","Amnoid_bead","CHIR_MEF","ActA_MEF","ActA_noMEF","FGF_noMEF","BMP_MEF","BMP_noMEF","SB43_MEF")
C1 <- cor( log(as.matrix(AvExp$integrated[,order2]) +1), log(as.matrix(AvExp$integrated[,order1]) +1) )
C2 <- cor( log(as.matrix(AvExp$integrated[order2])+1) )


CGlobal1 <- cor( as.matrix(C1), as.matrix(C2), method = "pearson" )
CGlobal2 <- cor( as.matrix(C1), as.matrix(C2), method = "spearman" )
write.csv(as.data.frame(CGlobal1), file=paste(saveext,"/Correlation_pearson.csv",sep=""))
write.csv(as.data.frame(CGlobal2), file=paste(saveext,"/Correlation_spearman.csv",sep=""))


redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), breaks = mat_breaks, gaps_row = c(2,4), gaps_col = c(3,5), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoBulkCrossCorrelation",".pdf",sep=""),width=3.7,height=4)

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.90, length.out = 60)
pheatmap(C2, color =  redblue1(60), breaks = mat_breaks, border_color = NA, gaps_row = c(2,4), gaps_col = c(2,4), display_numbers = round(C2,2), fontsize_number = 5,  cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/PseudoBulkInvitroCorrelation",".pdf",sep=""),width=4.2,height=4)

C1 <- cor( log(as.matrix(AvExp$RNA[,order2]) +1), log(as.matrix(AvExp$RNA[,order1]) +1) )
C2 <- cor( log(as.matrix(AvExp$RNA[order2])+1) )

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), breaks = mat_breaks, gaps_row = c(2,4), gaps_col = c(3,5), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation",".pdf",sep=""),width=3.7,height=4)

#First do some clustering to break down the groups?
mammal.combined$ID2 <- Idents(mammal.combined)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined$IDs <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Cl5","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Cl5","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Cl5","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Cl5","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.9)
mammal.combined$Cl9 <- Idents(mammal.combined)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type2_Cl9","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type2_Cl9","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

uID <- paste("Cl",mammal.combined$Cl5,"_",mammal.combined$IDs,sep="")
Idents(mammal.combined) <- uID

TypeCount <- table(Idents(mammal.combined))
TypeKeep <- rownames(TypeCount)[which(TypeCount>=5)]

DefaultAssay(mammal.combined) <- "RNA"
AvExp <- AverageExpression(mammal.combined)

C1 <- cor( log(as.matrix(AvExp$integrated[,TypeKeep]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=TRUE,cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterInt",".pdf",sep=""),width=8,height=8)

C1 <- cor( log(as.matrix(AvExp$RNA[,TypeKeep]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=TRUE,cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterRNA",".pdf",sep=""),width=8,height=8)

list1 <- c("Cl0_EmDisc_CS5",  
"Cl0_EmDisc_CS6", 
"Cl0_EmDisc_CS7",
"Cl3_EmDisc_CS6",  
"Cl3_EmDisc_CS7",  
"Cl3_Am_CS6",  
"Cl1_Am_CS6",    
"Cl1_Am_CS7" ,
"Cl2_Am_CS7")


list2 <- c("Cl0_EmD",   
"Cl0_EmDisc", 
"Cl0_ActA_MEF",
"Cl0_ActA_noMEF",
"Cl3_ActA_noMEF", 
"Cl3_FGF_noMEF",
"Cl3_CHIR_MEF",
"Cl3_Am", 
"Cl1_Amnoid_bead",
"Cl1_Am",
"Cl1_BMP_noMEF",  
"Cl1_BMP_MEF",
"Cl1_SB43_MEF",
"Cl2_Am",
"Cl2_BMP_MEF",
"Cl2_FGF_noMEF",        
"Cl2_SB43_MEF",
"Cl2_CHIR_MEF")

C1 <- cor( log(as.matrix(AvExp$RNA[,list1]) +1), log(as.matrix(AvExp$RNA[,list2]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60),gaps_row = c(3,6,8), gaps_col = c(4,8,13), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterRNACross",".pdf",sep=""),width=8,height=4)

C1 <- cor( log(as.matrix(AvExp$integrated[,list1]) +1), log(as.matrix(AvExp$integrated[,list2]) +1) )
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60),gaps_row = c(3,6,8), gaps_col = c(4,8,13), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation_ClusterIntCross",".pdf",sep=""),width=8,height=4)







#Add in preimplantation
marmoset_data_Pre <- CreateSeuratObject(counts = raw_countsPre[,which(KeyPre$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Pre) <- KeyPre$Primary.lineage[which(KeyPre$QC>0)]
marmoset_data_Pre$LOC <- KeyPre$Location[which(KeyPre$QC>0)]
marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Epi_CS3"))
marmoset_data_Pre <- NormalizeData(marmoset_data_Pre, verbose = FALSE)
marmoset_data_Pre$Stage <- "Pre"
marmoset_data_Pre$Dataset <- "InVivo"

marmoset_dataInVivo <- merge(marmoset_data_CS5, y = c(marmoset_data_Pre,marmoset_data_CS6,marmoset_data_CS7,marmoset_data_CRUK,marmoset_data_CRUK2,marmoset_data_CRUK3, marmoset_data_20307, marmoset_data_20308), project = "merged")
marmoset_dataInVivo$Dataset <- "2) Marmoset in vivo"


#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo, marmoset_data_InVitro), dims = 1:20, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- merge(marmoset_data, y = c(marmoset_data2), project = "merged")
#mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 20000)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#Now do dim red plots
uLab <- as.character(Idents(mammal.combined))
uLab[which(uLab=="Am_CS5_PGC")] <- "Am_CS5"
uLab[which(uLab=="Am_CS7_EmDisc")] <- "Am_CS7"
uLab[which(uLab=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
uLab[which(uLab=="PGC_CS5")] <- "Am_CS5"
uLab[which(uLab=="Am_CS6_EmDisc")] <- "Am_CS6"
uLab[which(uLab=="EmDisc_CS5_Am")] <- "EmDisc_CS5"
uLab[which(uLab=="EmDisc_CS5_Gast")] <- "EmDisc_CS5"
uLab[which(uLab=="EmDisc_CS6_Am")] <- "EmDisc_CS6"
uLab[which(uLab=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
uLab[which(uLab=="EmDisc_CS7_Gast")] <- "EmDisc_CS7"
uLab[which(uLab=="EmDisc_Gast_CS6")] <- "EmDisc_CS6"
uLab[which(uLab=="EmDisc_Gast_CS7")] <- "EmDisc_CS7"
uLab[which(uLab=="Stalk_CS7_PGC")] <- "PGC_CS7"

Idents(mammal.combined) <- uLab
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="Amnoid_bead")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Am")] <- 1
cluster_letters[which(Idents(mammal.combined)=="EmDisc")] <- 2
cluster_letters[which(Idents(mammal.combined)=="EmD")] <- 2
cluster_letters[which(Idents(mammal.combined)=="CHIR_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="SB43_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="ActA_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="BMP_MEF")] <- 2
cluster_letters[which(Idents(mammal.combined)=="FGF_noMEF")] <- 1
cluster_letters[which(Idents(mammal.combined)=="BMP_noMEF")] <- 1
cluster_letters[which(Idents(mammal.combined)=="ActA_noMEF")] <- 1

cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined, metadata = cluster_letters,col.name = 'cell.orig')

DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by = "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split_wPre.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split_wPre.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by = "cell.orig", reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge_wPre.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge_wPre.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


mammal.combined2 <- subset(mammal.combined,idents=c("CHIR_MEF","SB43_MEF","BMP_MEF","ActA_MEF","FGF_noMEF","BMP_noMEF","ActA_noMEF"),invert=TRUE)
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_split1_wPre.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_split1_wPre.pdf",sep=""),width = 16, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig",  reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_merge1_wPre.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined2, cols = coluse, pt.size = 4, shape.by =  "cell.orig", reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_merge1_wPre.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



#NOW DO VOLCANO
Ae <- AvExp$RNA
Ae$gene <- rownames(Ae)
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

#Now do the volcano plots
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS5"), ident.1 = c("Cl0_EmD"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_EmD_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_EmD_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS6"), ident.1 = c("Cl0_EmD"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_Em_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_EmD_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS7"), ident.1 = c("Cl0_EmD"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_Em_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_EmD_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS5"), ident.1 = c("Cl0_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_MEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_MEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS5"), ident.1 = c("Cl0_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_noMEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_noMEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS5"), ident.1 = c("Cl0_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_MEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl0_ActA_MEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS6"), ident.1 = c("Cl0_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_ActA_noMEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_ActA_noMEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS6"), ident.1 = c("Cl0_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_ActA_MEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl0_ActA_MEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS7"), ident.1 = c("Cl0_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_ActA_noMEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_ActA_noMEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl0_EmDisc_CS7"), ident.1 = c("Cl0_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_ActA_MEF_Cl0_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl0_ActA_MEF_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Now cluster3 ...
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS5"), ident.1 = c("Cl3_FGF_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_FGF_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_FGF_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS6"), ident.1 = c("Cl3_FGF_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_FGF_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_FGF_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS7"), ident.1 = c("Cl3_FGF_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_FGF_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_FGF_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS5"), ident.1 = c("Cl3_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_MEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_MEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS5"), ident.1 = c("Cl3_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS5"), ident.1 = c("Cl3_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_MEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS5_Cl3_ActA_MEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS6"), ident.1 = c("Cl3_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_ActA_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_ActA_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS6"), ident.1 = c("Cl3_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_ActA_MEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS6_Cl3_ActA_MEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS7"), ident.1 = c("Cl3_ActA_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_ActA_noMEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_ActA_noMEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl3_EmDisc_CS7"), ident.1 = c("Cl3_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_ActA_MEF_Cl3_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDiscCS7_Cl3_ActA_MEF_Cl3_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS6"), ident.1 = c("Cl1_Am"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_Am_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_Am_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS6"), ident.1 = c("Cl1_Am"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_Am_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_Am_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS7"), ident.1 = c("Cl1_Am"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_Am_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_Am_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS6"), ident.1 = c("Cl1_BMP_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_BMP_noMEF_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_BMP_noMEF_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS7"), ident.1 = c("Cl1_BMP_noMEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_BMP_noMEF_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_BMP_noMEF_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS6"), ident.1 = c("Cl1_BMP_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_BMP_MEF_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_BMP_MEF_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS7"), ident.1 = c("Cl1_BMP_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_BMP_MEF_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_BMP_MEF_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS6"), ident.1 = c("Cl1_SB43_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_SB43_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS6_Cl1_SB43_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cluster 1 ... amnion
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("Cl1_Am_CS7"), ident.1 = c("Cl1_SB43_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_SB43_Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"AmCS7_Cl1_SB43_Cl1_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#
#
#
#
#

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(mammal.combined, ident.2 = "Cl0_EmDisc_CS7", ident.1 = "Cl0_EmD",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(mammal.combined, ident.2 = "Cl1_Am_CS7", ident.1 = "Cl1_Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC7_Cl0Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC7_Cl0Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(mammal.combined, ident.2 = "Cl0_EmDisc_CS6", ident.1 = "Cl0_EmD",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(mammal.combined, ident.2 = "Cl1_Am_CS6", ident.1 = "Cl1_Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC6_Cl0Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC6_Cl0Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(mammal.combined, ident.2 = "Cl3_EmDisc_CS7", ident.1 = "Cl3_FGF_noMEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(mammal.combined, ident.2 = "Cl1_Am_CS7", ident.1 = "Cl1_Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC7_Cl3Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC7_Cl3Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(mammal.combined, ident.2 = "Cl3_EmDisc_CS6", ident.1 = "Cl3_FGF_noMEF",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(mammal.combined, ident.2 = "Cl1_Am_CS6", ident.1 = "Cl1_Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC6_Cl3Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_Am_EmD_CSC6_Cl3Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Next compaer to epiblast
saveRDS(Idents(mammal.combined),file=paste(saveext,"ClIDs.rds",sep=""))

KeyPre <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/CS1-3/CS1-3Key.csv",sep=",",header = T, row.names=1)
raw_countsPre <- read.table("/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Annotations_and_manifests/CS1-3/featurecountsCS1-3.csv",sep=",",header = T, row.names=1)


marmoset_data_Pre <- CreateSeuratObject(counts = raw_countsPre[,which(KeyPre$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Pre) <- KeyPre$Primary.Annotation[which(KeyPre$QC>0)]
marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Epi_CS3_"))
marmoset_data_Pre <- NormalizeData(marmoset_data_Pre, verbose = FALSE)

list1 <- WhichCells(mammal.combined,idents="Cl0_EmD")
list2 <- WhichCells(mammal.combined,idents="Cl0_EmDisc_CS6")
list3 <- WhichCells(mammal.combined,idents="Cl0_EmDisc_CS7")
list4 <- WhichCells(mammal.combined,idents="Cl1_Am_CS6")
list5 <- WhichCells(mammal.combined,idents="Cl1_Am_CS7")
list6 <- WhichCells(mammal.combined,idents="Cl3_EmDisc_CS6")
list7 <- WhichCells(mammal.combined,idents="Cl3_EmDisc_CS7")
list8 <- WhichCells(mammal.combined,idents="Cl0_ActA_MEF")
list9 <- WhichCells(mammal.combined,idents="Cl0_ActA_noMEF")
list10 <- WhichCells(mammal.combined,idents="Cl3_ActA_noMEF")
list11 <- WhichCells(mammal.combined,idents="Cl3_FGF_noMEF")
list12 <- WhichCells(mammal.combined,idents="Cl3_CHIR_MEF")
list13 <- WhichCells(mammal.combined,idents="Cl1_Am")
list14 <- WhichCells(mammal.combined,idents="Cl1_BMP_noMEF")
list15 <- WhichCells(mammal.combined,idents="Cl1_BMP_MEF")

marmoset_dataMerge2 <- merge(marmoset_data_CS5, y = c(marmoset_data_CS6,marmoset_data_CS7,marmoset_data_Pre,marmoset_data_InVitro), project = "merged")
Idents(object=marmoset_dataMerge2,cells=list1) <- "Cl0_EmD"
Idents(object=marmoset_dataMerge2,cells=list2) <- "Cl0_EmDisc_CS6"
Idents(object=marmoset_dataMerge2,cells=list3) <- "Cl0_EmDisc_CS7"
Idents(object=marmoset_dataMerge2,cells=list4) <- "Cl1_Am_CS6"
Idents(object=marmoset_dataMerge2,cells=list5) <- "Cl1_Am_CS7"
Idents(object=marmoset_dataMerge2,cells=list6) <- "Cl3_EmDisc_CS6"
Idents(object=marmoset_dataMerge2,cells=list7) <- "Cl3_EmDisc_CS7"
Idents(object=marmoset_dataMerge2,cells=list8) <- "Cl0_ActA_MEF"
Idents(object=marmoset_dataMerge2,cells=list9) <- "CCl0_ActA_noMEF"
Idents(object=marmoset_dataMerge2,cells=list10) <- "Cl3_ActA_noMEF"
Idents(object=marmoset_dataMerge2,cells=list11) <- "Cl3_FGF_noMEF"
Idents(object=marmoset_dataMerge2,cells=list12) <- "Cl3_CHIR_MEF"
Idents(object=marmoset_dataMerge2,cells=list13) <- "Cl1_Am"
Idents(object=marmoset_dataMerge2,cells=list14) <- "Cl1_BMP_noMEF"
Idents(object=marmoset_dataMerge2,cells=list15) <- "Cl1_BMP_MEF"

list16 <- WhichCells(mammal.combined,idents="Cl0_EmDisc_CS5")
Idents(object=marmoset_dataMerge2,cells=list16) <- "Cl0_EmDisc_CS5"
list17 <- WhichCells(mammal.combined,idents="Cl1_Am_CS5")
Idents(object=marmoset_dataMerge2,cells=list17) <- "Cl1_Am_CS5"


DefaultAssay(marmoset_dataMerge2) <- "RNA"
AvExp2 <- AverageExpression(marmoset_dataMerge2)

Ae <- AvExp2$RNA
Ae$gene <- rownames(Ae)

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl0_EmDisc_CS7", ident.1 = "Epi_CS3_",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl1_Am_CS7", ident.1 = "Epi_CS3_", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmDAm_CS7_Epi_Cl0Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmDAm_CS7_Epi_Cl0Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl0_EmDisc_CS5", ident.1 = "Epi_CS3_",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl1_Am_CS5", ident.1 = "Epi_CS3_", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmDAm_CS5_Epi_Cl0Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmDAm_CS5_Epi_Cl0Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Now cross comparisons EmD - EmDiscc_CS6+, Am - Am_CS5+
Cl1 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl0_EmD", ident.1 = "Epi_CS3_",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl1_Am", ident.1 = "Epi_CS3_", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmAm_Epi_Cl0Cl1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "", y = "")
ggsave(filename=paste(saveext,"Scatter_EmDAm_Epi_Cl0Cl1_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Next ...
Cl1 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl0_EmD", ident.1 = "Epi_CS3_",verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl2 <- FindMarkers(marmoset_dataMerge2, ident.2 = "Cl1_Am", ident.1 = "Epi_CS3_", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl1 <- FindMarkers(marmoset_dataMerge2, ident.2 = c("Cl0_EmDisc_CS7"), ident.1 = c("Cl3_ActA_MEF"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDAm_Cl0_Epi_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"EmDAm_Cl0_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





marmoset_dataInVivo2 <- subset(marmoset_dataInVivo,idents = c("Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS5","VE_CS5","VE_CS6","SYS_CS5","SYS_CS7","SYS_CS6","ExMes_CS5","ExMes_CS6","ExMes_CS7"),invert=TRUE)
#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo2, marmoset_dataInVitro2), dims = 1:20, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_Lab","_20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type3_Lab","20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type3_Lab","_20dim_5k_merge.pdf",sep=""),width = 16, height = 16)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_Lab","_20dim_5k_merge.pdf",sep=""),width = 16, height = 16)


DefaultAssay(mammal.combined) <- "RNA"
FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2C", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_TFAP2C.pdf",sep=""),width =16, height = 8)


DefaultAssay(mammal.combined) <- "RNA"
FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_TFAP2A.pdf",sep=""),width =16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2A", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_TFAP2A.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_VTCN1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "VTCN1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_VTCN1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_SOX2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_SOX2.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_POU5F1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "POU5F1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_POU5F1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_T.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "T", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_T.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_PDGFRA.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "PDGFRA", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_PDGFRA.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_SOX17.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX17", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_SOX17.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "HESX1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_HESX1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "HESX1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_HESX1.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "MIXL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_MIXL1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "MIXL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_MIXL1.pdf",sep=""),width = 16, height = 8)


FeaturePlot(mammal.combined,  reduction = "pca", features = "EOMES", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_EOMES.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "EOMES", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_EOMES.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "MESP1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_MESP2.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "MESP2", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_MESP2.pdf",sep=""),width = 16, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "ISL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_amnioids3_ISL1.pdf",sep=""),width = 16, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "ISL1", combine=TRUE, split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_amnioids3_ISL1.pdf",sep=""),width = 16, height = 8)


DefaultAssay(mammal.combined) <- "integrated"
mammal.combined$IDs <- Idents(mammal.combined)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_Cl5","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)
mammal.combined$Cl10 <- Idents(mammal.combined)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_Cl10","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
mammal.combined$Cl15 <- Idents(mammal.combined)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_Cl15","_20dim_5k_merge.pdf",sep=""),width = 16, height = 8)

ID11 <- as.character(mammal.combined$Dataset)
ID11[which(ID11=="1) Marmoset in vitro")] <- "InVitro"
ID11[which(ID11=="2) Marmoset in vivo")] <- "InVivo"
IDCl <- paste(ID11,mammal.combined$Cl10,sep="")
Idents(mammal.combined) <- IDCl



#0,5,1,8,3,2,4
#DefaultAssay(marmoset_dataDylanMax) <- "RNA"
Ae <- AverageExpression(marmoset_dataDylanMax, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)



SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)

SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1


AvExp <- AverageExpression(marmoset_dataInVitro)
DefaultAssay(marmoset_dataInVitro) <- "RNA"

Ae <- AvExp$RNA
Ae$gene <- rownames(Ae)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "EmDisc_CS5", ident.1 = "EmD", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmDisc_CS5,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS5_EmD_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmDisc_CS5,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS5_EmD_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "Am_CS6", ident.1 = "Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am_CS6,Am)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_CS6_Am_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am_CS6,Am)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_CS6_Am_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "EmDisc_CS6", ident.1 = "EmD", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmDisc_CS6,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS6_EmD_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmDisc_CS6,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS6_EmD_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


#Run emdisc vs emd
Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "EmDisc_CS6", ident.1 = "EmD", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval <- 1
Ae5$FC <- 0
Ae5$Indi <- 0
#Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "EmDisc_CS6", ident.1 = "EmD", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5[rownames(Cl1),"Pval"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC"] <- Cl1$avg_logFC
Ae5$Indi[which( abs(Ae5$FC)>log2(1.2) & Ae5$Pval>-log2(0.05) )] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
#Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-4
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5) & Cl1$p_val_adj<0.05)],TF)
genes.to.label2 = intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5) & Cl1$p_val_adj<0.05)],SIGNAL1)
p1 <- ggplot(Ae5, aes(FC,Pval)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'blue')
ggsave(filename=paste(saveext,"Volcano_EmDisc_CS6_EmD_ECM_TF.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro, ident.2 = "EmDisc_CS6", ident.1 = "EmD", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmDisc_CS6,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS6_EmD_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmDisc_CS6,EmD)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS6_EmD_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "Amnoid_bead", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,Amnoid_bead)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_Amnoid_bead_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,Amnoid_bead)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_Amnoid_bead_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "EmDisc", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,EmDisc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_EmDisc_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,EmDisc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_EmDisc_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)



Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "Am", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,Am)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_Am_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,Am)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_Am_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)



Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmDisc", ident.1 = "Amnoid_bead", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmDisc,Amnoid_bead)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_Amnoid_bead_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmDisc,Amnoid_bead)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_Amnoid_bead_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)



Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "ActA_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,ActA_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_ActA_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,ActA_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_ActA_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "ActA_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,ActA_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_ActA_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,ActA_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_ActA_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "BMP_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,BMP_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_BMP_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,BMP_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_BMP_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "BMP_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,BMP_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_BMP_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,BMP_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_BMP_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "CHIR_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,CHIR_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_CHIR_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,CHIR_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_CHIR_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "FGF_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,FGF_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_FGF_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,FGF_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_FGF_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "EmD", ident.1 = "SB43_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmD,SB43_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_SB43_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(EmD,SB43_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmD_SB43_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)



Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "ActA_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,ActA_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_ActA_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,ActA_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_ActA_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "ActA_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,ActA_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_ActA_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,ActA_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_ActA_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "BMP_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,BMP_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_BMP_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,BMP_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_BMP_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "BMP_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,BMP_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_BMP_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,BMP_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_BMP_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "CHIR_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,CHIR_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_CHIR_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,CHIR_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_CHIR_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "FGF_noMEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,FGF_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_FGF_noMEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,FGF_noMEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_FGF_noMEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataInVitro2, ident.2 = "Am", ident.1 = "SB43_MEF", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Am,SB43_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_SB43_MEF_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
genes.to.label = intersect(rownames(Cl1),SIGNAL3)
p1 <- ggplot(Ae5, aes(Am,SB43_MEF)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_SB43_MEF_ECM.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "VE_CS5", ident.1 = "HYPO", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(VE_CS5,HYPO)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"VE_HYPO.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Tb_CS3", ident.1 = "HYPO", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Tb_CS3,HYPO)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_HYPO.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Hyp_CS3", ident.1 = "HYPO", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Hyp_CS3,HYPO)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Hyp_HYPO.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "EmDisc_CS6", ident.1 = "HYPO", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(EmDisc_CS6,HYPO)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_CS6_HYPO.pdf",sep=""),width = 13, height = 13, plot = p1)



#Now Dylan comparison
Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "OKAEP5esc", ident.1 = "newTSP3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(OKAEP5esc,newTSP3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"OKAEP5esc_newTSP3.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Cl0", ident.1 = "newTSP3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl0,newTSP3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl0_newTSP3.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Cl1", ident.1 = "newTSP3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl1,newTSP3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl1_newTSP3.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Cl2", ident.1 = "newTSP3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl2,newTSP3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl2_newTSP3.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.2 = "Cl3", ident.1 = "newTSP3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl3,newTSP3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl3_newTSP3.pdf",sep=""),width = 13, height = 13, plot = p1)



Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.1 = "OKAEP5esc", ident.2 = "Cl0", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl0,OKAEP5esc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl0_OKAEP5esc.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.1 = "OKAEP5esc", ident.2 = "Cl1", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl1,OKAEP5esc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl1_OKAEP5esc.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.1 = "OKAEP5esc", ident.2 = "Cl2", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl2,OKAEP5esc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl2_OKAEP5esc.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(marmoset_dataDylanMax, ident.1 = "OKAEP5esc", ident.2 = "Cl3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(Cl3,OKAEP5esc)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Tb_Cl3_OKAEP5esc.pdf",sep=""),width = 13, height = 13, plot = p1)



DefaultAssay(mammal.combined) <- "RNA"
Ae <- AverageExpression(mammal.combined, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)
#

TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo2", ident.1 = "InVitro2", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo2,InVitro2)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl2.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo4", ident.1 = "InVitro4", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo4,InVitro4)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl4.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo0", ident.1 = "InVitro0", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo0,InVitro0)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl0.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo5", ident.1 = "InVitro5", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo5,InVitro5)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl5.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo1", ident.1 = "InVitro1", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo1,InVitro1)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl1.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo1", ident.1 = "InVitro1", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo1,InVitro1)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl1.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo8", ident.1 = "InVitro8", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo8,InVitro8)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl8.pdf",sep=""),width = 13, height = 13, plot = p1)


Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(mammal.combined, ident.2 = "InVivo3", ident.1 = "InVitro3", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(2))
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(InVivo3,InVitro3)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Cl3.pdf",sep=""),width = 13, height = 13, plot = p1)


IDs <- as.character(mammal.combined$IDs)
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$IDs=="Am_CS5")] <- "Am_CS5"
IDs[which(mammal.combined$IDs=="Am_CS6")] <- "Am_CS6"
IDs[which(mammal.combined$IDs=="Am_CS7")] <- "Am_CS7"
IDs[which(mammal.combined$IDs=="Am")] <- "Am"
IDs[which(mammal.combined$IDs=="Amnoid_bead")] <- "Amnoid_bead"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_amnion_only.pdf",sep=""),width = 16, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type3_amnion_only.pdf",sep=""),width = 16, height = 8)


IDs <- as.character(mammal.combined$IDs)
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$IDs=="Am_CS5")] <- "Am_CS5"
IDs[which(mammal.combined$IDs=="Am_CS6")] <- "Am_CS6"
IDs[which(mammal.combined$IDs=="Am_CS7")] <- "Am_CS7"
IDs[which(mammal.combined$IDs=="Am")] <- "Am"
IDs[which(mammal.combined$IDs=="Amnoid_bead")] <- "Amnoid_bead"
IDs[which(mammal.combined$IDs=="SB43_MEF")] <- "SB43_MEF"
IDs[which(mammal.combined$IDs=="BMP_MEF")] <- "BMP_MEF"
IDs[which(mammal.combined$IDs=="ActA_MEF")] <- "ActA_MEF"
IDs[which(mammal.combined$IDs=="FGF_noMEF")] <- "FGF_noMEF"
IDs[which(mammal.combined$IDs=="ActA_noMEF")] <- "ActA_noMEF"
IDs[which(mammal.combined$IDs=="BMP_noMEF")] <- "BMP_noMEF"

Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_amnion_only2.pdf",sep=""),width = 16, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type3_amnion_only2.pdf",sep=""),width = 16, height = 8)

IDs <- as.character(mammal.combined$IDs)
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$IDs=="EmDisc_CS5")] <- "EmDisc_CS5"
IDs[which(mammal.combined$IDs=="EmDisc_CS6")] <- "EmDisc_CS6"
IDs[which(mammal.combined$IDs=="EmDisc_CS7")] <- "EmDisc_CS7"
IDs[which(mammal.combined$IDs=="EmD")] <- "EmD"
IDs[which(mammal.combined$IDs=="EmDisc")] <- "EmDisc"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type3_emdisc_only.pdf",sep=""),width = 16, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, split.by = "Dataset",) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type3_emdisc_only.pdf",sep=""),width = 16, height = 8)


DefaultAssay(mammal.combined) <- "integrated"
AvExp <- AverageExpression(mammal.combined)

order1 <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS6","Am_CS7")
order2 <- c("EmD","EmDisc","Am","Amnoid_bead","CHIR_MEF","ActA_MEF","ActA_noMEF","FGF_noMEF","BMP_MEF","BMP_noMEF","SB43_MEF")
C1 <- cor( log(as.matrix(AvExp$integrated[,order2]) +1), log(as.matrix(AvExp$integrated[,order1]) +1) )
C2 <- cor( log(as.matrix(AvExp$integrated[order2])+1) )

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), breaks = mat_breaks, gaps_row = c(2,4), gaps_col = c(3), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoBulkCrossCorrelation3",".pdf",sep=""),width=3.7,height=4)

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.5, 0.95, length.out = 60)
pheatmap(C2, color =  redblue1(60), breaks = mat_breaks, border_color = NA, gaps_row = c(2,4), gaps_col = c(2,4), display_numbers = round(C2,2), fontsize_number = 5,  cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/PseudoBulkInvitroCorrelation3",".pdf",sep=""),width=4.2,height=4)

C1 <- cor( log(as.matrix(AvExp$RNA[,order2]) +1), log(as.matrix(AvExp$RNA[,order1]) +1) )
C2 <- cor( log(as.matrix(AvExp$RNA[order2])+1) )

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.6, 0.9, length.out = 60)
pheatmap(C1, color =  redblue1(60), breaks = mat_breaks, gaps_row = c(2,4), gaps_col = c(3), display_numbers = round(C1,2), fontsize_number = 5, border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/PseudoRNABulkCrossCorrelation3",".pdf",sep=""),width=3.7,height=4)





marmoset_dataInVivo2 <- subset(marmoset_dataInVivo,idents = c("Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS5","VE_CS5","VE_CS6","SYS_CS5","SYS_CS7","SYS_CS6","ExMes_CS5","ExMes_CS6","ExMes_CS7"),invert=TRUE)
#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo2, marmoset_dataInVitro2), dims = 1:20, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type4_Lab","_20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type4_Lab","20dim_5k_split.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type4_Lab","_20dim_5k_merge.pdf",sep=""),width = 16, height = 16)

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type4_Lab","_20dim_5k_merge.pdf",sep=""),width = 16, height = 16)



#Grey plots too


mammal.combined2 <- FindClusters(mammal.combined2, resolution = 0.5)
mammal.combined2$Cl1 <- Idents(mammal.combined2)









DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_splitembryo_genomescomp.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_splitembryo_genomescomp.pdf",sep=""),width = 16, height = 8)

DimPlot(mammal.combined, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_genomescomp.pdf",sep=""),width = 16, height = 16)




types <- Idents(mammal.combined)

Idents(mammal.combined) <- droplevels(Idents(mammal.combined))
cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_splitembryo.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined,cols = coluse,  pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_splitembryo.pdf",sep=""),width = 42, height = 8)



dsakldsalkdsandklsadanlsadandl

mammal.combined$uID2 <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$species
D0 <- subset(mammal.combined,idents="PreImp")
D1 <- subset(mammal.combined,idents="E15B")
D2 <- subset(mammal.combined,idents="E15C")
D3 <- subset(mammal.combined,idents="E25A")

#D3 <- subset(mammal.combined,idents=c("E25A","E25A","E25C"))

#Rename to keep variable convention for the join species modelling
mammal.anchors3 <- FindIntegrationAnchors(object.list = list(D1,D2,D3), dims = 1:20, anchor.features = 4000, k.filter=30)
mammal.combined3 <- IntegrateData(anchorset = mammal.anchors3, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined3) <- "integrated"
mammal.combined3 <- ScaleData(mammal.combined3, verbose = FALSE)
mammal.combined3 <- RunPCA(mammal.combined3, npcs = 20, verbose = FALSE)
mammal.combined3 <- RunUMAP(mammal.combined3, reduction = "pca", dims = 1:20)
mammal.combined3 <- RunTSNE(mammal.combined3, reduction = "pca", dims = 1:20)
mammal.combined3 <- FindNeighbors(mammal.combined3, reduction = "pca", dims = 1:20)

Idents(mammal.combined3) <- mammal.combined3$uID2

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined3, cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_splitembryoalign.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined3,cols = coluse,  pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_splitembryoalign.pdf",sep=""),width = 42, height = 8)


#Rename to keep variable convention for the join species modelling
mammal.anchors4 <- FindIntegrationAnchors(object.list = list(D0,D1,D2,D3), dims = 1:20, anchor.features = 4000, k.filter=30)
mammal.combined4 <- IntegrateData(anchorset = mammal.anchors4, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined4) <- "integrated"
mammal.combined4 <- ScaleData(mammal.combined4, verbose = FALSE)
mammal.combined4 <- RunPCA(mammal.combined4, npcs = 20, verbose = FALSE)
mammal.combined4 <- RunUMAP(mammal.combined4, reduction = "pca", dims = 1:20)
mammal.combined4 <- RunTSNE(mammal.combined4, reduction = "pca", dims = 1:20)
mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20)

Idents(mammal.combined4) <- mammal.combined4$uID2

colind <- integer( length( levels(Idents(mammal.combined4)) )  )
for (i in 1:length( levels(Idents(mammal.combined4)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined4))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined4, cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_splitembryoalign2.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined4,cols = coluse,  pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_splitembryoalign2.pdf",sep=""),width = 42, height = 8)

mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20, k.param = 10)
KNN_K10 <- mammal.combined4@graphs$integrated_nn
SNN_K10 <- mammal.combined4@graphs$integrated_snn
mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20, k.param = 5)
KNN_K5 <- mammal.combined4@graphs$integrated_nn
SNN_K5 <- mammal.combined4@graphs$integrated_snn
mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20, k.param = 20)
KNN_K20 <- mammal.combined4@graphs$integrated_nn
SNN_K20 <- mammal.combined4@graphs$integrated_snn



inds1 <- which(mammal.combined4$species=="PreImp")
inds2 <- which(mammal.combined4$species=="E15B")
inds3 <- which(mammal.combined4$species=="E15C")
inds4 <- which(mammal.combined4$species=="E25A")

mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20, k.param = 20)
KNN_K20 <- mammal.combined4@graphs$integrated_nn
SNN_K20 <- mammal.combined4@graphs$integrated_snn
L1 <- SNN_K20[inds2,inds1]
L2 <- SNN_K20[inds3,inds2]
L3 <- SNN_K20[inds4,inds3]
L1max <- apply(as.data.frame(L1), 1, max) 
L2max <- apply(as.data.frame(L2), 1, max) 
L3max <- apply(as.data.frame(L3), 1, max) 
L1max_ <- apply(as.data.frame(L1), 1, which.max)
L2max_ <- apply(as.data.frame(L2), 1, which.max)
L3max_ <- apply(as.data.frame(L3), 1, which.max)
IDs1 <- as.character(Idents(mammal.combined4)[inds1])
IDs2 <- as.character(Idents(mammal.combined4)[inds2])
IDs3 <- as.character(Idents(mammal.combined4)[inds3])
IDs4 <- as.character(Idents(mammal.combined4)[inds4])

from1 <- IDs1[L1max_]
to1 <- IDs2
from2 <- IDs2[L2max_]
to2 <- IDs3
from3 <- IDs3[L3max_]
to3 <- IDs4
F <- c(from1[which(L1max>0)],from2[which(L2max>0)],from3[which(L3max>0)])
T <- c(to1[which(L1max>0)],to2[which(L2max>0)],to3[which(L3max>0)])
IDs1_ <- colnames(mammal.combined4)[inds1]
IDs2_ <- colnames(mammal.combined4)[inds2]
IDs3_ <- colnames(mammal.combined4)[inds3]
IDs4_ <- colnames(mammal.combined4)[inds4]
from1 <- IDs1_[L1max_]
to1 <- IDs2_
from2 <- IDs2_[L2max_]
to2 <- IDs3_
from3 <- IDs3_[L3max_]
to3 <- IDs4_
F2 <- c(from1[which(L1max>0)],from2[which(L2max>0)],from3[which(L3max>0)])
T2 <- c(to1[which(L1max>0)],to2[which(L2max>0)],to3[which(L3max>0)])
FT <- data.frame(x=F,y=T,x2=F2,y2=T2)
write.csv(FT, file=paste(saveext,"/FT_20.csv",sep=""))



mammal.combined4 <- FindNeighbors(mammal.combined4, reduction = "pca", dims = 1:20, k.param = 50)
KNN_K20 <- mammal.combined4@graphs$integrated_nn
SNN_K20 <- mammal.combined4@graphs$integrated_snn
inds1 <- which(mammal.combined4$species=="PreImp")
inds2 <- which(mammal.combined4$species=="E15B")
inds3 <- which(mammal.combined4$species=="E15C")
inds4 <- which(mammal.combined4$species=="E25A")
L1 <- SNN_K20[inds2,inds1]
L2 <- SNN_K20[inds3,inds2]
L3 <- SNN_K20[inds4,inds3]
L1max <- apply(as.data.frame(L1), 1, max) 
L2max <- apply(as.data.frame(L2), 1, max) 
L3max <- apply(as.data.frame(L3), 1, max) 
L1max_ <- apply(as.data.frame(L1), 1, which.max)
L2max_ <- apply(as.data.frame(L2), 1, which.max)
L3max_ <- apply(as.data.frame(L3), 1, which.max)
IDs1 <- as.character(Idents(mammal.combined4)[inds1])
IDs2 <- as.character(Idents(mammal.combined4)[inds2])
IDs3 <- as.character(Idents(mammal.combined4)[inds3])
IDs4 <- as.character(Idents(mammal.combined4)[inds4])
from1 <- IDs1[L1max_]
to1 <- IDs2
from2 <- IDs2[L2max_]
to2 <- IDs3
from3 <- IDs3[L3max_]
to3 <- IDs4
F <- c(from1[which(L1max>0)],from2[which(L2max>0)],from3[which(L3max>0)])
T <- c(to1[which(L1max>0)],to2[which(L2max>0)],to3[which(L3max>0)])
IDs1_ <- colnames(mammal.combined4)[inds1]
IDs2_ <- colnames(mammal.combined4)[inds2]
IDs3_ <- colnames(mammal.combined4)[inds3]
IDs4_ <- colnames(mammal.combined4)[inds4]
from1 <- IDs1_[L1max_]
to1 <- IDs2_
from2 <- IDs2_[L2max_]
to2 <- IDs3_
from3 <- IDs3_[L3max_]
to3 <- IDs4_
F2 <- c(from1[which(L1max>0)],from2[which(L2max>0)],from3[which(L3max>0)])
T2 <- c(to1[which(L1max>0)],to2[which(L2max>0)],to3[which(L3max>0)])
FT <- data.frame(x=F,y=T,x2=F2,y2=T2)
write.csv(as.data.frame(FT), file=paste(saveext,"/FT_50.csv",sep=""))



#choose N nearest neightbours, there are x 
#phyper(x, m, n, k)
#x == overalp
#m	the number of white balls in the urn. == No. of objective sample
#n	the number of black balls in the urn. == No. match sample
#k	the number of balls drawn from the urn. == NNs

u1 <- unique(IDs1)
Len1 <- length(IDs1)
P1 <- matrix(0,dim( L1 )[1],4)
for (i in 1:dim( L1 )[1] ) {
  n1 <- length(which(IDs1=="ICM_CS3"))
  n2 <- length(which(IDs1=="Epi_CS3"))
  n3 <- length(which(IDs1=="Hyp_CS3"))
  n4 <- length(which(IDs1=="Tb_CS3"))  
  xx <- IDs1[which(L1[i,]>0)]
  x1 <- length(which(xx=="ICM_CS3"))
  x2 <- length(which(xx=="Epi_CS3"))
  x3 <- length(which(xx=="Hyp_CS3"))
  x4 <- length(which(xx=="Tb_CS3"))  
  P1[i,1] <- 1-phyper(x1-1, n1, Len1-n1, 50)
  P1[i,2] <- 1-phyper(x2-1, n2, Len1-n2, 50)
  P1[i,3] <- 1-phyper(x3-1, n3, Len1-n3, 50)
  P1[i,4] <- 1-phyper(x4-1, n4, Len1-n4, 50)  
}

u1 <- unique(IDs2)
Len1 <- length(IDs2)
P2 <- matrix(0,dim( L2 )[1],7)
for (i in 1:dim( L2 )[1] ) {
  n1 <- length(which(IDs2=="Tb_CS5"))
  n2 <- length(which(IDs2=="PGC_CS5"))
  n3 <- length(which(IDs2=="ExMes_CS5"))
  n4 <- length(which(IDs2=="VE_CS5"))  
  n5 <- length(which(IDs2=="EmDisc_CS5"))  
  n6 <- length(which(IDs2=="SYS_CS5"))  
  n7 <- length(which(IDs2=="Am_CS5"))    
  xx <- IDs2[which(L2[i,]>0)]
  x1 <- length(which(xx=="Tb_CS5"))
  x2 <- length(which(xx=="PGC_CS5"))
  x3 <- length(which(xx=="ExMes_CS5"))
  x4 <- length(which(xx=="VE_CS5"))  
  x5 <- length(which(xx=="EmDisc_CS5"))  
  x6 <- length(which(xx=="SYS_CS5"))  
  x7 <- length(which(xx=="Am_CS5"))    
  P2[i,1] <- 1-phyper(x1-1, n1, Len1-n1, 50)
  P2[i,2] <- 1-phyper(x2-1, n2, Len1-n2, 50)
  P2[i,3] <- 1-phyper(x3-1, n3, Len1-n3, 50)
  P2[i,4] <- 1-phyper(x4-1, n4, Len1-n4, 50)  
  P2[i,5] <- 1-phyper(x5-1, n5, Len1-n5, 50)  
  P2[i,6] <- 1-phyper(x6-1, n6, Len1-n6, 50)  
  P2[i,7] <- 1-phyper(x7-1, n7, Len1-n7, 50)    
}



u1 <- unique(IDs3)
Len1 <- length(IDs3)
P3 <- matrix(0,dim( L3 )[1],7)
for (i in 1:dim( L3 )[1] ) {
  n1 <- length(which(IDs3=="Tb_CS6"))
  n2 <- length(which(IDs3=="PGC_CS6"))
  n3 <- length(which(IDs3=="ExMes_CS6"))
  n4 <- length(which(IDs3=="VE_CS6"))  
  n5 <- length(which(IDs3=="EmDisc_CS6"))  
  n6 <- length(which(IDs3=="SYS_CS6"))  
  n7 <- length(which(IDs3=="Am_CS6"))    
  xx <- IDs3[which(L3[i,]>0)]
  x1 <- length(which(xx=="Tb_CS6"))
  x2 <- length(which(xx=="PGC_CS6"))
  x3 <- length(which(xx=="ExMes_CS6"))
  x4 <- length(which(xx=="VE_CS6"))  
  x5 <- length(which(xx=="EmDisc_CS6"))  
  x6 <- length(which(xx=="SYS_CS6"))  
  x7 <- length(which(xx=="Am_CS6"))    
  P3[i,1] <- 1-phyper(x1-1, n1, Len1-n1, 50)
  P3[i,2] <- 1-phyper(x2-1, n2, Len1-n2, 50)
  P3[i,3] <- 1-phyper(x3-1, n3, Len1-n3, 50)
  P3[i,4] <- 1-phyper(x4-1, n4, Len1-n4, 50)  
  P3[i,5] <- 1-phyper(x5-1, n5, Len1-n5, 50)  
  P3[i,6] <- 1-phyper(x6-1, n6, Len1-n6, 50)  
  P3[i,7] <- 1-phyper(x7-1, n7, Len1-n7, 50)    
}

aID <- as.character(mammal.combined$AltID)
aID[which(mammal.combined$species2=="1) Marmoset")] <- as.character(mammal.combined$CID[which(mammal.combined$species2=="1) Marmoset")])
Idents(mammal.combined) <- aID #mammal.combined$AltID

DimPlot(mammal.combined, pt.size = 4, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_altID.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_altID.pdf",sep=""),width = 42, height = 8)


DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_altIDmerge.pdf",sep=""),width = 25, height = 15)

DimPlot(mammal.combined, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_altIDmerge.pdf",sep=""),width = 25, height = 15)


mammal.combined$CID <- Idents(mammal.combined)



IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(Idents(mammal.combined)=="EmDisc_CS5")] <- "EmDisc_CS5"
IDs[which(Idents(mammal.combined)=="EmDisc_CS6")] <- "EmDisc_CS6"
IDs[which(Idents(mammal.combined)=="EmDisc_CS7")] <- "EmDisc_CS7"
IDs[which(Idents(mammal.combined)=="EmDisc_CS6/7")] <- "EmDisc_CS6/7"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS5")] <- "EmDiscPS_CS5"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS6")] <- "EmDiscPS_CS6"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS6/7")] <- "EmDiscPS_CS6/7"

Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_EmDisc.pdf",sep=""),width = 42, height = 8)




IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(mammal.combined$CID=="ExMes_CS5")] <- "ExMes_CS5"
IDs[which(mammal.combined$CID=="ExMes_CS6")] <- "ExMes_CS6"
IDs[which(mammal.combined$CID=="ExMes_CS7")] <- "ExMes_CS7"
IDs[which(mammal.combined$CID=="ExMes_CS6/7")] <- "ExMes_CS6/7"

Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_ExMes.pdf",sep=""),width = 42, height = 8)





IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(mammal.combined$CID=="VE_CS5")] <- "VE_CS5"
IDs[which(mammal.combined$CID=="VE_CS6")] <- "VE_CS6"
IDs[which(mammal.combined$CID=="VE_CS7")] <- "VE_CS7"
IDs[which(mammal.combined$CID=="VE_CS6/7")] <- "VE_CS6/7"

IDs[which(mammal.combined$CID=="SYS_CS5")] <- "SYS_CS5"
IDs[which(mammal.combined$CID=="SYS_CS6")] <- "SYS_CS6"
IDs[which(mammal.combined$CID=="SYS_CS7")] <- "SYS_CS7"
IDs[which(mammal.combined$CID=="SYS_CS6/7")] <- "SYS_CS6/7"


Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_VESYS.pdf",sep=""),width = 42, height = 8)




IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(mammal.combined$CID=="Tb_CS5")] <- "Tb_CS5"
IDs[which(mammal.combined$CID=="Tb_CS6")] <- "Tb_CS6"
IDs[which(mammal.combined$CID=="Tb_CS7")] <- "Tb_CS7"
IDs[which(mammal.combined$CID=="Tb_CS6/7")] <- "Tb_CS6/7"


Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_Tb.pdf",sep=""),width = 42, height = 8)




IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(mammal.combined$CID=="Am_CS5")] <- "Am_CS5"
IDs[which(mammal.combined$CID=="Am_CS6")] <- "Am_CS6"
IDs[which(mammal.combined$CID=="Am_CS7")] <- "Am_CS7"
IDs[which(mammal.combined$CID=="Am_CS6/7")] <- "Am_CS6/7"


Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_Am.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "ISL1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/ISL1.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "NANOS3", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/NANOS3.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX17", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/SOX17.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "VTCN1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VTCN1.pdf",sep=""),width = 42, height = 8)


FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2A", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/TFAP2A.pdf",sep=""),width = 42, height = 8)



FeaturePlot(mammal.combined,  reduction = "pca", features = "PDGFRA", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PDGFRA.pdf",sep=""),width = 42, height = 8)



FeaturePlot(mammal.combined,  reduction = "umap", features = "NANOS3", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/NANOS3_UMAP.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX17", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/SOX17_UMAP.pdf",sep=""),width = 42, height = 8)

IDD <- as.character(Idents(marmoset_data1))
IDD[which(IDD=="Tb_CS5")] <- "Tb"
IDD[which(IDD=="Tb_CS6")] <- "Tb"
IDD[which(IDD=="Tb_CS7")] <- "Tb"
IDD[which(IDD=="Am_CS5")] <- "Am"
IDD[which(IDD=="Am_CS6")] <- "Am"
IDD[which(IDD=="Am_CS7")] <- "Am"
IDD[which(IDD=="ExMes_CS5")] <- "ExMes"
IDD[which(IDD=="ExMes_CS6")] <- "ExMes"
IDD[which(IDD=="ExMes_CS7")] <- "ExMes"
Idents(marmoset_data1) <- IDD

for (i in 1:300) {
  VlnPlot(marmoset_data1,  idents=c("Am","Tb"), features = rownames(Amn1)[i], pt.size = 4)
  ggsave(filename=paste(saveext,"Markers/VlAmGenes_",rownames(Amn1)[i],".pdf",sep=""),width = 8, height = 8)
  
  VlnPlot(marmoset_data1,  idents=c("ExMes","Tb"), features = rownames(Mes1)[i], pt.size = 4)
  ggsave(filename=paste(saveext,"Markers/VlExMesGenes_",rownames(Mes1)[i],".pdf",sep=""),width = 8, height = 8)  
  
  
#  FeaturePlot(mammal.combined,  reduction = "pca", features = rownames(Amn1)[i], combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
#  ggsave(filename=paste(saveext,"Markers/AmGenes_",rownames(Amn1)[i],".pdf",sep=""),width = 42, height = 8)
  
#  FeaturePlot(mammal.combined,  reduction = "pca", features = rownames(Mes1)[i], combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
#  ggsave(filename=paste(saveext,"Markers/ExMesGenes_",rownames(Mes1)[i],".pdf",sep=""),width = 42, height = 8)  
  
}



FeaturePlot(mammal.combined,  reduction = "umap", features = "KDR", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/KDR_UMAP.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "KDR", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/KDR_PCA.pdf",sep=""),width = 42, height = 8)


#Subset for refined annotation

FeaturePlot(mammal.combined,  reduction = "umap", features = "HLA-G", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/HLAG_UMAP.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "HLA-G", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/HLAG_PCA.pdf",sep=""),width = 42, height = 8)


Tb1 <- FindMarkers(marmoset_data1,ident.2 = c("ExMes"),ident.1 = c("Tb"),only.pos = TRUE,test.use="MAST")
Tb2 <- FindMarkers(marmoset_data1,ident.2 = c("Am"),ident.1 = c("Tb"),only.pos = TRUE,test.use="MAST")

list1 <- intersect(rownames(Tb1),rownames(Tb2))
for (i in 1:length(list1)) {
  VlnPlot(marmoset_data1,  idents=c("Am","ExMes","Tb"), features = list1[i], pt.size = 4)
  ggsave(filename=paste(saveext,"Markers/VlTbGenes_",list1[i],".pdf",sep=""),width = 8, height = 8)
  
#  VlnPlot(marmoset_data1,  idents=c("ExMes","Tb"), features = list1[i], pt.size = 4)
#  ggsave(filename=paste(saveext,"Markers/VlExMesGenes_",list1[i],".pdf",sep=""),width = 8, height = 8)  
  
  
   FeaturePlot(mammal.combined,  reduction = "pca", features = list1[i], combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
    ggsave(filename=paste(saveext,"Markers/TbGenes_",list1[i],".pdf",sep=""),width = 42, height = 8)
  
  #  FeaturePlot(mammal.combined,  reduction = "pca", features = rownames(Mes1)[i], combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
  #  ggsave(filename=paste(saveext,"Markers/ExMesGenes_",rownames(Mes1)[i],".pdf",sep=""),width = 42, height = 8)  
  
}






#Rename to keep variable convention for the join species modelling
marmoset_data3 <- subset(marmoset_data,idents=c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","PGC_CS5","PGC_CS6"))
cynomolgous_data3 <- subset(cynomolgous_data,idents=c("EmDisc_CS5","EmDisc_CS6/7","EmDiscPS_CS6/7","ExMes_CS5","ExMes_CS6/7","PGC_CS5","PGC_CS6/7"))
human_dataA3 <- subset(human_dataA,idents=c("EmDisc_CS5","EmDisc_CS6","EmDiscPS_CS6"))
mammal.anchors2 <- FindIntegrationAnchors(object.list = list(marmoset_data3, cynomolgous_data3, human_dataA3), dims = 1:20, anchor.features = 4000, k.filter = 30)
mammal.combined2 <- IntegrateData(anchorset = mammal.anchors2, dims = 1:20)
DefaultAssay(mammal.combined2) <- "integrated"
mammal.combined2 <- ScaleData(mammal.combined2, verbose = FALSE)
mammal.combined2 <- RunPCA(mammal.combined2, npcs = 20, verbose = FALSE)
mammal.combined2 <- RunUMAP(mammal.combined2, reduction = "pca", dims = 1:20)
mammal.combined2 <- RunTSNE(mammal.combined2, reduction = "pca", dims = 1:20)
mammal.combined2 <- FindNeighbors(mammal.combined2, reduction = "pca", dims = 1:20)


cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined2,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_EmbryonicLineages.pdf",sep=""),width = 24, height = 8)

DimPlot(mammal.combined2,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4_EmbryonicLineages.pdf",sep=""),width = 24, height = 8)

FeaturePlot(mammal.combined2,  reduction = "umap", features = "VTCN1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"Markers/Embryonic_UMAP_VTCN1.pdf",sep=""),width = 42, height = 8)





FeaturePlot(mammal.combined2,  reduction = "pca", features = "ISL1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryonic_PCA_ISL1.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "NANOS3", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryonic_PCA_NANOS3.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "SOX17", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryonic_PCA_SOX17.pdf",sep=""),width = 42, height = 8)

FeaturePlot(mammal.combined2,  reduction = "pca", features = "VTCN1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryonic_PCA_VTCN1.pdf",sep=""),width = 42, height = 8)


FeaturePlot(mammal.combined2,  reduction = "pca", features = "TFAP2A", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species2", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryonic_PCA_TFAP2A.pdf",sep=""),width = 42, height = 8)






mammal.combined2 <- FindClusters(mammal.combined2, resolution = 0.5)
mammal.combined2$Cl1 <- Idents(mammal.combined2)

DimPlot(mammal.combined2, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/PCA_Type_Lab","_4_Cl.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined2, pt.size = 4, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_Type_Lab","_4_Cl.pdf",sep=""),width = 42, height = 8)


Cl <- as.character(Idents(mammal.combined2))
Cl[which(Cl=="2")] <- "EmDiscCl"
Cl[which(Cl=="3")] <- "EmDiscCl"
Cl[which(Cl=="5")] <- "EmDiscCl"

Cl[which(Cl=="1")] <- "AmCl"
Cl[which(Cl=="0")] <- "ExMesCl"
Cl[which(Cl=="6")] <- "PGCCl"
Cl[which(Cl=="4")] <- "PSCCl"


Sp <- mammal.combined2$species2
Sp[which(Sp=="1) Marmoset")] <- "m"
Sp[which(Sp=="2) Cynomolgous")] <- "c"
Sp[which(Sp=="3) Human (in vitro)")] <- "h"

Idents(mammal.combined2) <- paste(Sp,Cl,sep="")

genes <- intersect(intersect(rownames(human_dataA3),rownames(marmoset_data3)),rownames(cynomolgous_data3))

genes2 <- intersect(rownames(mammal.combined2),genes)
#avg.Cl0 <- log1p(AverageExpression(mammal.combined2, features = genes, verbose = FALSE)$RNA)
avg.Cl0 <- log1p(AverageExpression(mammal.combined2, features = genes, verbose = FALSE, assays = "RNA"))


#avg.Cl0 <- log1p(AverageExpression(mammal.combined2,  verbose = FALSE)$RNA)
avg.Cl0$gene <- rownames(avg.Cl0)
#avg.Cl0 <- avg.Cl0[genes,]


#87,21 testes
#74,82,ovaries

TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1


#Do all PGCs
Cl1 <- FindMarkers(mammal.combined2, ident.2 = "mEmDiscCl", ident.1 = "cEmDiscCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl1) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mEmDiscCl,cEmDiscCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_marm_cyno.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "mEmDiscCl", ident.1 = "hEmDiscCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mEmDiscCl,hEmDiscCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_marm_human.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "cEmDiscCl", ident.1 = "hEmDiscCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(cEmDiscCl,hEmDiscCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"EmDisc_cyno_human.pdf",sep=""),width = 20, height = 20, plot = p1)


#Do all PGCs
Cl1 <- FindMarkers(mammal.combined2, ident.2 = "mPSCCl", ident.1 = "cPSCCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl1) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mPSCCl,cPSCCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"PS_marm_cyno.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "mPSCCl", ident.1 = "hPSCCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mPSCCl,hPSCCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"PS_marm_human.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "cPSCCl", ident.1 = "hPSCCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(cPSCCl,hPSCCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"PS_cyno_human.pdf",sep=""),width = 20, height = 20, plot = p1)




#Do all PGCs
Cl1 <- FindMarkers(mammal.combined2, ident.2 = "mExMesCl", ident.1 = "cExMesCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl1) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mExMesCl,cExMesCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"ExMes_marm_cyno.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "mExMesCl", ident.1 = "hExMesCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mExMesCl,hExMesCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"ExMes_marm_human.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "cExMesCl", ident.1 = "hExMesCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(cExMesCl,hExMesCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"ExMes_cyno_human.pdf",sep=""),width = 20, height = 20, plot = p1)



#Do all PGCs
Cl1 <- FindMarkers(mammal.combined2, ident.2 = "mAmCl", ident.1 = "cAmCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl1) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mAmCl,cAmCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_marm_cyno.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "mAmCl", ident.1 = "hAmCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(mAmCl,hAmCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_marm_human.pdf",sep=""),width = 20, height = 20, plot = p1)

Cl2 <- FindMarkers(mammal.combined2, ident.2 = "cAmCl", ident.1 = "hAmCl", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
genes.to.label = rownames(Cl2) #intersect(rownames(Cl1),TF)
p1 <- ggplot(avg.Cl0, aes(cAmCl,hAmCl)) + geom_point(color="lightgrey") + theme_classic()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"Am_cyno_human.pdf",sep=""),width = 20, height = 20, plot = p1)




#Bored now
dasdljksadsalkdslkdsalkdsaldlksamdlksalkdaslkdnlasnd
















#And plot these
#DimPlot(mammal.combined, reduction = "pca", pt.shape = "cell.orig", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/PCA_Type",".pdf",sep=""),width = 26, height = 8)
#DimPlot(mammal.combined, reduction = "umap", pt.shape = "cell.orig", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 26, height = 8)
#DimPlot(mammal.combined, reduction = "tsne", pt.shape = "cell.orig", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 26, height = 8)

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))
dir.create(paste(saveext,"/Docs/",sep=""))
dir.create(paste(saveext,"/Monocle/",sep=""))
dir.create(paste(saveext,"/GO/",sep=""))

#56E600
#48BF00
#56E600	Zygote_CS1
#00BF30	8-cell_CS2
##48BF00	4-cell_CS2
#009926	Mor_CS3
#00E6E6	earlyICM_CS3
#00BFBF	Epi_CS3
#0039E6	EmDisc_CS5
#0233BF	EmDisc_CS7
#1A0873	Am_CS7
#4576E7	Am_CS5
#921FE6	Tb_CS5
#7813BF	Tb_CS7
#BF0489	Tb_CS3
#F04C04	VE_CS5
#E68600	SYS_protrusion_CS5
#BF7104	SYS_CS7
#E6B500	Hyp_CS3
#BF9600	ExMes_CS5
#967700	ExMes_CS7
#E6E600	PGC_CS5

#Load marmoset data key
TF<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Leaving_package/Dimensionality\ reduction\ techniques\ smart-seq2\ -\ Boroviaklab\ data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)

#List of cell cycle markers loaded in with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Load marmoset data key
cyBS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cynomolgous/cyKey.csv",sep=",", header = T, row.names=1)
BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/Keycorrect_CPAll2.csv",sep=",",header = T, row.names=1)
TF<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Leaving_package/Dimensionality\ reduction\ techniques\ smart-seq2\ -\ Boroviaklab\ data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)

#Load monkey
raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cynomolgous/featurecountsNakamura3.csv",sep=",",header = T, row.names=1)
cynomolgous_data <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
#cyisinvit <- cyBS$AllNoES
cylabs <- cyBS$LABEL_4 
#cylabs <- cylabs[which(cyisinvit>0)]
cylabs2 <- cyBS$LABEL_2
#cylabs2 <- cylabs2[which(cyisinvit>0)]
raw_counts1<- raw_counts1[,which(cyisinvit>0)]
cynomolgous_data <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
cynomolgous_data$species <- "2) Cynomolgous"
cynomolgous_data$divergence1 <- "Cyno"
cynomolgous_data <- subset(cynomolgous_data, subset = nFeature_RNA > 0)
cynomolgous_data$divergence2 <- "primate"
cynomolgous_data <- NormalizeData(cynomolgous_data, verbose = FALSE)
cynomolgous_data <- FindVariableFeatures(cynomolgous_data, selection.method = "vst", nfeatures = 20000)
cynomolgous_data <- ScaleData(cynomolgous_data, verbose = FALSE)
Idents(cynomolgous_data) <- cylabs
cynomolgous_data$AltID <- cylabs2

#Load the marmoset data ... to get ordering
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
isinvit <- BS$vsCynoMin * BS$QC #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$GSPostImpNoMat2 #* BS[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$E25
labs <- BS$Annotation2 #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$Location
labs <- labs[which(isinvit>0)]

#Load raw counts
raw_counts3<- raw_counts3[,which(isinvit>0)]
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "1) Marmoset"
marmoset_data$divergence1 <- "Marm"
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)
marmoset_data <- ScaleData(marmoset_data, verbose = FALSE)
Idents(marmoset_data) <- labs
marmoset_data$AltID <- labs

#raw_countsB<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/GSE130114_MF1453.csv",sep=",",header = T, row.names=1)
#cyno_dataB <- CreateSeuratObject(counts = na.omit(raw_countsB), assay = "RNA",min.cells = 0, min.features = 0)
cyivBS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cyinvitLab.csv",sep=",", header = T, row.names=1)

#cyno_dataB$species <- "3) Cynomologous (in vitro)"
#cyno_dataB$divergence1 <- "Cyno_invit"
#cyno_dataB <- subset(cyno_dataB, subset = nFeature_RNA > 0)
#cyno_dataB <- NormalizeData(cyno_dataB, verbose = FALSE)
#cyno_dataB <- FindVariableFeatures(cyno_dataB, selection.method = "vst", nfeatures = 20000)
#cyno_dataB <- ScaleData(cyno_dataB, verbose = FALSE)
#Idents(cyno_dataB) <- cyivBS$Type


hBS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/Human_invitro2.csv",sep=",",header = T, row.names=1)
#hisinvit <- which(hBS$Type5 != "MISSING" & hBS$Type5 != "Tb_CS4" &  hBS$Type5 !=  "VE_CS4" & hBS$Type5 !=  "Epi_CS4" )
hisinvit <- which(hBS$Type5 != "MISSING" ) # & hBS$Type5 != "Tb_CS4" &  hBS$Type5 !=  "VE_CS4" & hBS$Type5 !=  "Epi_CS4" )
hlabs <- hBS$Type5[hisinvit]
hlabs2 <- hBS$Type2[hisinvit]
raw_countsA<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/GSE136447_555-samples-fpkm.csv",sep=",",header = T, row.names=1)
human_dataA <- CreateSeuratObject(counts = na.omit(raw_countsA[,hisinvit]), assay = "RNA",min.cells = 0, min.features = 0)
human_dataA$species <- "4) Human (in vitro)"
human_dataA$divergence1 <- "Human_invit"
human_dataA <- subset(human_dataA, subset = nFeature_RNA > 0)
human_dataA <- NormalizeData(human_dataA, verbose = FALSE)
human_dataA <- FindVariableFeatures(human_dataA, selection.method = "vst", nfeatures = 20000)
#Idents(human_dataA) <- hBS$Type5
#human_dataA <- subset(human_dataA, idents = c("MISSING"), invert = TRUE)
human_dataA <- ScaleData(human_dataA, verbose = FALSE)
Idents(human_dataA) <- hlabs
human_dataA$AltID <- hlabs2


raw_countsB<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/GSE130114_MF1453.csv",sep=",",header = T, row.names=1)
isinvit <- cyivBS$vsCyno
labs <- cyivBS$Type
labs <- labs[which(isinvit>0)]
raw_countsB<- raw_countsB[,which(isinvit>0)]
Cyno_dataB <- CreateSeuratObject(counts = na.omit(raw_countsB), assay = "RNA",min.cells = 0, min.features = 0)
Cyno_dataB$species <- "3) Cynomologous (in vitro)"
Cyno_dataB$divergence1 <- "Cyno_invit"
Cyno_dataB <- subset(Cyno_dataB, subset = nFeature_RNA > 0)
Cyno_dataB <- NormalizeData(Cyno_dataB, verbose = FALSE)
Cyno_dataB <- FindVariableFeatures(Cyno_dataB, selection.method = "vst", nfeatures = 20000)
Cyno_dataB <- ScaleData(Cyno_dataB, verbose = FALSE)
Idents(Cyno_dataB) <- labs #cyivBS$Type 

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_data, cynomolgous_data, human_dataA), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20, perplexity = 40)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

types <- Idents(mammal.combined)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
#mammal.combined$Cl1 <- Idents(mammal.combined)
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.25)
#mammal.combined$Cl2 <- Idents(mammal.combined)
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.75)
#mammal.combined$Cl3 <- Idents(mammal.combined)
#Idents(mammal.combined) <- types
#mammal.combined$IDs <- Idents(mammal.combined)
         

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]



DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined, pt.size = 4, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4.pdf",sep=""),width = 42, height = 8)

mammal.combined$CID <- Idents(mammal.combined)
IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(Idents(mammal.combined)=="EmDisc_CS5")] <- "EmDisc_CS5"
IDs[which(Idents(mammal.combined)=="EmDisc_CS6")] <- "EmDisc_CS6"
IDs[which(Idents(mammal.combined)=="EmDisc_CS7")] <- "EmDisc_CS7"
IDs[which(Idents(mammal.combined)=="EmDisc_CS6/7")] <- "EmDisc_CS6/7"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS5")] <- "EmDiscPS_CS5"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS6")] <- "EmDiscPS_CS6"
IDs[which(Idents(mammal.combined)=="EmDiscPS_CS6/7")] <- "EmDiscPS_CS6/7"

Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_EmDisc.pdf",sep=""),width = 42, height = 8)




mammal.combined$CID <- Idents(mammal.combined)
IDs <- as.character(mammal.combined$CID)
IDs[1:length(IDs)] <- "other"
IDs[which(Idents(mammal.combined)=="ExMes_CS5")] <- "ExMes_CS5"
IDs[which(Idents(mammal.combined)=="ExMes_CS6")] <- "ExMes_CS6"
IDs[which(Idents(mammal.combined)=="ExMes_CS7")] <- "ExMes_CS7"
IDs[which(Idents(mammal.combined)=="ExMes_CS6/7")] <- "ExMes_CS6/7"

Idents(mammal.combined) <- IDs

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","PGC_CS5","PGC_CS6","PGC_CS6/7","PGC_E50","EmDiscPS_CS5","EmDiscPS_CS6","EmDiscPS_CS6/7","other")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#E6E600","#E6E600","#E6E600","#E6E600","#000000","#000000","#000000","#d3d3d3")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4_ExMes.pdf",sep=""),width = 42, height = 8)




#cols <- c("#1A0873","#13084D","#0039E6","#0233BF","#0233BF","#0233BF","#00BFBF","#00BFBF","#F04C04","#BF3C04","#BF3C04","#992F03","#E6B500","#00E6E6","#E6E600","#E6E600","#BF7104","#BF7104","#BF7104","#E68600","#BF0489","#BF0489","#921FE6","#7813BF","#7813BF","#7813BF","#E6B500","#BF9600")

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="VE_CS4")] <- 1
cluster_letters <- as.factor(cluster_letters)

#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')

#And plot these
#DimPlot(mammal.combined, reduction = "pca", pt.shape = "cell.orig", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/PCA_Type",".pdf",sep=""),width = 26, height = 8)
#DimPlot(mammal.combined, reduction = "umap", pt.shape = "cell.orig", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 26, height = 8)
#DimPlot(mammal.combined, reduction = "tsne", pt.shape = "cell.orig", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 26, height = 8)


ElbowPlot(mammal.combined)
ggsave(filename=paste(saveext,"/DimRed/Components",".pdf",sep=""),width = 4, height = 4)


DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_12.pdf",sep=""),width = 42, height = 8)


DimPlot(mammal.combined,, cols = coluse, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_8.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 


DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "pca", split.by = "species", label = FALSE, repel = TRUE) + theme(panel.background = element_rect(fill = 'black', colour = 'black'))
ggsave(filename=paste("~/Desktop/PCA_Type_Lab","_6_w.pdf",sep=""),width = 32, height = 8)


DimPlot(mammal.combined, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_TEST2.pdf",sep=""),width = 42, height = 8)


DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "pcap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_12.pdf",sep=""),width = 42, height = 8)


DimPlot(mammal.combined,  cols = coluse, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp2.pdf",sep=""),width = 42, height = 8, useDingbats=FALSE)

split <- SplitObject(mammal.combined, split.by = "species")
AllMarkers <- FindAllMarkers(split$`1) Marmoset`, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


##### ICM_CS3_Human_invit
IntG <- intersect( intersect(rownames(marmoset_data), rownames(cynomolgous_data)), rownames(human_dataA) )
dir.create(paste(saveext,"/EpiCS3/",sep=""))
uGenes1 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Epi_CS3")] )
for (i in 1:length(uGenes1)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes1[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/EpiCS3/Markerscatter_",as.character(uGene1s[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/HypCS3/",sep=""))
uGenes2 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Hyp_CS3")] )
for (i in 1:length(uGenes2)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes2[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/HypCS3/Markerscatter_",as.character(uGenes2[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/TbCS3/",sep=""))
uGenes3 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS3")] )
for (i in 1:length(uGenes3)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes3[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/TbCS3/Markerscatter_",as.character(uGenes3[i]),".pdf",sep=""),width = 35, height = 8)
}


dir.create(paste(saveext,"/EmDisc_CS5/",sep=""))
uGenes4 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="EmDisc_CS5")] )
for (i in 1:length(uGenes4)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes4[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/EmDisc_CS5/Markerscatter_",as.character(uGenes4[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/SYS_CS5/",sep=""))
uGenes5 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="SYS_CS5")] )
for (i in 1:length(uGenes5)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes5[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/SYS_CS5/Markerscatter_",as.character(uGenes5[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/Tb_CS5/",sep=""))
uGenes6 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS5")] )
for (i in 1:length(uGenes6)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes6[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/Tb_CS5/Markerscatter_",as.character(uGenes6[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/ExMes_CS5/",sep=""))
uGenes7 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="ExMes_CS5")] )
for (i in 1:length(uGenes7)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes7[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/ExMes_CS5/Markerscatter_",as.character(uGenes7[i]),".pdf",sep=""),width = 35, height = 8)
}


library("ggtern")


subs_data <-  subset(mammal.combined, idents = c("ICM_CS3"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 

ICM1 <- FindMarkers(subs_data,ident.1 = "ICM_CS3_Marm", ident.2 = "ICM_CS3_Cyno")
ICM2 <- FindMarkers(subs_data,ident.1 = "ICM_CS3_Marm", ident.2 = "ICM_CS3_Human_invit")
ICM3 <- FindMarkers(subs_data,ident.1 = "ICM_CS3_Cyno", ident.2 = "ICM_CS3_Human_invit")


#merge_Synth1 <- mammal.combined
#Cl.cells <- subset(merge_Synth1, idents = c(4))
#Idents(Cl.cells) <- "divergence1"
#avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
#avg.Cl.cells$gene <- rownames(avg.Cl.cells)
#genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
##p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
#ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl4.pdf",sep=""),width = 20, height = 20, plot = p1)



subs_data <-  subset(mammal.combined, idents = c("ICM_CS3"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 
ggplot(log1p(avexp$RNA), aes(ICM_CS3_Marm, ICM_CS3_Cyno)) + geom_point() + theme_classic()  #+ ggtitle("CD4 Naive T Cells")
ggsave(filename=paste(saveext,"/Markers/ICM_CS3","_MarmCyno.pdf",sep=""),width = 8, height = 8)
ggplot(log1p(avexp$RNA), aes(ICM_CS3_Marm, ICM_CS3_Human_invit)) + geom_point() + theme_classic()  #+ ggtitle("CD4 Naive T Cells")
ggsave(filename=paste(saveext,"/Markers/ICM_CS3","_MarmHuman.pdf",sep=""),width = 8, height = 8)

ICM1 <- FindMarkers(subs_data,ident.1 = "ICM_CS3_Marm", ident.2 = "ICM_CS3_Cyno")
ICM2 <- FindMarkers(subs_data,ident.1 = "ICM_CS3_Marm", ident.2 = "ICM_CS3_Human_invit")

#Epi_CS4 
#Epi_CS3 
#Epi
subs_data <-  subset(mammal.combined, idents = c("Epi_CS4","Epi_CS3"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 
Epi1 <- FindMarkers(subs_data,ident.1 = "Epi_CS3_Marm", ident.2 = "Epi_CS3_Cyno")
Epi2 <- FindMarkers(subs_data,ident.1 = "Epi_CS3_Marm", ident.2 = "Epi_CS4_Human_invit")
#Hyp
subs_data <-  subset(mammal.combined, idents = c("Hyp_CS3","VE_CS4"), invert = FALSE )
newLabs <- as.character(Idents(subs_data))
newLabs[which(newLabs=="VE_CS4")] <- "Hyp_CS3"
Idents(subs_data) <- newLabs
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 


Hyp1 <- FindMarkers(subs_data,ident.1 = "Hyp_CS3_Marm", ident.2 = "Hyp_CS3_Cyno")
Hyp2 <- FindMarkers(subs_data,ident.1 = "Hyp_CS3_Marm", ident.2 = "Hyp_CS3_Human_invit")
#Tb
subs_data <-  subset(mammal.combined, idents = c("Tb_CS3","Tb_CS4"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 
Tb1 <- FindMarkers(subs_data,ident.1 = "Tb_CS3_Marm", ident.2 = "Tb_CS3_Cyno")
Tb2 <- FindMarkers(subs_data,ident.1 = "Tb_CS3_Marm", ident.2 = "Tb_CS4_Human_invit")
#EmDisc_CS5
subs_data <-  subset(mammal.combined, idents = c("EmDisc_CS5"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) 
EmDisc1 <- FindMarkers(subs_data,ident.1 = "EmDisc_CS5_Marm", ident.2 = "EmDisc_CS5_Cyno")
EmDisc2 <- FindMarkers(subs_data,ident.1 = "EmDisc_CS5_Marm", ident.2 = "EmDisc_CS5_Human_invit")

IntG <- intersect( intersect(rownames(marmoset_data), rownames(cynomolgous_data)), rownames(human_dataA) )

allMarkers <- rbind(ICM1,ICM2,Epi1,Epi2,Hyp1,Hyp2,Tb1,Tb2,EmDisc1,EmDisc2)
allMarkers <- allMarkers[which(allMarkers$p_val_adj<0.05),]

uGenes <- intersect(IntG, rownames(allMarkers) )

#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_6.pdf",sep=""),width = 42, height = 8)
dir.create(paste(saveext,"/MarkersDiff/",sep=""))
#uGenes <- unique( rownames(allMarkers))
for (i in 1:length(uGenes)) {
  FeaturePlot(mammal.combined, features = as.character(uGenes[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 2)
  ggsave(filename=paste(saveext,"/MarkersDiff/Markerscatter_",as.character(uGenes[i]),".pdf",sep=""),width = 16, height = 8)
}


FeaturePlot(mammal.combined, features = "T", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste("~/Desktop/Markerscatter_T.pdf",sep=""),width = 35, height = 8)

FeaturePlot(mammal.combined, features = "POU5F1", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste("~/Desktop/Markerscatter_POU5F1.pdf",sep=""),width = 35, height = 8)

FeaturePlot(mammal.combined, features = "KDR", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_KDR.pdf",sep=""),width = 16, height = 4)

FeaturePlot(mammal.combined, features = "TTR", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_TTR.pdf",sep=""),width = 16, height = 4)


FeaturePlot(mammal.combined, features = "TTR", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste("~/Desktop/Markerscatter_TTR1.pdf",sep=""),width = 35, height = 8)

FeaturePlot(mammal.combined, features = "TTR", reduction = "pca", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste("~/Desktop/Markerscatter_TTR2.pdf",sep=""),width = 42, height = 8)




#####
FeaturePlot(mammal.combined, reduction = "pca", features = "TFAP2C", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 2)

ggsave(filename=paste(saveext,"/MarkersDiff/Markerscatter_TFAP2c.pdf",sep=""),width = 8, height = 4)

fulllist <- intersect(intersect(rownames(marmoset_data), rownames(cynomolgous_data)), rownames(human_dataA))
#data_A[data_A$index %in% data_B$index,]

subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Tb_CS3"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
avexp2 <- subset(avexp$RNA, rownames(avexp$RNA) %in% fulllist)

C1 <- cor(  log2(avexp2 +1),log2(avexp2 +1) )
#C1 <- cor(avexp2,avexp2)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 2), breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes.pdf",sep=""))

C1 <- cor( log2(avexp$integrated +1) ,log2(avexp$integrated +1 ) )
#C1 <- cor(avexp$integrated,avexp$integrated)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 2), breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_integratedgenes.pdf",sep=""))




subs_data2 <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Tb_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS6/7"), invert = FALSE )
subs_data2$celltype.cond <- paste(Idents(subs_data2), subs_data2$divergence1, sep = "_")
Idents(subs_data2) <- "celltype.cond"
avexp_  <- AverageExpression(object = subs_data2) #, use.scale = TRUE)
avexp2_ <- subset(avexp_$RNA, rownames(avexp_$RNA) %in% fulllist)

C1 <- cor(  log2(avexp2_ +1),log2(avexp2_ +1) )
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 2), fontsize_number = 3, breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended_allgenes.pdf",sep=""))

C1 <- cor( log2(avexp_$integrated +1) ,log2(avexp_$integrated +1 ) )
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended_inntegratedgenes.pdf",sep=""))



subs_data <-  subset(mammal.combined, idents = c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS6/7"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
avexp2 <- subset(avexp$RNA, rownames(avexp$RNA) %in% fulllist)

C1 <- cor(  log2(avexp2 +1),log2(avexp2 +1) )
#C1 <- cor(avexp2,avexp2)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 2), breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended2_allgenes.pdf",sep=""))


subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Tb_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS6/7","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7","SYS_CS6/7"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
C1 <- cor( log2(avexp_$integrated +1) ,log2(avexp_$integrated +1 ) )
#C1 <- cor(avexp$integrated,avexp$integrated)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended2_integratedgenes.pdf",sep=""))



#subs_data <-  mammal.combined 
#subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
#Idents(subs_data) <- "celltype.cond"
#avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
#C1 <- cor(avexp$integrated,avexp$integrated)
#mat_breaks <- seq(0.5, 1, length.out = 20)
#redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended2.pdf",sep=""))

#write.csv(as.data.frame(avexp$RNA), file=paste(saveext,"/DimRed/ExpressionAllDatasets.csv",sep=""))



#Now incorporate other datasets ...
#Rename to keep variable convention for the join species modelling
mammal.anchors1 <- FindIntegrationAnchors(object.list = list(marmoset_data, cynomolgous_data, human_dataA, Cyno_dataB), dims = 1:20, anchor.features = intersect( rownames(Cyno_dataB) , rownames(avexp$integrated)) )
mammal.combined <- IntegrateData(anchorset = mammal.anchors1, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

types <- Idents(mammal.combined)



cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


#cols <- c("#1A0873","#13084D","#0039E6","#0233BF","#0233BF","#0233BF","#00BFBF","#00BFBF","#F04C04","#BF3C04","#BF3C04","#992F03","#E6B500","#00E6E6","#E6E600","#E6E600","#BF7104","#BF7104","#BF7104","#E68600","#BF0489","#BF0489","#921FE6","#7813BF","#7813BF","#7813BF","#E6B500","#BF9600")

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="VE_CS4")] <- 1
cluster_letters <- as.factor(cluster_letters)

#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_4dataset_","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_4dataset_","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_4dataset_","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_4dataset_","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_4dataset_","_12.pdf",sep=""),width = 42, height = 8)



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab_4dataset_","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab_4dataset_","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab_4dataset_","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "pcap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab_4dataset_","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab_4dataset_","_12.pdf",sep=""),width = 42, height = 8)



#FeaturePlot(mammal.combined, features = "NODAL", split.by = "species")
#ggsave(filename=paste(saveext,"UMAP_Nodal_Split",".pdf",sep=""),width = 24, height = 8)

#Can do the same for violin plot representations
#VlnPlot(mammal.combined, features = "NODAL", split.by = "species")
#ggsave(filename=paste(saveext,"Cluster_Nodal_Split",".pdf",sep=""),width = 24, height = 24)



subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Tb_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS6/7"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
C1 <- cor(avexp$integrated,avexp$integrated)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended_4dataset.pdf",sep=""))

subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Tb_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
C1 <- cor(avexp$integrated,avexp$integrated)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_4dataset.pdf",sep=""))

subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7"), invert = FALSE )
subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
Idents(subs_data) <- "celltype.cond"
avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
C1 <- cor(avexp$integrated,avexp$integrated)
mat_breaks <- seq(0.5, 1, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_emb_4dataset.pdf",sep=""))





#subs_data <-  subset(mammal.combined, idents = c("ICM_CS3","Hyp_CS3","Epi_CS3","Epi_CS4","Tb_CS3","Tb_CS4","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","EmDisc_CS6/7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_CS6/7","VE_CS4","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7","SYS_CS6/7"), invert = FALSE )
#subs_data$celltype.cond <- paste(Idents(subs_data), subs_data$divergence1, sep = "_")
#subs_data$celltype <- Idents(subs_data)
#Idents(subs_data) <- "celltype.cond"
#avexp  <- AverageExpression(object = subs_data) #, use.scale = TRUE)
#C1 <- cor(avexp$integrated,avexp$integrated)
#mat_breaks <- seq(0.5, 1, length.out = 20)
#redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, display_numbers = round(C1, digits = 2), fontsize_number = 3, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_extended2_4dataset.pdf",sep=""))











#Idents(mammal.combined) <- mammal.combined$Cl1
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl1","_4.pdf",sep=""),width = 42, height = 8)
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl1","_6.pdf",sep=""),width = 42, height = 8)
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl1","_8.pdf",sep=""),width = 42, height = 8)
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl1","_10.pdf",sep=""),width = 42, height = 8)
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl1","_12.pdf",sep=""),width = 42, height = 8)

#Idents(mammal.combined) <- mammal.combined$Cl2
#DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl2","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl2","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl2","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl2","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl2","_12.pdf",sep=""),width = 42, height = 8)

Idents(mammal.combined) <- mammal.combined$Cl3
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl3","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl3","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl3","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl3","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Cl3","_12.pdf",sep=""),width = 42, height = 8)


#No look for differences in genes so we can do motif analysis ..
Idents(mammal.combined) <- mammal.combined$Cl1
merge_Synth1 <- mammal.combined
merge_Synth1$celltype.cond <- paste(Idents(merge_Synth1), merge_Synth1$divergence1, sep = "_")
Idents(merge_Synth1) <- "celltype.cond"

#1,4,7

#Tb_CS3
Cl_1_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = "1_Newworld", ident.2 = "1_Oldworld", verbose = FALSE)
Cl_1_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = "1_Newworld", ident.2 = "1_Ape", verbose = FALSE)
Cl_1_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = "1_Oldworld", ident.2 = "1_Ape", verbose = FALSE)

#ICM_CS3
Cl_4_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = "4_Newworld", ident.2 = "4_Oldworld", verbose = FALSE)
Cl_4_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = "4_Newworld", ident.2 = "4_Ape", verbose = FALSE)
Cl_4_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = "4_Oldworld", ident.2 = "4_Ape", verbose = FALSE)

#Hyp_CS3
Cl_7_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = "7_Newworld", ident.2 = "7_Oldworld", verbose = FALSE)
Cl_7_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = "7_Newworld", ident.2 = "7_Ape", verbose = FALSE)
Cl_7_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = "7_Oldworld", ident.2 = "7_Ape", verbose = FALSE)

#Epi_CS3
Cl_6_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = "6_Newworld", ident.2 = "6_Oldworld", verbose = FALSE)
Cl_6_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = "6_Newworld", ident.2 = "6_Ape", verbose = FALSE)
Cl_6_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = "6_Oldworld", ident.2 = "6_Ape", verbose = FALSE)

#VE_CS5
Cl_10_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = "10_Newworld", ident.2 = "10_Oldworld", verbose = FALSE)
Cl_10_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = "10_Newworld", ident.2 = "10_Ape", verbose = FALSE)
Cl_10_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = "10_Oldworld", ident.2 = "10_Ape", verbose = FALSE)

#Tb_CS5/7
Cl_02_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = c("0_Newworld","2_Newworld"), ident.2 = c("0_Oldworld","2_Oldworld"), verbose = FALSE)
Cl_02_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = c("0_Newworld","2_Newworld"), ident.2 = c("0_Ape","2_Ape"), verbose = FALSE)
Cl_02_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = c("0_Oldworld","2_Oldworld"), ident.2 = c("0_Ape","2_Ape"), verbose = FALSE)

#Epi_CS5/7
Cl_39_MarmCyno <- FindMarkers(merge_Synth1, ident.1 = c("3_Newworld","9_Newworld"), ident.2 = c("3_Oldworld","9_Oldworld"), verbose = FALSE)
Cl_39_MarmHuman <- FindMarkers(merge_Synth1, ident.1 = c("3_Newworld","9_Newworld"), ident.2 = c("3_Ape","9_Ape"), verbose = FALSE)
Cl_39_CynoHuman <- FindMarkers(merge_Synth1, ident.1 = c("3_Oldworld","9_Oldworld"), ident.2 = c("3_Ape","9_Ape"), verbose = FALSE)

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(1))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log2(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl1.pdf",sep=""),width = 20, height = 20, plot = p1)

library("ggtern")

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(4))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl4.pdf",sep=""),width = 20, height = 20, plot = p1)



merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(7))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl7.pdf",sep=""),width = 20, height = 20, plot = p1)

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(6))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl6.pdf",sep=""),width = 20, height = 20, plot = p1)

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(10))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl10.pdf",sep=""),width = 20, height = 20, plot = p1)
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl10_nolab.pdf",sep=""),width = 20, height = 20, plot = p1)

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(0,2))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl02.pdf",sep=""),width = 20, height = 20, plot = p1)
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl02_nolab.pdf",sep=""),width = 20, height = 20, plot = p1)

merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c(3,9))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindMarkers(Cl.cells, ident.1 = "Oldworld", ident.2 = "Newworld"))
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl39.pdf",sep=""),width = 20, height = 20, plot = p1)
p1 <- ggplot(avg.Cl.cells, aes(Newworld,Oldworld)) + geom_point() + ggtitle("Marm vs Cyno") + theme_minimal()
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/Marm_Cyno_Cl39_nolab.pdf",sep=""),width = 20, height = 20, plot = p1)




#Now whole compare labelled datasets
Idents(mammal.combined) <- mammal.combined$IDs
merge_Synth1 <- mammal.combined
Cl.cells <- subset(merge_Synth1, idents = c("ICM_CS3"))
Idents(Cl.cells) <- "divergence1"
avg.Cl.cells <- log1p(AverageExpression(Cl.cells, verbose = FALSE)$RNA)
avg.Cl.cells$gene <- rownames(avg.Cl.cells)
genes.to.label = rownames(FindAllMarkers(Cl.cells))
D1 <- as.data.frame(avg.Cl.cells)
D1$colsi <- 0
D1[rownames(FindAllMarkers(Cl.cells)),5] <- 1
D1$colsi <- as.factor(D1$colsi)
D1 <-D1[order(D1$colsi),]
D2 <- D1[which(D1$colsi==1),]
p1 <- ggtern(data=D1, aes(x=Newworld,y=Oldworld, z=Ape)) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("grey", "black")) 
p1 <- p1 +  annotate(geom  = 'text',x     = D2$Newworld/(D2$Newworld+D2$Oldworld+D2$Ape),y= D2$Oldworld/(D2$Newworld+D2$Oldworld+D2$Ape),z= D2$Ape/(D2$Newworld+D2$Oldworld+D2$Ape),label = rownames(D2),color = c("red")) 
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/ICM_CS3.pdf",sep=""),width = 40, height = 40, plot = p1)


p1 <- ggtern(data=D1, aes(x=Newworld,y=Oldworld, z=Ape)) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("grey", "black")) 
p1 <- p1 +  annotate(geom  = 'text',x     = D2$Newworld/(D2$Newworld+D2$Oldworld+D2$Ape),y= D2$Oldworld/(D2$Newworld+D2$Oldworld+D2$Ape),z= D2$Ape/(D2$Newworld+D2$Oldworld+D2$Ape),label = rownames(D2),color = c("red")) 
ggsave(filename=paste("~/Desktop/Thorsten/FINAL/CCAFigForPaper/ICM_CS3.pdf",sep=""),width = 40, height = 40, plot = p1)



#split$`1) Marmoset`
#split$`2) Cynomolgous`
#split$`4) Human (in vitro)`

#Align to other annotations
Idents(mammal.combined) <- types

colKey1 <- levels(Idents(mammal.combined))
ID1 <- factor(Idents(mammal.combined)) #factor(c(as.character(labs), as.character(cylabs),  as.character(hlabs)   ) )

Idents(mammal.combined) <- factor(as.character(mammal.combined$AltID)) #factor(c(as.character(labs), as.character(cylabs2),  as.character(hlabs2)   ) )
colKey2 <- levels(Idents(mammal.combined))
ID2 <- Idents(mammal.combined) #factor(c(as.character(labs), as.character(cylabs2),  as.character(hlabs2)   ) )

ind1 <- integer(length(colKey2))
for (i in 1:length(colKey2)) {
  ind1[i] <- which(colKey1==ID1[which(ID2==colKey2[i])[1]])
}

necols <- coluse[ind1]



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_4.pdf",sep=""),width = 42, height = 8, Dingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_6.pdf",sep=""),width = 42, height = 8, Dingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_8.pdf",sep=""),width = 42, height = 8, Dingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_10.pdf",sep=""),width = 42, height = 8, Dingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_12.pdf",sep=""),width = 42, height = 8, Dingbats = FALSE)



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_4.pdf",sep=""),width = 42, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_6.pdf",sep=""),width = 42, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_8.pdf",sep=""),width = 42, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 10, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_10.pdf",sep=""),width = 42, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 12, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_12.pdf",sep=""),width = 42, height = 8, useDingbats = FALSE)


#Check the clusters




#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_data, cynomolgous_data, cyno_dataB), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)

Cls <- Idents(mammal.combined)
Idents(mammal.combined) <- factor(c(as.character(labs), as.character(cylabs),  as.character(cyivBS$CS)   ) )



cols <- c("#1A0873","#13084D","#0039E6","#1A0873","#921FE6","#0039E6","#0233BF","#0233BF","#0233BF","#00BFBF","#00BFBF","#F04C04","#BF3C04","#BF3C04","#992F03","#E6B500","#00E6E6","#E6E600","#E6E600","#BF7104","#BF7104","#BF7104","#E68600","#BF0489","#BF0489","#921FE6","#7813BF","#7813BF","#7813BF","#E6B500","#BF9600")

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="VE_CS4")] <- 1
cluster_letters <- as.factor(cluster_letters)

#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cyno","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cyno","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cyno","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cyno","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cyno","_12.pdf",sep=""),width = 42, height = 8)


Idents(mammal.combined) <- Cls

DimPlot(mammal.combined,  shape.by = 'cell.orig', pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_cynoCl","_4.pdf",sep=""),width = 42, height = 8)

mammal.combined$CellType <- factor(c(as.character(labs), as.character(cylabs),  as.character(cyivBS$CS)   ) )

#Write everythign out to csv
write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsCl_cls.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNE_cls.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddings_cls.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey_cls.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCA_cls.csv",sep=""))



cyivBS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cyinvitLab.csv",sep=",", header = T, row.names=1)



#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_data, cynomolgous_data, human_dataA, Cyno_dataB), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)


Idents(mammal.combined) <- factor(c(as.character(labs), as.character(cylabs),  as.character(hlabs), as.character(cyivBS$Type)   ) )

cols = c("#1A0873","#13084D","#13084D","#0039E6","#0233BF","#0233BF","#0233BF","#00BFBF","#00BFBF","#F04C04","#F04C04","#BF3C04","#BF3C04","#992F03","#E6B500","#00E6E6","#E6E600","#E6E600","#E6E600","#BF7104","#BF7104","#BF7104","#E68600","#BF0489","#BF0489","#921FE6","#921FE6","#7813BF","#7813BF","#BF9600","#BF9600")



cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="VE_CS4")] <- 1
cluster_letters <- as.factor(cluster_letters)

#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')


DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_all","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_all","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_all","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_all","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = cols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_all","_12.pdf",sep=""),width = 42, height = 8)



