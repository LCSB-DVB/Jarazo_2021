library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(scales)
library(plyr)
library(tidyverse)
library(cowplot)
library(ggpubr)

###############
## Load data ##
###############
N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N1/Analysis_20200904_165644/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N2/Analysis_20200904_165936/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N3/Analysis_20200904_170046/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N1/Analysis_20200904_170151/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N2/Analysis_20200904_170307/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N3/Analysis_20200904_170423/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N7=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N1/Analysis_20200904_170535/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N8=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N2/Analysis_20200904_170654/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N9=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N3/Analysis_20200904_170814/ObjectsTuj1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N10=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N1/Analysis_20200904_165644/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N11=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N2/Analysis_20200904_165936/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N12=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d7_N3/Analysis_20200904_170046/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N13=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N1/Analysis_20200904_170151/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N14=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N2/Analysis_20200904_170307/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N15=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d14_N3/Analysis_20200904_170423/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N16=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N1/Analysis_20200904_170535/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N17=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N2/Analysis_20200904_170654/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N18=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Neuromorph/JJ_20171219_20X_Diff_Apop_d21_N3/Analysis_20200904_170814/ObjectsTH.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

Table_7_Tuj1 = rbind(N1,N2,N3)
Table_14_Tuj1 = rbind(N4,N5,N6)
Table_21_Tuj1 = rbind (N7,N8,N9)
Table_7_TH= rbind(N10,N11,N12)
Table_14_TH= rbind(N13,N14,N15)
Table_21_TH= rbind (N16,N17,N18)

Table_7 = full_join(Table_7_Tuj1,Table_7_TH, by=c("File","Barcode","AreaName", "ROW", "COL","Field"))
Table_14 = full_join(Table_14_Tuj1,Table_14_TH, by=c("File","Barcode","AreaName", "ROW", "COL","Field"))
Table_21 = full_join(Table_21_Tuj1,Table_21_TH, by=c("File","Barcode","AreaName", "ROW", "COL","Field"))

Table_7$Timepoint = 7
Table_14$Timepoint = 14
Table_21$Timepoint = 21

Table_all = rbind(Table_7,Table_14,Table_21)

Controls = c("T129","K7","KOLF")
Patients = c("2122","2124","947","4C1")
Table_all$Category =  ifelse(Table_all$AreaName %in% Patients, 'Patient', ifelse(Table_all$AreaName %in% Controls, 'Controls', 'Corrected'))

Table_all %>% 
  dplyr::filter(NucVol >20000)->Table_all

###########################
## Saving_tables_heatmap ##
###########################
Table_all %>% 
  dplyr::filter(AreaName!="2124") %>% #line removed for missing GC. Not included in the paper
  group_by(Barcode,Timepoint,Category) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE) ->toto2
#write.csv(toto2, "P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_neuromorph_2.csv",row.names=FALSE)

interesting_features = c("NucVol","Tuj1Area",  
                         "MinTuj1Area","MaxTuj1Area","MeanTuj1Area",   "StdTuj1Area",   
                         "MedTuj1Area","MadTuj1Area","TotTuj1Volume",  "CountTuj1", 
                         "MedMajorAxisLenghtTuj1","MedMinorAxisLenghtTuj1","MedIntensityTuj1",  "Tuj1SkelPixels",
                         "Tuj1PerimPixels","Tuj1BodyPixels", "Tuj1ShapeBySurface","Tuj1Fragmentation",
                         "TotalNodeCountTuj1","TotalLinkCountTuj1","NodesPerTuj1Mean",  "LinksPerTuj1Mean", 
                         "AverageNodeDegreeTuj1","MedianNodeDegreeTuj1","StdNodeDegreeTuj1", "MadNodeDegreeTuj1",
                         "THArea","MinTHArea",  "MaxTHArea",  "MeanTHArea",
                         "StdTHArea",  "MedTHArea",  "MadTHArea",  "TotTHVolume",   
                         "CountTH","MedMajorAxisLenghtTH"   ,"MedMinorAxisLenghtTH","MedIntensityTH",
                         "THSkelPixels",   "THPerimSomaPixels", "THPerimNeuritePixels","THBodySomaPixels", 
                         "THBodyNeuritePixels","THFragmentation","TotalNodeCountTH",  "TotalLinkCountTH", 
                         "NodesPerTHMean", "LinksPerTHMean", "AverageNodeDegreeTH","MedianNodeDegreeTH",
                         "StdNodeDegreeTH","MadNodeDegreeTH")

metadata =c("files","Row","Column","field",                 
            "AreaName","Barcode","Timepoint","Category")

Table_all %>% 
  dplyr::select(metadata,interesting_features) %>% 
  pivot_longer(-metadata, names_to = "Feature", values_to = "measure") -> Table_all_long

Table_all_long %>% 
  group_by(Barcode,AreaName,Feature) %>% 
  dplyr::summarise(q75=(quantile(measure,probs=0.75,na.rm = TRUE))) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("Barcode","AreaName","Feature")) %>% 
  group_by(Feature) %>% 
  dplyr::summarise(limmax = pretty(max(q75))[2]) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("Feature"))-> All_tidy

saveRDS(All_tidy,"All_tidy_neuromorph.rdata")