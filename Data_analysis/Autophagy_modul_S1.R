library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(scales)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(dunn.test)

###############
## Load data ##
###############

N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200617_40X_Autophagy_Modulation_Chloro_S1_R1/Analysis_20200911_100737/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200617_40X_Autophagy_Modulation_Chloro_S1_R2/Analysis_20200911_101157/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200617_40X_Autophagy_Modulation_Chloro_S1_R3/Analysis_20200911_101621/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R1/Analysis_20200911_104427/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R2/Analysis_20200911_104816/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table('S:/HCS_Platform/Data/JavierJarazo/Autophagy_modulation/S1/JJ_20200717_40X_Autophagy_Modulation_Chloro_S1_C2_R3/Analysis_20200911_105224/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

Table_Chloro = rbind(N1,N2,N3,N4,N5,N6)

Table_Chloro$AreaName = str_replace(Table_Chloro$AreaName,"_Un_","_Chloro_")

Table_all = Table_Chloro
tmp = do.call(rbind,strsplit(Table_all$AreaName,"_"))
Table_all$Line = tmp[,1]
Table_all$Treatment = tmp[,2]
Table_all$Concentration = tmp[,3]

Treatments = c("Baf","Chloro","Rapa")
Controls = c("T12","K7","KOLF","68")
Patients = c("2122","825","4c1")
Table_all$Condition = ifelse(Table_all$Treatment %in% Treatments, 'Treated', 'Untreated')
Table_all$Category =  ifelse(Table_all$Line %in% Patients, 'Patient', ifelse(Table_all$Line %in% Controls, 'Controls', 'Corrected'))
Table_all$THTotalArea = Table_all$THNeuriteArea + Table_all$THSomaArea
Table_all$THTotalAreaNorm = Table_all$THTotalArea/Table_all$NucArea
Table_all$THSomaAreaNorm = Table_all$THSomaArea/Table_all$NucArea
Table_all$THNeuriteAreaNorm = Table_all$THNeuriteArea/Table_all$NucArea
Table_all$ColoLLTSAreaNormNeu = Table_all$ColoLLTSArea/Table_all$THTotalArea
Table_all$CountColoLLTSNormNeu = Table_all$CountColoLLTS/Table_all$THTotalArea
Table_all$ColoLLTNAreaNormNeu = Table_all$ColoLLTNArea/Table_all$THTotalArea
Table_all$CountColoLLTNNormNeu = Table_all$CountColoLLTN/Table_all$THTotalArea
Table_all$ColoLCTSAreaNormNeu = Table_all$ColoLCTSArea/Table_all$THTotalArea
Table_all$CountColoLCTSNormNeu = Table_all$CountColoLCTS/Table_all$THTotalArea
Table_all$ColoLCTNAreaNormNeu = Table_all$ColoLCTNArea/Table_all$THTotalArea
Table_all$CountColoLCTNNormNeu = Table_all$CountColoLCTN/Table_all$THTotalArea
Table_all$ColoLATSAreaNormNeu = Table_all$ColoLATSArea/Table_all$THTotalArea
Table_all$CountColoLATSNormNeu = Table_all$CountColoLATS/Table_all$THTotalArea
Table_all$ColoLATNAreaNormNeu = Table_all$ColoLATNArea/Table_all$THTotalArea
Table_all$CountColoLATNNormNeu = Table_all$CountColoLATN/Table_all$THTotalArea
Table_all$ColoLLTHNorm = (Table_all$ColoLLTSArea + Table_all$ColoLLTNArea)/Table_all$NucArea

Table_all %>% 
  dplyr::filter(NucArea %between% c(10000,300000))->Table_all

saveRDS(Table_all,"Table_all_Autophagy_modul_S1.rdata")

interesting_features = c("NucArea","NucDeadArea","THTotalArea","THTotalAreaNorm",
                         "THSomaArea","THNeuriteArea","THSomaAreaNorm","THNeuriteAreaNorm",  
                         "LC3Area","LC3AreaNorm","CountLC3","CountLC3Norm",
                         "LAMPArea","LAMPAreaNorm","CountLAMP","CountLAMPNorm",                  
                         "ColoLLArea","ColoLLAreaNorm","CountColoLL","CountColoLLNorm",
                         "ColoLLTSArea","ColoLLTSAreaNorm","CountColoLLTS","CountColoLLTSNorm","ColoLLTSAreaNormNeu","CountColoLLTSNormNeu",
                         "ColoLLTNArea","ColoLLTNAreaNorm","CountColoLLTN","CountColoLLTNNorm","ColoLLTNAreaNormNeu","CountColoLLTNNormNeu",
                         "ColoLCTSArea","ColoLCTSAreaNorm","CountColoLCTS","CountColoLCTSNorm","ColoLCTSAreaNormNeu","CountColoLCTSNormNeu",                
                         "ColoLCTNArea","ColoLCTNAreaNorm","CountColoLCTN","CountColoLCTNNorm","ColoLCTNAreaNormNeu","CountColoLCTNNormNeu",               
                         "ColoLATSArea","ColoLATSAreaNorm","CountColoLATS","CountColoLATSNorm","ColoLATSAreaNormNeu","CountColoLATSNormNeu",
                         "ColoLATNArea","ColoLATNAreaNorm","CountColoLATN","CountColoLATNNorm","ColoLATNAreaNormNeu","CountColoLATNNormNeu",
                         "ColoLLTHNorm")

metadata =c("File","Barcode",                          
            "AreaName","ROW","COL","Field",
            "Line","Treatment","Concentration",                   
            "Condition","Category")

Table_all %>% 
  dplyr::select(metadata,interesting_features) %>% 
  pivot_longer(-metadata, names_to = "Feature", values_to = "measure") -> Table_all_long

Table_all_long %>% 
  group_by(Barcode,AreaName,Feature) %>% 
  dplyr::summarise(q75=(quantile(measure,probs=0.75,na.rm = TRUE))) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("Barcode","AreaName","Feature")) %>% 
  group_by(Treatment,Feature) %>% 
  dplyr::summarise(limmax = pretty(max(q75))[2]) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("Treatment","Feature"))-> All_tidy

saveRDS(All_tidy,"All_tidy.rdata")

###########################
## GRAPHS for figure size##
###########################
###################################################
## Violin with line through mnedian of timepoints #
###################################################
Table_all %>%
  dplyr::filter(Line!="68") %>% 
  dplyr::filter(Category!="Corrected") %>% 
  dplyr::filter(Treatment=="Chloro") ->Table_all

Table_all$Concentration[Table_all$Concentration==5]<-"50"
Table_all$Concentration[Table_all$Concentration==0]<-"Un"
Table_all$Concentration[Table_all$Concentration==1]<-"500"
Table_all$Concentration[Table_all$Concentration==2]<-"1"
Table_all$Concentration[Table_all$Concentration==3]<-"5"
Table_all$Concentration[Table_all$Concentration==4]<-"10"

Features = c("LC3AreaNorm","CountLC3Norm","LAMPAreaNorm","CountLAMPNorm",
             "ColoLCTNAreaNorm","ColoLCTNAreaNormNeu","ColoLCTSAreaNorm","ColoLCTSAreaNormNeu",
             "ColoLATNAreaNorm","ColoLATNAreaNormNeu","ColoLATSAreaNorm","ColoLATSAreaNormNeu",
             "ColoLLTNAreaNorm","ColoLLTNAreaNormNeu","ColoLLTSAreaNorm","ColoLLTSAreaNormNeu",
             "NucArea","THTotalArea","ColoLLAreaNorm","ColoLLTHNorm")

levels_x_axis = c("Un","500","1","5","10","50")
Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), X2= c(1,0.2,0.1,0.01,0.1,0.2,0.2,0.2,0.005,0.01,0.025,0.025,0.001,0.002,0.005,0.005,3e5,2e5,0.03,0.005))   

graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(Table_all,aes_string(y=feat, x="Concentration", fill="Category")) + 
    geom_violin(aes(group=interaction(Category,Concentration)),alpha=0.4,position=position_dodge(width = 0.6),lwd=0.1) + 
    geom_boxplot(aes(group=interaction(Category,Concentration)), fill="white",position=position_dodge(width = 0.6),outlier.size = 0.1,lwd=0.1,width=0.3) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("cornflowerblue","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 4),expand = c(0,0)) +
    scale_x_discrete(limits = levels_x_axis) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=8),axis.text.x= element_text(size=5),axis.text.y= element_text(size=5),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(Table_all,aes(y=LC3AreaNorm, x=Concentration, fill=Category)) + 
  geom_violin(aes(group=interaction(Category,Concentration)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(aes(group=interaction(Category,Concentration)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
  stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("cornflowerblue","firebrick1")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .45, width = .15, height = .15) +
  draw_plot(graph[[2]], x = .15, y = .45, width = .15, height = .15) +
  draw_plot(graph[[3]], x = .3, y = .45, width = .15, height = .15) +
  draw_plot(graph[[4]], x = .45, y = .45, width = .15, height = .15) +
  draw_plot(graph[[5]], x = 0, y = .3, width = .15, height = .15) +
  draw_plot(graph[[6]], x = .15, y = .3, width = .15, height = .15) +
  draw_plot(graph[[7]], x = .3, y = .3, width = .15, height = .15) +
  draw_plot(graph[[8]], x = .45, y = .3, width = .15, height = .15) +
  draw_plot(graph[[9]], x = .0, y = .15, width = .15, height = .15) +
  draw_plot(graph[[10]], x = .15, y = .15, width = .15, height = .15) +
  draw_plot(graph[[11]], x = .3, y = .15, width = .15, height = .15) +
  draw_plot(graph[[12]], x = .45, y = .15, width = .15, height = .15) +
  draw_plot(graph[[13]], x = .0, y = .0, width = .15, height = .15) +
  draw_plot(graph[[14]], x = .15, y = .0, width = .15, height = .15) +
  draw_plot(graph[[15]], x = .3, y = .0, width = .15, height = .15) +
  draw_plot(graph[[16]], x = .45, y = .0, width = .15, height = .15) +
  draw_plot(graph[[17]], x = .0, y = .6, width = .15, height = .15) +
  draw_plot(graph[[18]], x = .15, y = .6, width = .15, height = .15) +
  draw_plot(graph[[19]], x = .3, y = .6, width = .15, height = .15) +
  draw_plot(graph[[20]], x = .45, y = .6, width = .15, height = .15) +
  draw_plot(my_legend, x = .93, y = .2, width = .00000001, height = .00000001)

ggsave("Autophagy_modul_S1.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


#########
##STATS##
######### 
assign_significance = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 'ns'
  significance[pvals <0.05 & pvals >0.01] = '*' 
  significance[pvals <=0.01 & pvals >0.001] = '**'
  significance[pvals <=0.001 & pvals >0.0001] = '***'
  significance[pvals <=0.0001] = '****'
  return(significance)
}

assign_SM = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 0
  significance[pvals <0.05 & pvals >0.01] = 1 
  significance[pvals <=0.01 & pvals >0.001] = 2
  significance[pvals <=0.001 & pvals >0.0001] = 3
  significance[pvals <=0.0001] = 4
  return(significance)
}

assign_Comp = function(Comparison){
  significance = rep(NA,length(Comparison))
  significance[Comparison =="Controls_Un - Patient_Un"] = 1
  significance[Comparison == "Controls_500 - Patient_500"] = 2 
  significance[Comparison =="Controls_1 - Patient_1"] = 3
  significance[Comparison =="Controls_5 - Patient_5"] = 4
  significance[Comparison =="Controls_10 - Patient_10"] = 5
  significance[Comparison =="Controls_50 - Patient_50"] = 6
  significance[Comparison =="Controls_500 - Controls_Un"] = 1
  significance[Comparison =="Controls_1 - Controls_500"] = 2
  significance[Comparison =="Controls_1 - Controls_5"] = 3
  significance[Comparison =="Controls_10 - Controls_5"] = 4
  significance[Comparison =="Controls_10 - Controls_50"] = 5
  significance[Comparison =="Patient_500 - Patient_Un"] = 1
  significance[Comparison =="Patient_1 - Patient_500"] = 2
  significance[Comparison =="Patient_1 - Patient_5"] = 3
  significance[Comparison =="Patient_10 - Patient_5"] = 4
  significance[Comparison =="Patient_10 - Patient_50"] = 5
  return(significance)
}

###################
## Kruskal-Wallis #
###################
#All features wanted comparisons
Table_all$Handle = paste(Table_all$Category, Table_all$Concentration, sep="_")
wanted =c("Controls_Un - Patient_Un", "Controls_500 - Patient_500","Controls_1 - Patient_1",
          "Controls_5 - Patient_5","Controls_10 - Patient_10","Controls_50 - Patient_50",
          "Controls_500 - Controls_Un","Controls_1 - Controls_500","Controls_1 - Controls_5",
          "Controls_10 - Controls_5","Controls_10 - Controls_50",
          "Patient_500 - Patient_Un","Patient_1 - Patient_500","Patient_1 - Patient_5",
          "Patient_10 - Patient_5","Patient_10 - Patient_50")
results = list()
for(feat in Features){
  kpvals = kruskal.test(Table_all[[feat]],as.factor(Table_all$Handle))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(interesting_features))
  tmp = dunn.test(Table_all[[feat]],as.factor(Table_all$Handle), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(interesting_features))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = as.data.table(do.call(rbind,results))
table_results$significancekw = assign_significance(as.double(table_results$kpvalsadj))
table_results$significancedt = assign_significance(as.double(table_results$padjusted))
table_results$significanceSM = assign_SM(as.double(table_results$padjusted))
table_results$order = assign_Comp(table_results$wanted)
table_results = dplyr::select(table_results,"wanted", "feat", "kpvals","kpvalsadj","significancekw","pvals","padjusted","significancedt","significanceSM","order")

write.table(table_results, "KW_wanted_comparisons_All_features_Automodul_S1.csv", sep = ',', row.names = FALSE)

tmp = do.call(rbind,strsplit(table_results$wanted," - "))
table_results$Handle1 = tmp[,1]
table_results$Handle2 = tmp[,2]

All_tidy %>% 
  dplyr::filter(Treatment=="Chloro") %>% 
  dplyr::filter(Feature %in% Features) %>%
  rename(feat = Feature) %>% 
  dplyr::group_by(Category,Concentration,feat) %>% 
  summarise(mymedian = median(measure,na.rm = TRUE)) %>% 
  ungroup()->All_medians

All_medians$Concentration[All_medians$Concentration==5]<-"50"
All_medians$Concentration[All_medians$Concentration==0]<-"Un"
All_medians$Concentration[All_medians$Concentration==1]<-"500"
All_medians$Concentration[All_medians$Concentration==2]<-"1"
All_medians$Concentration[All_medians$Concentration==3]<-"5"
All_medians$Concentration[All_medians$Concentration==4]<-"10"

All_medians$Handle1= paste0(All_medians$Category,"_",All_medians$Concentration)

All_medians %>% 
  rename(median1= mymedian) %>% 
  dplyr::select(Handle1,feat,median1)->All_medians1

toto = merge(table_results,All_medians1,by = c("Handle1","feat"))

All_medians1 %>% 
  rename(median2= median1,Handle2=Handle1) %>% 
  dplyr::select(Handle2,feat,median2)->All_medians2

All_merged = merge(toto,All_medians2,by = c("Handle2","feat"))

All_merged$Concentration1 = do.call(rbind,strsplit(All_merged$Handle1,"_"))[,2]
All_merged$Concentration2 = do.call(rbind,strsplit(All_merged$Handle2,"_"))[,2]
All_merged$HandleCat = paste0(do.call(rbind,strsplit(All_merged$Handle1,"_"))[,1],"_",do.call(rbind,strsplit(All_merged$Handle2,"_"))[,1])

All_merged$Concentration1[All_merged$Concentration1=="Un"]<-"0"
All_merged$Concentration1[All_merged$Concentration1=="500"]<-"0.5"
All_merged$Concentration2[All_merged$Concentration2=="Un"]<-"0"
All_merged$Concentration2[All_merged$Concentration2=="500"]<-"0.5"

All_merged$Multiplier1 = ifelse(All_merged$Concentration1>= All_merged$Concentration2,1,-1)

# Ramp creation -> for the comparisons of the different concentrations within category
# the higher the ramp value the higher the median in the higher concentration
# In the comparison Control Patient per concentration the higher the ramp value
# the higher the median in the Control
# This is since we are making if median1>median2 *1
All_merged$Multiplier2 = ifelse(All_merged$median1>= All_merged$median2,1,-1)
All_merged$SM = All_merged$Multiplier1*All_merged$Multiplier2*All_merged$significanceSM
write.table(All_merged, "KW_wanted_comparisons_All_features_Automodul_S1_SM.csv", sep = ',', row.names = FALSE)

###########################
##Selecting median images #
###########################
Features = c("LC3AreaNorm","CountLC3Norm","LAMPAreaNorm","CountLAMPNorm",
             "ColoLCTNAreaNorm","ColoLCTNAreaNormNeu","ColoLCTSAreaNorm","ColoLCTSAreaNormNeu",
             "ColoLATNAreaNorm","ColoLATNAreaNormNeu","ColoLATSAreaNorm","ColoLATSAreaNormNeu",
             "ColoLLTNAreaNorm","ColoLLTNAreaNormNeu","ColoLLTSAreaNorm","ColoLLTSAreaNormNeu",
             "NucArea","THTotalArea","ColoLLAreaNorm","ColoLLTHNorm")

Table_all %>% 
  dplyr::filter(Category!="Corrected") %>%
  dplyr::select(metadata,ColoGFAPS100bVolByNucVol) %>% 
  dplyr::filter(Timepoint==21)->Table_filtered 

Table_filtered %>% 
  group_by(Category) %>% 
  summarise(median= median(ColoGFAPS100bVolByNucVol)) %>% 
  ungroup() %>% 
  full_join(Table_filtered, by="Category") %>% 
  dplyr::filter(ColoGFAPS100bVolByNucVol>=median*0.9 & ColoGFAPS100bVolByNucVol<=median*1.1 ) ->Table_median

