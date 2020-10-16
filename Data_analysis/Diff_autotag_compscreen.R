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
N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C1_R1/Analysis_20200828_175649/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C1_R2/Analysis_20200828_180226/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C1_R3/Analysis_20200828_180503/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C2_R1/Analysis_20200828_180742/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C2_R2/Analysis_20200828_181021/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S1_C2_R3/Analysis_20200828_181301/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N7=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C1_R1/Analysis_20200828_181536/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N8=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C1_R2/Analysis_20200828_181815/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N9=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C1_R3/Analysis_20200828_182056/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N10=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C2_R1/Analysis_20200828_182334/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N11=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C2_R2/Analysis_20200828_184353/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N12=read.table('S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/JJ_20200810_60X_Auto_S2_C2_R3/Analysis_20200828_184708/ObjectsGrouped.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N13=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C1_R1/Analysis_20200907_190701/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N14=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C1_R2/Analysis_20200907_190945/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N15=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C1_R3/Analysis_20200907_191022/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N16=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C2_R1/Analysis_20200907_191058/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N17=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C2_R2/Analysis_20200907_191133/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N18=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S1_C2_R3/Analysis_20200907_191209/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N19=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C1_R1/Analysis_20200907_191245/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N20=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C1_R2/Analysis_20200907_191321/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N21=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C1_R3/Analysis_20200907_191357/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N22=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C2_R1/Analysis_20200907_191432/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N23=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C2_R2/Analysis_20200907_191507/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N24=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Autotag_Compscreen/LC3/Revision/KOLF/JJ_20200811_60X_Mito_S2_C2_R3/Analysis_20200907_191542/ObjectsGrouped.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)

#There was a mislabelling in the concetrations in the lay file in the Opera for the C2 plates
#It should have been rows 2 and 3 5um, rows 4 and 5 10um, rows 6 and 7 50um
Table_c2 = rbind(N4,N5,N6,N10,N11,N12,N16,N17,N18,N22,N23,N24)
tmp = do.call(rbind,strsplit(Table_c2$AreaName,"_"))
Table_c2$Concentration = tmp[,3]

Table_c2$Concentration[Table_c2$Concentration=="50"]<-"5"
Table_c2$Concentration[Table_c2$Concentration=="10"]<-"50"
Table_c2$Concentration[Table_c2$Concentration=="500"]<-"10"

Table_c1 = rbind(N1,N2,N3,N7,N8,N9,N13,N14,N15,N19,N20,N21)
tmp = do.call(rbind,strsplit(Table_c1$AreaName,"_"))
Table_c1$Concentration = tmp[,3]

Table_all = rbind(Table_c1,Table_c2)

Table_all$GroupCountNorm = Table_all$GroupCount / Table_all$Areaofcells
Table_all$sum_AreaNorm = Table_all$sum_Area / Table_all$Areaofcells

tmp = do.call(rbind,strsplit(Table_all$AreaName,"_"))
Table_all$Line = tmp[,1]
Table_all$Tag = tmp[,2]


Controls = c("68","KOLF")
Patients = c("2122","825","826","4C1")
Table_all$Category =  ifelse(Table_all$Line %in% Patients, 'Patient', ifelse(Table_all$Line %in% Controls, 'Control', 'Corrected'))

Table_all %>% 
  dplyr::filter(Areaofcells %between% c(5000,400000))->Table_all

saveRDS(Table_all,"Diff_autotag_comp_screen.rdata")

interesting_features = c("Areaofcells","GroupCount", "sum_Area",   "mean_Area",  "sem_Area",   "std_Area",  
                         "median_Area","GroupCountNorm","sum_AreaNorm" )

metadata =c("Barcode","AreaName","Column","Row","Field","TimePoint",
            "Line","Tag","Concentration","Category")

Table_all %>% 
  dplyr::select(metadata,interesting_features) %>% 
  pivot_longer(-metadata, names_to = "Feature", values_to = "measure") -> Table_all_long

Table_all_long %>% 
  group_by(Barcode,AreaName,TimePoint,Concentration,Feature) %>% 
  dplyr::summarise(q75=(quantile(measure,probs=0.75,na.rm = TRUE))) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("Barcode","AreaName","TimePoint","Concentration","Feature")) %>% 
  group_by(TimePoint,Feature) %>% 
  dplyr::summarise(limmax = pretty(max(q75))[2]) %>% 
  ungroup() %>% 
  full_join(Table_all_long, by =c("TimePoint","Feature"))-> All_tidy

saveRDS(All_tidy,"All_tidy_diff_autotag_compscreen.rdata")


###########################
## GRAPHS for figure size##
###########################
###################################################
## Violin with line through mnedian of timepoints #
###################################################
Table_all %>% 
  dplyr::filter(TimePoint==5) ->Table_all

Features = c("GroupCount","sum_Area","GroupCountNorm","sum_AreaNorm")

levels_x_axis = c("Un","500","1","5","10","50")
Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0), X2= c(600,15000,0.005,0.1))   

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
legend_plot = ggplot(Table_all,aes(y=Areaofcells, x=Concentration, fill=Category)) + 
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
  draw_plot(graph[[3]], x = 0, y = .3, width = .15, height = .15) +
  draw_plot(graph[[4]], x = .15, y = .3, width = .15, height = .15) +
  draw_plot(my_legend, x = .63, y = .2, width = .00000001, height = .00000001)

ggsave("Diff_autotag_compscreen.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


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
  significance[Comparison =="Control_Un - Patient_Un"] = 1
  significance[Comparison == "Control_500 - Patient_500"] = 2 
  significance[Comparison =="Control_1 - Patient_1"] = 3
  significance[Comparison =="Control_5 - Patient_5"] = 4
  significance[Comparison =="Control_10 - Patient_10"] = 5
  significance[Comparison =="Control_50 - Patient_50"] = 6
  significance[Comparison =="Control_500 - Control_Un"] = 1
  significance[Comparison =="Control_1 - Control_500"] = 2
  significance[Comparison =="Control_1 - Control_5"] = 3
  significance[Comparison =="Control_10 - Control_5"] = 4
  significance[Comparison =="Control_10 - Control_50"] = 5
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

wanted =c("Control_Un - Patient_Un", "Control_500 - Patient_500","Control_1 - Patient_1",
          "Control_5 - Patient_5","Control_10 - Patient_10","Control_50 - Patient_50",
          "Control_500 - Control_Un","Control_1 - Control_500","Control_1 - Control_5",
          "Control_10 - Control_5","Control_10 - Control_50",
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

write.table(table_results, "KW_wanted_comparisons_All_features_diff_autotag.csv", sep = ',', row.names = FALSE)

tmp = do.call(rbind,strsplit(table_results$wanted," - "))
table_results$Handle1 = tmp[,1]
table_results$Handle2 = tmp[,2]

All_tidy %>% 
  dplyr::filter(TimePoint==5) %>% 
  dplyr::filter(Feature %in% Features) %>%
  rename(feat = Feature) %>% 
  dplyr::group_by(Category,Concentration,feat) %>% 
  summarise(mymedian = median(measure,na.rm = TRUE)) %>% 
  ungroup()->All_medians

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
write.table(All_merged, "KW_wanted_comparisons_All_features_diff_autotag_SM.csv", sep = ',', row.names = FALSE)
