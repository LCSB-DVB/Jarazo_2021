library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(gganimate)

##LOADING OF DATA AND PREPARING TABLE

N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Organoids/Cluster_Output/PINK1_N/JJ_20190723_PINK1_N_A3A4_20190723_124757/data_all.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Organoids/Cluster_Output/PINK1_N/JJ_20190725_PINK1_N_A1-A4_20190725_093410/data_all.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Organoids/Cluster_Output/PINK1_N/JJ_20190729_PINK1_N_A1-A4_20190729_111528/data_all.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Organoids/Cluster_Output/PINK1_N/JJ_20190730_PINK1_N_A4_20190730_105423/data_all.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

TablesAll = rbind(N1,N2,N3,N4)

#Normalization
TablesAll$SkelNorm = TablesAll$SkelTH/TablesAll$THMaskSum
TablesAll$NodesNorm = TablesAll$Nodes/TablesAll$SkelTH
TablesAll$LinksNorm = TablesAll$Links/TablesAll$SkelTH
TablesAll$NucMaskDeadNorm = TablesAll$NucMaskDead/TablesAll$NucMaskSum

tmp = do.call(rbind,strsplit(TablesAll$AreaName,"_"))
tmp[,4][tmp[,4]=="4\t"]<- "4"
tmp[,4][tmp[,4]=="4\t\t"]<- "4"
TablesAll$Line = tmp[,1]
TablesAll$Condition = tmp[,2]
TablesAll$Experiment = tmp[,4]
Patients = c("4C1", "2122")
Controls = c("K7", "163")
TablesAll$Category = ifelse(TablesAll$Line %in% Patients, 'Patient', 'Control')
TablesAll = TablesAll[TablesAll$Line!="163",] #Line presented a chromosomal aberration in ch1
TablesAll$Handle = paste(TablesAll$Category, TablesAll$Condition, sep="_")

###########
## Plots ##
###########
Features = c("NucMaskSum", "Tuj1MaskSum","THMaskSum","MAP2MaskSum","MAP2ByTuj1","Tuj1ByNuc","Tuj1ByNucAlive","MAP2ByNuc","MAP2ByNucAlive","THByNuc","THByNucAlive","THByTuj1","THByMAP2","NucMaskDead","NucMaskAlive","THFragmentation","SkelTH","Nodes","Links","THPercentofAllNuc","Tuj1PercentofAllNuc","SkelNorm", "NodesNorm","LinksNorm","NucMaskDeadNorm") 
levels_x_axis = c("Un","500","1","5","10")

##################
## per Category ##
##################
Limit_y_axis = data.table(Features = Features, X1= c(0,4000000,0,4000000,0.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,80,0.02,0,0.04,0), X2= c(20000000,16000000,6000000,10000000,1.2,5,5,4,4,1,1,0.6,0.6,750,15000000,0.9,300000,15000,30000,80,110,0.08,0.08,0.14,0.001))

graph= list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TablesAll, aes_string(x='Condition', y=feat,fill = 'Category')) + 
    geom_boxplot(alpha=0.3,outlier.color = 'white',outlier.alpha = 0,lwd=0.25) + 
    geom_point(size=0.1, alpha=0.5,  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), aes(group=Category)) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    scale_x_discrete(limits=levels_x_axis) +
    scale_y_continuous (limits = c(lmmin,lmmax), expand = c(0,0)) +
    ggtitle(feat)+
    theme_classic() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +  
    theme(axis.text.x = element_text( size = 6), axis.title.x = element_blank(), axis.title.y=element_blank()) +
    theme(plot.title = element_text(size =12, hjust = 0.5))
}
Reduce(`+`, graph) * theme(legend.position = "none") 
ggsave("Organoids_PINK1_N_per_category.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")


###########
## STATS ##
###########
#All features all comparisons

results = list()
for(feat in Features){
  kpvals = kruskal.test(TablesAll[[feat]],as.factor(TablesAll$Handle))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
  tmp = dunn.test(TablesAll[[feat]],as.factor(TablesAll$Handle), list=TRUE)
  pvals = tmp$P
  wanted = tmp$comparisons
  padjusted =p.adjust(pvals,method="BH",n=length(Features)*66)
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)

write.table(table_results, "All_comparisons_All_features.csv", sep = ',', row.names = FALSE)

#wanted comparisons all features

wanted =c("Control_Un - Patient_Un", "Patient_500 - Patient_Un","Patient_1 - Patient_Un","Patient_5 - Patient_Un","Patient_10 - Patient_Un" )
results = list()
for(feat in Features){
  kpvals = kruskal.test(TablesAll[[feat]],as.factor(TablesAll$Handle))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
  tmp = dunn.test(TablesAll[[feat]],as.factor(TablesAll$Handle), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(Features))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)

write.table(table_results, "wanted_comparisons_All_features.csv", sep = ',', row.names = FALSE)

#wanted comparisons selected features
TableSelected =TablesAll
Selected_Features = c("THPercentofAllNuc") 
wanted =c("Control_Un - Patient_Un", "Patient_500 - Patient_Un","Patient_1 - Patient_Un","Patient_5 - Patient_Un","Patient_10 - Patient_Un")
results = list()
for(feat in Selected_Features){
  kpvals = kruskal.test(TableSelected[[feat]],as.factor(TableSelected$Handle))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Selected_Features))
  tmp = dunn.test(TableSelected[[feat]],as.factor(TableSelected$Handle), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(Selected_Features))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)


write.table(table_results, "wanted_comparisons_selected_features_THPercentofAllNuc.csv", sep = ',', row.names = FALSE)

#wanted comparisons selected features
TableSelected =TablesAll
Selected_Features = c("Tuj1PercentofAllNuc") 
wanted =c("Control_Un - Patient_Un", "Patient_500 - Patient_Un","Patient_1 - Patient_Un","Patient_5 - Patient_Un","Patient_10 - Patient_Un")
results = list()
for(feat in Selected_Features){
  kpvals = kruskal.test(TableSelected[[feat]],as.factor(TableSelected$Handle))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Selected_Features))
  tmp = dunn.test(TableSelected[[feat]],as.factor(TableSelected$Handle), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(Selected_Features))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)


write.table(table_results, "wanted_comparisons_selected_features_Tuj1PercentofAllNuc.csv", sep = ',', row.names = FALSE)

#Filtering for selecting images

TablesFiltered = TablesAll[TablesAll$Line=="K7",]
wantedcolumns = c(1,2,3,23,30)
TablesFiltered = TablesFiltered[wantedcolumns]
TablesFiltered = TablesFiltered[TablesFiltered$Condition=="Un",]

#Filtering for calculating number fields per condition

TablesFiltered = TablesAll[TablesAll$Category=="Patient",]
TablesFiltered = TablesFiltered[TablesFiltered$Condition=="10",]

wantedcolumns = c(1,2,3,23,30)
TablesFiltered = TablesFiltered[wantedcolumns]

