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
library(cowplot)
library(ggpubr)
library(scales)

##LOADING OF DATA AND PREPARING TABLE
N1=read.table('S:/HCS_Platform/Data/JavierJarazo/TFEB/Neurons/JJ_20190802_20X_Neurons_TFEB_1/Analysis_20190803_104559/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/TFEB/Neurons/JJ_20190802_20X_Neurons_TFEB_2/Analysis_20190803_104807/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/TFEB/Neurons/JJ_20190802_20X_Neurons_TFEB_3/Analysis_20190803_105010/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

TablesAll = rbind(N1,N2,N3)
Features = c("NucArea","NucAliveNorm","NucDeadNorm","TFEBArea","TFEBAreaNorm","CountTFEBObj","CountTFEBObjNorm","TFEBIntenNorm","ColoNucMaskArea","ColoNucMaskNorm","CountColoNucObj","CountColoNucObjNorm","ColoNucIntenNorm","ColoTHMaskArea","ColoTHMaskNorm","CountColoTHObj","CountColoTHObjNorm","ColoTHIntenNorm")
Coordinates= read.table('P:/Experiments/Extra_data_PINK1_paper/TFEB/96_target_coordinates.csv',header=TRUE, sep=',', stringsAsFactors = FALSE)
TablesAll = merge (TablesAll, Coordinates, by=c("ROW", "COL"))
Patients = c("2122","825","4C1")
Controls = c("K7","T129","KOLF")
TablesAll$Category = ifelse(TablesAll$Line %in% Controls, 'Control', 'Patient')
TablesAll$Handle = paste(TablesAll$Category, TablesAll$Condition, sep="_")

###########
## Plots ##
###########
levels_x_axis = c("Un","500","1","5","10")
Features = c("CountColoNucObjNorm")
Limit_y_axis = data.table(Features = Features, X1= c(0), X2= c(1))
TableSelected = TablesAll[TablesAll$ColoTHMaskArea != 0,]

graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="Condition", fill="Category")) + 
    geom_violin(aes(group=interaction(Category,Condition)),alpha=0.4,position=position_dodge(width = 0.8)) + 
    geom_boxplot(aes(group=interaction(Category,Condition)), fill="white",position=position_dodge(width = 0.8),outlier.size = 1,width=0.2) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("cornflowerblue","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    scale_x_discrete(limits=levels_x_axis) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=CountColoNucObjNorm, x=Condition, fill=Category)) + 
  geom_violin(aes(group=interaction(Category,Condition)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(aes(group=interaction(Category,Condition)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
  stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("cornflowerblue","firebrick1")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = 0, width = .25, height = .25) +
  draw_plot(my_legend, x = .45, y = .2, width = .00000001, height = .00000001)

ggsave("TFEB_neurons_per_category.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

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
TableSelected = TablesAll[TablesAll$ColoTHMaskArea != 0,]
Selected_Features = c("CountColoNucObjNorm","CountColoTHObjNorm") 
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


write.table(table_results, "wanted_comparisons_selected_features_TFEB_Neurons.csv", sep = ',', row.names = FALSE)

#Filtering for selecting images
TablesFiltered = TablesAll[TablesAll$Line=="2122",]
TablesFiltered = TablesFiltered[TablesFiltered$ColoTHMaskArea != 0,]
wantedcolumns = c(1:4,6,26,27,23,18)
TablesFiltered = TablesFiltered[wantedcolumns]
TablesFiltered = TablesFiltered[TablesFiltered$Condition=="Un",]

#Filtering for calculating number fileds per condition
TablesFiltered = TablesAll[TablesAll$Category=="Patient",]
TablesFiltered = TablesFiltered[TablesFiltered$Condition== "500",]
wantedcolumns = c(1:4,6,26,27,23,18)
TablesFiltered = TablesFiltered[wantedcolumns]
TablesFiltered = TablesFiltered[TablesFiltered$Condition=="Un",]


###########################
##Selecting median images #
###########################
Features = c("CountColoNucObjNorm")
metadata = c("ROW", "COL", "File","Barcode","AreaName","Field","COORDINATES",        
             "Line","Condition","Category","Handle")
TablesAll %>% 
  dplyr::select(metadata,CountColoNucObjNorm) %>% 
  dplyr::filter(Condition=="5")->Table_filtered 

Table_filtered %>% 
  group_by(Category) %>% 
  summarise(median= median(CountColoNucObjNorm)) %>% 
  ungroup() %>% 
  full_join(Table_filtered, by="Category") %>% 
  dplyr::filter(CountColoNucObjNorm>=median*0.9 & CountColoNucObjNorm<=median*1.1 ) ->Table_median