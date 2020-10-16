library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(kSamples)
library(tidyverse)
library(cowplot)
##LOADING OF DATA AND PREPARING TABLE
#Open the .csv and save as Text 

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_ECAR_All_measurements_per_line_N1_1column.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N2=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_ECAR_All_measurements_per_line_N2_1column.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N3=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_ECAR_All_measurements_per_line_N3_1column.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N1$Replicate = 1 
N2$Replicate = 2
N3$Replicate = 3 
TablesAll = rbind(N1,N2,N3)
TablesAll = TablesAll[!TablesAll=="",]
TablesAll = TablesAll[complete.cases(TablesAll),]
TablesAll = TablesAll[TablesAll$ECAR > 0,]
Lines = unique(TablesAll$Line)
Controls = c("34769 c1","K7 L","K7 WT","KOLF","T12.9 N"," T12.9 N")
Patients = c("2122","2124","825","826")
Isogenic = c("2122 GC")
TablesAll$Category = ifelse(TablesAll$Line %in% Patients, 'Patient', ifelse(TablesAll$Line %in% Isogenic, 'Isogenic', 'Control'))

Treatments = unique(TablesAll$Treatment)
##Preparation of table for calculating ECAR CHANGE Olygomycin-Basal
Table_diff = TablesAll %>%
  filter(Treatment %in% c("Basal","Oligomycin")) %>%
  group_by(Line, Well, Replicate,Category) %>%
  do(data.frame(diffOB = .$ECAR[.$Treatment=="Oligomycin"]- .$ECAR[.$Treatment=="Basal"]))

##########
## GRAPHS
##########
Limit_y_axis = data.table(Treatments = Treatments, X1= c(0,0,0,0), X2= c(50,60,30,40))
for(tr in Treatments){
  TableSelected = TablesAll[TablesAll$Treatment==tr,]
  lmmin=Limit_y_axis[Treatments==tr,X1]
  lmmax=Limit_y_axis[Treatments==tr,X2]
  plot=ggplot(TableSelected, aes_string(y='ECAR', x="Category", fill="Category")) + 
    geom_boxplot(alpha=0.3, outlier.color = 'grey',outlier.alpha = 0.5,outlier.size =0.5,lwd=0.25) + 
    #geom_point(size=0.1, alpha=0.5,  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), aes(group=Category)) +
    scale_x_discrete(limits=levels_x_axis) +
    scale_fill_manual(values=c("cornflowerblue","lightgoldenrod1","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax), expand = c(0,0)) +
    ggtitle(tr)+
    theme_classic() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +  
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_blank(), axis.title.y=element_blank()) +
    theme(plot.title = element_text(size =12, hjust = 0.5))
  ggdraw () +
    draw_plot(plot, x = 0, y = 0, width = 0.2, height = 0.2)
  ggsave(sprintf("Seahorse_glyco_%s.pdf",tr), plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")
  
}

Limit_y_axis = data.table(Features_OB = Features_OB, X1= c(0,0,0,0), X2= c(50,60,30,40))
Features_OB =c("diffOB")
for(feat in Features_OB){
  # lmmin=Limit_y_axis[Features_OB==feat,X1]
  # lmmax=Limit_y_axis[Features_OB==feat,X2]
  plot=ggplot(Table_diff, aes_string(y='diffOB', x="Category", fill="Category")) + 
    geom_boxplot(alpha=0.3, outlier.color = 'grey',outlier.alpha = 0.5,outlier.size =0.5,lwd=0.25) + 
    scale_x_discrete(limits=levels_x_axis) +
    scale_fill_manual(values=c("cornflowerblue","lightgoldenrod1","firebrick1")) +
    scale_y_continuous (limits = c(0,20), expand = c(0,0)) +
    ggtitle(feat)+
    theme_classic() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +  
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_blank(), axis.title.y=element_blank()) +
    theme(plot.title = element_text(size =12, hjust = 0.5))
  ggdraw () +
    draw_plot(plot, x = 0, y = 0, width = 0.2, height = 0.2)
  ggsave(sprintf("Seahorse_glyco_%s.pdf",feat), plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")
  
}

###########
## STATS ##
###########
#All features all comparisons
results = list()
for(tr in Treatments){
  TableSelected = TablesAll[TablesAll$Treatment==tr,]
  kpvals = kruskal.test(TableSelected$ECAR,as.factor(TableSelected$Category))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Treatments))
  tmp = dunn.test(TableSelected$ECAR,as.factor(TableSelected$Category), list=TRUE)
  pvals = tmp$P
  wanted = tmp$comparisons
  padjusted =p.adjust(pvals,method="BH",n=length(Treatments)*length(wanted))
  results[[tr]] = cbind(wanted,tr,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)

#wanted comparisons selected features
Selected_Treatments = c("Basal") 
wanted =c("Control - Patient", "Isogenic - Patient")
results = list()
for(tr in Selected_Treatments){
  TableSelected = TablesAll[TablesAll$Treatment==tr,]
  kpvals = kruskal.test(TableSelected$ECAR,as.factor(TableSelected$Category))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Selected_Treatments))
  tmp = dunn.test(TableSelected$ECAR,as.factor(TableSelected$Category), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(Selected_Treatments)*length(wanted))
  results[[tr]] = cbind(wanted,tr,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)
write.table(table_results, "wanted_comparisons_selected_features_Seahorse_glyco.csv", sep = ',', row.names = FALSE)

#wanted comparisons selected features
TableSelected = Table_diff
Selected_Features = c("diffOB") 
wanted =c("Control - Patient", "Isogenic - Patient")
results = list()
for(feat in Selected_Features){
  kpvals = kruskal.test(TableSelected[[feat]],as.factor(TableSelected$Category))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Selected_Features))
  tmp = dunn.test(TableSelected[[feat]],as.factor(TableSelected$Category), list=TRUE)
  pvals = tmp$P[match(wanted,tmp$comparisons)]
  padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(Selected_Features))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)
write.table(table_results, "wanted_comparisons_selected_features_Seahorse_glyco_diffOB.csv", sep = ',', row.names = FALSE)
