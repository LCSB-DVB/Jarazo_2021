library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(kSamples)
library(cowplot)

##LOADING OF DATA AND PREPARING TABLE
#Open the .csv and save as Text 

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_All_measurements_per_line_N1.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N2=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_All_measurements_per_line_N2.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N3=read.table("S:/HCS_Platform/Data/JavierJarazo/Seahorse/R_All_measurements_per_line_N3.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
TablesAll = rbind(N1,N2,N3)
TablesAll = TablesAll[complete.cases(TablesAll),]
Lines = unique(TablesAll$Group.Name)
Controls = c("34769 c1","K7 L","K7 WT","KOLF","T12.9 N")
Patients = c("2122","2124","825","826")
Isogenic = c("2122 GC")
TablesAll$Category = ifelse(TablesAll$Group.Name %in% Patients, 'Patient', ifelse(TablesAll$Group.Name %in% Isogenic, 'Isogenic', 'Control'))

#without isogenics
TablesAll= TablesAll[TablesAll$Group.Name!='2122 GC' & TablesAll$Group.Name!='34769 c1' ,]
Features = colnames(TablesAll [,3:10])

##########
## GRAPHS
##########

Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0,0,-10,0.5,0.2), X2= c(30,60,60,20,40,20,2,1.2))
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  plot=ggplot(TablesAll, aes_string(x='Category', y=feat,fill = 'Category')) + 
    geom_boxplot(alpha=0.3, outlier.color = 'grey',outlier.alpha = 0.5,outlier.size =0.5,lwd=0.25) + 
    scale_x_discrete(limits=levels_x_axis) +
    scale_fill_manual(values=c("cornflowerblue","lightgoldenrod1","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax), expand = c(0,0)) +
    ggtitle(feat)+
    theme_classic() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +  
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_blank(), axis.title.y=element_blank()) +
    theme(plot.title = element_text(size =12, hjust = 0.5))
  ggdraw () +
    draw_plot(plot, x = 0, y = 0, width = 0.2, height = 0.2)
  ggsave(sprintf("Seahorse_mitochondria_%s.pdf",feat), plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")
  
}

###########
## STATS ##
###########
#All features all comparisons
results = list()
for(feat in Features){
  kpvals = kruskal.test(TablesAll[[feat]],as.factor(TablesAll$Category))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
  tmp = dunn.test(TablesAll[[feat]],as.factor(TablesAll$Category), list=TRUE)
  pvals = tmp$P
  wanted = tmp$comparisons
  padjusted =p.adjust(pvals,method="BH",n=length(Features)*length(wanted))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)

#wanted comparisons selected features
TableSelected = TablesAll
Selected_Features = c("ATP.Production","Coupling.Efficiency") 
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
write.table(table_results, "wanted_comparisons_selected_features_Seahorse_mitochondria.csv", sep = ',', row.names = FALSE)
