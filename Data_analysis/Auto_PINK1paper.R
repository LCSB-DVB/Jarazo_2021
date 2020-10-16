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
library(ggpubr)
library(scales)
##LOADING OF DATA AND PREPARING TABLE
#Open the .csv and save as Text 

TablesAll=read.table("S:/HCS_Platform/Data/JavierJarazo/Rosella/LC3/PINK1_paper/Summary_All.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
TablesAll_2=read.table("S:/HCS_Platform/Data/JavierJarazo/Rosella/LC3/PINK1_paper/ObjectsGrouped_1.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)

# barcodes to keep TablesAll
# "JA20161203_AllBasal_AutoMito_LC3-VPS35-PINK1" -> PINK1 -> "PINK1_basal"
# "JA_20161127_WTcontrol_AutoSensor" -> WT ->"basal_control"
# "JA_20161127_PINK1_AutoSensor"-> PINK1 -> "basal_control"
# "JA_20161203_Basal_WT_LRRK2" -> WT -> "WT_basal"
# Only barcode in TablesAll_2 -> "JA_20170225_autolines_003" -> AreaName -> "PINK1_basal_control" "WT_basal_control"

#Pre-processing TablesAll
Barcodes_PINK1 = c("JA20161203_AllBasal_AutoMito_LC3-VPS35-PINK1","JA_20161127_PINK1_AutoSensor")
Barcodes_WT = c("JA_20161127_WTcontrol_AutoSensor","JA_20161203_Basal_WT_LRRK2")
TablesAll_PINK1 = TablesAll[TablesAll$Barcode==Barcodes_PINK1,]
TablesAll_WT = TablesAll[TablesAll$Barcode==Barcodes_WT,]
AreaName_PINK1 = c("PINK1_basal", "basal_control" )
AreaName_WT = c("WT_basal", "basal_control" )
TablesAll_PINK1 = TablesAll_PINK1[TablesAll_PINK1$AreaName==AreaName_PINK1,]
TablesAll_WT = TablesAll_WT[TablesAll_WT$AreaName==AreaName_WT,]
TablesAll_PINK1$Category = "PINK1"
TablesAll_WT$Category = "WT"

#Pre-processing TablesAll_2
TablesAll_2 = TablesAll_2[TablesAll_2$AreaName=="PINK1_basal_control" | TablesAll_2$AreaName=="WT_basal_control",]
tmp = do.call(rbind,strsplit(TablesAll_2$AreaName,"_"))
TablesAll_2$Category = tmp[,1]
TableFinal = rbind(TablesAll_PINK1,TablesAll_WT,TablesAll_2)

###########
## Plots ##
###########
Features = c("sum_AutophagosomeDecision","sum_PhagophoreDecision","sum_LateAutolysosomeDecision","sum_EarlyAutolysosomeDecision")
levels_x_axis = c("WT","PINK1")

###################
# Boxplot dotplot # 
###################
##################
## per Category ##
##################
Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0), X2= c(12,750,750,75))
graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableFinal,aes_string(y=feat, x="Category", fill="Category")) + 
    geom_violin(alpha=0.4,position=position_dodge(width=0.6)) + 
    geom_boxplot(fill="white",outlier.size = 0.3,width=0.2,position=position_dodge(width=0.6)) +
    scale_x_discrete(limits=levels_x_axis) +
    scale_fill_manual(values=c("firebrick1","cornflowerblue")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 4),expand = c(0,0),labels = comma) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.text.x= element_text(size=8),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableFinal,aes(y=sum_AutophagosomeDecision, x=Category, fill=Category)) + 
  geom_violin(alpha=0.4,position=position_dodge(width = 1),lwd=0.1) + 
  geom_boxplot(fill="white",position=position_dodge(width = 1),outlier.size = 0.3,lwd=0.3,width=0.2) +
  scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[2]], x = 0, y = .2, width = .2, height = .2) +
  draw_plot(graph[[1]], x = .2, y = .2, width = .2, height = .2) +
  draw_plot(graph[[4]], x = 0, y = 0, width = .2, height = .2) +
  draw_plot(graph[[3]], x = .2, y = 0, width = .2, height = .2) +
  draw_plot(my_legend, x = .64, y = .2, width = .00000001, height = .00000001)

ggsave("Autophagy.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

###########
## STATS ##
###########
#All features all comparisons
results = list()
for(feat in Features){
  kpvals = kruskal.test(TableFinal[[feat]],as.factor(TableFinal$Category))$p.value
  kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
  tmp = dunn.test(TableFinal[[feat]],as.factor(TableFinal$Category), list=TRUE)
  pvals = tmp$P
  wanted = tmp$comparisons
  padjusted =p.adjust(pvals,method="BH",n=length(Features)*length(wanted))
  results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
}

table_results = do.call(rbind,results)

write.table(table_results, "All_comparisons_All_features_LC3_PINK1paper.csv", sep = ',', row.names = FALSE)
