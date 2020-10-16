library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(kSamples)
library(ggpubr)
library(scales)
library(cowplot)

##LOADING OF DATA AND PREPARING TABLE
N1=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d7_N1/Analysis_20171229_124627/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d7_N2/Analysis_20171229_124738/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d7_N3/Analysis_20171229_124752/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d14_N1/Analysis_20171229_124806/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d14_N2/Analysis_20171229_124820/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d14_N3/Analysis_20171229_124834/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N7=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d21_N1/Analysis_20171229_124848/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N8=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d21_N2/Analysis_20171229_124902/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N9=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d21_N3/Analysis_20171229_124916/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)
N10=read.table("S:/HCS_Platform/Data/JavierJarazo/Staging/Apop_JJ_20171219_20X_Diff_Apop_d42/Analysis_20171229_124929/Objects.txt",header=TRUE, sep=',',stringsAsFactors = FALSE)

N1$TimePoint = 7
N2$TimePoint= 7
N3$TimePoint = 7
N4$TimePoint = 14
N5$TimePoint = 14
N6$TimePoint = 14
N7$TimePoint = 21
N8$TimePoint = 21
N9$TimePoint = 21
N10$TimePoint = 42 #Not included because one replicate only
TablesAll = rbind(N1,N2,N3,N4,N5,N6,N7,N8,N9)
TablesAll = TablesAll[TablesAll$NucArea != 0,]
TablesAll = setDT(TablesAll)[!(NucArea %between% c(0,20000)),]
Patients = c("4C1","2122","2124","947")
Controls = c("K7","T129","KOLF")
TablesAll$Category = ifelse(TablesAll$AreaName %in% Patients, 'Patient', 'Control')
Features =c("PARPAreaNormTHArea")

###########################
## Saving_tables_heatmap ##
###########################
TablesAll %>% 
  dplyr::filter(AreaName!="2124") %>% #line removed for missing GC. Not included in the paper
  group_by(Barcode,TimePoint,Category) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE) ->toto2
#write.csv(toto2, "P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_apop_2.csv",row.names=FALSE)

Limit_y_axis = data.table(Features = Features, X1= c(0), X2= c(0.06))
graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TablesAll,aes_string(y=feat, x="TimePoint", fill="Category")) + 
    geom_violin(aes(group=interaction(Category,TimePoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
    geom_boxplot(aes(group=interaction(Category,TimePoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("cornflowerblue","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 4),expand = c(0,0)) +
    scale_x_discrete(limits = c(7,14,21)) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}
p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .2, width = .2, height = .2)
ggsave("Timeseries_Apoptosis.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

##STATS##  
 
  assign_significance = function(pvals){
    significance = rep(NA,length(pvals))
    significance[pvals >= 0.05] = 'ns'
    significance[pvals <0.05 & pvals >0.01] = '*' 
    significance[pvals <=0.01 & pvals >0.001] = '**'
    significance[pvals <=0.001 & pvals >0.0001] = '***'
    significance[pvals <=0.0001] = '****'
    return(significance)
  }

  ###################
  ## Kruskal-Wallis #
  ###################
  Features =c("PARPAreaNormTHArea")
  #All features wanted comparisons
  TablesAll %>% 
    dplyr::filter(AreaName!="2124") ->Table_all #line removed for missing GC. Not included in the paper
  Table_all$Handle = paste(Table_all$Category, Table_all$TimePoint, sep="_")
  
  wanted =c("Control_7 - Patient_7", "Control_14 - Patient_14","Control_21 - Patient_21")
  results = list()
  for(feat in Features){
    kpvals = kruskal.test(Table_all[[feat]],as.factor(Table_all$Handle))$p.value
    kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
    tmp = dunn.test(Table_all[[feat]],as.factor(Table_all$Handle), list=TRUE)
    pvals = tmp$P[match(wanted,tmp$comparisons)]
    padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(Features))
    results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
  }
  
  table_results = as.data.table(do.call(rbind,results))
  table_results$significancekw = assign_significance(as.double(table_results$kpvalsadj))
  table_results$significancedt = assign_significance(as.double(table_results$padjusted))
  table_results = dplyr::select(table_results,"wanted", "feat", "kpvals","kpvalsadj","significancekw","pvals","padjusted","significancedt")
  
  write.table(table_results, "KW_wanted_comparisons_Diff_staging_apop.csv", sep = ',', row.names = FALSE)
  