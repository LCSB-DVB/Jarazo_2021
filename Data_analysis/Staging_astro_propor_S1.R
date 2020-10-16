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
N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_7_S1_R1/Analysis_20200702_173603/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_7_S1_R2/Analysis_20200702_181011/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_7_S1_R3/Analysis_20200702_183629/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_14_S1_R1/Analysis_20200702_191136/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_14_S1_R2/Analysis_20200702_220058/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200608_10X_Diff_Staging_Day_14_S1_R3/Analysis_20200702_224118/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N7=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200611_10X_Diff_Staging_Day_21_S1_R1/Analysis_20200702_232455/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N8=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200611_10X_Diff_Staging_Day_21_S1_R2/Analysis_20200703_000812/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N9=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_propor/JJ_20200611_10X_Diff_Staging_Day_21_S1_R3/Analysis_20200703_004859/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

Table_7 = rbind(N1,N2,N3)
Table_14 = rbind(N4,N5,N6)
Table_21 = rbind (N7,N8,N9)

Table_7$Timepoint = 7
Table_14$Timepoint = 14
Table_21$Timepoint = 21

Table_all = rbind(Table_7,Table_14,Table_21)

Controls = c("T129","K7","KOLF","68")
Patients = c("2122","825","826","4C1")

Table_all$Category =  ifelse(Table_all$AreaName %in% Patients, 'Patient', ifelse(Table_all$AreaName %in% Controls, 'Controls', 'Corrected'))

Table_all %>% 
  dplyr::filter(NucVol >10000)->Table_all

###########################
## Saving_tables_heatmap ##
###########################
Table_all %>% 
  dplyr::filter(AreaName!="826") %>% 
  group_by(Barcode,Timepoint,Category) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE) ->toto2
#write.csv(toto2, "P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_diff_percent_new_2.csv",row.names=FALSE)

interesting_features = c("THbyNucVol","THbyTHVol",
                         "Tuj1byTuj1Vol","Tuj1byNucVol", "GFAPbyGFAPVol","GFAPbyNucVol",
                         "THVolByTuj1Vol","Tuj1VolByNucVol","THVolByNucVol","GFAPVolByNucVol",       
                         "GFAPVolByTHVol","GFAPVolByTuj1Vol","ColoTHTuj1VolByNucVol","ColoTHTuj1VolByTuj1Vol",
                         "Tuj1Vol",   "THVol","NucVol",    "GFAPVol",  
                         "ColoTHTuj1Vol")

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

saveRDS(All_tidy,"All_tidy_Staging_astro_propor.rdata")

###########################
## GRAPHS for figure size##
###########################
###################################################
## Violin with line through mnedian of timepoints #
###################################################
  Table_all= Table_all[Table_all$AreaName!='826_GC' & Table_all$AreaName!='826',]
  Table_all %>% 
    filter(Timepoint==14)->Table_all
  
  Features = c("THVolByTuj1Vol","Tuj1VolByNucVol","THVolByNucVol","Tuj1Vol","THVol","NucVol")
  
  Limit_y_axis = data.table(Features = Features, X1= c(0,1,0,20000,0,25000), X2= c(1.5,6,4,200000,200000,100000))  
  graph = list()
  for(feat in Features){
    lmmin=Limit_y_axis[Features==feat,X1]
    lmmax=Limit_y_axis[Features==feat,X2]
    graph[[feat]]=ggplot(Table_all,aes_string(y=feat, x="Category", fill="Category")) + 
      geom_violin(alpha=0.4,position=position_dodge(width = 1),lwd=0.1) + 
      geom_boxplot(fill="white",position=position_dodge(width = 1),outlier.size = 0.3,lwd=0.3,width=0.2) +
      scale_fill_manual(values=c("cornflowerblue","lightgoldenrod1","firebrick1")) +
      # scale_y_continuous (breaks = extended_breaks(lmmin,lmmax,m = 4,only.loose = FALSE), expand = c(0,0)) +
      scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0),labels = comma) +
      theme_classic() +
      theme(plot.title= element_text(hjust=0.5, size=10),axis.text.x= element_text(size=5),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
      ggtitle(feat)
  }
  
  #Get the legend as a separate plot
  legend_plot = ggplot(Table_all,aes(y=NucVol, x=Category, fill=Category)) + 
    geom_violin(alpha=0.4,position=position_dodge(width = 1),lwd=0.1) + 
    geom_boxplot(fill="white",position=position_dodge(width = 1),outlier.size = 0.3,lwd=0.3,width=0.2) +
    scale_fill_manual(values=c("cornflowerblue","lightgoldenrod1","firebrick1")) +
    theme(legend.background=element_rect(fill = "transparent", colour = NA)) 
  
  my_legend = get_legend(legend_plot)
  
  p=ggdraw() +
    draw_plot(graph[[1]], x = 0, y = .2, width = .15, height = .2) +
    draw_plot(graph[[2]], x = .15, y = .2, width = .15, height = .2) +
    draw_plot(graph[[3]], x = .3, y = .2, width = .15, height = .2) +
    draw_plot(graph[[4]], x = .15, y = 0, width = .15, height = .2) +
    draw_plot(graph[[5]], x = 0, y = 0, width = .15, height = .2) +
    draw_plot(graph[[6]], x = .3, y = 0, width = .15, height = .2) +
    draw_plot(my_legend, x = .63, y = .2, width = .00000001, height = .00000001)
  
  ggsave("Isogenics.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)
  
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
  
  Table_all = data.table(Table_all)
  
  ###################
  ## Kruskal-Wallis #
  ###################
  #All features wanted comparisons
  Table_all$Handle = paste(Table_all$Category, Table_all$Timepoint, sep="_")
  
  wanted =c("Controls_14 - Patient_14", "Controls_14 - Corrected_14","Corrected_14 - Patient_14")
  results = list()
  for(feat in Features){
    kpvals = kruskal.test(Table_all[[feat]],as.factor(Table_all$Handle))$p.value
    kpvalsadj = p.adjust(kpvals,method="BH",n=length(Features))
    tmp = dunn.test(Table_all[[feat]],as.factor(Table_all$Handle), list=TRUE)
    pvals = tmp$P[match(wanted,tmp$comparisons)]
    padjusted =p.adjust(pvals,method="BH",n=length(wanted)*length(interesting_features))
    results[[feat]] = cbind(wanted,feat,kpvals,kpvalsadj,pvals,padjusted)
  }
  
  table_results = as.data.table(do.call(rbind,results))
  table_results$significancekw = assign_significance(as.double(table_results$kpvalsadj))
  table_results$significancedt = assign_significance(as.double(table_results$padjusted))
  table_results = dplyr::select(table_results,"wanted", "feat", "kpvals","kpvalsadj","significancekw","pvals","padjusted","significancedt")
  
  write.table(table_results, "KW_wanted_comparisons_All_features_astro_propor.csv", sep = ',', row.names = FALSE)
  