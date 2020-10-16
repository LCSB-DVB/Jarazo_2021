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
library(dunn.test)

###############
## Load data ##
###############
N1=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_7_S2_R1/Analysis_20200703_024436/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N2=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_7_S2_R2/Analysis_20200703_031140/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N3=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_7_S2_R3/Analysis_20200703_033746/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N4=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_14_S2_R1/Analysis_20200703_040510/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N5=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_14_S2_R2/Analysis_20200703_043204/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N6=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200619_10X_Diff_Staging_Day_14_S2_R3/Analysis_20200703_045217/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N7=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200624_10X_Diff_Staging_Day_21_S2_R1/Analysis_20200703_050915/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N8=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200624_10X_Diff_Staging_Day_21_S2_R2/Analysis_20200703_052552/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
N9=read.table('S:/HCS_Platform/Data/JavierJarazo/Staging/Revision/Astro_activ/JJ_20200624_10X_Diff_Staging_Day_21_S2_R3/Analysis_20200703_054229/Objects.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

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
  dplyr::filter(NucVol >20000)->Table_all

###########################
## Saving_tables_heatmap ##
###########################
Table_all %>% 
  dplyr::filter(AreaName!="826") %>%
  group_by(Barcode,Timepoint,Category) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE) ->toto2
#write.csv(toto2, "P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/Table_all_astro_percent_2.csv",row.names=FALSE)

interesting_features = c("THbyNucVol","THbyTHVol", "S100bbyS100bVol",
                         "S100bbyNucVol","GFAPbyGFAPVol","GFAPbyNucVol",
                         "THVolByS100bVol", "S100bVolByNucVol","THVolByNucVol",  
                         "GFAPVolByNucVol", "GFAPVolByTHVol",  "GFAPVolByS100bVol",         
                         "ColoGFAPS100bVolByNucVol"  , "ColoGFAPS100bVolByS100bVol","S100bVol", 
                         "THVol",     "NucVol",    "GFAPVol",  
                         "ColoGFAPS100bVol")

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

saveRDS(All_tidy,"All_tidy_Staging_astro_active.rdata")


  ###########################
  ## GRAPHS for figure size##
  ###########################
  ###################################################
  ## Violin with line through mnedian of timepoints #
  ###################################################
  Table_all %>% 
    dplyr::filter(Category!="Corrected") %>% 
    dplyr::filter(AreaName!="826" | AreaName!="826_GC")->Table_all
  
  Features = c("GFAPVolByNucVol","S100bVolByNucVol","ColoGFAPS100bVolByNucVol","GFAPVol","S100bVol","NucVol")
  
  Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0,0,20000), X2= c(0.5,4,0.2,40000,150000,100000))   
  graph = list()
  for(feat in Features){
    lmmin=Limit_y_axis[Features==feat,X1]
    lmmax=Limit_y_axis[Features==feat,X2]
    graph[[feat]]=ggplot(Table_all,aes_string(y=feat, x="Timepoint", fill="Category")) + 
      geom_violin(aes(group=interaction(Category,Timepoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
      geom_boxplot(aes(group=interaction(Category,Timepoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.1,lwd=0.1,width=0.8) +
      scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
      stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
      scale_color_manual(values=c("cornflowerblue","firebrick1")) +
      scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 4),expand = c(0,0)) +
      scale_x_discrete(limits = c(7,14,21)) +
      theme_classic() +
      theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
      ggtitle(feat)
  }
  
  #Get the legend as a separate plot
  legend_plot = ggplot(Table_all,aes(y=NucVol, x=Timepoint, fill=Category)) + 
    geom_violin(aes(group=interaction(Category,Timepoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
    geom_boxplot(aes(group=interaction(Category,Timepoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("cornflowerblue","firebrick1")) +
    theme(legend.background=element_rect(fill = "transparent", colour = NA)) 
  
  my_legend = get_legend(legend_plot)
  
  p=ggdraw() +
    draw_plot(graph[[1]], x = 0, y = .2, width = .15, height = .2) +
    draw_plot(graph[[2]], x = .15, y = .2, width = .15, height = .2) +
    draw_plot(graph[[3]], x = .3, y = .2, width = .15, height = .2) +
    draw_plot(graph[[5]], x = .15, y = 0, width = .15, height = .2) +
    draw_plot(graph[[4]], x = 0, y = 0, width = .15, height = .2) +
    draw_plot(graph[[6]], x = .3, y = 0, width = .15, height = .2) +
    draw_plot(my_legend, x = .63, y = .2, width = .00000001, height = .00000001)
  
  ggsave("Timeseries_astro_activ.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

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
  
  #Filter the values removing outliers
  Table_all %>% 
    dplyr::filter(GFAPByNucVol<0.5) %>% 
    dplyr::filter(S100bVolByNucVol<4) %>% 
    dplyr::filter(ColoGFAPS100bVolByNucVol<0.2) %>%
    dplyr::filter(GFAPVol<40000) %>%
    dplyr::filter(S100bVol<150000)->Table_all

  Table_all = data.table(Table_all)
  
  ###################
  ## Kruskal-Wallis #
  ###################
  #All features wanted comparisons
  Table_all$Handle = paste(Table_all$Category, Table_all$Timepoint, sep="_")
  
  wanted =c("Controls_7 - Patient_7", "Controls_14 - Patient_14","Controls_21 - Patient_21")
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
  
  write.table(table_results, "KW_wanted_comparisons_All_features_astro_active.csv", sep = ',', row.names = FALSE)
  
  
###########################
##Selecting median images #
###########################
  Features = c("GFAPVolByNucVol","S100bVolByNucVol","ColoGFAPS100bVolByNucVol","GFAPVol","S100bVol","NucVol")
  
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
  