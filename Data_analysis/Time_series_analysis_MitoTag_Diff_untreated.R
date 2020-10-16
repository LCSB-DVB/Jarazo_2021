library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(kSamples)
library(npIntFactRep)
library(nparLD)
library(tidyverse)
library(foreach)
library(cowplot)
library(ggpubr)
library(scales)

##LOADING OF DATA AND PREPARING TABLE
#Open the .csv and save as Text 

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Mitotag_Compscreen/ATP5C1/JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_1/Analysis_20171227_195919/Objects.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N2=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Mitotag_Compscreen/ATP5C1/JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_2/Analysis_20171227_201349/Objects.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N3=read.table("S:/HCS_Platform/Data/JavierJarazo/Diff_Mitotag_Compscreen/ATP5C1/JJ_20171206_60X_NESC_Mito_Lyso_tags_DIFF_3/Analysis_20171227_202837/Objects.txt",header=TRUE, sep='/t',stringsAsFactors = FALSE)

TablesAll = rbind(N1,N2,N3)
TablesAll$All_Events = TablesAll$GroupCount_ObjectsMitochondria + TablesAll$GroupCount_ObjectsMitophagy
TablesAll$MitophagyEventsByMitochondriaEvents = TablesAll$GroupCount_ObjectsMitophagy/TablesAll$GroupCount_ObjectsMitochondria
TablesAll$MitophagyAreaByMitochondriaArea = TablesAll$sum_Area_ObjectsMitophagy/TablesAll$sum_Area_ObjectsMitochondria
TablesAll$MitophagyAreaMeanByMitochondriaAreaMean = TablesAll$mean_Area_ObjectsMitophagy/TablesAll$mean_Area_ObjectsMitochondria
tmp = do.call(rbind,strsplit(TablesAll$AreaName,"_"))
TablesAll$line = tmp[,1]
TablesAll$Treatment = tmp[,2]
Patients = c("4C1")
Controls = c("EPI")
Category = c("Control", "Patient")
TablesAll$Category = ifelse(TablesAll$line %in% Patients, 'Patient', 'Control')
TablesAll$Condition = ifelse(TablesAll$Treatment=='DMSO', 'Untreated', 'Treated') 
TablesAll$Combined= paste(TablesAll$Treatment, TablesAll$Category, sep = "_")
Treatments = unique(TablesAll$Treatment)
TreatedConditions = setdiff(Treatments,"DMSO")
TimePoints = unique(TablesAll$TimePoint)
Features = colnames(TablesAll [,7:20])

TablesAll = data.frame(setDT(TablesAll)[!(sum_Area_ObjectsMitochondria %between% c(0,2000)),])

#########
## GRAPHS 
#########
#################################################################
##BOX PLOT + DOTPLOT + LINE THROUGH MEDIAN ONLY P_DMSO V C_DMSO##
#################################################################
Features =c('mean_Area_ObjectsMitochondria','GroupCount_ObjectsMitophagy','sum_Area_ObjectsMitophagy','mean_Area_ObjectsMitophagy')
Limit_y_axis = data.table(Features = Features, X1= c(40,0,0,10), X2= c(200,200,10000,80))
graph = list()
TableSelected = TablesAll[TablesAll$Treatment=='DMSO',]
TableSelected$TimePoint= as.factor(TableSelected$TimePoint)
TableSelected = data.frame(TableSelected)

for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="TimePoint", fill="Category")) + 
    geom_violin(aes(group=interaction(Category,TimePoint)),position=position_dodge(width = 0.85),alpha=0.4,lwd=0.1) + 
    geom_boxplot(aes(group=interaction(Category,TimePoint)), position=position_dodge(width = 0.85),fill="transparent",outlier.size = 0.15,lwd=0.3,width=0.4) +
    scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
    stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("cornflowerblue","firebrick1")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}
#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=mean_Area_ObjectsMitochondria, x=TimePoint, fill=Category)) + 
  geom_violin(aes(group=interaction(Category,TimePoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(aes(group=interaction(Category,TimePoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
  stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("cornflowerblue","firebrick1")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .3, width = .3, height = .3) +
  draw_plot(graph[[2]], x = .3, y = .3, width = .3, height = .3) +
  draw_plot(graph[[3]], x = 0, y = 0, width = .3, height = .3) +
  draw_plot(graph[[4]], x = .3, y = 0, width = .3, height = .3) +
  draw_plot(my_legend, x = .64, y = .3, width = .00000001, height = .00000001)

ggsave("Diff_Mitotag_DMSO.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


#################################################################
##BOX PLOT + DOTPLOT + LINE THROUGH MEDIAN  P_DMSO V P_tr##
#################################################################
Features =c('mean_Area_ObjectsMitochondria','GroupCount_ObjectsMitophagy','sum_Area_ObjectsMitophagy','mean_Area_ObjectsMitophagy')
Limit_y_axis = data.table(Features = Features, X1= c(40,0,0,10), X2= c(200,200,10000,80))
TablesAll$TimePoint= as.factor(TablesAll$TimePoint)
TablesAll=data.frame(TablesAll)

graph = list()
TableCondi = TablesAll[TablesAll$Treatment=='Cy' | TablesAll$Treatment=='DMSO',]
TableSelected = TableCondi [TableCondi$Category=='Patient',]
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  TableSelected$Condition = factor(TableSelected$Condition, levels = c('Untreated','Treated'),ordered = TRUE)
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="TimePoint", fill="Condition")) + 
    geom_violin(aes(group=interaction(Condition,TimePoint)),position=position_dodge(width = 0.85),alpha=0.4,lwd=0.1) + 
    geom_boxplot(aes(group=interaction(Condition,TimePoint)), position=position_dodge(width = 0.85),fill="transparent",outlier.size = 0.15,lwd=0.3,width=0.4) +
    scale_fill_manual(values=c("firebrick1", "#006E82")) +
    stat_summary(fun.y=median, geom="line", aes(color=Condition, group=factor(Condition)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("firebrick1", "#006E82")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=mean_Area_ObjectsMitochondria, x=TimePoint, fill=Condition)) + 
  geom_violin(aes(group=interaction(Condition,TimePoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(aes(group=interaction(Condition,TimePoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("firebrick1", "#006E82")) +
  stat_summary(fun.y=median, geom="line", aes(color=Condition, group=factor(Condition)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("firebrick1", "#006E82")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .25, width = .25, height = .2) +
  draw_plot(graph[[2]], x = .25, y = .25, width = .25, height = .2) +
  draw_plot(graph[[3]], x = 0, y = 0, width = .25, height = .2) +
  draw_plot(graph[[4]], x = .25, y = 0, width = .25, height = .2) +
  draw_plot(my_legend, x = .64, y = .3, width = .00000001, height = .00000001)

ggsave("Diff_Mitotag_Cy.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

############################################################################
##BOX PLOT + DOTPLOT + LINE THROUGH MEDIAN  P_DMSO V P_tr + LINE of C_DMSO##
############################################################################
Features =c('mean_Area_ObjectsMitochondria','GroupCount_ObjectsMitophagy','sum_Area_ObjectsMitophagy','mean_Area_ObjectsMitophagy')
Limit_y_axis = data.table(Features = Features, X1= c(40,0,0,10), X2= c(200,200,10000,80))
TablesAll$TimePoint= as.factor(TablesAll$TimePoint)
TablesAll=data.frame(TablesAll)

graph = list()
TableCondi = TablesAll[TablesAll$Treatment=='Cy' | TablesAll$Treatment=='DMSO',]
TableSelected = TableCondi [TableCondi$Category=='Patient' |TableCondi$Category=='Control',]
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="TimePoint", fill="AreaName")) + 
    geom_violin(data = . %>% filter(AreaName %in% c('4C1_Cy','4C1_DMSO')),aes(group=interaction(AreaName,TimePoint)),position=position_dodge(width = 0.85),alpha=0.4,lwd=0.1) + 
    geom_boxplot(data = . %>% filter(AreaName %in% c('4C1_Cy','4C1_DMSO')),aes(group=interaction(AreaName,TimePoint)), position=position_dodge(width = 0.85),fill="transparent",outlier.size = 0.15,lwd=0.3,width=0.4) +
    scale_fill_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
    stat_summary(fun.y=median, geom="line", aes(color=AreaName, group=factor(AreaName)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=mean_Area_ObjectsMitochondria, x=TimePoint, fill=AreaName)) + 
  geom_violin(data = . %>% filter(AreaName %in% c('4C1_Cy','4C1_DMSO')),aes(group=interaction(AreaName,TimePoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(data = . %>% filter(AreaName %in% c('4C1_Cy','4C1_DMSO')),aes(group=interaction(AreaName,TimePoint)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
  stat_summary(fun.y=median, geom="line", aes(color=AreaName, group=factor(AreaName)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .25, width = .25, height = .2) +
  draw_plot(graph[[2]], x = .25, y = .25, width = .25, height = .2) +
  draw_plot(graph[[3]], x = 0, y = 0, width = .25, height = .2) +
  draw_plot(graph[[4]], x = .25, y = 0, width = .25, height = .2) +
  draw_plot(my_legend, x = .64, y = .3, width = .00000001, height = .00000001)

ggsave("Diff_Mitotag_Cy.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)  

############################################################################
##BOX PLOT + DOTPLOT + LINE THROUGH MEDIAN  P_DMSO V P_tr + C_DMSO V C_tr ##
############################################################################
Features =c('mean_Area_ObjectsMitochondria','GroupCount_ObjectsMitophagy','sum_Area_ObjectsMitophagy','mean_Area_ObjectsMitophagy')
Limit_y_axis = data.table(Features = Features, X1= c(40,0,0,10), X2= c(200,200,10000,80))
TablesAll$TimePoint= as.factor(TablesAll$TimePoint)
TablesAll=data.frame(TablesAll)

graph = list()
TableCondi = TablesAll[TablesAll$Treatment=='Cy' | TablesAll$Treatment=='DMSO',]
TableSelected = TableCondi [TableCondi$Category=='Patient' |TableCondi$Category=='Control',]
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="TimePoint", fill="AreaName")) + 
    geom_violin(aes(group=interaction(AreaName,TimePoint)),position=position_dodge(width = 0.85),alpha=0.4,lwd=0.1) + 
    geom_boxplot(aes(group=interaction(AreaName,TimePoint)), position=position_dodge(width = 0.85),fill="transparent",outlier.size=0.01,outlier.alpha = 0.2,lwd=0.1,width=0.4) +
    scale_fill_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
    stat_summary(fun.y=median, geom="line", aes(color=AreaName, group=factor(AreaName)), size=0.5, linetype="dashed")+
    scale_color_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=mean_Area_ObjectsMitochondria, x=TimePoint, fill=AreaName)) + 
  geom_violin(aes(group=interaction(AreaName,TimePoint)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(aes(group=interaction(AreaName,TimePoint)), fill="white",position=position_dodge(width = 6),outlier.size=0.01,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
  stat_summary(fun.y=median, geom="line", aes(color=AreaName, group=factor(AreaName)), size=0.5, linetype="dashed")+
  scale_color_manual(values=c("#006E82","firebrick1","yellow","cornflowerblue")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .25, width = .25, height = .2) +
  draw_plot(graph[[2]], x = .25, y = .25, width = .25, height = .2) +
  draw_plot(graph[[3]], x = 0, y = 0, width = .25, height = .2) +
  draw_plot(graph[[4]], x = .25, y = 0, width = .25, height = .2) +
  draw_plot(my_legend, x = .64, y = .3, width = .00000001, height = .00000001)

ggsave("Diff_Mitotag_Cy_all_groups.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)  


#############
####STATS#### 
#############
# Preparation of tables for analysis
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
TablesAll$subj = sprintf("%s_%s_%s_%s",TablesAll$Column,TablesAll$Row,TablesAll$Field,substrRight(TablesAll$Barcode,1))
dic = TablesAll %>% group_by(subj) %>% summarize(treat = unique(Treatment),cat = unique(Category), con = unique(Condition), comb = unique(Combined))
TablesComplete = TablesAll %>% complete(subj,TimePoint)

nas = which(is.na(TablesComplete$Treatment))
for(n in nas){
  TablesComplete$Treatment[n] = as.character(dic[dic$subj ==TablesComplete$subj[n],"treat"])
  TablesComplete$Category[n] = as.character(dic[dic$subj ==TablesComplete$subj[n],"cat"])
  TablesComplete$Condition[n] = as.character(dic[dic$subj ==TablesComplete$subj[n],"con"])
  TablesComplete$Combined[n] = as.character(dic[dic$subj ==TablesComplete$subj[n],"comb"])
}
assign_significance = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 'ns'
  significance[pvals <0.05 & pvals >0.01] = '*' 
  significance[pvals <=0.01 & pvals >0.001] = '**'
  significance[pvals <=0.001 & pvals >0.0001] = '***'
  significance[pvals <=0.0001] = '****'
  return(significance)
}

#### Means and medians calculation per feature per treatment per category of line
#TablesAll = TablesAll[complete.cases(TablesAll),] # For calculating the means and medians use this run not the TablesComplete
means_treated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
means_untreated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
medians_treated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
medians_untreated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))

means_treated = list()
means_untreated = list()
medians_treated = list()
medians_untreated = list()

for(cat in Category){
  TableCat = TablesAll[TablesAll$Category==cat,]
  for (tm in TimePoints){
    TableTm = TableCat[TableCat$TimePoint==tm,]
   for(tr in TreatedConditions){
    TableSelected = TableTm[TableTm$Treatment=='DMSO'| TableTm$Treatment==tr,]
    for(feat in Features){
      tmp = aggregate(TableSelected[[feat]],list(TableSelected$Condition),mean)
      tmp2 = aggregate(TableSelected[[feat]],list(TableSelected$Condition),median)
      means_untreated [feat,tr,tm] = tmp$x[tmp$Group.1=="Untreated"]
      means_treated[[feat]][[tr]][[tm]] = tmp$x[tmp$Group.1=="Treated"]
      medians_untreated[[feat]][[tr]][[tm]] = tmp2$x[tmp2$Group.1=="Untreated"]
      medians_treated[[feat]][[tr]][[tm]] = tmp2$x[tmp2$Group.1=="Treated"]
     }
   }
  }
  means_untreated = do.call(rbind,unlist(means_untreated[[feat]][[tr]][[tm]],recursive=FALSE))
  write.csv(means_untreated, sprintf("%s_means_untreated_All_lines.csv",cat))
  write.csv(means_treated, sprintf("%s_means_treated_All_lines.csv",cat))
  write.csv(medians_untreated, sprintf("%s_medians_untreated_All_lines.csv",cat))
  write.csv(medians_treated, sprintf("%s_medians_treated_All_lines.csv",cat))
  }

###nparLD###
# Comparison DMSO_CONTROL V DMSO_PATIENT
TableSelected = TablesComplete[TablesComplete$Treatment=='DMSO',]

res_list = list()
for(feat in Features){
  res = f1.ld.f1(TableSelected[[feat]],time=TableSelected$TimePoint,group=TableSelected$Category,subject = TableSelected$subj,group.name = "Category")
  res_list[[feat]] = res$ANOVA.test.mod.Box
}
tab = as.data.frame(do.call(rbind,res_list))

tab$sig = assign_significance(as.double(tab$"p-value"))
tab$Feature = Features

write.csv(tab, "nparLD_DMSO.csv")

# Comparison TREATMENT_PATIENT V DMSO_PATIENT
TableSelected = TablesComplete[TablesComplete$Category=='Patient',]
statis_tr=list()
for(tr in TreatedConditions){
  
  for(feat in Features){
    res = f1.ld.f1(TableSelected[[feat]],time=TableSelected$TimePoint,group=TableSelected$Condition,subject = TableSelected$subj,group.name = "Condition")
    statis_tr[[tr]][[feat]] = res$ANOVA.test.mod.Box
  }
  
  tab1 = as.data.frame(do.call(rbind,statis_tr[[tr]]))
  
  tab1$sig = assign_significance(as.double(tab1$"p-value"))
  tab1$Feature = Features
  
  write.csv(tab1, sprintf("nparLD_P-t_vs_P-un_%s.csv",tr))
}

# Comparison TREATMENT_PATIENT V DMSO_CONTROL 
for(tr in TreatedConditions){
  TableSelected = TablesComplete[TablesComplete$Combined %in% c('DMSO_Control',sprintf('%s_Patient',tr)),]
  res_list = list()
  for(feat in Features){
    res = f1.ld.f1(TableSelected[[feat]],time=TableSelected$TimePoint,group=TableSelected$Combined,subject = TableSelected$subj,group.name = "Combined")
    res_list[[feat]] = res$ANOVA.test.mod.Box
  }
  tab2 = as.data.frame(do.call(rbind,res_list))
  
  tab2$sig = assign_significance(as.double(tab2$"p-value"))
  tab2$Feature = Features
  
  write.csv(tab2, sprintf("nparLD_P-t_vs_C-un_%s.csv",tr))
}

# Comparison of all the treatments in a DMSO_Control v DMSO_Patient v tr_Control v tr_Patient manner 

for(tr in TreatedConditions){
  out= list()
  for(feat in Features){
    TableSelected = TablesComplete[TablesComplete$Treatment=='DMSO'| TablesComplete$Treatment==tr,]
    res = f2.ld.f1(TableSelected[[feat]],time=TableSelected$TimePoint,group1=TableSelected$Category,group2 = TableSelected$Treatment,subject = TableSelected$subj,group1.name = "Category",group2.name = "Treatment")
    tmp = res$ANOVA.test.mod.Box
    out[[feat]] = cbind(tmp,sig = assign_significance(tmp$"p-value"),feature=feat,treatment=tr)
    
  }
  sta= do.call(rbind,out)
  write.csv(sta, sprintf("nparLD_All_compar%s.csv",tr))
}

##Kolmogrov Smirnov##
####Kolmogorov-Smirnov for comparing the different distribution of the selected pairs combinations for each feature. The several tests performed were corrected
####by using p adjustment of Benjamini-Hochberg
ksPvals = list()
TablesAll= data.table(TablesAll)
for(tr in TreatedConditions) {
  for (tm in TimePoints){
    TableSelected = TablesAll[(TablesAll$Treatment=='DMSO'| TablesAll$Treatment==tr) & TablesAll$TimePoint == tm,]
    for(feat in Features){
      p1 = ks.test(as.matrix(TableSelected[Category == "Control" & Condition == "Untreated",..feat]), as.matrix(TableSelected[Category == "Patient" & Condition == "Untreated",..feat]), alternative = "two.sided")$p.value # "Un_Control" "Un_Patient"
      p2 = ks.test(as.matrix(TableSelected[Category == "Patient" & Condition == "Untreated",..feat]), as.matrix(TableSelected[Category == "Patient" & Condition == "Treated",..feat]), alternative = "two.sided")$p.value #"Un_Patient" "TU_Patient"
      p3 = ks.test(as.matrix(TableSelected[Category == "Control" & Condition == "Untreated",..feat]), as.matrix(TableSelected[Category == "Patient" & Condition == "Treated",..feat]), alternative = "two.sided")$p.value #"Un_Control" "TU_Patient"

      n1 = paste0(unique(TableSelected[Category == "Control" & Condition == "Untreated", Combined]), " - ", unique(TableSelected[Category == "Patient" & Condition == "Untreated", Combined]))
      n2 = paste0(unique(TableSelected[Category == "Patient" & Condition == "Untreated", Combined]), " - ", unique(TableSelected[Category == "Patient" & Condition == "Treated", Combined]))
      n3 = paste0(unique(TableSelected[Category == "Control" & Condition == "Untreated", Combined]), " - ", unique(TableSelected[Category == "Patient" & Condition == "Treated", Combined]))

      ksPvals[[tr]][[feat]][[tm]] = data.frame(feature=feat, treatment=tr, day=tm, names = c(n1,n2,n3), pvals3=p.adjust(c(p1,p2,p3),method="BH",n=3),pvals42=p.adjust(c(p1,p2,p3),method="BH",n=3*length(TimePoints)))
    }
  }
  
}
for(tr in TreatedConditions) {
  x = do.call(rbind,unlist(ksPvals[[tr]],recursive=FALSE))
  x$pvals = x$pvals42
  x$pvals42 = NULL
  x$sig = assign_significance(x$pvals)
  write.csv(x,sprintf("time_point_%s_adj42.csv",tr),row.names=FALSE)
  
  y = unlist(ksPvals[[tr]],recursive=FALSE)
  correct = sapply(y, function(a){
    all((a$pvals42<0.05) == c(TRUE,TRUE,FALSE))
  })
  write.table(cbind(names(y)[correct]),sprintf("significant_days_%s_adj42.csv",tr),row.names=FALSE,col.names=FALSE)
}
for(tr in TreatedConditions) {
  x = do.call(rbind,unlist(ksPvals[[tr]],recursive=FALSE))
  x$pvals = x$pvals3
  x$pvals3 = NULL
  x$sig = assign_significance(x$pvals)
  write.csv(x,sprintf("time_point_%s_adj3.csv",tr),row.names=FALSE)
  
  y = unlist(ksPvals[[tr]],recursive=FALSE)
  correct = sapply(y, function(a){
    all((a$pvals3<0.05) == c(TRUE,TRUE,FALSE))
  })
  write.table(cbind(names(y)[correct]),sprintf("significant_days_%s_adj3.csv",tr),row.names=FALSE,col.names=FALSE)
}
