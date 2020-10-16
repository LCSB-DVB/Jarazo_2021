library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)
library(ggsignif)
library(ggpubr)
library(scales)
library(cowplot)

##LOADING OF DATA AND PREPARING TABLE
#Open the .csv and save as Text 
N1=read.table("S:/HCS_Platform/Data/JavierJarazo/Differentiation/JJ_20171215_10X_Diff_Parkin_N1/Analysis_20171226_120759/Objects.csv",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N2=read.table("S:/HCS_Platform/Data/JavierJarazo/Differentiation/JJ_20171215_10X_Diff_Parkin_N2/Analysis_20171226_121003/Objects.csv",header=TRUE, sep='/t',stringsAsFactors = FALSE)
N3=read.table("S:/HCS_Platform/Data/JavierJarazo/Differentiation/JJ_20171215_10X_Diff_Parkin_N3/Analysis_20171226_121054/Objects.csv",header=TRUE, sep='/t',stringsAsFactors = FALSE)

TablesAll=rbind(N1,N2,N3)
tmp = do.call(rbind,strsplit(TablesAll$AreaName,"_"))
TablesAll$line = tmp[,1]
TablesAll$Treatment = tmp[,2]
Lines=unique(TablesAll$line) 
Patients = c("29369C9","29369C16")
Controls = c("K7","T129N","KOLF","3C1","K7L")
TablesAll$Category = ifelse(TablesAll$line %in% Patients, 'Patient', 'Control')

TablesAll = TablesAll[TablesAll$NucVol != 0 & TablesAll$THbyTHVol != "NaN" ,]
TablesAll = setDT(TablesAll)[!(NucVol %between% c(0,10000)),]
TablesAll = setDT(TablesAll)[!(Tuj1Vol %between% c(0,50000)),]

Features = colnames(TablesAll)[c(7:16)]
Treatments=unique(TablesAll$Treatment)
Category = unique(TablesAll$Category)

TreatedConditions = c("TU","Cy100","F","C","J")
TablesAll$Condition = ifelse(TablesAll$Treatment=='DMSO', 'Untreated', 'Treated') 
TablesAll$Combined= paste(TablesAll$Treatment, TablesAll$Category, sep = "_")

TablesAll = TablesAll[TablesAll$Treatment=='DMSO'| TablesAll$Treatment=='Cy100',]
TableSelected = TablesAll[TablesAll$line=='K7'| TablesAll$line=='K7L'|TablesAll$line=='29369C9',]
TableSelected = TableSelected[TableSelected$Category !='Control',]
Features = c("THVolByTuj1Vol","Tuj1VolByNucVol","THVolByNucVol","Tuj1Vol","THVol","NucVol")

Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0,0,0), X2= c(1,4,1.5,150000,80000,100000)) #to have the same y axis as the median_q
graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  graph[[feat]]=ggplot(TableSelected,aes_string(y=feat, x="Combined", fill="Combined")) + 
    geom_violin(alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
    geom_boxplot(fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.3) +
    scale_fill_manual(values=c("#006E82","indianred2")) +
    scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 5),expand = c(0,0)) +
    scale_x_discrete(limits = c("DMSO_Patient","Cy100_Patient")) +
    theme_classic() +
    theme(plot.title= element_text(hjust=0.5, size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
    ggtitle(feat)
}

#Get the legend as a separate plot
legend_plot = ggplot(TableSelected,aes(y=NucVol, x=Combined, fill=Combined)) + 
  geom_violin(alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
  geom_boxplot(fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
  scale_fill_manual(values=c("#006E82","indianred2")) +
  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .2, width = .2, height = .2) +
  draw_plot(graph[[2]], x = .2, y = .2, width = .2, height = .2) +
  draw_plot(graph[[3]], x = .4, y = .2, width = .2, height = .2) +
  draw_plot(graph[[4]], x = .2, y = 0, width = .2, height = .2) +
  draw_plot(graph[[5]], x = 0, y = 0, width = .2, height = .2) +
  draw_plot(graph[[6]], x = .4, y = 0, width = .2, height = .2) +
  draw_plot(my_legend, x = .7, y = .2, width = .00000001, height = .00000001)

ggsave("Parkin_differentiation.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


#### Means and medians calculation per feature per treatment per category of line
means_treated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
means_untreated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
medians_treated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
medians_untreated = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))

for(cat in Category){
  TableCat = TablesAll[TablesAll$Category==cat,]
  for(tr in TreatedConditions){
    TableSelected = TableCat[TableCat$Treatment=='DMSO'| TableCat$Treatment==tr,]
    for(feat in Features){
      tmp = aggregate(TableSelected[[feat]],list(TableSelected$Condition),mean)
      tmp2 = aggregate(TableSelected[[feat]],list(TableSelected$Condition),median)
      means_untreated[feat,tr] = tmp$x[tmp$Group.1=="Untreated"]
      means_treated[feat,tr] = tmp$x[tmp$Group.1=="Treated"]
      medians_untreated[feat,tr] = tmp2$x[tmp2$Group.1=="Untreated"]
      medians_treated[feat,tr] = tmp2$x[tmp2$Group.1=="Treated"]
      write.csv(means_untreated, sprintf("%s_means_untreated_Parkin.csv",cat))
      write.csv(means_treated, sprintf("%s_means_treated_Parkin.csv",cat))
      write.csv(medians_untreated, sprintf("%s_medians_untreated_Parkin.csv",cat))
      write.csv(medians_treated, sprintf("%s_medians_treated_Parkin.csv",cat))
    }
    
  }
}

####Kruskal Wallis for comparing the different features of Treated Patients against Untreated Patients

pvalsCompoundEval = matrix(NA,nrow=length(Features),ncol=length(TreatedConditions),dimnames=list(Features,TreatedConditions))
postHocPvals = list()

for(tr in TreatedConditions){
  TableSelected = TablesAll[TablesAll$Treatment=='DMSO'| TablesAll$Treatment==tr,]
  for(feat in Features){
    pvalsCompoundEval[feat,tr] = kruskal.test(TableSelected[[feat]],as.factor(TableSelected$Combined))$p.value
    postHocPvals[[feat]][[tr]] = dunn.test(TableSelected[[feat]],as.factor(TableSelected$Combined),method="bh", altp=TRUE)
  }
  
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

all_sigs = list()
for(tr in TreatedConditions){
  for(feat in Features){
    x = postHocPvals[[feat]][[tr]]
    y = data.frame(comparisons = x$comparisons,altP.adjusted = x$altP.adjusted)
    y$dec = assign_significance(y$altP.adjusted)
    write.csv(y, sprintf("%s_%s_DunnPvalues_Parkin.csv",feat,tr),row.names=FALSE)
    all_sigs[[feat]][[tr]] = y
    
  }
  
}
