library(RColorBrewer)
library(grid)
library(ggplot2)
library(patchwork)
library(data.table)
library(FSA)
library(dunn.test)

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/SysMedPD/WP1/1.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)

assign_significance = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 'ns'
  significance[pvals <0.05 & pvals >0.01] = '*' 
  significance[pvals <=0.01 & pvals >0.001] = '**'
  significance[pvals <=0.001 & pvals >0.0001] = '***'
  significance[pvals <=0.0001] = '****'
  return(significance)
}
Patients = c(2122,825,2124,826)
N1$Category = ifelse(N1$AreaName %in% Patients,'Patient','Control')
p1 = wilcox.test(THMaskNorm ~ Category, data=N1)$p.value
c1 = "Control - Patient"
mytable = data.frame(Comparison = c(c1), Feature=c("THMaskNorm"), pvalue = c(p1))
mytable$significance = assign_significance(mytable$pvalue)
write.table(mytable, "MW_OrganoPlate.csv", sep = ',', row.names = FALSE)
TablesAll =N1

#########
## GRAPHS
#########

Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0,0,0,0,0), X2= c(3000000,4000000,1,200000000,2000000000,1,1.5,3))
graph = list()
for(feat in Features){
  lmmin=Limit_y_axis[Features==feat,X1]
  lmmax=Limit_y_axis[Features==feat,X2]
  sig= all_sigs[[feat]]
  label = sig$dec 
  graph[[feat]]=ggplot(TablesAll, aes_string(y=feat, x="Category")) + 
    geom_point(aes(color=Category), size=4, alpha=0.25,  position = position_jitter(w = 0.2, h = 0)) +
    scale_color_manual(values=c("cornflowerblue","firebrick1")) + 
    stat_summary(fun.y=median,geom="point", shape= 95, size=5) +
    stat_summary(fun.y=quantile, fun.args=list(probs=c(0.25, 0.75)),geom="line") + 
    scale_y_continuous (limits = c(lmmin,lmmax), expand = c(0,0)) +
    theme_classic() + 
    theme(plot.title= element_text(hjust=0.5, size=10),plot.subtitle =element_text(hjust=0.5, size=10), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
    ggtitle(feat,subtitle=format(label))
}

Reduce(`+`, graph)
ggsave("Differentiation_OrganoPlate.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")
