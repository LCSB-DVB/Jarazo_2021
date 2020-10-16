create_plots = function(Table_all){ 

	###################################################
	## Violin with line through mnedian of timepoints #
	###################################################
	Table_all %>% 
	  dplyr::filter(TimePoint==5) ->Table_all

	Features = c("GroupCount","sum_Area","GroupCountNorm","sum_AreaNorm")

	levels_x_axis = c("Un","500","1","5","10","50")
	Limit_y_axis = data.table(Features = Features, X1= c(0,0,0,0), X2= c(600,15000,0.005,0.1))   

	graph = list()
	for(feat in Features){
	  lmmin=Limit_y_axis[Features==feat,X1]
	  lmmax=Limit_y_axis[Features==feat,X2]
	  graph[[feat]]=ggplot(Table_all,aes_string(y=feat, x="Concentration", fill="Category")) + 
		geom_violin(aes(group=interaction(Category,Concentration)),alpha=0.4,position=position_dodge(width = 0.6),lwd=0.1) + 
		geom_boxplot(aes(group=interaction(Category,Concentration)), fill="white",position=position_dodge(width = 0.6),outlier.size = 0.1,lwd=0.1,width=0.3) +
		scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
		stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
		scale_color_manual(values=c("cornflowerblue","firebrick1")) +
		# scale_y_continuous (breaks = extended_breaks(lmmin,lmmax,m = 4,only.loose = FALSE), expand = c(0,0)) +
		scale_y_continuous (limits = c(lmmin,lmmax),breaks = breaks_extended(n = 4),expand = c(0,0)) +
		scale_x_discrete(limits = levels_x_axis) +
		theme_classic() +
		theme(plot.title= element_text(hjust=0.5, size=8),axis.text.x= element_text(size=5),axis.text.y= element_text(size=5),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
		ggtitle(feat)
	}

	#Get the legend as a separate plot
	legend_plot = ggplot(Table_all,aes(y=Areaofcells, x=Concentration, fill=Category)) + 
	  geom_violin(aes(group=interaction(Category,Concentration)),alpha=0.4,position=position_dodge(width = 6),lwd=0.1) + 
	  geom_boxplot(aes(group=interaction(Category,Concentration)), fill="white",position=position_dodge(width = 6),outlier.size = 0.3,lwd=0.3,width=0.8) +
	  scale_fill_manual(values=c("cornflowerblue","firebrick1")) +
	  stat_summary(fun.y=median, geom="line", aes(color=Category, group=factor(Category)), size=0.5, linetype="dashed")+
	  scale_color_manual(values=c("cornflowerblue","firebrick1")) +
	  theme(legend.background=element_rect(fill = "transparent", colour = NA)) 

	my_legend = get_legend(legend_plot)

  return(list(graphs=graph, legend=my_legend))
}