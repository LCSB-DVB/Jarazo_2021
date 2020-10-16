### Take plots already created and add the significance matrix (already calculated) at the bottom.
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(scales)
library(grid)

# color palettes
control_pat_col_fun = colorRamp2(c(+4, 0, -4), c("cornflowerblue", "gray90", "firebrick1"))
controls_col_fun = colorRamp2(c(-4, 0, +4), c("grey28", "gray90", "cornflowerblue"))
patients_col_fun = colorRamp2(c(-4, 0, +4), c("grey28", "gray90", "firebrick1"))
treat_control_pat_col_fun = colorRamp2(c(+4, 0, -4), c("yellow", "gray90", "#006E82"))
treat_controls_col_fun = colorRamp2(c(-4, 0, +4), c("grey28", "gray90", "yellow"))
treat_patients_col_fun = colorRamp2(c(-4, 0, +4), c("grey28", "gray90", "#006E82"))

control_pat_col_lgd = Legend(col_fun = control_pat_col_fun, title = "Controls_Patient")
controls_col_lgd = Legend(col_fun = controls_col_fun, title = "Controls_Controls")
patients_col_lgd = Legend(col_fun = patients_col_fun, title = "Patient_Patient")
treat_control_pat_col_lgd = Legend(col_fun = treat_control_pat_col_fun, title = "CD_Controls_Patient")
treat_controls_col_lgd = Legend(col_fun = treat_controls_col_fun, title = "CD_Controls_Controls")
treat_patients_col_lgd = Legend(col_fun = treat_patients_col_fun, title = "CD_Patient_Patient")
pd1 = packLegend(control_pat_col_lgd, controls_col_lgd, patients_col_lgd,
                direction = "horizontal")
pd2 = packLegend(treat_control_pat_col_lgd,treat_controls_col_lgd, treat_patients_col_lgd,
                 direction = "horizontal")

# save legends
pdf("significance_legends.pdf")
pushViewport(viewport(width = 1, height = 1))
grid.rect()
draw(pd1, x = unit(1, "cm"), y = unit(5, "cm"), just = c("left", "bottom"))
draw(pd2, x = unit(1, "cm"), y = unit(3, "cm"), just = c("left", "top"))
popViewport()
dev.off()


# split a tibble into list of tibbles based on group_by
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

# create for each matrix the associated annotation
tibble_to_annotation = function(tibble_df,type) {
  df = tibble_df %>% 
    pivot_wider(HandleCat,names_from=order,values_from=SM) %>% 
    select(sort(colnames(.))) %>% 
    column_to_rownames("HandleCat")
  if(type=="untreated"){
    ha = columnAnnotation(
      Controls_Patient = as.numeric(df["Controls_Patient",]),
      Controls_Controls = as.numeric(df["Controls_Controls",]),
      Patient_Patient = as.numeric(df["Patient_Patient",]),
      col = list(Controls_Patient = control_pat_col_fun,
                 Controls_Controls = controls_col_fun,
                 Patient_Patient = patients_col_fun),
      na_col = "white",
      annotation_name_side ="left",
      annotation_name_gp = gpar(fontsize=4),
      gp = gpar(col = "white"),
      width = unit(1, "cm"),
      annotation_height = unit(0.1, "cm"))
  } else {
    ha = columnAnnotation(
      Controls_Patient = as.numeric(df["Controls_Patient",]),
      Controls_Controls = as.numeric(df["Controls_Controls",]),
      Patient_Patient = as.numeric(df["Patient_Patient",]),
      col = list(Controls_Patient = treat_control_pat_col_fun,
                 Controls_Controls = treat_controls_col_fun,
                 Patient_Patient = treat_patients_col_fun),
      na_col = "white",
      annotation_name_side ="left",
      annotation_name_gp = gpar(fontsize=4),
      gp = gpar(col = "white"),
      width = unit(1, "cm"),
      annotation_height = unit(0.1, "cm"))
  }
  return(ha)
}

# create annotation for the significance matrix
create_anno = function(comparison_file,type){
  sig_matrix = read.csv(comparison_file)
  
  minimal_sig_matrix = sig_matrix %>%
    select(feat,HandleCat,order,SM) 
  
  matrices_list = named_group_split(minimal_sig_matrix,feat)
  anno_list = lapply(matrices_list,tibble_to_annotation,type)
  return(anno_list)
}

# create plots with Javi's code wrapped in a function
create_data_plots = function(rdata_file,function_file){
  Table_all = readRDS(rdata_file)
  source(function_file)
  plots = create_plots(Table_all)
}

# plot the graphs and the significance matrix together
plot_together = function(rdata_file,function_file,comparison_file,type){
  plots = create_data_plots(rdata_file,function_file)
  graphs = plots$graphs
  leg = plots$legend
  
  anno = create_anno(comparison_file,type)
  
  i=1
  x_start=x=0
  y=.8
  xmax=.8
  y_padding = .02
  xmin=ymin=0
  plot_wide=plot_high=anno_wide=.15
  anno_high=0.05
  p=ggdraw()
  
  while (i<=length(graphs)){
    this_feature = names(graphs)[i]
    this_plot = graphs[[this_feature]]
    this_anno = anno[[this_feature]]
    
    p = p +
      draw_plot(this_plot, x = x, y = y, width = plot_wide, height = plot_high) +
      draw_plot(grid.grabExpr(draw(this_anno)), x = x, y = y-anno_high-y_padding, width = anno_wide, height = anno_high)
    
    x = x+plot_wide
      
    if (x>=xmax){
      x=x_start
      y = y-plot_high-anno_high-y_padding
    }
    i=i+1
  }
  
  p1= p+ draw_plot(leg, x = .9, y = .1, width = .00000001, height = .00000001)
  
  out_file = gsub(".rdata","_annotated.pdf",gsub("Table_all_","",rdata_file))
  ggsave(out_file, plot=p1, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)
  
}


### RUN ###
plot_together("Table_all_Autophagy_modul_S1.rdata",
              "Autophagy_modul_S1_function.R",
              "KW_wanted_comparisons_All_features_Automodul_S1_SM.csv",
              "untreated")
plot_together("Table_all_Autophagy_modul_S3.rdata",
              "Autophagy_modul_S3_function.R",
              "KW_wanted_comparisons_All_features_Automodul_S3_SM.csv",
              "untreated")

plot_together("Diff_autotag_comp_screen.rdata",
              "Diff_autotag_compscreen_function.R",
              "KW_wanted_comparisons_All_features_diff_autotag_SM_fixed.csv",
              "treated")

plot_together("Table_all_CD_untag_S1.rdata",
              "CD_untag_S1_function.R",
              "KW_wanted_comparisons_All_features_CD_untag_S1_SM.csv",
              "treated")
plot_together("Table_all_CD_untag_S2.rdata",
              "CD_untag_S2_function.R",
              "KW_wanted_comparisons_All_features_CD_untag_S2_SM.csv",
              "treated")
plot_together("Table_all_CD_untag_S3.rdata",
              "CD_untag_S3_function.R",
              "KW_wanted_comparisons_All_features_CD_untag_S3_SM.csv",
              "treated")

