############################################################
# For PhyloAcc
# Boxplots
# Han Yan and Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(ape)
#library(tidyr)
library(ggbeeswarm)
#library(ggsignif)
#library(pROC)
#library(PRROC)
library(here)
source(here("scripts", "lib", "design.r"))
source(here("scripts", "lib", "read_results.r"))

## KEY
## nameLabel -> node_labels
## birdname -> label_map
## namePrint -> new_node_labels
## caseLabel -> scenarios
## elem_Acc -> st_elem
## elem_GT -> gt_elem
## res_backup -> beast_targets
## res_backup2 -> beast_non_targets

############################################################

ratite_tree_file = here("data", "ratite.tree")

label_map_file = here("data", "SimulatedData", "Name_Label_Mapping.txt")

node_labels = c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","cryCin","tinGut","eudEle", "notPer","anoDid","strCam",
                "aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme","cryCin.tinGut",
                "eudEle.notPer","cryCin.eudEle","cryCin.anoDid", "aptHaa.cryCin","aptHaa.strCam")

save_fig = F

############################################################

cat(as.character(Sys.time()), " | Fig10: Adjusting node labels...\n")

new_node_labels = adjustLabels(label_map_file, node_labels)

## Re-label the nodes based on the simulation letters
############################################################

cat(as.character(Sys.time()), " | Fig10: Assigning target lineages...\n")

# assign_targets_result = assignTargets(new_node_labels)
# target = assign_targets_result[[1]]
# nontarget = assign_targets_result[[2]]
# target_B = assign_targets_result[[3]]
# nontarget_B = assign_targets_result[[4]]

case_targets = assignTargets(node_labels, new_node_labels)

## Assign target lineages
############################################################

# cat(as.character(Sys.time()), " | Fig10: Reading likelihood files...\n")
# 
# st_elem_lik_file = here("results", "simulations", "Acc300elem_lik.txt")
# gt_elem_lik_file = here("results", "simulations", "GT300elem_lik.txt")
# # Results files from 300 loci simulate
# 
# st_elem_lik = readElemLik(st_elem_lik_file, subset_cols=TRUE)
# # Reading the species tree method results and getting adding the model with
# # the max probability
# 
# gt_elem_lik = readElemLik(gt_elem_lik_file)
# # Reading the gene tree method results and getting adding the model with
# # the max probability

############################################################

cat(as.character(Sys.time()), " | Fig10: Reading files and generating plots...\n")

case_labels = c("A", "B", "C")
## CASE LABELS DO NOT CORRESPOND TO PANEL LABELS
## 2 = case 2 = single paraphyletic = panel B
## 3 = case 3 = monophyletic = panel A
## 4 = case 4 = multi-paraphyletic = panel C

scenarios = c("8", "5", "6")
gt_extra = c("", "", "")
# This is how the files for the different cases are labeled...

method_order = c("PhyloAcc-GT", "PhyloAcc", "*BEAST2")
# Ordering for the legend


titles = c(expression(paste("Single acceleration, 3x", theta)), 
           expression(paste("Single acceleration, 6x", theta)),
           expression(paste("Single acceleration, 10x", theta)))

targ_p_list = list()
non_targ_p_list = list()
p_list = list()
# Initialize lists of plots

c = 2
#target_ind = 1

for(theta_scale in c(3,6,10)){
  
  cur_case_targets = case_targets[[c+1]]
  
  st_dir = here("results", "simulations", "Theta_Scaling", paste0("Theta", theta_scale), paste0("Result_Case", case_labels[c]), "PhyloAcc")
  gt_dir = here("results", "simulations", "Theta_Scaling", paste0("Theta", theta_scale), paste0("Result_Case", case_labels[c]))
  beast_dir = here("results", "simulations", "Theta_Scaling", paste0("Theta", theta_scale), paste0("Result_Case", case_labels[c]), "starBEAST2")
  # Input directories for the different run types
  
  st_elem_lik_file = paste0(st_dir, "/simu_200_100_2-", scenarios[c], "_elem_lik.txt")
  gt_elem_lik_file = paste0(gt_dir, "/simu_100_2-", scenarios[c], gt_extra[c], "_elem_lik.txt")
  # The elem_lik files for both PhyloAcc methods for the current scenario
  
  st_m1_postz_file = paste0(st_dir, "/", "simu_200_100_2-", scenarios[c], "_rate_postZ_M1.txt")
  st_m2_postz_file = paste0(st_dir, "/", "simu_200_100_2-", scenarios[c], "_rate_postZ_M2.txt")
  # Rate file paths for st
  
  gt_m1_postz_file = paste0(gt_dir, "/", "simu_100_2-", scenarios[c], gt_extra[c], "_rate_postZ_M1.txt")
  gt_m2_postz_file = paste0(gt_dir, "/", "simu_100_2-", scenarios[c], gt_extra[c], "_rate_postZ_M1.txt")
  # Rate file paths for gt
  
  beast_file = paste0(beast_dir, "/", "result1_2-", scenarios[c], ".RData", sep="")
  # File names based on the current scenario
  
  ## Some setup
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Reading likelihood files\n")
  
  print(st_elem_lik_file)
  st_elem_lik = readElemLik(st_elem_lik_file)
  # Reading the species tree method results and getting adding the model with
  # the max probability
  
  print(gt_elem_lik_file)
  gt_elem_lik = readElemLik(gt_elem_lik_file)
  # Reading the gene tree method results and getting adding the model with
  # the max probability
  
  ## Read likelihood files
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c],"- Theta", theta_scale, "X : Reading ST probs\n")
  
  #st_elem_lik_scenario = subset(st_elem_lik, No. == scenarios[c])
  st_opt_postz = readPostZ(st_m1_postz_file, st_m2_postz_file, st_elem_lik)
  
  # Read st results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Reading GT probs\n")

  #gt_elem_lik_scenario = subset(gt_elem_lik, No. == scenarios[c])
  print(gt_m1_postz_file)
  print(gt_m2_postz_file)
  gt_opt_postz = readPostZ(gt_m1_postz_file, gt_m2_postz_file, gt_elem_lik)
  
  # Read gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Adjusting columns\n")
  
  adjust_cols_result = adjustCols(st_opt_postz, gt_opt_postz, node_labels, new_node_labels)
  st_opt_postz = adjust_cols_result[[1]]
  gt_opt_postz = adjust_cols_result[[2]]
  #targ_ind = adjust_cols_result[[3]]

  # Organize column names
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Reading BEAST results: ", beast_file, "\n")
  
  #beast_file_path = here("results", "simulations", beast_dir, beast_file)
  
  read_beast_results = readBeast(beast_file, cur_case_targets, node_labels, new_node_labels, theta=TRUE)
  beast_targets = read_beast_results[[1]]
  beast_non_targets = read_beast_results[[2]]

  # Load BEAST results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Converting target lineages to long data\n")
  
  st_targ = select(st_opt_postz, cur_case_targets[["targets"]]$new.labels)
  st_targ_long = reshape2::melt(st_targ)
  st_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Metho
  
  gt_targ = select(gt_opt_postz, cur_case_targets[["targets"]]$new.labels)
  gt_targ_long = reshape2::melt(gt_targ)
  gt_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Metho
  
  beast_targ_long = reshape2::melt(beast_targets)
  names(beast_targ_long)[1] = "variable"
  beast_targ_long$Method = "*BEAST2"
  
  target_results = rbind(st_targ_long, gt_targ_long, beast_targ_long)
  names(target_results)[1:2] = c("species", "acc")
  # Combine st and gt results
  
  # Select target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Converting non-target linages to long data\n")
  
  st_non_targ_long = select(st_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  st_non_targ_long = reshape2::melt(st_non_targ_long)
  st_non_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Metho
  
  gt_non_targ_long = select(gt_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  gt_non_targ_long = reshape2::melt(gt_non_targ_long)
  gt_non_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Metho
  
  beast_non_targ_long = reshape2::melt(beast_non_targets)
  names(beast_non_targ_long)[1] = "variable"
  beast_non_targ_long$Method = "*BEAST2"
  
  non_target_results = rbind(st_non_targ_long, gt_non_targ_long, beast_non_targ_long)
  names(non_target_results)[1:2] = c("species", "acc")
  # Combine st and gt results
  
  # Select non-target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Target boxplot\n")
  
  target_results$Method = factor(target_results$Method, levels=method_order)
  
  target_p = ggplot(target_results, aes(x=species, y=acc, color=Method, fill=Method)) +
    geom_quasirandom(size=1, width=0.1, alpha=0.1, dodge.width=0.8) +
    geom_boxplot(outlier.shape=NA, alpha=0.1) +
    xlab("Nodes simulated with acceleration") +
    ylab("P(acceleration)") +
    scale_color_manual(values=c("#CC79A7","#56B4E9","#F0E442")) +
    bartheme() +
    theme(legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.text.x=element_text(angle=40, hjust=1, size=8),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=12),
          plot.margin=margin(1,0.1,0,0.1, unit="cm"))
  #panel.grid.major.x = element_line(color="#d3d3d3", size=0.25))
  
  if(theta_scale == 3){
    target_p = target_p + theme(legend.position="bottom")
    fig_legend = get_legend(target_p)
    target_p = target_p + theme(legend.position="none")
  }
  
  #targ_p_list[[c]] = p
  
  # Targets
  #####
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], "- Theta", theta_scale, "X : Non-target boxplot\n")
  
  non_target_results$Method = factor(non_target_results$Method, levels=method_order)
  
  non_target_p = ggplot(non_target_results, aes(x=species, y=acc, color=Method, fill=Method)) +
    geom_quasirandom(size=1, width=0.1, alpha=0.1, dodge.width=0.8) +
    geom_boxplot(outlier.shape=NA, fill="transparent") +
    xlab("Nodes simulated without acceleration") +
    ylab("") +
    scale_color_manual(values=c("#CC79A7","#56B4E9","#F0E442")) +
    bartheme() +
    theme(legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.text.x=element_text(angle=40, hjust=1, size=8),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          plot.margin=margin(1,0.1,0,0.1, unit="cm"))
  
  #non_targ_p_list[[c]] = p
  
  
  row_label = "C"
  panel_labels = c("E","F")
  title = titles[3]
  h_adj = -1.5
  if(theta_scale == 3){
    row_label = "A"
    panel_labels = c("A", "B")
    title = titles[1]
    h_adj = -1.6
  }else if(theta_scale == 6){
    row_label = "B"
    panel_labels = c("C","D")
    title = titles[2]
    h_adj = -1.6
  }
  
  p_title = ggdraw() + 
    draw_label(
      title,
      x = 0,
      vjust = 1,
      hjust = h_adj
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  p_combo = plot_grid(target_p, non_target_p, ncol=2, labels=panel_labels, label_y=1)
  p_panel = plot_grid(p_combo, p_title, nrow=2, rel_heights=c(1,0.2))
  print(p_panel)
  p_list[[row_label]] = p_panel
  # Combine the target and non-target plots
  
  
  
  
  # Non-targets
  #####
  # Generate boxplots
  ##########
}

cat(as.character(Sys.time()), " | Fig10: Combining plots...\n")

fig_main = plot_grid(plotlist=p_list[c("A", "B", "C")], nrow=3)

#fig_left = plot_grid(plotlist=targ_p_list, nrow=3, labels=c("A", "B", "C"))
#fig_right = plot_grid(plotlist=non_targ_p_list, nrow=3)
#fig_main = plot_grid(fig_left, fig_right, ncol=2)
fig = plot_grid(fig_main, fig_legend, nrow=2, rel_heights=c(1, 0.1))

print(fig)

# Combine the plots
######################

if(save_fig){
  figfile = "../figs/fig10.pdf"
  cat(as.character(Sys.time()), " | Fig10: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=8, height=8, units="in")
}

# Save the figure
######################








