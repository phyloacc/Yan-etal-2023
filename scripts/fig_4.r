############################################################
# For PhyloAcc
# Boxplots
# Han Yan and Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

############################################################

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(ggbeeswarm)
library(stringr)
library(here)
source(here("scripts", "lib", "design.r"))
source(here("scripts", "lib", "read_results.r"))

## KEY from old variable names
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

node_labels = c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "cryCin", "tinGut", "eudEle", "notPer", "anoDid", "strCam",
                "aptHaa.aptOwe", "aptHaa.aptRow", "casCas.droNov", "aptHaa.casCas", "rheAme.rhePen", "aptHaa.rheAme", "cryCin.tinGut",
                "eudEle.notPer", "cryCin.eudEle", "cryCin.anoDid", "aptHaa.cryCin", "aptHaa.strCam")

save_fig = F

# Options and inputs
############################################################

cat(as.character(Sys.time()), " | Fig4: Adjusting node labels...\n")

new_node_labels = adjustLabels(label_map_file, node_labels)

## Re-label the nodes based on the simulation letters
############################################################

cat(as.character(Sys.time()), " | Fig4: Assigning target lineages...\n")

case_targets = assignTargets(node_labels, new_node_labels)

## Assign target lineages
############################################################

cat(as.character(Sys.time()), " | Fig4: Reading likelihood files...\n")

st_elem_lik_file = here("results", "simulations", "Acc300elem_lik.txt")
gt_elem_lik_file = here("results", "simulations", "GT300elem_lik.txt")
# Results files from 300 loci simulate

st_elem_lik = readElemLik(st_elem_lik_file, subset_cols=TRUE)
# Reading the species tree method results and getting adding the model with
# the max probability

gt_elem_lik = readElemLik(gt_elem_lik_file)
# Reading the gene tree method results and getting adding the model with
# the max probability

############################################################

cat(as.character(Sys.time()), " | Fig4: Reading postz files and generating plots...\n")

case_labels = c("A", "B", "C")
## CASE LABELS DO NOT CORRESPOND TO PANEL LABELS
# Case 2 / scenario 8 / label B / two independent accelerations
# Case 3 / scenario 5 / label A / single acceleration
# Case 4 / scenario 6 / label C / three independent accelerations

scenarios = c("8", "5", "6")
# This is how the files for the different cases are labeled...

method_order = c("PhyloAcc-GT", "PhyloAcc", "*BEAST2")
# Ordering for the legend

titles = c("Single acceleration", "Two independent accelerations", "Three independent accelerations")
# Titles for the rows of the figure

targ_p_list = list()
non_targ_p_list = list()
p_list = list()
# Initialize lists of plots

for(c in 1:3){
  
  cur_case_targets = case_targets[[c+1]]
  # Look up the current set of target and non-target lineages

  st_dir = paste("PhyloAcc_result_Case", case_labels[c], sep="")
  gt_dir = paste("PhyloAcc-GT_result_Case", case_labels[c], sep="")
  beast_dir = paste("starBeast2_Result_Case", case_labels[c], sep="")
  # Input directories for the different run types
  
  m1_postz_file = paste("simu_200_100_2-", scenarios[c], "_rate_postZ_M1.txt", sep="")
  m2_postz_file = paste("simu_200_100_2-", scenarios[c], "_rate_postZ_M2.txt", sep="")
  beast_file = paste("result1_2-", scenarios[c], ".RData", sep="")
  # File names based on the current scenario
  
  ## Some setup
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Reading ST probs\n")
  
  st_m1_postz_file = here("results", "simulations", st_dir, m1_postz_file)
  st_m2_postz_file = here("results", "simulations", st_dir, m2_postz_file)
  # Rate file paths for st
  
  st_elem_lik_scenario = subset(st_elem_lik, No. == scenarios[c])
  st_opt_postz = readPostZ(st_m1_postz_file, st_m2_postz_file, st_elem_lik_scenario)

  # Read st results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Reading GT probs\n")
  
  gt_m1_postz_file = here("results", "simulations", gt_dir, m1_postz_file)
  gt_m2_postz_file = here("results", "simulations", gt_dir, m2_postz_file)
  # Rate file paths for gt
  
  gt_elem_lik_scenario = subset(gt_elem_lik, No. == scenarios[c])
  gt_opt_postz = readPostZ(gt_m1_postz_file, gt_m2_postz_file, gt_elem_lik_scenario)

  # Read gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Adjusting columns\n")

  adjust_cols_result = adjustCols(st_opt_postz, gt_opt_postz, node_labels, new_node_labels)
  st_opt_postz = adjust_cols_result[[1]]
  gt_opt_postz = adjust_cols_result[[2]]
  #targ_ind = adjust_cols_result[[3]]

  # Organize column names
  ##########
  
  beast_file_path = here("results", "simulations", beast_dir, beast_file)
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Reading BEAST results:", beast_file_path, "\n")

  read_beast_results = readBeast(beast_file_path, cur_case_targets, node_labels, new_node_labels)
  beast_targets = read_beast_results[[1]]
  beast_non_targets = read_beast_results[[2]]

  # Load BEAST results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Converting target lineages to long data\n")

  st_targ = select(st_opt_postz, cur_case_targets[["targets"]]$new.labels)
  st_targ_long = reshape2::melt(st_targ)
  st_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Method

  gt_targ = select(gt_opt_postz, cur_case_targets[["targets"]]$new.labels)
  gt_targ_long = reshape2::melt(gt_targ)
  gt_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Method

  beast_targ_long = reshape2::melt(beast_targets)
  names(beast_targ_long)[1] = "variable"
  beast_targ_long$Method = "*BEAST2"

  target_results = rbind(st_targ_long, gt_targ_long, beast_targ_long)
  names(target_results)[1:2] = c("species", "acc")
  # Combine st and gt results

  # Select target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Converting non-target linages to long data\n")
  
  st_non_targ_long = select(st_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  st_non_targ_long = reshape2::melt(st_non_targ_long)
  st_non_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Method
  
  gt_non_targ_long = select(gt_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  gt_non_targ_long = reshape2::melt(gt_non_targ_long)
  gt_non_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Method
  
  beast_non_targ_long = reshape2::melt(beast_non_targets)
  names(beast_non_targ_long)[1] = "variable"
  beast_non_targ_long$Method = "*BEAST2"
  
  non_target_results = rbind(st_non_targ_long, gt_non_targ_long, beast_non_targ_long)
  names(non_target_results)[1:2] = c("species", "acc")
  # Combine st and gt results
  
  # Select non-target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Target boxplot\n")
  
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
  
  if(c == 1){
    target_p = target_p + theme(legend.position="bottom")
    fig_legend = get_legend(target_p)
    target_p = target_p + theme(legend.position="none")
  }
  # Get the legend from the first plot
  
  #targ_p_list[[c]] = p
  
  # Targets
  #####
  
  cat(as.character(Sys.time()), " | ------->", case_labels[c], ": Non-target boxplot\n")
  
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
  
  
  target_results$grouping = "Accelerated"
  non_target_results$grouping = "Non-accelerated"
  all_results = rbind(target_results, non_target_results)
  all_results$Method = factor(all_results$Method, levels=method_order)
  
  num_accel = nrow(cur_case_targets$targets)
  num_non_accel = nrow(cur_case_targets$non.targets)
  
  red = rep("red", num_accel)
  black = rep("black", num_non_accel)
  
  all_cols = c(red, black)
  
  line_y = -0.7
  label_y = line_y - 0.1
  
  all_p = ggplot(all_results, aes(x=species, y=acc, color=Method, fill=Method)) +
    geom_quasirandom(size=1, width=0.1, alpha=0.1, dodge.width=0.8) +
    geom_boxplot(outlier.shape=NA, fill="transparent") +
    #geom_text(position = position_dodge(width = 1), aes(x=grouping, y=0), label=grouping) +
    #facet_wrap(~grouping, strip.position = "bottom", scales = "free_x") +
    xlab("") +
    ylab("") +
    scale_color_manual(values=c("#CC79A7","#56B4E9","#F0E442")) +
    bartheme() +
    theme(legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.text.x=element_text(angle=40, hjust=1, size=8, color=all_cols),
          axis.title.y=element_text(size=10),
          plot.margin=margin(1,0.1,0,0.1, unit="cm"),
          panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(),
          strip.placement = "outside") +
    coord_cartesian(xlim=c(1,25), ylim=c(0,1), clip="off") +
    annotate("segment", x = 1, xend = num_accel, y = line_y, yend = line_y) +
    annotate("text", x = num_accel/2, y = label_y, label="Accelerated") +
    annotate("segment", x = num_accel+1, xend = num_accel+num_non_accel, y = line_y, yend = line_y) +
    annotate("text", x = ((num_accel+1)+(num_accel+num_non_accel))/2, y = label_y, label="Non-Accelerated")
  
  print(all_p)
  
  row_label = "C"
  panel_labels = c("E","F")
  title = titles[3]
  h_adj = -1
  if(c == 1){
    row_label = "B"
    panel_labels = c("C","D")
    title = titles[2]
    h_adj = -1.05
  }else if(c == 2){
    row_label = "A"
    panel_labels = c("A", "B")
    title = titles[1]
    h_adj = -2
  }
  # Title and labels depending on panel
  
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
  # Title for current row
  
  p_combo = plot_grid(target_p, non_target_p, ncol=2, labels=panel_labels, label_y=1)
  p_panel = plot_grid(p_combo, p_title, nrow=2, rel_heights=c(1,0.2))
  #p_panel = plot_grid(p_title, all_p, nrow=2, rel_heights=c(0.2,1))
  print(p_panel)
  p_list[[row_label]] = p_panel
  # Combine the target and non-target plots
  
  # Non-targets
  #####
  # Generate boxplots
  ##########
}

cat(as.character(Sys.time()), " | Fig4: Combining plots...\n")

fig_main = plot_grid(plotlist=p_list[c("A", "B", "C")], nrow=3)
fig = plot_grid(fig_main, fig_legend, nrow=2, rel_heights=c(1, 0.1))
print(fig)

# Combine the plots
######################

if(save_fig){
  figfile = "../figs/fig4.pdf"
  cat(as.character(Sys.time()), " | Fig4: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=8, height=8, units="in")
}

# Save the figure
######################








