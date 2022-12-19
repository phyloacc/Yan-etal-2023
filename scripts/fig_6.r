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
library(ggbeeswarm)
library(stringr)
library(here)
source(here("scripts", "lib", "design.r"))
source(here("scripts", "lib", "read_results.r"))

############################################################

ratite_tree_file = here("data", "ratite.tree")

label_map_file = here("data", "SimulatedData", "Name_Label_Mapping.txt")

node_labels = c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","cryCin","tinGut","eudEle", "notPer","anoDid","strCam",
                "aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme","cryCin.tinGut",
                "eudEle.notPer","cryCin.eudEle","cryCin.anoDid", "aptHaa.cryCin","aptHaa.strCam")

save_fig = F

# Setup
############################################################

cat(as.character(Sys.time()), " | Fig4: Adjsuting node labels...\n")

new_node_labels = adjustLabels(label_map_file, node_labels)

## Re-label the nodes based on the simulation letters
############################################################

cat(as.character(Sys.time()), " | Fig4: Assigning target lineages...\n")

case_targets = assignTargets(node_labels, new_node_labels)

## Assign target lineages
############################################################

cat(as.character(Sys.time()), " | Fig6: Generating boxplots for 3 scenarios...\n")

scenarios = c(1, 2, 3)
########scenario1: input 2 true 3: M2 correct. 2 sets intersects

target_sets = c(3, 4, 3)
file_labels = c("NC_true5spec8", "trueC4inputC2", "NC_true5spec4")
panel_labels = list(c("A", "B"), c("C", "D"), c("E", "F"))
row_labels = c("A", "B", "C")
# Initialize scenario variables and names...


input_targets = list(c("A1", "A2", "A3", "(A1,A2)", "(A1,A3)", "C1", "C2", "(C1,C2)"),
                     c("A1", "A2", "A3", "C1", "C2", "(A1,A2)", "(A1,A3)", "(C1,C2)"),
                     c("C1", "C2", "D1", "F1", "(C1,C2)"))

target_x_col = list(c("red", "red", "red", "black", "black", "red", "red", "black", "black"),
                    c("red", "red", "red", "black", "black", "red", "red", "black", "black", "red", "red", "black", "black", "red", "black"),
                    c("black", "black", "black", "black", "black", "black", "black", "black", "black"))

non_target_x_col = list(c("red", "red", "black", "black", "black", "black", "black", "black", "red", "black", "black", "black", "black", "black", "black", "black"),
                        c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
                        c("red", "red", "black", "black", "black", "black", "red", "red", "red", "black", "black", "black", "black", "black", "black", "black"))
# Setup lists for the color of the x-axis labels to indicate the input target lineages... I couldn't find a better way to do
# this with the data I have...

titles = c("Scenario 1", "Scenario 2", "Scenario 3")
# Titles for plots

p_list = list()
# A list in which to save the generated boxplots

for(s in scenarios){
  
  cur_case_targets = case_targets[[target_sets[s]]]
  # Look up the current set of target and non-target lineages
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Reading likelihood files\n")
  
  st_dir = here("results", "simulations", "MisspecifiedInput", paste0("Scenario", s), "Acc")
  gt_dir = here("results", "simulations", "MisspecifiedInput", paste0("Scenario", s))
  # Directories for st and gt files
  
  st_elem_lik_file = here(st_dir, paste0("simu_", file_labels[s], "_elem_lik.txt"))
  gt_elem_lik_file = here(gt_dir, paste0("simu_", file_labels[s], "_elem_lik.txt"))
  # Results files from the current scenario
  
  st_elem_lik = readElemLik(st_elem_lik_file)
  # Reading the species tree method results and getting adding the model with
  # the max probability
  
  gt_elem_lik = readElemLik(gt_elem_lik_file)
  # Reading the gene tree method results and getting adding the model with
  # the max probability
  
  
  ############################################################
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Reading probability files...\n")
  
  m1_postz_file = paste0("simu_", file_labels[s], "_rate_postZ_M1.txt")
  m2_postz_file = paste0("simu_", file_labels[s], "_rate_postZ_M2.txt")
  #beast_file = paste("result1_2-", scenarios[c], ".RData", sep="")
  # File names based on the current scenario
  
  ##########
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Reading ST probs\n")
  
  st_m1_postz_file = here(st_dir, m1_postz_file)
  st_m2_postz_file = here(st_dir, m2_postz_file)
  # Prob file paths for st
  
  st_opt_postz = readPostZ(st_m1_postz_file, st_m2_postz_file, st_elem_lik)
  
  # Read st results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Reading GT probs\n")
  
  gt_m1_postz_file = here(gt_dir, m1_postz_file)
  gt_m2_postz_file = here(gt_dir, m2_postz_file)
  # Prob file paths for gt
  
  gt_opt_postz = readPostZ(gt_m1_postz_file, gt_m2_postz_file, gt_elem_lik)
  
  # Read gt results
  ############################################################
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Adjusting columns\n")
  
  adjust_cols_result = adjustCols(st_opt_postz, gt_opt_postz, node_labels, new_node_labels)
  st_opt_postz = adjust_cols_result[[1]]
  gt_opt_postz = adjust_cols_result[[2]]
  
  # Organize column names
  ############################################################

  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Converting target lineages to long data\n")
  
  st_targ = select(st_opt_postz, cur_case_targets[["targets"]]$new.labels)
  st_targ_long = reshape2::melt(st_targ)
  st_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Metho
  
  gt_targ = select(gt_opt_postz, cur_case_targets[["targets"]]$new.labels)
  gt_targ_long = reshape2::melt(gt_targ)
  gt_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Metho
  
  target_results = rbind(st_targ_long, gt_targ_long)
  names(target_results)[1:2] = c("species", "acc")
  # Combine st and gt results
  
  # Select target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Converting non-target lineages to long data\n")
  
  st_non_targ_long = select(st_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  st_non_targ_long = reshape2::melt(st_non_targ_long)
  st_non_targ_long$Method = "PhyloAcc"
  # Convert from wide to long for st target branches and add a label for Metho
  
  gt_non_targ_long = select(gt_opt_postz, cur_case_targets[["non.targets"]]$new.labels)
  gt_non_targ_long = reshape2::melt(gt_non_targ_long)
  gt_non_targ_long$Method = "PhyloAcc-GT"
  # Convert from wide to long for gt target branches and add a label for Metho
  
  non_target_results = rbind(st_non_targ_long, gt_non_targ_long)
  names(non_target_results)[1:2] = c("species", "acc")
  # Combine st and gt results
  
  # Select non-target branches and combine st and gt results
  ##########
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Target boxplot\n")
  
  method_order = c("PhyloAcc-GT", "PhyloAcc")
  # Ordering for the legend
  
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
          axis.text.x=element_text(angle=40, hjust=1, size=8, color=target_x_col[[s]]),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=12),
          plot.margin=margin(1,0.1,0,0.1, unit="cm"))
  
  if(s == 1){
    target_p = target_p + theme(legend.position="bottom")
    fig_legend = get_legend(target_p)
    target_p = target_p + theme(legend.position="none")
  }
  # Get the legend from the first plot
  
  #print(target_p)
  # Targets
  #####
  
  cat(as.character(Sys.time()), " | ------->", scenarios[s], ": Non-target boxplot\n")
  
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
          axis.text.x=element_text(angle=40, hjust=1, size=8, color=non_target_x_col[[s]]),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=12),
          plot.margin=margin(1,0.1,0,0.1, unit="cm"))
  
  #print(non_target_p)
  # Non-targets
  #####
  
  p_title = ggdraw() + 
    draw_label(
      titles[s],
      x = 0,
      vjust = 1,
      hjust = -4
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  # Current row title
  
  p_combo = plot_grid(target_p, non_target_p, ncol=2, labels=panel_labels[[s]])
  p_panel = plot_grid(p_combo, p_title, nrow=2, rel_heights=c(1,0.2))
  p_list[[row_labels[s]]] = p_panel
  # Combine the target and non-target plots
  
  # Generate boxplots
  ##########
}

cat(as.character(Sys.time()), " | Fig6: Combining plots...\n")

fig_main = plot_grid(plotlist=p_list[c("A", "B", "C")], nrow=3)
fig = plot_grid(fig_main, fig_legend, nrow=2, rel_heights=c(1, 0.1))
print(fig)
# Combine the panels and the legend

# Combine the plots
######################

if(save_fig){
  figfile = "../figs/fig6.png"
  cat(as.character(Sys.time()), " | Fig6: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=8, height=8, units="in")
}

# Save the figure
######################


