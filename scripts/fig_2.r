############################################################
# For PhyloAcc
# Figure 2
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(viridis)
library(ggbeeswarm)
library(ggsignif)
library(cowplot)
library(here)
source("lib/design.r")
source("lib/get_tree_info.r")

############################################################

save_int_figs = F
save_final_fig = F
outdir = here("figs")

# Options
############################################################

cat(as.character(Sys.time()), " | Fig2: Reading trees...\n")

full_treefile = here("data", "trees", "ratite-subset-generic.tree")
full_sim_tree = read.tree(full_treefile)

treefile = here("data", "trees", "ratite-subset-generic-reduced.tree")
sim_tree = read.tree(treefile)
tree_to_df_list = treeToDF(sim_tree)
sim_tree_info = tree_to_df_list[["info"]]
# Read the tree and parse with treetoDF

sim_tree_info = sim_tree_info[order(sim_tree_info$node), ]
# Re-sort the data frame by R node order after the merge so the trees still work

#tree_info_file = "../data/ratite-sim.csv"
#write.csv(file=tree_info_file, sim_tree_info, row.names=F)
# Write the tree info to a csv

sim_info_file = here("data", "ratite-sim-cases.csv")
sim_info = read.csv(sim_info_file, header=T)
# Read the tree info with simulation cases

cols = c("N"="#333333", "Y"=corecol(pal="wilke", numcol=1))

# Input files
############################################################

cat(as.character(Sys.time()), " | Fig2: Generating panel A (full tree)...\n")

full_p = ggtree(full_sim_tree, size=1, ladderize=F) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        legend.title=element_text(size=8))
print(full_p)

# Full tree
############################################################

cat(as.character(Sys.time()), " | Fig2: Generating panel B (single acceleration)...\n")

sim_info$Case.3 = as.factor(sim_info$Case.3)

mono_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.3), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  xlim(0,10) +
  scale_color_manual(name='Lineages simulated with accelerated sequence evolution:', values=cols, breaks=c("Y")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  ggtitle(paste("Single acceleration\n")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        legend.title=element_text(size=8),
        plot.title=element_text(hjust=0.5, size=12),
        plot.margin=margin(0.8,0.1,0,0.1, unit="cm")) +
  #guides(colour=guide_legend(nrow=2)) +
  coord_flip()
print(mono_p)

# Single acceleration (Case 3)
############################################################

cat(as.character(Sys.time()), " | Fig2: Generating panel C (two accelerations)...\n")

conv_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.2), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  #xlim(0,10) +
  scale_color_manual(name='Lineages simulated with accelerated sequence evolution:', values=cols, breaks=c("Y","N")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  ggtitle(paste("Two independent\naccelerations")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        plot.title=element_text(hjust=0.5, size=12),
        plot.margin=margin(0.8,0.1,0,0.1, unit="cm")) +
  coord_flip()
print(conv_p)

# Two accelerations (Case 2)
############################################################

cat(as.character(Sys.time()), " | Fig2: Generating panel D (three accelerations)...\n")

para_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.4), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  xlim(0,10) +
  scale_color_manual(name='Lineages simulated with accelerated sequence evolution:', values=cols, breaks=c("Y","N")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  ggtitle(paste("Three independent\naccelerations")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        plot.title=element_text(hjust=0.5, size=12),
        plot.margin=margin(0.8,0.1,0,0.1, unit="cm")) +
  coord_flip()
print(para_p)

# Three accelerations (Case 4)
############################################################

cat(as.character(Sys.time()), " | Fig2: Combining panels...\n")

leg = get_legend(mono_p)
# Get the legend

#panels = list(mono_p + theme(legend.position="none"),
#              conv_p + theme(legend.position="none"),
#              para_p + theme(legend.position="none"))
# Combine the 3 simulation trees and remove legens

fig_full = plot_grid(NULL, full_p, NULL, ncol=3, rel_widths=c(0.25, 0.5, 0.25), labels=c("", "A", ""), label_size=14)
# First panel with full tree and NULL margins

fig_panel = plot_grid(mono_p + theme(legend.position="none"),
                      conv_p + theme(legend.position="none"),
                      para_p + theme(legend.position="none"),
                      nrow=1, labels=c("B","C","D"), label_size=14)
# Combine the 3 simulation trees and remove legends

fig_main = plot_grid(fig_full, fig_panel, nrow=2, rel_heights=c(1, 0.7))
# Combine top and bottom

fig = plot_grid(fig_main, leg, nrow=2, rel_heights=c(1,0.1))
print(fig)
# Add legend

# Combine panels
############################################################

if(save_final_fig){
  figfile = here(outdir, "fig2.pdf")
  cat(as.character(Sys.time()), " | Fig3: Saving figure:", figfile, "\n")
  ggsave(figfile, fig, width=6.5, height=6, unit="in")
}

# Save the final figure
############################################################




