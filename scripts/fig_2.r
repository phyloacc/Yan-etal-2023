################################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

################################################################

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

################################################################

save_int_figs = F
save_final_fig = F
outdir = here("figs")

# Options
################################################################

treefile = here("data", "ratite-sim.tre")
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
################################################################
# 
# plotCladogram<-function(tree){
#   foo<-function(tree,x){
#     n<-1:tree$Nnode+length(tree$tip.label)
#     setNames(sapply(n,function(n,x,t) mean(x[Descendants(t,n,
#                                                          "tips")[[1]]]),x=x,t=tree),n)
#   }
#   tree<-reorder(tree,"cladewise")
#   x<-setNames(1:length(tree$tip.label),tree$tip.label)
#   phenogram(compute.brlen(tree),c(x,foo(tree,x)),ylab="")
# }
# 
# plotCladogram(sim_tree)

sim_info$Case.3 = as.factor(sim_info$Case.3)


mono_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.3), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  xlim(0,10) +
  scale_color_manual(name='Accelerated lineages:', values=cols, breaks=c("Y")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  #ggtitle(paste("Case", case)) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        legend.title=element_text(size=8)) +
  #guides(colour=guide_legend(nrow=2)) +
  coord_flip()
print(mono_p)

# Monophyletic acceleration (Case 3)
################################################################

para_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.4), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  xlim(0,10) +
  scale_color_manual(name='Accelerated?', values=cols, breaks=c("Y","N")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  #ggtitle(paste("Case", case)) +
  scale_x_reverse() +
  theme(legend.position="bottom") +
  coord_flip()
print(para_p)

# Paraphyletic acceleration (Case 4)
################################################################

conv_p = ggtree(sim_tree, size=2, ladderize=F, aes(color=sim_info$Case.2), position=position_nudge(x=0.4, y = 0.3)) +
  #tree_p = ggtree(sim_tree, size=2, ladderize=F, aes_string(data=sim_info, color=case_str)) +
  #xlim(0,10) +
  scale_color_manual(name='Accelerated?', values=cols, breaks=c("Y","N")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2) +
  #ggtitle(paste("Case", case)) +
  scale_x_reverse() +
  theme(legend.position="bottom") +
  coord_flip()
print(conv_p)

# Convergent acceleration (Case 2)
################################################################

leg = get_legend(mono_p)

panels = list(mono_p + theme(legend.position="none"),
              para_p + theme(legend.position="none"),
              conv_p + theme(legend.position="none"))

fig_panel = plot_grid(mono_p + theme(legend.position="none"),
                      para_p + theme(legend.position="none"),
                      conv_p + theme(legend.position="none"),
                      nrow=1, labels=c("A","B","C"), label_size=14)
fig = plot_grid(fig_panel, leg, nrow=2, rel_heights=c(1,0.1))

if(save_final_fig){
  outfilename = paste(outdir, "fig2.png", sep="")
  ggsave(outfilename, fig, width=6.5, height=2, unit="in")
}

# Save the final figure
################################################################




