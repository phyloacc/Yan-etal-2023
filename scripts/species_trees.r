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

write_ratite_info = F
save_ratite_fig = F
write_mammal_info = F
save_mammal_fig = T

outdir = here("figs")

# Options
################################################################

ratite_treefile = here("data", "trees", "ratite-full.tree")
ratite_info_file = here("data", "trees", "ratites.csv")

mammal_treefile = here("data", "trees", "mammal-full.tree")
mammal_info_file = here("data", "trees", "mammals.csv")

# Input files
################################################################

cat(as.character(Sys.time()), "  | Reading ratite tree\n", sep="")

ratite_tree = read.tree(ratite_treefile)
ratite_info = read.csv(ratite_info_file, header=T)
tree_to_df_list = treeToDF(ratite_tree)
tree_info = tree_to_df_list[["info"]]
ratite_info = merge(ratite_info, tree_info, by="label", all=T)
# Read the tree and parse with treetoDF

ratite_info = ratite_info[order(ratite_info$node), ]
# Re-sort the data frame by R node order after the merge so the trees still work

if(write_ratite_info){
  cat(as.character(Sys.time()), "  | Writing ratite tree info\n", sep="")
  ratite_tree_info_out = here("data", "trees", "ratite-tree.csv")
  write.csv(ratite_info, ratite_tree_info_out, row.names=F)
}
# Write the full table if specified

cat(as.character(Sys.time()), "  | Re-labeling ratite tree tips\n", sep="")
ratite_tree$tip.label = ratite_info[ratite_info$node.type=="tip",]$s.name

cat(as.character(Sys.time()), "  | Plotting ratite tree\n", sep="")
ratite_p = ggtree(ratite_tree, size=1, ladderize=F) +
  xlim(0,1) +
  #scale_color_manual(name='Accelerated lineages:', values=cols, breaks=c("Y")) +
  geom_tiplab(color="#333333", size=4, hjust=-0.2, fontface="italic") +
  #ggtitle(paste("Case", case)) +
  #scale_x_reverse() +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        legend.title=element_text(size=8))
#guides(colour=guide_legend(nrow=2)) +
#coord_flip()

print(ratite_p)

if(save_ratite_fig){
  figfile = "../figs/full-ratite-tree.pdf"
  cat(as.character(Sys.time()), " | Saving ratite tree figure:", figfile, "\n")
  ggsave(filename=figfile, ratite_p, width=6, height=7, units="in")
}

# Ratites
################################################################

cat(as.character(Sys.time()), "  | Reading mammal tree\n", sep="")

mammal_tree = read.tree(mammal_treefile)
mammal_info = read.csv(mammal_info_file, header=T)
tree_to_df_list = treeToDF(mammal_tree)
tree_info = tree_to_df_list[["info"]]
mammal_info = merge(mammal_info, tree_info, by="label", all=T)
# Read the tree and parse with treetoDF

mammal_info = mammal_info[order(mammal_info$node), ]
# Re-sort the data frame by R node order after the merge so the trees still work

if(write_mammal_info){
  cat(as.character(Sys.time()), "  | Writing mammal tree info\n", sep="")
  mammal_tree_info_out = here("data", "trees", "mammal-tree.csv")
  write.csv(mammal_info, mammal_tree_info_out, row.names=F)
}
# Write the full table if specified

cat(as.character(Sys.time()), "  | Re-labeling mammal tree tips\n", sep="")
mammal_tree$tip.label = mammal_info[mammal_info$node.type=="tip",]$s.name

cat(as.character(Sys.time()), "  | Plotting mammal tree\n", sep="")
mammal_p = ggtree(mammal_tree, size=0.75, ladderize=F) +
  xlim(0,1) +
  #scale_color_manual(name='Accelerated lineages:', values=cols, breaks=c("Y")) +
  geom_tiplab(color="#333333", size=3, hjust=-0.2, fontface="italic") +
  #ggtitle(paste("Case", case)) +
  #scale_x_reverse() +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        legend.title=element_text(size=8))
#guides(colour=guide_legend(nrow=2)) +
#coord_flip()

print(mammal_p)

if(save_mammal_fig){
  figfile = "../figs/full-mammal-tree.pdf"
  cat(as.character(Sys.time()), " | Saving mammal tree figure:", figfile, "\n")
  ggsave(filename=figfile, mammal_p, width=6, height=7, units="in")
}

# Ratites
################################################################

