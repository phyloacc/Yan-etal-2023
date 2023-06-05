################################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggbeeswarm)
library(ggsignif)
library(cowplot)
source("lib/design.r")

################################################################

save_int_figs = F
save_final_fig = T
outdir = "../figs/"

# Options
################################################################

cat(as.character(Sys.time()), " | Fig11: Reading data...\n")

ratite_aln_stats_file = "../data/ratite-scf-tests-2/ratite-output-all-st-0.1-0.1/phyloacc-aln-stats.csv"
ratite_aln_stats = read.csv(ratite_aln_stats_file, header=T)
# Ratite locus data

sim_data_file = "../data/sim-benchmarks.csv"
sim_data = read.csv(sim_data_file, header=T, comment="#", quote='"')
sim_data = head(sim_data, -1)
# Read benchmark data

sim_data$maxmem = substr(sim_data$maxmem,1,nchar(sim_data$maxmem)-1)
sim_data$maxmem = as.numeric(sim_data$maxmem)
# Get max mem used

sim_data$sim = factor(sim_data$sim, levels=order(sim_data$sim))
# Factor the sims

sim_data_avg = sim_data %>% group_by(sim, type) %>% summarize(spec=mean(spec), length=mean(length), theta=mean(theta),
                                                              avg.runtime=(sum(runtime) / sum(num.loci)) / 60,
                                                              avg.cputime=(sum(cputime) / sum(num.loci)) / 60,
                                                              avg.maxmem=(sum(maxmem) / n()) / 1024)
# Calculate averages

sim_data_avg$spec = as.character(sim_data_avg$spec)
sim_data_avg$spec = factor(sim_data_avg$spec, levels=c("9", "13", "17"))
# Simulation data

cols = corecol(numcol=3)
cols = setNames(cols, levels(sim_data_avg$spec))
# Colors for some plots

options(scipen = 999)
# wat

# Input files
################################################################

cat(as.character(Sys.time()), " | Fig11: Plotting panel A...\n")

cputime_p = ggplot(sim_data_avg, aes(x=length, y=avg.cputime, color=as.character(spec), linetype=type)) +
  geom_line(size=1, alpha=0.5) +
  geom_point(size=3, alpha=0.5) +
  #scale_color_viridis(name='# of species in tree', option = "C", discrete=T) +
  scale_color_manual(name='# of species in tree: ', values=cols) +
  scale_linetype_manual(name='PhyloAcc model', labels=c("Gene\ntree", "Species\ntree"), values=c("solid", 22)) +
  xlab("Element length (bp)") +
  ylab("Avg. CPU time\nper element (minutes)") +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.75,0.15,0.5), "cm")) +
  guides(color="none", linetype=guide_legend(override.aes=list(alpha=1), title="PhyloAcc model: "))
print(cputime_p)
# The CPU time plot with the color legend

######

cputime_color_p = ggplot(sim_data_avg, aes(x=length, y=avg.cputime, color=as.character(spec), linetype=type)) +
  geom_line(size=1, alpha=0.5) +
  geom_point(size=3, alpha=0.5) +
  #scale_color_viridis(name='# of species in tree', option = "C", discrete=T) +
  scale_color_manual(name='# of species in tree: ', breaks=c("9", "13", "17"), values=cols) +
  scale_linetype_manual(name='PhyloAcc method', labels=c("Gene\ntree", "Species\ntree"), values=c("solid", 22)) +
  xlab("Element length (bp)") +
  ylab("Avg. CPU time\nper element (minutes)") +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")) +
  guides(linetype="none")
print(cputime_color_p)
# The same plot with the method legend...

spec_legend = get_legend(cputime_color_p)
type_legend = get_legend(cputime_p)
a_legend = plot_grid(spec_legend, type_legend, nrow=2)
# Get and combine both legends
###

panel_a = plot_grid(cputime_p + theme(legend.position="none"), a_legend, nrow=2, rel_heights=c(1,0.3))
print(panel_a)

if(save_int_figs){
  outfilename = paste(outdir, "sim-cputimes", sep="")
  outfilename = paste(outfilename, ".png", sep="")
  cat(as.character(Sys.time()), " | Fig11: Saving panel A: ", outfilename, "\n")
  ggsave(outfilename, cputime_p, width=5, height=4, unit="in")
}
# Save the intermediate figure

# CPU times
################################################################

cat(as.character(Sys.time()), " | Fig11: Plotting panel B...\n")

avg_ratite_len = mean(ratite_aln_stats$length)
med_ratite_len = median(ratite_aln_stats$length)
# Calculate average locus lengths

ratite_len_hist = ggplot(subset(ratite_aln_stats, length<500), aes(x=length)) +
  geom_histogram(fill=corecol(pal="wilke", numcol=1, offset=4), bins=50, color="#999999", size=0.1) +
  geom_vline(xintercept=med_ratite_len, size=0.5, linetype="11", color="#999999") +
  annotate("text", x=med_ratite_len+125, y=31000, label=paste("Median length: ", signif(med_ratite_len,4), "bp", sep=""), size=3, color="#999999") +
  scale_y_continuous(expand=c(0,0)) +
  xlab("locus length (bp)") +
  ylab("# loci") +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))
print(ratite_len_hist)

# Ratite aln length
################################################################

cat(as.character(Sys.time()), " | Fig11: Plotting panel C...\n")

ratite_scf_dist = ggplot(ratite_aln_stats, aes(x=node.scf.avg)) +
  geom_histogram(fill=corecol(pal="wilke", numcol=1, offset=4), bins=50, color="#999999", size=0.1) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Avg. sCF") +
  ylab("# loci") +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))
print(ratite_scf_dist)

if(save_int_figs){
  outfilename = paste(outdir, "ratite-scf-avgs", sep="")
  outfilename = paste(outfilename, ".png", sep="")
  cat(as.character(Sys.time()), " | Fig11: Saving panel C: ", outfilename, "\n")
  ggsave(outfilename, ratite_scf_dist, width=5, height=4, unit="in")
}
# Save the figure

# Ratite sCF avgs. per locus
################################################################

cat(as.character(Sys.time()), " | Fig11: Calculating runtimes based on sCF cutoffs...\n")

scf_avg_counts = data.frame()
cutoffs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
num_threads = c(14,8,16,32,64,128)
# Initialize parameters for runtime estimation

for(cutoff in cutoffs){
  num_gt = nrow(subset(ratite_aln_stats, node.scf.avg <= cutoff))
  # Get the number of loci with sCF below the current cutoff
  
  for(t in num_threads){
    #runtime = num_gt / t
    runtime = ((num_gt * 0.5) / t)
    # Assume half an hour per locus with GT method
    
    scf_avg_counts = rbind(scf_avg_counts, data.frame("cutoff"=cutoff, "num.loci"=num_gt, "threads"=t, "runtime"=runtime))
    # Add data for current cutoff and runtime to df
  }
}

scf_avg_counts$threads = factor(scf_avg_counts$threads)
# Factor the threads column

# Get expected number of loci and runtimes based on sCF cutoffs
################################################################

cat(as.character(Sys.time()), " | Fig11: Plotting panel D...\n")

scf_avg_loci_p = ggplot(scf_avg_counts, aes(x=cutoff, y=num.loci, color=threads)) +
  geom_point(size=3, color="#333333") +
  geom_line(size=1, color="#333333") +
  xlab("Avg. sCF cutoff") +
  ylab("Number of loci\nfor gene tree model") +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))
print(scf_avg_loci_p)

if(save_int_figs){
  outfilename = paste(outdir, "ratite-scf-avgs-numloci", sep="")
  outfilename = paste(outfilename, ".png", sep="")
  cat(as.character(Sys.time()), " | Fig11: Saving panel D: ", outfilename, "\n")
  ggsave(outfilename, scf_avg_loci_p, width=5, height=5, unit="in")
}
# Save the figure

# Ratite loci for gene tree model with varying sCF cutoffs
################################################################

cat(as.character(Sys.time()), " | Fig11: Plotting panel E...\n")

scf_avg_runtime_p = ggplot(scf_avg_counts, aes(x=cutoff, y=runtime, group=threads, color=threads)) +
  geom_point(size=3) +
  geom_line(size=1) +
  geom_hline(yintercept=8766, linetype="dashed", size=1, color="#333333") +
  #geom_hline(yintercept=720, linetype="dashed", size=1, color="#333333") +
  #geom_hline(yintercept=24, linetype="dashed", size=1, color="#333333") +
  xlab("Avg. sCF cutoff") +
  ylab("Runtime (hours)") +
  #scale_color_manual(values=(corecol(numcol=6, pal="wilke"))) +
  scale_color_viridis(name='# Threads', option = "D", discrete=T) +
  bartheme() +
  theme(legend.position="bottom",
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.box="vertical",
        plot.margin=unit(c(0.5,0.5,0.5,0), "cm")) +
  guides(colour = guide_legend(title.position="top"),
         linetype = guide_legend(override.aes = list(size = 2)))

print(scf_avg_runtime_p)

thread_legend = get_legend(scf_avg_runtime_p)

if(save_int_figs){
  outfilename = paste(outdir, "ratite-scf-avgs-runtime", sep="")
  outfilename = paste(outfilename, ".png", sep="")
  cat(as.character(Sys.time()), " | Fig11: Saving panel E: ", outfilename, "\n")
  ggsave(outfilename, scf_avg_runtime_p, width=5, height=5, unit="in")
}
# Save the figure

################################################################

cat(as.character(Sys.time()), " | Fig11: Combining panels...\n")

top = plot_grid(panel_a, ratite_len_hist, rel_widths=c(0.8,1), ncol=2, labels=c("A", "B"))
bottom = plot_grid(ratite_scf_dist, scf_avg_loci_p, scf_avg_runtime_p + theme(legend.position="none"), ncol=3, labels=c("C", "D", "E"))
last_leg = plot_grid(NULL, thread_legend, NULL, ncol=3, rel_widths=c(0.75, 0.2, 0.05))
fig = plot_grid(top, bottom, last_leg, nrow=3, rel_heights=c(1,1,0.2))

print(fig)

if(save_final_fig){
  figfile = paste(outdir, "fig11.pdf", sep="")
  cat(as.character(Sys.time()), " | Fig11: Saving figure: ", figfile, "\n")
  ggsave(figfile, fig, width=6.5, height=6, unit="in")
}


################################################################



