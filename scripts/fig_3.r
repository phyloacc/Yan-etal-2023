############################################################
# For rodent genomes
# Figure 5
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
#library(tidyr)
#library(ggbeeswarm)
#library(ggsignif)
#library(pROC)
library(PRROC)

library(here)
source(here("scripts", "lib", "design.r"))

############################################################

#calcPR <- function(chrome, dists, window_size, max_dist_mb){


############################################################


results = NULL

# read in data

simdir = here("results", "simulations", "Simu1")
# The directory with the simulation files named according to
# case

for(f in list.files(simdir, pattern="elemLik"))
{
  # List all elemLik files in the simulation directory
  
  p = here(simdir, f)
  # Append the file name to the path of the simulation directory
  
  dat = read.table(p, header = T)
  dat = dat[, c("ID", "logBF1", "logBF2")]
  # Read the current file and get only the columns of interest
  
  ff = strsplit(f, "_")[[1]]
  ff[4] = gsub(".txt","", ff[4])
  # Split the filename to get the current simulation case
  
  dat$case = ff[3]
  dat$input = ff[4]
  dat$method = ff[2]
  
  if(ff[3] == "trueM0"){
    dat$outcome = 0
  }else {
    dat$outcome = 1
  }
  # Add columns denoting the simulation case and method based
  # on the filename
  
  if(is.null(results))
  {
    results = dat
  }else{
    results = rbind(results, dat)
  }
  # Add the current data to the overall results data frame
  
}

results$ID2 = results$ID %% 100
# Assigns identical IDs to results with same method, case, and input but different models

# Read the input data
######################

m0 = subset(results, case == 'trueM0')
m0 = select(m0, -c("case"))

m1 = subset(results, case != 'trueM0')
full = merge(m1, m0, by=c("method", "input", "ID2"))
full = data.table(full)

# Label the truem0 cases and merge by ID2
######################

vars = list("method"=levels(as.factor(full$method)), "case"=levels(as.factor(full$case)), "input"=levels(as.factor(full$input)))
var_combos = expand.grid(vars)
# Get combinations of methods, cases, and inputs

ratios = c(1,10,20,50,70,100)
# The ratios to test for AUPRC

prs = data.frame("method"=c(), "case"=c(), "input"=c(), "ratio"=c(), "p"=c(), "r"=c(), "auprc"=c())
# A data frame to gather results for each combination of method, case, input, and ratio

for(i in 1:nrow(var_combos)){
  cur_method = as.character(var_combos[[1]][i])
  cur_case = as.character(var_combos[[2]][i])
  cur_input = as.character(var_combos[[3]][i])
  # Loop over all combinations of method, case, and input

  
  cur_data = subset(full, method==cur_method & case==cur_case & input==cur_input)
  # Subset the data for the current method, case, and input
  
  if(nrow(cur_data) > 0){
  # Some combinations of method, case, and input don't exist, so only calulate AUPRC for those that do
    
    xx = unique(cur_data$logBF1.x)
    lx = length(xx)
    # The number of unique BF1 values for loci with M1 is true
    
    ly = length(cur_data$logBF1.y)
    # The number of loci where M0 is true    

    for(ratio in ratios){
      if(ratio <= ly/lx){
        cur_result = pr.curve(xx, cur_data$logBF1.y[1:(ratio * lx)], curve=T)
      }else{
        cur_result = pr.curve(xx, rep(cur_data$logBF1.y, ratio * lx/ly), curve=T)
      }
      cur_pr = cur_result$curve
      cur_auprc = cur_result$auc.davis.goadrich
      # Calculate AUPRC for the current ratio
      
      prs = rbind(prs, data.frame("method"=cur_method, "case"=cur_case, "input"=cur_input, "ratio"=ratio, "p"=cur_pr[,2], "r"=cur_pr[,1], "auprc"=cur_auprc))
      # Add the result to the prs df
    }
  }
  
}


plot_cases = data.frame("case"=c('trueCase2', 'trueCase2', 'trueCase3', 'trueCase4'), "input"=c('inputCase2', 'inputCase4', 'inputCase3', 'inputCase4'))

auprc_list = list()
prc_list = list()

pr_ratio = 50

for(i in 1:nrow(plot_cases)){
  cur_pr = subset(prs, case == plot_cases[i,]$case & input == plot_cases[i,]$input)
  
  p = ggplot(cur_pr, aes(x=ratio, y=auprc, color=method)) +
    geom_point(size=2) + 
    geom_line(size=1) +
    scale_y_continuous(limits=c(0.7,1,by=0.1)) +
    xlab("") +
    ylab("") +
    scale_color_manual("Method", labels=c("Species tree", "Gene tree"), values=corecol(pal="wilke", numcol=2)) +
    bartheme() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=40, hjust=1, size=10),
          axis.title.y=element_text(size=10),
          plot.margin=margin(0.5,0.1,0,0.1, unit="cm"))
  
  if(i==1){
    p = p + ylab("AUPRC")
  }
  
  auprc_list[[i]] = p
  
  cur_pr = subset(cur_pr, ratio == pr_ratio)
  
  p = ggplot(cur_pr, aes(x=r, y=p, color=method)) +
    #geom_point(size=2) + 
    geom_line(size=1) +
    xlab("") +
    ylab("") +
    scale_color_manual("Method", labels=c("Species tree", "Gene tree"), values=corecol(pal="wilke", numcol=2)) +
    bartheme() +
    theme(legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.text.x=element_text(angle=40, hjust=1, size=10),
          axis.title.y=element_text(size=10),
          plot.margin=margin(0.5,0.1,0,0.1, unit="cm"))
  
  if(i==1){
    p = p + ylab("Precision")
    p = p + theme(legend.position="bottom")
    fig_leg = get_legend(p)
    p = p + theme(legend.position="none")
  }
  
  prc_list[[i]] = p
}
#print(p)

fig_top = plot_grid(plotlist=auprc_list, ncol=4, labels=c("A", "B", "C", "D")) +
  draw_label("Ratio of M0:M1 loci", x=0.5, y=0.05, vjust=-0.5, angle=0, size=10)

fig_bot = plot_grid(plotlist=prc_list, ncol=4, labels=c("", "", "", "")) +
  draw_label("Recall", x=0.5, y=0.05, vjust=-0.5, angle=0, size=10)

fig_main = plot_grid(fig_top, fig_bot, nrow=2)

fig = plot_grid(fig_main, fig_leg, nrow=2, rel_heights=c(1,0.1))

save_fig = T
if(save_fig){
  figfile = "../figs/fig3.png"
  cat(as.character(Sys.time()), " | Fig3: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=7.5, height=4, units="in")
}
# Save figure

# plot AUPRC
prs = prs[ (prs$case == 'trueCase2' & prs$input == 'inputCase2') | (prs$case == 'trueCase2' & prs$input == 'inputCase4') | 
             (prs$case == 'trueCase3' & prs$input == 'inputCase3') |  (prs$case == 'trueCase4' & prs$input == 'inputCase4') ,]


p = ggplot() + geom_line(data = prs2, aes(x = ratio, y = AUPRC, color = method)) +
  facet_wrap(vars(input,case), nrow=2, ncol=2,  scales="fixed") + 
  #facet_grid(rows=vars(case), cols = vars(input),  space="fixed")  + 
  bartheme()

print(p)







# prs2 = full[, {
#   list("AUPRC" = sapply(c(1,10,20,50,70,100), function(x) {
#     print(paste(method, case, input))
#     print(length(logBF1.x))
#     print(x)
#     xx = unique(logBF1.x)
#     ly = length(logBF1.y)
#     lx = length(xx)
#     
#     #print(ly)
#     #print(lx)
#     if(x <= ly/lx)
#     {
#       pr.curve(xx, logBF1.y[1:(x * lx)])$auc.davis.goadrich
#     }else{
#       pr.curve(xx, rep(logBF1.y, x * lx/ly))$auc.davis.goadrich
#     }
#   }),
#   "ratio" = c(1,10,20,50,70,100))
# }, by = c('method', 'case', 'input')]







