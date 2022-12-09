############################################################
# For PhyloAcc
# Figure 3
# Han Yan and Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

############################################################

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(PRROC)
library(here)
source(here("scripts", "lib", "design.r"))

############################################################

simdir = here("results", "simulations", "Simu1")
# The directory with the simulation files named according to
# case

save_fig = F
# Whether or not to save the final figure

# Options
############################################################

cat(as.character(Sys.time()), " | Fig7: Reading simulation data from ", simdir, "\n")

results = NULL
# Initialize the results df

for(f in list.files(simdir, pattern="elemLik"))
{
  # List all elemLik files in the simulation directory
  
  p = here(simdir, f)
  # Append the file name to the path of the simulation directory
  
  dat = read.table(p, header = T)
  dat = dat[, c("ID", "logBF1", "logBF2")]
  dat$logBF3 = dat$logBF1 - dat$logBF2
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

cat(as.character(Sys.time()), " | Fig7: Labeling simulation cases\n")

m0 = subset(results, case == 'trueM0')
m0 = select(m0, -c("case"))

m1 = subset(results, case != 'trueM0')
full = merge(m1, m0, by=c("method", "input", "ID2"))
full = data.table(full)

# Label the true m0 cases and merge by ID2
######################

cat(as.character(Sys.time()), " | Fig7: Calculating PR and AUPRC\n")

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
  
  cat(as.character(Sys.time()), " | -------> ", cur_method, " - ", cur_case, " - ", cur_input, "\n")
  
  cur_data = subset(full, method==cur_method & case==cur_case & input==cur_input)
  # Subset the data for the current method, case, and input
  
  if(nrow(cur_data) > 0){
    # Some combinations of method, case, and input don't exist, so only calulate AUPRC for those that do
    
    xx = unique(cur_data$logBF3.x)
    lx = length(xx)
    # The number of unique BF1 values for loci with M1 is true
    
    ly = length(cur_data$logBF3.y)
    # The number of loci where M0 is true    
    
    for(ratio in ratios){
      if(ratio <= ly/lx){
        cur_result = pr.curve(xx, cur_data$logBF3.y[1:(ratio * lx)], curve=T)
      }else{
        cur_result = pr.curve(xx, rep(cur_data$logBF3.y, ratio * lx/ly), curve=T)
      }
      cur_pr = cur_result$curve
      cur_auprc = cur_result$auc.davis.goadrich
      # Calculate AUPRC for the current ratio
      
      prs = rbind(prs, data.frame("method"=cur_method, "case"=cur_case, "input"=cur_input, "ratio"=ratio, "p"=cur_pr[,2], "r"=cur_pr[,1], "auprc"=cur_auprc))
      # Add the result to the prs df
      
    }
  }
}

# Calculate PR curves and AUPRC
######################

cat(as.character(Sys.time()), " | Fig7: Generating plots\n")

plot_cases = data.frame("case"=c('trueCase3', 'trueCase2', 'trueCase4'),
                        "input"=c('inputCase3', 'inputCase2', 'inputCase4'))
# The cases to plot

titles = c("Single acceleration\n\n", "Two independent\naccelerations\n", "Three independent\naccelerations\n", "Two independent\naccelerations with extra\ntarget lineages specified")

auprc_list = list()
prc_list = list()
# Lists to add figure grobs to

pr_ratio = 50
# The ratio at which to plot PR curves

for(i in 1:nrow(plot_cases)){
  cat(as.character(Sys.time()), " | -------> ", plot_cases[i,]$case, " - ", plot_cases[i,]$input, "\n")
  
  cur_pr = subset(prs, case == plot_cases[i,]$case & input == plot_cases[i,]$input)
  # For every scenario to plot, subset the pr data based on case and input
  
  p = ggplot(cur_pr, aes(x=ratio, y=auprc, color=method)) +
    geom_point(size=2) + 
    geom_line(size=1) +
    scale_y_continuous(limits=c(0,1,by=0.1)) +
    ggtitle(titles[i]) +
    xlab("") +
    ylab("") +
    scale_color_manual("Method", labels=c("Species tree", "Gene tree"), values=corecol(pal="wilke", numcol=2)) +
    bartheme() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=40, hjust=1, size=10),
          axis.title.y=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=8),
          plot.margin=margin(0.5,0.1,0,0.1, unit="cm"))
  # Generate the AUPRC
  
  if(i==1){
    p = p + ylab("AUPRC")
  }
  # For the first plot, add a title to the Y-axis
  
  if(i==4){
    p = p + theme(plot.title=element_text(vjust=3))
  }
  
  auprc_list[[i]] = p
  # Add the plot to the list of figures
  
  # Generate the AUPRC plot for the current scenario
  ##########
  
  cur_pr = subset(cur_pr, ratio == pr_ratio)
  
  p = ggplot(cur_pr, aes(x=r, y=p, color=method)) +
    #geom_point(size=2) + 
    geom_line(size=1) +
    #ggtitle("\n\n") +
    xlab("") +
    ylab("") +
    scale_color_manual("Method", labels=c("Species tree", "Gene tree"), values=corecol(pal="wilke", numcol=2)) +
    bartheme() +
    theme(legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.text.x=element_text(angle=40, hjust=1, size=10),
          #axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=8),
          plot.margin=margin(0.5,0.1,0,0, unit="cm"))
  # Render the PR curve
  
  if(i==1){
    p = p + ylab("Precision")
    p = p + theme(legend.position="bottom")
    fig_leg = get_legend(p)
    p = p + theme(legend.position="none")
  }
  # For the first plot, add a title to the y-axis, and get the legend to add back in later
  
  prc_list[[i]] = p
  # Add the plot to the list of figures
  
  # Generate the PR plot for the current scenario
  ##########
}

# Generate the plots for specified scenarios
######################

cat(as.character(Sys.time()), " | Fig7: Combining plots\n")

fig_top = plot_grid(plotlist=auprc_list, ncol=3, labels=c("A", "B", "C")) +
  draw_label("Ratio of non-accelerated to accelerated loci", x=0.5, y=0.05, vjust=-0.5, angle=0, size=10)
# Combine AUPRC plots with a single, shared x-axis title

fig_bot = plot_grid(plotlist=prc_list, ncol=3, labels=c("", "", "")) +
  draw_label("Recall", x=0.5, y=0.05, vjust=-0.5, angle=0, size=10)
# Combine PR plots with a single, shared x-axis title

fig_main = plot_grid(fig_top, fig_bot, nrow=2, rel_heights=c(1,0.87), align='vh')
# Combine plots

fig = plot_grid(fig_main, fig_leg, nrow=2, rel_heights=c(1,0.1))
# Add legend

print(fig)

# Combine the plots into the final figure
######################

if(save_fig){
  figfile = "../figs/fig7.png"
  cat(as.character(Sys.time()), " | Fig7: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=7.5, height=5.5, units="in")
}

# Save the figure
######################









