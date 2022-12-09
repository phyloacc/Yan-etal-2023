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

simdir = here("results", "simulations", "ThetaNull_elemLik")
# The directory with the simulation files named according to
# case

save_fig = F
# Whether or not to save the final figure

# Options
############################################################

cat(as.character(Sys.time()), " | Fig9: Reading simulation data from ", simdir, "\n")

results = NULL
# Initialize the results df

for(theta_scale in c(3,6,10)){
  theta_dir = paste0(simdir, "/theta", theta_scale)
  
  cat(as.character(Sys.time()), " | ----> ", theta_scale, " - ", theta_dir, "\n")
  
  for(f in list.files(theta_dir, pattern="elemLik"))
  {
    # List all elemLik files in the simulation directory
    
    p = here(theta_dir, f)
    # Append the file name to the path of the simulation directory
    
    cat(as.character(Sys.time()), " | -------> ", theta_scale, " - ", p, "\n")
    
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
    
    dat$theta = theta_scale
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
  
}

results$ID2 = results$ID %% 100
# Assigns identical IDs to results with same method, case, and input but different models

# Read the input data
######################

cat(as.character(Sys.time()), " | Fig9: Labeling simulation cases\n")

m0 = subset(results, case == 'trueM0')
m0 = select(m0, -c("case"))

m1 = subset(results, case != 'trueM0')
full = merge(m1, m0, by=c("method", "input", "theta", "ID2"))
full = data.table(full)

# Label the truem0 cases and merge by ID2
######################

cat(as.character(Sys.time()), " | Fig9: Calculating PR and AUPRC\n")

vars = list("method"=levels(as.factor(full$method)), "case"=levels(as.factor(full$case)), "input"=levels(as.factor(full$input)), "theta"=levels(as.factor(full$theta)))
var_combos = expand.grid(vars)
# Get combinations of methods, cases, and inputs

ratios = c(1,10,20,50,70,100)
# The ratios to test for AUPRC

prs = data.frame("method"=c(), "case"=c(), "input"=c(), "theta"=c(), "ratio"=c(), "p"=c(), "r"=c(), "auprc"=c())
# A data frame to gather results for each combination of method, case, input, and ratio

for(i in 1:nrow(var_combos)){
  cur_method = as.character(var_combos[[1]][i])
  cur_case = as.character(var_combos[[2]][i])
  cur_input = as.character(var_combos[[3]][i])
  cur_theta = as.character(var_combos[[4]][i])
  # Loop over all combinations of method, case, and input
  
  cat(as.character(Sys.time()), " | -------> ", cur_method, " - ", cur_case, " - ", cur_input, " - ", cur_theta, "\n")
  
  cur_data = subset(full, method==cur_method & case==cur_case & input==cur_input & theta==cur_theta)
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
      
      prs = rbind(prs, data.frame("method"=cur_method, "case"=cur_case, "input"=cur_input, "theta"=cur_theta, "ratio"=ratio, "p"=cur_pr[,2], "r"=cur_pr[,1], "auprc"=cur_auprc))
      # Add the result to the prs df
    }
  }
}

# Calculate PR curves and AUPRC
######################

cat(as.character(Sys.time()), " | Fig3: Generating plots\n")

plot_cases = data.frame("case"=c('trueCase3', 'trueCase2', 'trueCase4', 'trueCase3', 'trueCase2', 'trueCase4', 'trueCase3', 'trueCase2', 'trueCase4'),
                        "input"=c('inputCase3', 'inputCase2', 'inputCase4', 'inputCase3', 'inputCase2', 'inputCase4', 'inputCase3', 'inputCase2', 'inputCase4'),
                        "theta"=c(3, 3, 3, 6, 6, 6, 10, 10, 10))
# The cases to plot

titles = c("Single acceleration\n\n", "Two independent\naccelerations\n", "Three independent\naccelerations\n", "\n\n", "\n\n", "\n\n", "\n\n", "\n\n", "\n\n")
# Plot titles

auprc_list = list()
prc_list = list()
# Lists to add figure grobs to

pr_ratio = 50
# The ratio at which to plot PR curves

for(i in 1:nrow(plot_cases)){
  cat(as.character(Sys.time()), " | -------> ", plot_cases[i,]$case, " - ", plot_cases[i,]$input, " - ", plot_cases[i,]$theta, "\n")
  
  cur_pr = subset(prs, case == plot_cases[i,]$case & input == plot_cases[i,]$input & theta == plot_cases[i,]$theta)
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
          plot.margin=margin(0.1,0.1,-1,0.1, unit="cm"))
  # Generate the AUPRC
  
  # if(plot_cases[i,]$input == "inputCase4"){
  #   t = as.character(plot_cases[i,]$theta)
  #   p = p + ylab(bquote(paste(.(t),"x", theta))) +
  #     scale_y_continuous(limits=c(0,1,by=0.1), position = 'right', sec.axis = dup_axis()) + 
  #     theme(axis.title.y.left = element_blank(),
  #           axis.title.y.right = element_text(angle=0, vjust=0.5),
  #           axis.ticks.y.right = element_blank(),
  #           axis.text.y.right = element_blank(),
  #           axis.line.y.right = element_blank())
  # }
  # If we want the theta values on the right (messes up scaling)
  
  if(i==1 || i==4 || i==7){
    t = as.character(plot_cases[i,]$theta)
    p = p + ylab(bquote(paste("AUPRC ", .(t),"x", theta)))
  }
  # For the first plot, add a title to the Y-axis
  
  #if(i==4){
  #  p = p + theme(plot.title=element_text(vjust=3))
  #}
  # Not sure what this did...
  
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
  
  # Generate the PR plot for the current scenario (not plotted for Fig 9)
  ##########
}

# Generate the plots for specified scenarios
######################

cat(as.character(Sys.time()), " | Fig9: Combining plots\n")


fig_main = plot_grid(plotlist=auprc_list, ncol=3, nrow=3)
# Combine panels 

fig = plot_grid(fig_main, fig_leg, nrow=2, rel_heights=c(1,0.1)) +
  draw_label("Ratio of non-accelerated to accelerated loci", x=0.5, y=0.075, angle=0, size=10)
# Add legend and x-axis title

print(fig)

# Combine the plots into the final figure
######################

if(save_fig){
  figfile = "../figs/fig9.png"
  cat(as.character(Sys.time()), " | Fig9: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig, width=6, height=6 ,units="in")
}

# Save the figure
######################









