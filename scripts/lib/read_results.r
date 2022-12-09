############################################################
# For PhyloAcc
# Scripts to read PhyloAcc and Beast results
# Gregg Thomas
############################################################

readTree <- function(tree_file) {
# In case we need to read the species tree in nexus format (unused)

  ratite_tree = read.tree(file=ratite_tree_file)

  tipname_nexus = c('Moa', 'Kiwi1', 'Kiwi2', 'Kiwi3', 'Cassowary','Tinamou1','Emu', 'Tinamou3','Chicken','Eagle', 'Ibis','Tinamou4','Rhea1','Rhea2','Ostrich','Tinamou2')

  name_nexus = c(tipname_nexus[1],tipname_nexus[c(6,16)],'(Tinamou1,Tinamou2)',tipname_nexus[c(8,12)],'(Tinamou3,Tinamou4)','(Tinamou1,Tinamou3)','(Tinamou1,Moa)',
                 tipname_nexus[c(2,3)],'(Kiwi1,Kiwi2)',tipname_nexus[4],'(Kiwi1,Kiwi3)',tipname_nexus[c(5,7)],'(Cassowary,Emu)','(Kiwi1,Cassowary)',
                 tipname_nexus[c(13,14)],'(Rhea1,Rhea2)','(Kiwi1,Rhea1)','(Kiwi1,Tinamou1)',tipname_nexus[15],'(Kiwi1,Ostrich)',tipname_nexus[9],tipname_nexus[c(10,11)],
                 'taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')

  return(ratite_tree)
}

#########################

adjustLabels <- function(mapping_file, orig_node_labels) {
# Reads label map file to convert old ratite node labels to new generic node labels for simulation cases

  label_map = read.csv(mapping_file, header=T, sep="\t")
  # Read the label map file

  new_node_labels = orig_node_labels
  # Initialize the new labels

  for(i in 1:13){
    new_node_labels[i]=label_map[which(label_map[,1]==orig_node_labels[i]),2]
  }
  # Adjust the tips

  for(i in 14:length(orig_node_labels)){
    t1=label_map[which(label_map[,1]==substr(orig_node_labels[i],1,6)),2][1]
    t2=label_map[which(label_map[,1]==substr(orig_node_labels[i],8,8+6)),2]
    new_node_labels[i]=paste0("(", t1, ",", t2, ")")
  }
  # Adjust the internal nodes

  return(new_node_labels)
}

#########################

assignTargets <- function(old_node_labels, new_node_labels){
# Compiles lists and label mappings of target branches for each simulation scenario

  target_list = list()
  target_list[[2]] = c("aptHaa","aptOwe","aptRow","rheAme","rhePen","aptHaa.aptOwe", "aptHaa.aptRow","rheAme.rhePen")
  target_list[[3]] = c("aptHaa","aptOwe", "aptRow","casCas","droNov","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas")
  target_list[[4]] = c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","anoDid","strCam","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme")
  # Lists of target branches with old node labels for:
  # Case 2 / scenario 8 / label B / two independent accelerations
  # Case 3 / scenario 5 / label A / single acceleration
  # Case 4 / scenario 6 / label C / three independent accelerations

  case_targets = list()
  # Initialize new case_targets list that stores targets and non-targets for each case with old and new node labels

  for(c in c(2,3,4)){
    case_targets[[c]] = list()
    # The current case will have its own list of targets and non-targets

    cur_targets = target_list[[c]]
    cur_target_inds = match(cur_targets, old_node_labels)
    cur_targets_new = new_node_labels[cur_target_inds]
    # Get the target labels from the current case and match with new labels

    targets = data.frame("old.labels"=cur_targets, "inds"=cur_target_inds, "new.labels"=cur_targets_new)
    # Compile the old and new labels into a data frame

    case_targets[[c]][["targets"]] = targets
    # Add the targets to the current case_targets list

    cur_non_targets = old_node_labels[!old_node_labels %in% cur_targets]
    cur_non_target_inds = match(cur_non_targets, old_node_labels)
    cur_non_targets_new = new_node_labels[cur_non_target_inds]
    # Get the non-target labels from the current case and match with new labels

    non_targets = data.frame("old.labels"=cur_non_targets, "inds"=cur_non_target_inds, "new.labels"=cur_non_targets_new)
    # Compile the old and new labels into a data frame

    case_targets[[c]][["non.targets"]] = non_targets
    # Add the non-targets to the current case_targets list
  }

  return(case_targets)
}

#########################

readElemLik <- function(elem_lik_file, subset_cols=FALSE) {
# Reads an elem_lik file from PhyloAcc and determines model with max probability

    elem_lik = as.matrix(read.table(elem_lik_file, header=T))

    if(subset_cols){
        elem_lik = elem_lik[,1:7]
    }

    optM = apply(elem_lik[,3:5], 1, function(x) which(x==max(x))[1]-1)
    elem_lik = cbind(elem_lik, optM)
    elem_lik = as.data.frame(elem_lik)
    # Reading the results and getting adding the model with
    # the max probability

    return(elem_lik)
}

#########################

readPostZ <- function(m1_postz_file, m2_postz_file, cur_elem_lik, use_m2=FALSE){
# Reads a postZ files from PhyloAcc and combines them in a df with the probabilities
# from the model that has the highest likelihood

  m1_postz = read.table(m1_postz_file, header=T, row.names=1)
  m1_postz = m1_postz[order(as.numeric(rownames(m1_postz))),]
  # Read the st M1 rate file and order by row name
  
  m2_postz = read.table(m2_postz_file, header=T, row.names=1)
  m2_postz = m2_postz[order(as.numeric(rownames(m2_postz))),]
  # Read the st M2 rate file and order by row name
    
  opt_postz = m1_postz
  # Initialize the main st rate df
  
  m0_ind = which(cur_elem_lik$optM == 0)
  m2_ind = which(cur_elem_lik$optM == 2)
  # Get the row numbers of the loci that have optimal model 0 or 2
  
  if(use_m2){
    opt_postz = m2_postz
    # For some instances comparing M2 performance, just use all rows from
    # the m2 file
  }else{
    if(length(m0_ind) > 0){
        opt_postz[m0_ind,] = 0
    }
    # Replace rows in the main df that have optimal model 0 with 0
    
    if(length(m2_ind > 0)){
        opt_postz[m2_ind,] = m2_postz[m2_ind,]
    }
    # Replace the rows in the main df that have optimal model 2 witih
    # the probabilities from model 2
  }

  return(opt_postz)
}

#########################

adjustCols <- function(st_opt_postz, gt_opt_postz, old_colnames, new_colnames) {
# Selects probabilities for each branch in state 2 and removes certain columns

    st_opt_postz = st_opt_postz[,-seq(1,17)]
    st_opt_postz = st_opt_postz[,grepl("_3", colnames(st_opt_postz))]

    gt_opt_postz = gt_opt_postz[,-seq(1,17)]
    gt_opt_postz = gt_opt_postz[,grepl("_3", colnames(gt_opt_postz))]
    # Get only the probabilites for z=2 for non-outgroup branches

    colnames(st_opt_postz) = colnames(gt_opt_postz) = substr(colnames(gt_opt_postz), 1, str_length(colnames(gt_opt_postz))-2)
    # Remove the _3 label from the column names so columns match node labels

    outRelated=c('taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
    st_opt_postz = select(st_opt_postz, -outRelated)
    gt_opt_postz = select(gt_opt_postz, -outRelated)
    # Remove remaining branches leading to the outgroup

    colnames(gt_opt_postz) = colnames(st_opt_postz) = new_colnames
    # Rename columns with new node labels

    return(list(st_opt_postz, gt_opt_postz))
}

#########################

readBeast <- function(beast_file, cur_case_targets, old_colnames, new_colnames, theta=FALSE) {
# Reads output from BEAST runs already processed into an R data format

  load(beast_file)
  # Load beast results as "res"

  beast_targets = as.data.frame(res$node)
  beast_non_targets = as.data.frame(res$`non-node`)
  # Get the results for the target and non-target lineages

  if(!theta){
    beast_targets = select(beast_targets, cur_case_targets[["targets"]]$old.labels)
    names(beast_targets) = cur_case_targets[["targets"]]$new.labels
    # Select target lineages and rename with new labels

    beast_non_targets = select(beast_non_targets, cur_case_targets[["non.targets"]]$old.labels)
    names(beast_non_targets) = cur_case_targets[["non.targets"]]$new.labels
    # Select non-target lineages and rename with new labels
  }else{
    # For the theta figures (10), the beast data is loaded slightly differently

    beast_all = cbind(beast_targets, beast_non_targets)
    beast_all = beast_all[, old_colnames]
    # Combine the target and non-target lineages and select the correct rows

    names(beast_all) = new_colnames
    # Rename columns with new labels

    beast_targets = select(beast_all, cur_case_targets[["targets"]]$new.labels)
    beast_non_targets = select(beast_all, cur_case_targets[["non.targets"]]$new.labels)
    # Separate target and non-target lineages
  }

  rownames(beast_targets) = rownames(beast_non_targets) = as.character(seq(0,99)) 
  # Add common row names

  return(list(beast_targets, beast_non_targets))
}

#########################