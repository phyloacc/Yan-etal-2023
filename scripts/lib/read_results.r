############################################################
# For PhyloAcc
# Scripts to read PhyloAcc and Beast results
# Gregg Thomas
############################################################

readTree <- function(tree_file) {
#   cat(as.character(Sys.time()), " | Fig4: Reading tree:", ratite_tree_file, "\n")

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

assignTargets <- function(new_node_labels){
  target = nontarget = list() #(1-1,2-8,2-5,2-6)
  target[[1]] = c()
  nontarget[[1]] = new_node_labels
  target[[2]] = c("aptHaa","aptOwe","aptRow","rheAme","rhePen","aptHaa.aptOwe", "aptHaa.aptRow","rheAme.rhePen")
  target[[3]] = c("aptHaa","aptOwe", "aptRow","casCas","droNov","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas")
  target[[4]] = c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","anoDid","strCam","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme")

  ## Target linaeges for each simulation
  ## 2 = case 2 = single paraphyletic = panel C
  ## 3 = case 3 = monophyletic = panel A
  ## 4 = case 4 = multi-paraphyletic = panel B
  ######################

  case = c("1-1","2-8","2-5","2-6")
  target_B = list()
  nontarget_B = list()
  target_B[[1]] = nontarget_B[[1]]=c()
  target_B[[2]] = c(10:14,19:21)
  target_B[[3]] = c(10:18)
  target_B[[4]] = c(1,10:22,24)
  for(i in 2:4){
    nontarget_B[[i]] = seq(1,31)
    nontarget_B[[i]] = nontarget_B[[i]][-target_B[[i]]]
  }
  ## Targets for the BEAST runs
  ######################  
  
  return(list(target, nontarget, target_B, nontarget_B))
}

#########################

readElemLik <- function(elem_lik_file, subset_cols=FALSE) {
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

readPostZ <- function(m1_postz_file, m2_postz_file, cur_elem_lik){

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
  
  if(length(m0_ind) > 0){
    opt_postz[m0_ind,] = 0
  }
  # Replace rows in the main st rate df that have optimal model 0 with 0
  
  if(length(m2_ind > 0)){
    opt_postz[m2_ind,] = m2_postz[m2_ind,]
  }

  return(opt_postz)

}

#########################

adjustCols <- function(st_opt_postz, gt_opt_postz, targets, old_colnames, new_colnames) {
    st_opt_postz = st_opt_postz[,-seq(1,17)]
    st_opt_postz = st_opt_postz[,grepl("_3", colnames(st_opt_postz))]

    gt_opt_postz = gt_opt_postz[,-seq(1,17)]
    gt_opt_postz = gt_opt_postz[,grepl("_3", colnames(gt_opt_postz))]
    # Get only the accelerated rates for non-outgroup branches

    colnames(st_opt_postz) = colnames(gt_opt_postz) = substr(colnames(gt_opt_postz), 1, str_length(colnames(gt_opt_postz))-2)
    # Remove the _3 label from the column names so columns match node labels

    outRelated=c('taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
    for(i in 1:length(outRelated)){
        ind = which(colnames(gt_opt_postz) == outRelated[i])
        gt_opt_postz = gt_opt_postz[,-ind]
        st_opt_postz = st_opt_postz[,-ind]
    }
    # Remove remaining branches leading to the outgroup

    target_indices = rep(NA, length(targets))
    for(i in 1:length(target_indices)){
        target_indices[i] = which(old_colnames == targets[i])
    }
    # Get the indices of the target branches

    colnames(gt_opt_postz) = colnames(st_opt_postz) = new_colnames
    # Rename columns with new node labels

    return(list(st_opt_postz, gt_opt_postz, target_indices))

}

#########################

readBeast <- function(beast_file, target_indices, old_colnames, new_colnames) {
  
  load(beast_file)
  # Load beast results as "res"

  beast_targets = res$node
  beast_non_targets = res$`non-node`
  # Get the results from the target and non-target branches

  for(i in 1:ncol(beast_targets)) {
    tmp = colnames(beast_targets)[i]

    if(str_length(tmp) > 6){
      colnames(beast_targets)[i] = paste0("(", new_colnames[which(old_colnames==substr(tmp,1,6))], ",", new_colnames[which(old_colnames==substr(tmp,8,8+6))], ")")
    }else{
      colnames(beast_targets)[i] = new_colnames[which(old_colnames==tmp)[1]]
    }
  }
  # Convert target names to new node labels
  
  for(i in 1:ncol(beast_non_targets)){
    tmp = colnames(beast_non_targets)[i]
    if(str_length(tmp) > 6){
      colnames(beast_non_targets)[i] = paste0("(", new_colnames[which(old_colnames==substr(tmp,1,6))], ",", new_colnames[which(old_colnames==substr(tmp,8,8+6))], ")")
    }else{
      colnames(beast_non_targets)[i] = new_colnames[which(old_colnames==tmp)[1]]
    }
  }
  # Convert non-target names to new node labels
  
  beast_targets = beast_targets[,colnames(gt_opt_postz)[target_indices]]
  if(c == 3){
    beast_non_targets = beast_non_targets[,seq(1,10)] #2-6
  }else{
    beast_non_targets = beast_non_targets[,seq(1,17)] #8 & 5 
  }
    
  beast_non_targets = beast_non_targets[,colnames(gt_opt_postz)[-target_indices]]
  rownames(beast_targets) = rownames(beast_non_targets) = as.character(seq(0,99)) 
  # Subset the data based on the current scenario

  return(list(beast_targets, beast_non_targets))
}

#########################