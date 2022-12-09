####same directory as the one for Figure 3 (well-specified case)
## same methodology. Just add a BF3, replacing BF2

library(ggplot2)
#library(pROC)
library(PRROC)
library(data.table)
# read in data
###########################AUPRC theta, fig 9

# prefix11="C:/bin/phyloacc/Yan-etal-2022/results/simulations/ThetaNull_elemLik/"
# setwd(prefix11)
# 
# # read in data
# results = NULL
# for(th in c(3,6,10)){
#   dir1=paste0(prefix11,'theta',th)
#   setwd(dir1)
#   for(f in list.files(pattern="elemLik"))
#   {
#     dat = read.table(f, header = T)
#     dat = dat[, c("ID", "logBF1", "logBF2")]
#     ff = strsplit(f, "_")[[1]]
#     
#     ff[4] = gsub(".txt","", ff[4])
#     dat$case = ff[3]
#     dat$input = ff[4]
#     if(ff[2]=="Acc"){
#       dat$method = "PhyloAcc"
#     }else{
#       dat$method = "PhyloAcc-GT"
#     }
#     
#     dat$outcome = 1
#     dat$theta=paste0(th,'theta')
#     
#     if(is.null(results))
#     {
#       results = dat
#     }else{
#       results = rbind(results, dat)
#     }
#   }
# }
# 
# results$ID2 = results$ID %% 100
# 
# prs_list=list()
# prs_list[[1]]=NA
# ca2 = c('A','B','C')
# for(ca in 2:4){
#   ## test logBF1, null: M0, alternative: case4, case2; input: Case4
#   dat = data.table(results)
#   dat[dat$case == 'trueM0', 'outcome'] = 0
#   dat_null = dat[case == 'trueM0' & input==paste0('inputCase',ca)]
#   dat_null[,c('case') := NULL]
#   dat_alter = dat[case == paste0('trueCase',ca)]
#   dat_alter = merge(dat_alter, dat_null, by = c('method', 'input','theta', 'ID2'), allow.cartesian = T)
#   
#   prs = dat_alter[, {
#     list("AUPRC" = sapply(c(1,10,20,50,70,100), function(x) {
#       xx = unique(logBF1.x)
#       ly = length(logBF1.y)
#       lx = length(xx)
#       
#       #print(ly)
#       #print(lx)
#       if(x <= ly/lx)
#       {
#         pr.curve(xx, logBF1.y[1:(x * lx)])$auc.davis.goadrich
#       }else{
#         pr.curve(xx, rep(logBF1.y, x * lx/ly))$auc.davis.goadrich
#       }
#     }),
#     "ratio" = c(1,10,20,50,70,100))
#   }, by = c('method', 'case','input', 'theta')]
#   
#   prs_list[[ca]] = prs[ (case == paste0('trueCase',ca) & theta == '3theta') | (case == paste0('trueCase',ca)  & theta == '6theta') | (case == paste0('trueCase',ca)  & theta == '10theta') ]
#   
# }
# 
# prs2 = prs_list
# for(j in 2:4){
#   for(i in 1:nrow(prs2[[j]])){
#     if(prs2[[j]][['case']][i]=='trueCase2'){
#       prs2[[j]][['case']][i]='trueCaseA'
#     }else if(prs2[[j]][['case']][i]=='trueCase3'){
#       prs2[[j]][['case']][i]='trueCaseB'
#     }else if(prs2[[j]][['case']][i]=='trueCase4'){
#       prs2[[j]][['case']][i]='trueCaseC'
#     }
#     
#     if(prs2[[j]][['input']][i]=='inputCase2'){
#       prs2[[j]][['input']][i]='inputCaseA'
#     }else if(prs2[[j]][['input']][i]=='inputCase3'){
#       prs2[[j]][['input']][i]='inputCaseB'
#     }else if(prs2[[j]][['input']][i]=='inputCase4'){
#       prs2[[j]][['input']][i]='inputCaseC'
#     }
#   }
#   names(prs2[[j]])[5:6]=c('AreaUnderPrecisionRecallCurve','No.AllConservedElement_No.AcceleratedElements_Ratio')
# }
# 
# plots=list()
# plots[[1]]=NA
# for(ca in 2:4){
#   plots[[ca]]=ggplot() + geom_line(data = prs2[[ca]], aes(x = No.AllConservedElement_No.AcceleratedElements_Ratio, y = AreaUnderPrecisionRecallCurve, color = method)) +
#     facet_wrap(vars(factor(theta, levels=c('3theta','6theta','10theta'))), nrow=1, ncol=3) + 
#     theme_minimal(base_size=14) +ggtitle(paste0('case ',ca2[ca-1]))+
#     theme(plot.title = element_text(size = 14)) + ylim(0,1)# ,strip.text.x=element_blank()
# }
# 
# figure=ggarrange(plots[[2]]+rremove('xlab')+rremove('ylab'),plots[[3]]+rremove('xlab'),plots[[4]]+rremove('ylab'),common.legend=T,legend="bottom",nrow=3)
# 
# annotate_figure(figure)#,top=text_grob("AUPRC at different theta scales and cases "))
# print(figure)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# stop("-------")
# 
# 
# 
# 
# 
# 
# setwd("C:/bin/phyloacc/Yan-etal-2022/results/simulations/Theta_Scaling")
# #get GT accuracy:
# prefix="C:/bin/phyloacc/Yan-etal-2022/results/simulations/Theta_Scaling/Theta"
# 
# ##setup
# caseLabel=c(8,5,6)
# nameLabel=c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","cryCin","tinGut","eudEle", "notPer","anoDid","strCam",
#             "aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme","cryCin.tinGut",
#             "eudEle.notPer","cryCin.eudEle","cryCin.anoDid", "aptHaa.cryCin","aptHaa.strCam")
# birdname=read.table('C:/bin/phyloacc/Yan-etal-2022/data/SimulatedData/Name_Label_Mapping.txt',sep='\t',stringsAsFactors = F)
# 
# namePrint=nameLabel
# for(i in 1:13){
#   namePrint[i]=birdname[which(birdname[,1]==nameLabel[i]),2]
# }
# for(i in 14:length(nameLabel)){
#   t1=birdname[which(birdname[,1]==substr(nameLabel[i],1,6)),2][1]
#   t2=birdname[which(birdname[,1]==substr(nameLabel[i],8,8+6)),2]
#   namePrint[i]=paste0("(",t1,",",t2,")")
# }
# 
# target=nontarget=list() #(1-1,2-8,2-5,2-6)
# target[[1]]=c()
# nontarget[[1]]=namePrint
# target[[2]]=c("aptHaa","aptOwe","aptRow","rheAme","rhePen","aptHaa.aptOwe", "aptHaa.aptRow","rheAme.rhePen")
# target[[3]]=c("aptHaa","aptOwe", "aptRow","casCas","droNov","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas")
# target[[4]]=c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","anoDid","strCam","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme")
# 
# #########BEAST
# case=c("1-1","2-8","2-5","2-6")
# target_B=list()
# nontarget_B=list()
# target_B[[1]]=nontarget_B[[1]]=c()
# target_B[[2]]=c(10:14,19:21)
# target_B[[3]]=c(10:18)
# target_B[[4]]=c(1,10:22,24)
# for(i in 2:4){
#   nontarget_B[[i]]=seq(1,31)
#   nontarget_B[[i]]=nontarget_B[[i]][-target_B[[i]]]
# }
# 
# spt="(((halLeu:0.05373562#0,nipNip:0.05176648#0)taeGut-aptFor:0.04155687#0.04578784,galGal:0.1660553#0)taeGut-galGal:0.0406841#0.04578784,((((((aptHaa:0.00138725#0,aptOwe:0.0016341#0)aptHaa-aptOwe:0.00305054#0.006321079,aptRow:0.00410494#0)aptHaa-aptRow:0.0277496#0.019876606,(casCas:0.0115461#0,droNov:0.0137332#0)casCas-droNov:0.0273773#0.029694415)aptHaa-casCas:0.0028599#0.086448824,(rheAme:0.00469588#0,rhePen:0.00533574#0)rheAme-rhePen:0.0566382#0.027801219)aptHaa-rheAme:0.00185668#0.058734321,(((cryCin:0.0470926#0,tinGut:0.0388556#0)cryCin-tinGut:0.0172068#0.037018559,(eudEle:0.0655012#0,notPer:0.073059#0)eudEle-notPer:0.0079941#0.057911894)cryCin-eudEle:0.0672093#0.095580579,anoDid:0.0492722#0)cryCin-anoDid:0.0253604#0.053495448)aptHaa-cryCin:0.0118742#0.064943475,strCam:0.051388#0)aptHaa-strCam:0.0406969#0.011627686)taeGut-aptHaa:0#0.04578784;"
# sptree=ape::read.tree(text=spt)
# #tipname_nexus2=c('anoDid','aptHaa', 'aptOwe', 'aptRow', 'casCas',   'cryCin','droNov','eudEle', 'galGal','halLeu','nipNip','notPer','rheAme','rhePen','strCam','tinGut')
# tipname_nexus= c('Moa',    'Kiwi1', 'Kiwi2', 'Kiwi3',  'Cassowary','Tinamou1','Emu', 'Tinamou3','Chicken','Eagle', 'Ibis','Tinamou4','Rhea1','Rhea2','Ostrich','Tinamou2')
# name_nexus=c(tipname_nexus[1],tipname_nexus[c(6,16)],'(Tinamou1,Tinamou2)',tipname_nexus[c(8,12)],'(Tinamou3,Tinamou4)','(Tinamou1,Tinamou3)','(Tinamou1,Moa)',
#              tipname_nexus[c(2,3)],'(Kiwi1,Kiwi2)',tipname_nexus[4],'(Kiwi1,Kiwi3)',tipname_nexus[c(5,7)],'(Cassowary,Emu)','(Kiwi1,Cassowary)',
#              tipname_nexus[c(13,14)],'(Rhea1,Rhea2)','(Kiwi1,Rhea1)','(Kiwi1,Tinamou1)',tipname_nexus[15],'(Kiwi1,Ostrich)',tipname_nexus[9],tipname_nexus[c(10,11)],
#              'taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
# sptree_name_all=c(sptree$tip.label,sptree$node.label)
# 
# 
# case = c("2-8" ,"2-5" ,"2-6")
# case2 = c('A','B','C')
# 
# c=2 #change to 1:3
# plot_list_tar=c()
# plot_list_nonTar=c()
# for(t in c(3,6,10)){
#   noGT_M1=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_rate_postZ_M1.txt"), header=T, row.names = 1)
#   noGT_M2=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_rate_postZ_M2.txt"), header=T, row.names = 1)
#   noGT=noGT_M1
#   elem_Acc=as.matrix(read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_elem_lik.txt"),header=T))
#   optM= apply(elem_Acc[,3:5],1,function(x) which(x==max(x))[1])
#   elem_Acc=cbind(elem_Acc, optM)
#   ind=which(elem_Acc[,'optM']==0)
#   if(length(ind)>0) noGT[ind,]=0
#   ind=which(elem_Acc[,'optM']==2)
#   if(length(ind)>0) noGT[ind,]=noGT_M2[ind,]
#   
#   GT_M1=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_rate_postZ_M1.txt"), header=T, row.names = 1)
#   GT_M2=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_rate_postZ_M2.txt"), header=T, row.names = 1)
#   GT=GT_M1
#   elem_GT=as.matrix(read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_elem_lik.txt"),header=T))
#   optM=apply(elem_GT[,3:5],1,function(x) which(x==max(x))[1])
#   elemGT =cbind(elem_GT,optM)
#   ind=which(elemGT[,'optM']==0)
#   if(length(ind)>0) GT[ind,]=0
#   ind=which(elemGT[,'optM']==2)
#   if(length(ind)>0) GT[ind,]=GT_M2[ind,]
#   
#   noGT=noGT[,-seq(1,17)]
#   GT=GT[,-seq(1,17)]
#   noGT=noGT[,grepl("_3", colnames(noGT))]
#   GT=GT[,grepl("_3", colnames(GT))]
#   colnames(noGT)=colnames(GT)=substr(colnames(GT),1,str_length(colnames(GT))-2) #checked noGT, GT share col names
#   #remove br leading to out
#   outRelated=c('taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
#   for(i in 1:length(outRelated)){
#     ind=which(colnames(GT)==outRelated[i])
#     GT=GT[,-ind]
#     noGT=noGT[,-ind]
#   }
#   
#   targ_ind=rep(NA,length(target[[c+1]]))
#   for(i in 1:length(targ_ind)){
#     targ_ind[i]=which(nameLabel==target[[c+1]][i])
#   }
#   colnames(GT)=colnames(noGT)=namePrint
#   
#   #BEAST
#   load(paste0(prefix,t,'/Result_Case',case2[c],'/BEAST/result1_',case[c],'.RData')) #Result8/
#   res2=cbind(res$node,res$`non-node`)
#   #edit names to be consistent with GT names
#   for(i in 1:ncol(res2)){
#     tmp=colnames(res2)[i]
#     if(str_length(tmp)>6){
#       colnames(res2)[i]=paste0("(",namePrint[which(nameLabel==substr(tmp,1,6))],",",namePrint[which(nameLabel==substr(tmp,8,8+6))],")")
#     }else{
#       colnames(res2)[i]=namePrint[which(nameLabel==tmp)[1]]
#     }
#   }
#   res_backup=res2[,colnames(GT)[targ_ind]]
#   res_backup2=res2[,colnames(GT)[-targ_ind]]
#   rownames(res_backup)=rownames(res_backup2)=as.character(seq(0,99))
#   
#   #######draw boxplot
#   GTZ=GT[,targ_ind]
#   noGTZ=noGT[,targ_ind]
#   dat = cbind(reshape2::melt(GTZ), "Method" = rep("PhyloAcc-GT", ncol(GTZ)*nrow(GTZ))) #2-8: 8  #2-5: 8
#   dat = rbind(dat, cbind(reshape2::melt(noGTZ), "Method" = rep("PhyloAcc", ncol(noGTZ)*nrow(noGTZ))))
#   tmp=cbind(reshape2::melt(res_backup),"Method"=rep("*BEAST2", ncol(res_backup)*nrow(res_backup)))[,-1]
#   colnames(tmp)=colnames(dat)
#   dat = rbind(dat,tmp)
#   colnames(dat)[1:2] = c("species", "acc")
#   #dat$Method = factor(dat$Method)
#   dat$Method = factor(dat$Method, levels=c("PhyloAcc-GT", "PhyloAcc", "*BEAST2"))
#   plot_list_tar[[as.character(t)]]=ggplot(dat, aes(x = species, y = acc, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 12) +xlab("") + ylab("P(acceleration)") +theme(axis.text.x = element_text(angle = 30,hjust = 0.8)) + 
#     scale_fill_manual(values=c("#CC79A7","#56B4E9","#F0E442"))+ggtitle(paste0("accelerated branches, Theta=",t,"x")) #"#E69F00", 
#   
#   GTZ=GT[,-targ_ind]
#   noGTZ=noGT[,-targ_ind]
#   dat = cbind(reshape2::melt(GTZ), "Method" = rep("PhyloAcc-GT", ncol(GTZ)*nrow(GTZ))) #2-8: 8  #2-5: 8
#   dat = rbind(dat, cbind(reshape2::melt(noGTZ), "Method" = rep("PhyloAcc", ncol(noGTZ)*nrow(noGTZ))))
#   tmp=cbind(reshape2::melt(res_backup2),"Method"=rep("*BEAST2", ncol(res_backup2)*nrow(res_backup2)))[,-1]
#   colnames(tmp)=colnames(dat)
#   dat = rbind(dat,tmp)
#   colnames(dat)[1:2] = c("species", "acc")
#   #dat$Method = factor(dat$Method)
#   dat$Method = factor(dat$Method, levels=c("PhyloAcc-GT", "PhyloAcc", "*BEAST2"))
#   plot_list_nonTar[[as.character(t)]]=ggplot(dat, aes(x = species, y = acc, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 12) +xlab("")+ ylab("P(acceleration)") +theme(axis.text.x = element_text(angle = 30,hjust = 0.8)) + 
#     scale_fill_manual(values=c("#CC79A7","#56B4E9","#F0E442"))+ggtitle(paste0("non-accelerated branches, Theta=",t,"x")) #"#E69F00", 
# }
# figure=ggarrange(plot_list_tar[["3"]],plot_list_nonTar[["3"]],
#                  plot_list_tar[["6"]],plot_list_nonTar[["6"]],
#                  plot_list_tar[["10"]],plot_list_nonTar[["10"]],common.legend=T,legend="bottom",nrow=3,ncol=2)
# 
# annotate_figure(figure) #,top=text_grob(paste0("Posterior Probability of Acceleration, Case ",c+1))
# print(figure)
# #dev.off()












setwd("C:/bin/phyloacc/Yan-etal-2022/results/simulations/Theta_Scaling")
#get GT accuracy:
prefix="C:/bin/phyloacc/Yan-etal-2022/results/simulations/Theta_Scaling/Theta"

#setup
caseLabel=c(8,5,6)
nameLabel=c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","cryCin","tinGut","eudEle", "notPer","anoDid","strCam",
            "aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme","cryCin.tinGut",
            "eudEle.notPer","cryCin.eudEle","cryCin.anoDid", "aptHaa.cryCin","aptHaa.strCam")
birdname=read.table('C:/bin/phyloacc/Yan-etal-2022/data/SimulatedData/Name_Label_Mapping.txt',sep='\t',stringsAsFactors = F)

namePrint=nameLabel
for(i in 1:13){
  namePrint[i]=birdname[which(birdname[,1]==nameLabel[i]),2]
}
for(i in 14:length(nameLabel)){
  t1=birdname[which(birdname[,1]==substr(nameLabel[i],1,6)),2][1]
  t2=birdname[which(birdname[,1]==substr(nameLabel[i],8,8+6)),2]
  namePrint[i]=paste0("(",t1,",",t2,")")
}

target=nontarget=list() #(1-1,2-8,2-5,2-6)
target[[1]]=c()
nontarget[[1]]=namePrint
target[[2]]=c("aptHaa","aptOwe","aptRow","rheAme","rhePen","aptHaa.aptOwe", "aptHaa.aptRow","rheAme.rhePen")
target[[3]]=c("aptHaa","aptOwe", "aptRow","casCas","droNov","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas")
target[[4]]=c("aptHaa","aptOwe", "aptRow","casCas","droNov","rheAme","rhePen","anoDid","strCam","aptHaa.aptOwe", "aptHaa.aptRow","casCas.droNov","aptHaa.casCas", "rheAme.rhePen","aptHaa.rheAme")

#########BEAST
case=c("1-1","2-8","2-5","2-6")
target_B=list()
nontarget_B=list()
target_B[[1]]=nontarget_B[[1]]=c()
target_B[[2]]=c(10:14,19:21)
target_B[[3]]=c(10:18)
target_B[[4]]=c(1,10:22,24)
for(i in 2:4){
  nontarget_B[[i]]=seq(1,31)
  nontarget_B[[i]]=nontarget_B[[i]][-target_B[[i]]]
}

spt="(((halLeu:0.05373562#0,nipNip:0.05176648#0)taeGut-aptFor:0.04155687#0.04578784,galGal:0.1660553#0)taeGut-galGal:0.0406841#0.04578784,((((((aptHaa:0.00138725#0,aptOwe:0.0016341#0)aptHaa-aptOwe:0.00305054#0.006321079,aptRow:0.00410494#0)aptHaa-aptRow:0.0277496#0.019876606,(casCas:0.0115461#0,droNov:0.0137332#0)casCas-droNov:0.0273773#0.029694415)aptHaa-casCas:0.0028599#0.086448824,(rheAme:0.00469588#0,rhePen:0.00533574#0)rheAme-rhePen:0.0566382#0.027801219)aptHaa-rheAme:0.00185668#0.058734321,(((cryCin:0.0470926#0,tinGut:0.0388556#0)cryCin-tinGut:0.0172068#0.037018559,(eudEle:0.0655012#0,notPer:0.073059#0)eudEle-notPer:0.0079941#0.057911894)cryCin-eudEle:0.0672093#0.095580579,anoDid:0.0492722#0)cryCin-anoDid:0.0253604#0.053495448)aptHaa-cryCin:0.0118742#0.064943475,strCam:0.051388#0)aptHaa-strCam:0.0406969#0.011627686)taeGut-aptHaa:0#0.04578784;"
sptree=ape::read.tree(text=spt)
#tipname_nexus2=c('anoDid','aptHaa', 'aptOwe', 'aptRow', 'casCas',   'cryCin','droNov','eudEle', 'galGal','halLeu','nipNip','notPer','rheAme','rhePen','strCam','tinGut')
tipname_nexus= c('Moa',    'Kiwi1', 'Kiwi2', 'Kiwi3',  'Cassowary','Tinamou1','Emu', 'Tinamou3','Chicken','Eagle', 'Ibis','Tinamou4','Rhea1','Rhea2','Ostrich','Tinamou2')
name_nexus=c(tipname_nexus[1],tipname_nexus[c(6,16)],'(Tinamou1,Tinamou2)',tipname_nexus[c(8,12)],'(Tinamou3,Tinamou4)','(Tinamou1,Tinamou3)','(Tinamou1,Moa)',
             tipname_nexus[c(2,3)],'(Kiwi1,Kiwi2)',tipname_nexus[4],'(Kiwi1,Kiwi3)',tipname_nexus[c(5,7)],'(Cassowary,Emu)','(Kiwi1,Cassowary)',
             tipname_nexus[c(13,14)],'(Rhea1,Rhea2)','(Kiwi1,Rhea1)','(Kiwi1,Tinamou1)',tipname_nexus[15],'(Kiwi1,Ostrich)',tipname_nexus[9],tipname_nexus[c(10,11)],
             'taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
sptree_name_all=c(sptree$tip.label,sptree$node.label)


case = c("2-8" ,"2-5" ,"2-6")
case2 = c('A','B','C')

c=2 #change to 1:3
plot_list_tar=c() 
plot_list_nonTar=c()
for(t in c(3,6,10)){
  noGT_M1=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_rate_postZ_M1.txt"), header=T,row.names=1) ########
  noGT_M2=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_rate_postZ_M2.txt"), header=T,row.names=1)
  ##order rows of P(Z=2|Y) according to element ID
  noGT_M1= noGT_M1[order(as.numeric(rownames(noGT_M1))),] 
  noGT_M2= noGT_M2[order(as.numeric(rownames(noGT_M2))),]
  #store P(Z=2|Y) for all elements. Initialize to P(Z=2|Y,M1) for all elements
  noGT=noGT_M1 
  #read in elem_lik.txt file to determine optimal model for each elements.
  elem_Acc=as.matrix(read.table(paste0(prefix,t,"/Result_Case",case2[c],"/PhyloAcc/simu_200_100_",case[c],"_elem_lik.txt"),header=T))
  optM= apply(elem_Acc[,3:5],1,function(x) which(x==max(x))[1]-1) ##########
  elem_Acc=cbind(elem_Acc, optM)
  #if M0 is the best model, then all P(Z=2|Y) = 0 for this element
  ind=which(elem_Acc[,'optM']==0)
  if(length(ind)>0) noGT[ind,]=0
  #if M2 is the best model, get P(Z=2|Y, M2) for this element
  ind=which(elem_Acc[,'optM']==2)
  if(length(ind)>0) noGT[ind,]=noGT_M2[ind,]
  
  GT_M1=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_rate_postZ_M1.txt"), header=T, row.names=1)
  GT_M1 = GT_M1[order(as.numeric(rownames(GT_M1))),]
  GT_M2=read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_rate_postZ_M2.txt"), header=T, row.names=1)
  GT_M2 = GT_M2[order(as.numeric(rownames(GT_M1))),]
  GT=GT_M1
  elem_GT=as.matrix(read.table(paste0(prefix,t,"/Result_Case",case2[c],"/simu_100_",case[c],"_elem_lik.txt"),header=T))
  optM=apply(elem_GT[,3:5],1,function(x) which(x==max(x))[1]-1) ###########
  elem_GT =cbind(elem_GT,optM) #######
  ind=which(elem_GT[,'optM']==0)  ######
  if(length(ind)>0) GT[ind,]=0
  ind=which(elem_GT[,'optM']==2)
  if(length(ind)>0) GT[ind,]=GT_M2[ind,]
  
  #remove first 17 columns: not involving rates in non-outgroup branches.
  noGT=noGT[,-seq(1,17)] 
  GT=GT[,-seq(1,17)]
  # speciesLabel_3 is the estimated P(Z=2|Y) per species/branch.
  noGT=noGT[,grepl("_3", colnames(noGT))]
  GT=GT[,grepl("_3", colnames(GT))]
  colnames(noGT)=colnames(GT)=substr(colnames(GT),1,str_length(colnames(GT))-2) #checked noGT, GT share col names
  #remove br leading to out
  outRelated=c('taeGut.aptFor','taeGut.galGal','taeGut.aptHaa')
  for(i in 1:length(outRelated)){
    ind=which(colnames(GT)==outRelated[i])
    GT=GT[,-ind]
    noGT=noGT[,-ind]
  }
  
  #renaming branches to paper's naming convention
  targ_ind=rep(NA,length(target[[c+1]]))
  for(i in 1:length(targ_ind)){
    targ_ind[i]=which(nameLabel==target[[c+1]][i])
  }
  colnames(GT)=colnames(noGT)=namePrint
  
  #BEAST
  load(paste0(prefix,t,'/Result_Case',case2[c],'/BEAST/result1_',case[c],'.RData')) #Result8/
  res2=cbind(res$node,res$`non-node`)
  #edit names to be consistent with GT names
  for(i in 1:ncol(res2)){
    tmp=colnames(res2)[i]
    if(str_length(tmp)>6){
      colnames(res2)[i]=paste0("(",namePrint[which(nameLabel==substr(tmp,1,6))],",",namePrint[which(nameLabel==substr(tmp,8,8+6))],")")
    }else{
      colnames(res2)[i]=namePrint[which(nameLabel==tmp)[1]]
    }
  }
  res_backup=res2[,colnames(GT)[targ_ind]] #store rate estimates for target lineages
  res_backup2=res2[,colnames(GT)[-targ_ind]] #store rate estimates for non-target lineages
  rownames(res_backup)=rownames(res_backup2)=as.character(seq(0,99)) #add row index to be consistent with PhyloAcc
  
  #######draw boxplot for target lineages
  GTZ=GT[,targ_ind]
  noGTZ=noGT[,targ_ind]
  #convert to dataframe for boxplots
  dat = cbind(reshape2::melt(GTZ), "Method" = rep("PhyloAcc-GT", ncol(GTZ)*nrow(GTZ))) #2-8: 8  #2-5: 8
  dat = rbind(dat, cbind(reshape2::melt(noGTZ), "Method" = rep("PhyloAcc", ncol(noGTZ)*nrow(noGTZ))))
  tmp=cbind(reshape2::melt(res_backup),"Method"=rep("*BEAST2", ncol(res_backup)*nrow(res_backup)))[,-1]
  colnames(tmp)=colnames(dat)
  dat = rbind(dat,tmp)
  colnames(dat)[1:2] = c("species", "acc")
  dat$Method = factor(dat$Method, levels=c("PhyloAcc-GT", "PhyloAcc", "*BEAST2"))
  plot_list_tar[[as.character(t)]]=ggplot(dat, aes(x = species, y = acc, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 12) +xlab("") + ylab("P(acceleration)") +theme(axis.text.x = element_text(angle = 30,hjust = 0.8)) + 
    scale_fill_manual(values=c("#CC79A7","#56B4E9","#F0E442"))+ggtitle(paste0("accelerated branches, Theta=",t,"x")) #"#E69F00", 
  
  #draw boxplot for non-target lineages
  GTZ=GT[,-targ_ind]
  noGTZ=noGT[,-targ_ind]
  dat = cbind(reshape2::melt(GTZ), "Method" = rep("PhyloAcc-GT", ncol(GTZ)*nrow(GTZ))) #2-8: 8  #2-5: 8
  dat = rbind(dat, cbind(reshape2::melt(noGTZ), "Method" = rep("PhyloAcc", ncol(noGTZ)*nrow(noGTZ))))
  tmp=cbind(reshape2::melt(res_backup2),"Method"=rep("*BEAST2", ncol(res_backup2)*nrow(res_backup2)))[,-1]
  colnames(tmp)=colnames(dat)
  dat = rbind(dat,tmp)
  colnames(dat)[1:2] = c("species", "acc")
  dat$Method = factor(dat$Method, levels=c("PhyloAcc-GT", "PhyloAcc", "*BEAST2"))
  plot_list_nonTar[[as.character(t)]]=ggplot(dat, aes(x = species, y = acc, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 12) +xlab("")+ ylab("P(acceleration)") +theme(axis.text.x = element_text(angle = 30,hjust = 0.8)) + 
    scale_fill_manual(values=c("#CC79A7","#56B4E9","#F0E442"))+ggtitle(paste0("non-accelerated branches, Theta=",t,"x")) #"#E69F00", 
}
figure=ggarrange(plot_list_tar[["3"]],plot_list_nonTar[["3"]]+rremove('ylab'),
                 plot_list_tar[["6"]],plot_list_nonTar[["6"]]+rremove('ylab'),
                 plot_list_tar[["10"]],plot_list_nonTar[["10"]]+rremove('ylab'),common.legend=T,legend="bottom",nrow=3,ncol=2)

annotate_figure(figure) #,top=text_grob(paste0("Posterior Probability of Acceleration, Case ",c+1))

print(figure)
#dev.off()