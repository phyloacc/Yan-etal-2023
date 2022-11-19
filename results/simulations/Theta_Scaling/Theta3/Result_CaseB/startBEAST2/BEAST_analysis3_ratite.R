library(ape)
library(seqinr)
library(stringi)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

accPat="2-5"
#accPat=args[1]
cat(args[1],"\n")
cat(accPat," ",typeof(accPat),"\n")


getRateCompMatrix<-function (edgeRel,RateVec){
  matRes=matrix(0,nrow=31,ncol=4) #ch-rate, pa-rate, ch point acc, ch status
  for(i in 1:30){
    matRes[i,1:2]=RateVec[edgeRel[i,]]
    if(matRes[i,1]>matRes[i,2]){
      matRes[i,3]=1
    }else if(matRes[i,1]<matRes[i,2]){
      matRes[i,3]=-1
    }
  }
  matRes[31,]=c(1,1,0,0)
  
  for(i in seq(30,1)){
    if(matRes[i,3]==-1){
      next
    }
    if(matRes[i,3]==1 | matRes[edgeRel[i,2],4]==1){ #if child has acceleration, or if pa acc and child remain no chg, then ch status is ACC
      matRes[i,4]=1
    }
  }
  row.names(matRes)=row.names(edgeRel)
  return(matRes)
}

prefix1="species-simu_200_100_"
prefix2=".trees"

case=c("1-1","2-8","2-5","2-6")
target=list()
nontarget=list()
target[[1]]=nontarget[[1]]=c()
target[[2]]=c(10:14,19:21)
target[[3]]=c(10:18)
target[[4]]=c(1,10:22,24)
for(i in 2:4){
  nontarget[[i]]=seq(1,31)
  nontarget[[i]]=nontarget[[i]][-target[[i]]]
}
case1=which(case==accPat)
target1=target[[case1]]
Lt1=length(target1)
nontarget1=nontarget[[case1]]

spt="(((halLeu:0.05373562#0,nipNip:0.05176648#0)taeGut-aptFor:0.04155687#0.04578784,galGal:0.1660553#0)taeGut-galGal:0.0406841#0.04578784,((((((aptHaa:0.00138725#0,aptOwe:0.0016341#0)aptHaa-aptOwe:0.00305054#0.006321079,aptRow:0.00410494#0)aptHaa-aptRow:0.0277496#0.019876606,(casCas:0.0115461#0,droNov:0.0137332#0)casCas-droNov:0.0273773#0.029694415)aptHaa-casCas:0.0028599#0.086448824,(rheAme:0.00469588#0,rhePen:0.00533574#0)rheAme-rhePen:0.0566382#0.027801219)aptHaa-rheAme:0.00185668#0.058734321,(((cryCin:0.0470926#0,tinGut:0.0388556#0)cryCin-tinGut:0.0172068#0.037018559,(eudEle:0.0655012#0,notPer:0.073059#0)eudEle-notPer:0.0079941#0.057911894)cryCin-eudEle:0.0672093#0.095580579,anoDid:0.0492722#0)cryCin-anoDid:0.0253604#0.053495448)aptHaa-cryCin:0.0118742#0.064943475,strCam:0.051388#0)aptHaa-strCam:0.0406969#0.011627686)taeGut-aptHaa:0#0.04578784;"
sptree=ape::read.tree(text=spt)
tipname_nexus=c('anoDid', 'aptHaa', 'aptOwe', 'aptRow', 'casCas','cryCin','droNov','eudEle','galGal','halLeu','nipNip','notPer','rheAme','rhePen','strCam','tinGut')
name_nexus=c(tipname_nexus[1],tipname_nexus[c(6,16)],'cryCin-tinGut',tipname_nexus[c(8,12)],'eudEle-notPer','cryCin-eudEle','cryCin-anoDid',
             tipname_nexus[c(2,3)],'aptHaa-aptOwe',tipname_nexus[4],'aptHaa-aptRow',tipname_nexus[c(5,7)],'casCas-droNov','aptHaa-casCas',
             tipname_nexus[c(13,14)],'rheAme-rhePen','aptHaa-rheAme','aptHaa-cryCin',tipname_nexus[15],'aptHaa-strCam',tipname_nexus[9],tipname_nexus[c(10,11)],
             'taeGut-aptFor','taeGut-galGal','taeGut-aptHaa')
sptree_name_all=c(sptree$tip.label,sptree$node.label)

edgeRelation=matrix(NA,nrow=31,ncol=6) 
#row: edge. #col: child species ID (correpond to tipname_nexus); pa sp ID; child rate; pa rate; rate change: 0=no, 1=acc, -1=conserved; status
sptree_name_all=c(sptree$tip.label,sptree$node.label)
for(i in 1:30){
  edgeRelation[i,1]=i
  temp_name=name_nexus[i]
  sp_ID=which(sptree_name_all==temp_name)
  pa_sp_ID=sptree$edge[which(sptree$edge[,2]==sp_ID),1]
  pa_name=sptree_name_all[pa_sp_ID]
  pa_nexus_ID=which(name_nexus==pa_name)
  edgeRelation[i,2]=pa_nexus_ID
}
edgeRelation[31,1:2]=c(31,-1)
row.names(edgeRelation)=name_nexus

tree_template="((((1[&rate=0.9672256053231673]:0.049272200000000016,((6[&rate=0.9672256053231673]:0.04709260000000001,16[&rate=0.9672256053231673]:0.03885559999999999)[&rate=0.9672256053231673]:0.017206799999999994,(8[&rate=0.9672256053231673]:0.06550119999999998,12[&rate=0.9672256053231673]:0.07305899999999999)[&rate=0.9672256053231673]:0.007994100000000004)[&rate=0.9672256053231673]:0.06720930000000003)[&rate=0.9672256053231673]:0.025360399999999977,((((2[&rate=1.3301527822820187]:0.0013872500000000065,3[&rate=1.3301527822820187]:0.0016341000000000272)[&rate=1.3301527822820187]:0.0030505399999999905,4[&rate=1.3301527822820187]:0.004104940000000001)[&rate=1.3301527822820187]:0.027749599999999985,(5[&rate=1.3301527822820187]:0.011546100000000004,7[&rate=1.3301527822820187]:0.013733200000000001)[&rate=1.3301527822820187]:0.027377299999999993)[&rate=1.3301527822820187]:0.0028598999999999986,(13[&rate=0.9672256053231673]:0.0046958799999999995,14[&rate=0.9672256053231673]:0.0053357400000000055)[&rate=0.9672256053231673]:0.0566382)[&rate=0.9672256053231673]:0.0018566799999999994)[&rate=0.9672256053231673]:0.011874200000000001,15[&rate=0.9672256053231673]:0.05138799999999999)[&rate=0.9672256053231673]:0.04069690000000001,(9[&rate=0.9650521888254431]:0.1660553,(10[&rate=0.9650521888254431]:0.05373562000000001,11[&rate=0.9650521888254431]:0.051766480000000004)[&rate=0.9650521888254431]:0.041556869999999996)[&rate=0.9650521888254431]:0.0406841)[&rate=1.0]:0.0;"
template_tree=ape::read.tree(text=tree_template)

tree_rate=rep(NA,length(name_nexus))
names(tree_rate)=name_nexus

M=100
modelSel_all=matrix(NA,nrow=M,ncol=3) #model selection prob
trueSel_all=matrix(NA,nrow=M,ncol=Lt1) #target species selection prob
estRate_all=matrix(NA,nrow=M,ncol=2) #estiamted conserved and acc rate (mean. since can have diff acc(con) rates)
misMatch_tree=rep(0,M)
modelAcc=rep(0,M)

modelSelpt_all=rep(NA,M) #correpond to max Post
Ratept_all=matrix(NA,nrow=M,ncol=2)
trueSelpt_all=matrix(NA,nrow=M,ncol=Lt1)
for(m in 1:M){
  logfile=read.table(paste0("starbeast-simu_200_100_",accPat,"_",m-1,".log"),header=TRUE)
  mcmcL=length(logfile[,"posterior"])
  L=ceiling((mcmcL-1)*0.25)
  logfile=logfile[-(1:(L+1)),]
  indMax=which(logfile[,"posterior"]==max(logfile[,"posterior"]))
  indMax=indMax[length(indMax)]
  
  fileName=paste0(prefix1,accPat,"_",m-1,prefix2) 
  stringT=scan(fileName, what = "", sep = "\n", quiet = TRUE) #read in BEAST mcmc trees
  strL=length(stringT) 
  stringT=stringT[42:(strL-1)]
  # strL=length(stringT)
  # L=ceiling(strL*0.25)
  # stringT=stringT[L:strL]
  stringT=stringT[-(1:(L+1))]
  strL=length(stringT)
  ex=0
  
  modelSel=rep(0,strL)
  Ptarget1=rep(0,Lt1) #sampling freq of acc of target group species.
  estRate=matrix(NA,nrow=strL,ncol=2) #acc rate, con rate per t. later will get average and pass to estRate_all[m,]
  for(t in 1:strL){
    tree=ape::read.tree(text=stringT[t])
    if(sum(abs(template_tree$edge-tree$edge))!=0){
      ex=ex+1
      next
    }
    
    temp1=str_match_all(stringT[t],"rate=\\s*(.*?)\\s*\\]")[[1]][,2]
    L1=length(temp1)
    tree_rate=rep(0,L1)
    for(i in 1:L1){
      tree_rate[i]=as.numeric(str_split(temp1[i],",")[[1]][1])
    }
    #tree_rate=sapply(temp,as.numeric,USE.NAMES = FALSE)
    #tree_rate[31]=max(tree_rate[c(30,31)]) #outgrp cannot get acc, so take the ref rate to be the larger of 1 vs MCRA of outgrp.
    #But this will distort the actual rates.
    # indM_tmp=which(tree_rate[26:31]==max(tree_rate[26:31]))
    # indM_tmp=indM_tmp[length(indM_tmp)]
    # if(indM_tmp!=6){ #outgrp got accelerated at some point
    #   tmp_scale=tree_rate[25+indM_tmp]/tree_rate[31]
    #   tree_rate=tree_rate*tmp_scale
    # }
    ind_tmp=which(tree_rate[26:30]>tree_rate[31]) #manually adjust acc outgroup branch to non-acc
    if(length(ind_tmp)>0) tree_rate[25+ind_tmp]=tree_rate[31]
    edgeRelation[,3:6]=getRateCompMatrix(edgeRelation[,1:2],tree_rate)
    
    indA=which(edgeRelation[,6]==1)
    if(length(indA)>0){
      estRate[t,1]=mean(tree_rate[indA])
      if(t==indMax) Ratept_all[m,1]=estRate[t,1]
    }
    indC=which(edgeRelation[,6]==0)
    if(length(indC) !=0){
      ind=which(tree_rate[indC]!=1)
      if(length(ind)>0){
        estRate[t,2]=mean(tree_rate[indC[ind]])
        if(t==indMax) Ratept_all[m,2]=estRate[t,2]
      }
    }
    
    accBr=which(edgeRelation[,5]==1)
    accL=length(accBr)
    if(accL==0){ #no acceleration
      modelSel[t]=0
      if(t==indMax) modelSelpt_all[m]=0
    }else{
      if(length(intersect(accBr,nontarget[case1]))>0){
        modelSel[t]=2
        if(t==indMax) modelSelpt_all[m]=2
      }else{
        modelSel[t]=1
        if(t==indMax) modelSelpt_all[m]=1
      }
    }
    if(case1 !=1 ){
      Ptarget1=Ptarget1+edgeRelation[target1,6]
      if(t==indMax) trueSelpt_all[m,]=edgeRelation[target1,6]
    }
  }
  modelSel_all[m,]=c(sum(modelSel==0),sum(modelSel==1),sum(modelSel==2))/strL
  if(case1 !=1){
    if(which(modelSel_all[m,]==max(modelSel_all[m,]))==2) modelAcc[m]=1
    trueSel_all[m,]=Ptarget1/strL
  }else{
    if(which(modelSel_all[m,]==max(modelSel_all[m,]))==1) modelAcc[m]=1
  }
  #if(case1 != 1) trueSel_all[m,]=Ptarget1/strL
  estRate_all[m,]=apply(estRate,2,function(x) mean(x,na.rm = TRUE))
  misMatch_tree[m]=ex
}


res=list("model"=modelSel_all,"node"=trueSel_all,"rate"=estRate_all,"missedTree"=misMatch_tree, "modelAcc"=modelAcc, "AccP"=mean(modelAcc),
         "modelpt"=modelSelpt_all,"nodept"=trueSelpt_all,"ratept"=Ratept_all)
save(res,file=paste0("result1_",accPat,".RData"))

#modelSel_all Mx3: for each element (m): outoff 1500MCMC samples (t), how many times/proportion M0,M1,M2 selected. 
#trueSel_all Mx|target|:  for target branches: out off MCMC samples, how frequently its accelerated. i.e. posterior mean of Z=2.
#estRate_all Mx2: estimated accR, conR = mean(estRate): estRate: in each MCMC sample: get average of all accR (random local clock). Then average over all MCMC samples
#missedTree: if any tree topology diff from template
#modelAcc M: if M1 is preferred (proportion largest in modelSel_all)
#AccP: on average how many times M1 is correctly selected
#modelSelpt_all; trueSelpt_all; Ratept_all: point estimate correspond to max Posterior
