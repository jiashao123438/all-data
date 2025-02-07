

if (T) {
  dir.create("data")
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("PDFs")
  dir.create("01_origin_datas/GEO",recursive = T)
  dir.create("01_origin_datas/TCGA")
}
library(colorspace)
library(ggpubr)
library(ggsci)
library(reshape2)
library(ggpubr)
library(ggsci)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(corrplot)
library(survminer)
library(survival)
options(stringsAsFactors = F)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}
mg_violin <- function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL){
  library(ggplot2)
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  if(!is.null(ylim)){
    data_m=data_m[data_m[,2]<=max(ylim),]
    ylim[2]=1.2*ylim[2]
  }
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  if(ct<=10){
    p1=p1+ggsci::scale_fill_npg(name=leg.title)
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10)
                ,ggsci::pal_d3("category20", alpha = 0.6)(20)
                ,ggsci::pal_nejm("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, #
              axis.text.y=element_text(family="Times",face="plain"), #
              axis.title.y=element_text(family="Times",face="plain"), #
              #panel.border = element_blank(),axis.line = element_line(colour = "black"), #
              legend.text=element_text(face="plain", family="Times", colour="black"  #
              ),
              legend.title=element_text(face="plain", family="Times", colour="black" #
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
              #,panel.grid.major = element_blank(),   #
              #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
sig_boxplot_t<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="t.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
sig_boxplot_w<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}

my_boxplot=function(dat,group,group_cols,test_method='kruskal.test',fill= "Group",
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=SiaTs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+
    #scale_color_manual(values = sub.col) +  #
    scale_fill_manual(values = group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.format", method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}


my_mutiboxplot=function(dat,group,group_cols,test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],bw=T,
                    xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot()+
    #scale_color_manual(values = sub.col) +  #
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}

my_mutiboxplot_seg=function(dat,group,group_cols,test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],fill='group',
                        xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,legend.position='top',nrow,ncol){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+facet_wrap(~type,scales = 'free',nrow = nrow,ncol = ncol)+
    #scale_color_manual(values = sub.col) +  #
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill = fill,title =title) +
    #theme_light()+
    #theme_classic()+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}


mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  # if(jitter){
  #   if(is.null(point_size)){
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
  #   }else{
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
  #   }
  # }
  # 
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE 
  #              ,showP = T 
  #              ,droplines = F#
  #,colors = mg_colors[1:3] #
  #,rank="decreasing") #
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

#####data.pre#########
genecode=read.delim
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
#lncrna_genecode=genecode[which(genecode$TYPE=='lncRNA'),]
#####TCGA#######
#########
tcga_cli=read.delim('00_origin_datas/TCGA/Merge_CESC_clinical.txt',sep = '\t',header=T)
colnames(tcga_cli)
table(tcga_cli$A2_Event)
tcga_cli=data.frame(Samples=tcga_cli$A0_Samples,
                    Age=tcga_cli$A17_Age,
                    #Gender=tcga_cli$A18_Sex,
                    T.stage=tcga_cli$A3_T,
                    N.stage=tcga_cli$A4_N,
                    M.stage=tcga_cli$A5_M,
                    Stage=tcga_cli$A6_Stage,
                    Grade=tcga_cli$A7_Grade,
                    Status=tcga_cli$A2_Event,
                    OS.time=tcga_cli$A1_OS)
tcga_cli$Samples=paste0(tcga_cli$Samples,'-01')
rownames(tcga_cli)=tcga_cli$Samples
tcga_cli$OS.time
tcga_cli=tcga_cli %>% drop_na(OS.time)
tcga_cli=tcga_cli[which(tcga_cli$OS.time>0),]
table(tcga_cli$Status)
tcga_cli$OS=ifelse(tcga_cli$Status=='Alive',0,1)

fivenum(tcga_cli$Age)
tcga_cli$Age1=ifelse(tcga_cli$Age>46,'>46','<=46')

table(tcga_cli$T.stage)
tcga_cli$T.stage=substr(tcga_cli$T.stage,1,2)
tcga_cli$T.stage[which(tcga_cli$T.stage=='Ti'|tcga_cli$T.stage=='TX'|tcga_cli$T.stage=='')]=NA

table(tcga_cli$N.stage)
tcga_cli$N.stage[which(tcga_cli$N.stage==''|tcga_cli$N.stage=='NX')]=NA

table(tcga_cli$M.stage)
tcga_cli$M.stage[which(tcga_cli$M.stage==''|tcga_cli$M.stage=='MX')]=NA

table(tcga_cli$Stage)
tcga_cli$Stage=gsub('[ABC12]','',tcga_cli$Stage)
tcga_cli$Stage=gsub('Stage ','',tcga_cli$Stage)
tcga_cli$Stage[which(tcga_cli$Stage=='')]=NA

table(tcga_cli$Grade)
tcga_cli$Grade[which(tcga_cli$Grade=='GX'|tcga_cli$Grade=='Not Available')]=NA

head(tcga_cli)

#########
tcga_data=read.delim('00_origin_datas/TCGA/CESC_TPM.txt',sep = '\t',row.names = 1,check.names = F)
tcga_data[1:5,1:5]

range(tcga_data)
table(substr(colnames(tcga_data),14,15))

tcga_tpm_log_T=log2(tcga_data[,intersect(tcga_cli$Samples,colnames(tcga_data))]+1)
range(tcga_tpm_log_T)
dim(tcga_tpm_log_T)
#58938   291
tcga_cli=tcga_cli[intersect(tcga_cli$Samples,colnames(tcga_tpm_log_T)),]
dim(tcga_cli)
#291  11

sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#
length(sample_N)
length(sample_T)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 3    304 

tcga_mrna=tcga_tpm_log_T[intersect(mrna_genecode$SYMBOL,rownames(tcga_tpm_log_T)),]
dim(tcga_mrna)
#19503   291
tcga_exp=tcga_mrna[intersect(rownames(GSE44001_exp),rownames(tcga_mrna)),]
dim(tcga_exp)
#16033   291
############GSE44001#########
GSE44001=getGEOExpData('GSE44001')
save(GSE44001,file='00_origin_datas/GEO/GSE44001.RData')
load('00_origin_datas/GEO/GSE44001.RData')
GSE44001_cli=GSE44001$Sample
head(GSE44001_cli)
table(GSE44001_cli$status_of_dfs)
GSE44001_cli=data.frame(Samples=GSE44001_cli$Acc,
                        OS.time=GSE44001_cli$`disease_free_survival_(dfs)_(months)`*30,
                        OS=GSE44001_cli$status_of_dfs)
rownames(GSE44001_cli)=GSE44001_cli$Samples
GSE44001_cli$OS.time
table(GSE44001_cli$OS)

GSE44001_dat=GSE44001$Exp$GPL14951_29377_Data_col1
GSE44001_dat[1:5,1:5]
GSE44001_exp=exp_probe2symbol_v2(datExpr = GSE44001_dat,anno = GSE44001$Anno$GPL14951[c(1,12)],method = 'mean')
save(GSE44001_exp,file = '00_origin_datas/GEO/GSE44001_exp.RData')
load('00_origin_datas/GEO/GSE44001_exp.RData')
range(GSE44001_exp)
boxplot(GSE44001_exp[,1:5])
GSE44001_exp=GSE44001_exp[,intersect(GSE44001_cli$Samples,colnames(GSE44001_exp))]
GSE44001_cli=GSE44001_cli[intersect(GSE44001_cli$Samples,colnames(GSE44001_exp)),]
dim(GSE44001_exp)
dim(GSE44001_cli)

GSE44001_exp <- GSE44001_exp[intersect(rownames(GSE44001_exp),rownames(tcga_mrna)),]
##01.SiaTs#############
dir.create('01_SiaTs_gene')
SiaTs.gc=read.table("01_SiaTs_gene/SiaTs_PMID_37968767.txt",sep = "\t",quote ="",header = F)
SiaTs.gc=SiaTs.gc$V1
length(SiaTs.gc)
length(intersect(SiaTs.gc,rownames(tcga_mrna)))
setdiff(SiaTs.gc,rownames(tcga_mrna))

##1.1计算SiaTs score####
SiaTs=SiaTs.gc
SiaTs.score=t(ssGSEAScore_by_genes(gene.exp = tcga_mrna,genes = SiaTs))
#SiaTs.score=data.frame(Samples=rownames(SiaTs.score),score=as.numeric(SiaTs.score[,1]))
colnames(SiaTs.score)='SiaTs.score'
writeMatrix(SiaTs.score,'01_SiaTs_gene/SiaTs.score.txt',row = T)
#####1.2SiaTs ############
tcga.maf=getTCGAMAFByCode('CESC')
tcga.type.use=tcga_type
table(tcga.type.use$type)
colnames(tcga.type.use)[1]='Tumor_Sample_Barcode'
tcga.type.use$Tumor_Sample_Barcode=substr(tcga.type.use$Tumor_Sample_Barcode,1,12)
tcga.type.use.Normal=tcga.type.use[which(tcga.type.use$type=='Normal'),]
tcga.type.use.Tumor=tcga.type.use[which(tcga.type.use$type=='Tumor'),]
write.table(tcga.type.use.Normal,file='01_SiaTs_gene/tcga.type.Normal.txt')
write.table(tcga.type.use.Tumor ,file='01_SiaTs_gene/tcga.type.Tumor.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.type.use.Normal$Tumor_Sample_Barcode))
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = '01_SiaTs_gene/tcga.type.Normal.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.type.use.Tumor $Tumor_Sample_Barcode))
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = '01_SiaTs_gene/tcga.type.Tumor.txt')
tcga.maf2@clinical.data


type.color=c('#071952','#37B7C3')
# pdf('01_SiaTs_gene/normal_mut.pdf',height = 8,width =6)
# oncoplot(maf=tcga.maf1,
#          #clinicalFeatures = 'type',
#          genes = SiaTs.gc,
#          sortByAnnotation = T,
#          #annotationColor = list(type=c(Normal=as.character(type.color[1])))
# )
# dev.off()
pdf('01_SiaTs_gene/tumor_mut.pdf',height = 8,width =6)
oncoplot(maf=tcga.maf2,
         #clinicalFeatures = 'type',
         genes = SiaTs.gc,
         sortByAnnotation = T,
         #annotationColor = list(type=c(Tumor =as.character(type.color[2])))
)
dev.off()
###1.3####
SiaTs.score.cli=data.frame(SiaTs.score=SiaTs.score[tcga_cli$Samples,],tcga_cli[,c('OS','OS.time')])
SiaTs.score.cli$group=ifelse(SiaTs.score.cli$SiaTs.score>median(SiaTs.score.cli$SiaTs.score),'High','Low')

SiaTs.score.roc=ggplotTimeROC(SiaTs.score.cli$OS.time,
                              SiaTs.score.cli$OS,
                              SiaTs.score.cli$SiaTs.score,mks = c(1,3,5))
SiaTs.score.roc
SiaTs.score.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                         data = SiaTs.score.cli),
                             data=SiaTs.score.cli,
                             conf.int = T,pval = T,risk.table = T,
                             fun = "pct",size = 1,surv.median.line = 'hv',
                             title='TCGA-CESC',legend.title='SiaTs score',
                             legend.labs = c('High','Low'),
                             linetype = c("solid", "dashed","strata")[1],
                             #palette = risktype.col,
                             ylab='Overall Survival(OS)',
                             legend=c(0.85,0.8),#
                             ggtheme = theme_bw(base_size = 12))
SiaTs.score.km.OS=mg_merge_plot(SiaTs.score.km.OS$plot,SiaTs.score.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

SiaTs.score.km.OS
########1.4#####
tcga.mcp=immu_MCPcounter(tcga_exp)
tcga.cli.immu2=cbind.data.frame(`SiaTs score`=SiaTs.score.cli[tcga_cli$Samples,'SiaTs.score'],tcga.mcp[tcga_cli$Samples,])
colnames(tcga.cli.immu2)
cor_res <- Hmisc::rcorr(as.matrix(tcga.cli.immu2),type = 'pearson')
cor_res$P[is.na(cor_res$P)] <- 0
pdf('results/01.Necr_genes/Fig1b.pdf',height = 6,width = 6)
corrplot(corr = cor_res$r,
         p.mat = cor_res$P,
         mar = c(0,0,1,1),
         col=colorRampPalette(c('blue', 'white','red'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()

tcga.est=immu_estimate(tcga_exp)
save(tcga.est,file = '01_SiaTs_gene/tcga.est.RData')
#load('results/tcga.est.RData')
tcga.cli.immu=cbind.data.frame(`SiaTs score`=tcga_cli_use[tcga_cli$Samples,'SiaTs.score'],tcga.est[tcga_cli$Samples,])

ggcorplot <- function(plot_df,a,b,method="spearman"){
  corr_eqn <- function(x,y,digits=3) {
    test <- cor.test(x,y,method=method,exact=FALSE)
    paste(#paste0("Method = ",method),
      paste0("r = ",round(test$estimate,digits)),
      paste0("p.value= ",round(test$p.value,digits)),
      sep = ", ")
  }
  plot_df <- plot_df[,c(a,b)]
  names(plot_df) <- c("var1","var2")
  require(ggplot2)
  p=ggplot(plot_df,aes(var1,var2))+
    geom_point(col="black")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="red")+
    geom_rug(col="#006fbc")+
    theme_minimal()+
    xlab(a)+
    ylab(b)+
    ## 
    labs(title = corr_eqn(plot_df$var1,plot_df$var2),subtitle =paste0("Method = ",method) )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  return(p)
}

fig1c=list()
for (i in 1:length(colnames(tcga.cli.immu))) {
  print(i)
  fig1c[[i]] = ggcorplot(plot_df = tcga.cli.immu,a = 'SiaTs score', b = colnames(tcga.cli.immu)[i], method = 'spearman')
}
length(fig1c)
fig1c=mg_merge_plot(fig1c[2:length(fig1c)],ncol = 3,nrow = 1)
fig1c
savePDF('01_SiaTs_gene/Fig1c.pdf',fig1c,height = 5,width = 12)



#02.WGCNA####
dir.create('02_WGCNA')
library(WGCNA)
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 

my_mad <- function(x){mad(x,na.rm = TRUE)} #
wgcna_exp=t(tcga_exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga_exp)
#  19503   291
#
tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]
# dim(tpm_T2)
#  60 13836
#tpm_T2=tcga_exp[which(apply(tcga_exp,1,sd)>0.5),]
#tpm_T2=(2^tpm_T2-1)
range(tpm_T2)

pdf('02_WGCNA/1.pdf',width = 10,height = 10)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

# minModuleSize = 30,    ##
# mergeCutHeight = 0.25, ##
tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.3,
                                 minModuleSize=60)


table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))
#17
write.csv(tpm_T2.module$Modules,file = "02_WGCNA/WGCNA_Modules.csv",row.names = T)
pdf('02_WGCNA/2.pdf',height = 5,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#
MODULE.gene.num <- as.data.frame(table(tpm_T2.module$Modules[,2]))
write.table(MODULE.gene.num,file = "02_WGCNA/MODULE.gene.num.txt",sep = "\t",quote = F,row.names =F)
pdf('02_WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02_WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()



SiaTs.score
##choose 
tcga_cli_use <-cbind.data.frame(tcga_cli[tcga_cli$Samples,],SiaTs.score=SiaTs.score[tcga_cli$Samples,])
head(tcga_cli_use)
tcga_cli_use=tcga_cli_use[,-1]
colnames(tcga_cli_use)
tcga_cli_use.part=data.frame(SiaTs.score=tcga_cli_use[tcga_cli$Samples,c(11)])
str(tcga_cli_use.part)
tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))


spms=tcga_cli_use.part
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02_WGCNA/5.pdf',width = 4,height =6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "black"
column = match(module, modNames)
column


moduleGenes<- (tpm_T2.module$Modules[,'mergedColors']==module)

tcga.wgcna.gene=c(names(which(moduleGenes)))
length(tcga.wgcna.gene)
# 1220
write.table(tcga.wgcna.gene,file = "02_WGCNA/tcga.wgcna.gene.csv",sep = "\t",quote = F,row.names = F)
pdf('02_WGCNA/S_model_GREEN.pdf',width = 6,height = 6)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
##3.####
dir.create("03_gsea")


##3.2####

cox.gene.enrichment=mg_clusterProfiler(genes =tcga.wgcna.gene )
fig3a=list()
fig3a[[1]]=enrichplot::dotplot(cox.gene.enrichment$KEGG)+ggtitle('KEGG')
fig3a[[2]]=enrichplot::dotplot(cox.gene.enrichment$GO_BP)+ggtitle('Biological Process')
fig3a[[3]]=enrichplot::dotplot(cox.gene.enrichment$GO_CC)+ggtitle('Cellular Component')
fig3a[[4]]=enrichplot::dotplot(cox.gene.enrichment$GO_MF)+ggtitle('Molecular Function')

fig3=mg_merge_plot(mg_merge_plot(fig3a[[1]],fig3a[[2]],labels = c('A','B'),ncol=2),
                   mg_merge_plot(fig3a[[3]],fig3a[[4]],ncol=2,nrow=1,labels = LETTERS[3:4],heights = c(1,1),widths = c(1,2)),
                   nrow=2,heights = c(1,1))
savePDF('03_gsea/Fig3.pdf',fig3,height = 12,width = 18)


#04.#######
dir.create('04_model')
##4.0####
tcga.SiaTs.score.gene <- data.frame(SiaTs_SCORE=SiaTs.score[rownames(tcga_cli),1],t(tcga_exp[tcga.wgcna.gene,rownames(tcga_cli)]))
outTab <- data.frame()
SiaTs_score <- "SiaTs_SCORE"                                  
j=tcga.SiaTs.score.gene[,SiaTs_score]
y=j
for (gene in colnames(tcga.SiaTs.score.gene)) {
  i=tcga.SiaTs.score.gene[,gene]
  x=i
  corT=cor.test(x,y)
  
  z=lm(y~x)
  cor=corT$estimate
  cor=round(cor,4)
  pvalue=corT$p.value
  pval=signif(pvalue,4)
  pval=format(pval, scientific = TRUE)
  
  outTab=rbind(outTab,
               cbind(SiaTs_score=SiaTs_score,gene=gene,cor=cor,
                     pvalue=pval))
}

corMatrix <- outTab[order(outTab$cor,decreasing = F),]
corMatrix$cor <- as.numeric(corMatrix$cor )
corMatrix$pvalue <- as.numeric(corMatrix$pvalue )
cor=0.3
cor.p=0.05
corMatrix <-corMatrix[which(abs(corMatrix$cor)>cor&corMatrix$pvalue<cor.p),]
dim(corMatrix)
write.table(corMatrix,file = "04_model/corTable-SiaTs_SCORE-gene.xls",quote = F,sep = "\t",row.names = F)
#######
dir.create('files/model_select')
select_gene_zscore <- function(dat1, dat2 = NULL, dat3 = NULL, 
                               a = 1, n = 100, 
                               ratio = 0.5, cut_p = 0.05, 
                               years = c(1,3,5)){
  library(timeROC)
  library(survival)
  library(glmnet)
  # library(mosaic)
  sample.index <- data.frame()
  SampleingTime <- c()
  SamplingSeed <- c()
  tra.auc <- c()
  test.auc <- c()
  geo1.auc <- c()
  geo2.auc <- c()
  all.auc <- c()
  tra.p <- c()
  test.p <- c()
  geo1.p <- c()
  geo2.p <- c()
  all.p <- c()
  tcga.dat <- dat1
  geo1.dat <- dat2
  geo2.dat <- dat3
  gene.list <- c()
  GeneNums <- c()
  for (i in a:n) {
    # i=2
    set.seed(i)
    par(mfrow = c(2, 3))
    myd.index <- seq(1,nrow(tcga.dat),1)
    tra.index <- sample(myd.index,size = round(nrow(tcga.dat)*ratio))
    tra.dat <- tcga.dat[tra.index,]
    test.dat <- tcga.dat[-tra.index,]
    write.table(rownames(tra.dat), file = paste0('files/model_select/tra.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    write.table(rownames(test.dat), file = paste0('files/model_select/test.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    tra.cox <- t(apply(tra.dat[,3:c(ncol(tra.dat))],2,function(x){
      
      vl=as.numeric(x)
      tm=tra.dat$OS.time
      ev=tra.dat$OS
      #ev=ifelse(ev=='Alive',0,1)
      dat=data.frame(tm,ev,vl)[which(tm > 0 & !is.na(vl)),]
      return(coxFun(dat))
    }))
    colnames(tra.cox) <- c('p.value','HR','Low 95%CI','High 95%CI')
    tra.cox <- as.data.frame(tra.cox)
    write.table(tra.cox, file = paste0('files/model_select/sig_cox_zscore_', i, '.txt'), sep = '\t', quote = F)
    cut.cox=tra.cox[which(tra.cox[,1]<cut_p),] 
    if (nrow(cut.cox) > 1) {
      print(paste("Processing....: ",i," resampling",sep=""))
      flush.console()
      geneid=rownames(cut.cox)
      geneid=geneid[which(!is.na(geneid))]
      if (length(geneid)>1) {
        sample.index<-c(sample.index,data.frame(tra.index))
        set.seed(i)
        cv.fit <- cv.glmnet(as.matrix(tra.dat[,geneid]), cbind(time=tra.dat$OS.time, 
                                                               status=tra.dat$OS)
                            ,family="cox", nlambda=100, alpha=1) 
        sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
        write.table(sig.coef, file = paste0('files/model_select/sig.coef_zscore_', i, '.txt'), sep = '\t', quote = F)
        if (length(names(sig.coef)>1)) {
          dat1 <- cbind(time=tra.dat$OS.time,
                        status=tra.dat$OS,
                        tra.dat[,match(names(sig.coef),
                                       colnames(tra.dat))])
          fmla <- as.formula(paste0("Surv(time, status) ~"
                                    ,paste0(names(sig.coef),collapse = '+')))
          cox <- coxph(fmla, data =as.data.frame(dat1))
          # lan=coef(cox)
          # print(lan)
          cox1=step(cox, trace = 0)
          lan=coef(cox1)
          write.table(lan, file = paste0('files/model_select/lan_zscore_', i, '.txt'), sep = '\t', quote = F)
          # lan=sig.coef
          final_gene=names(cox1$coefficients)
          # final_gene=names(cox$coefficients)
          GeneNums <-c(GeneNums, length(final_gene))
          
          risk.tra=as.numeric(lan%*%as.matrix(t(tra.dat[,final_gene])))
          # risk.tra=zscore(risk.tra)
          ROC.DSST <- timeROC(T=tra.dat$OS.time,
                              delta=tra.dat$OS,
                              marker=risk.tra,
                              cause=1,
                              weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),
                              iid=TRUE)
          tra.auc=c(tra.auc,max(ROC.DSST$AUC))
          trap=plotKMCox(data.frame(tra.dat$OS.time,tra.dat$OS,ifelse(risk.tra>=median(risk.tra),'H','L')))
          tra.p=c(tra.p,trap)
          risk.test=as.numeric(lan%*%as.matrix(t(test.dat[,final_gene])))
          # risk.test=zscore(risk.test)
          ROC.DSST1=timeROC(T=test.dat$OS.time,delta=test.dat$OS,marker=risk.test,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          test.auc=c(test.auc,max(ROC.DSST1$AUC))
          testp=plotKMCox(data.frame(test.dat$OS.time,test.dat$OS,ifelse(risk.test>=median(risk.test),'H','L')))
          test.p=c(test.p,testp)
          risk.all=as.numeric(lan%*%as.matrix(t(tcga.dat[,final_gene])))
          # risk.all=zscore(risk.all)
          ROC.DSST3=timeROC(T=tcga.dat$OS.time,delta=tcga.dat$OS,marker=risk.all,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          all.auc=c(all.auc,max(ROC.DSST3$AUC))
          allp=plotKMCox(data.frame(tcga.dat$OS.time,tcga.dat$OS,ifelse(risk.all>=median(risk.all),'H','L')))
          all.p=c(all.p,allp)
          # final_gene1=as.character(gene.type[final_gene,1])
          
          if (length(intersect(final_gene,colnames(geo1.dat)))==length(final_gene)) {
            risk.geo1=as.numeric(lan%*%as.matrix(t(geo1.dat[,intersect(final_gene,colnames(geo1.dat))])))
            # risk.geo1=zscore(risk.geo1)
            ROC.DSST2=timeROC(T=geo1.dat$OS.time,delta=geo1.dat$OS,marker=risk.geo1,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo1.auc=c(geo1.auc,max(ROC.DSST2$AUC))
            geop=plotKMCox(data.frame(geo1.dat$OS.time,geo1.dat$OS,ifelse(risk.geo1>=median(risk.geo1),'H','L')))
            geo1.p=c(geo1.p,geop)
          }else{
            geo1.auc=c(geo1.auc,NA)
            geo1.p=c(geo1.p,NA)
          }
          
          if (length(intersect(final_gene,colnames(geo2.dat)))==length(final_gene)) {
            risk.geo2=as.numeric(lan%*%as.matrix(t(geo2.dat[,intersect(final_gene,colnames(geo2.dat))])))
            # risk.geo2=zscore(risk.geo2)
            ROC.DSST4=timeROC(T=geo2.dat$OS.time,delta=geo2.dat$OS,marker=risk.geo2,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo2.auc=c(geo2.auc,max(ROC.DSST4$AUC))
            geop=plotKMCox(data.frame(geo2.dat$OS.time,geo2.dat$OS,ifelse(risk.geo2>=median(risk.geo2),'H','L')))
            geo2.p=c(geo2.p,geop)
          }else{
            geo2.auc=c(geo2.auc,NA)
            geo2.p=c(geo2.p,NA)
          }
          print(c('AUC',max(ROC.DSST$AUC),max(ROC.DSST1$AUC),max(ROC.DSST3$AUC)))
          print(c('P',trap,testp,allp))
          print(c('num',length(final_gene)))
          print("........")
          SampleingTime=c(SampleingTime,i)
          gene.list=c(gene.list,as.character(final_gene))
          SamplingSeed=c(SamplingSeed,rep(i,length(final_gene)))
        }
      }
    }
  }
  myd.clustering.df=data.frame("SamplingTime"=SampleingTime,
                               "TrainRiskP"=tra.p,"TrainRiskAUC"=tra.auc,
                               "TestRiskP"=test.p,"TestRiskAUC"=test.auc,
                               "TCGARiskP"=all.p,"TCGARiskAUC"=all.auc,
                               "GEO1RiskP"=geo1.p,"GEO1RiskAUC"=geo1.auc,
                               "GEO2RiskP"=geo2.p,"GEO2RiskAUC"=geo2.auc,
                               "GeneNums"=GeneNums)
  sample.index=as.data.frame(sample.index)
  colnames(sample.index)=paste("seed",SampleingTime,sep="")
  myd.clustering.genes=data.frame("SamplingSeed"=SamplingSeed,"Genes"=gene.list)
  return(list(myd.clustering.df,sample.index,myd.clustering.genes))
}

tcga_model_data <- cbind(tcga_cli[, c("OS.time", "OS")],
                         t(tcga_exp[intersect(corMatrix$gene,rownames(tcga_exp)), tcga_cli$Samples]))
dim(tcga_model_data)
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
tcga_model_data <- crbind2DataFrame(tcga_model_data)
myd_exp_resampling <- select_gene_zscore(dat1 = tcga_model_data,
                                         # dat2 = silu_model_data,
                                         # dat3 = GSE84437_model_data,
                                         a = 1,
                                         n = 5000,
                                         ratio = 0.7,
                                         cut_p = 0.05,
                                         years = c(1:5))
myd_exp_resampling[[1]]

write.csv(myd_exp_resampling[[1]],file = "04_model/myd_exp_resampling.csv")
###############
SamplingTime <- read.table("04_model/select_SamplingTime.txt",header = T,sep = "\t")
SamplingTime<- SamplingTime$SamplingTime
length(SamplingTime)
num <- SamplingTime[18]#OK  18：ok
tra.samples <- rownames(read.delim(paste0('files/model_select1/tra.dat_zscore_',num,'.txt'), 
                                   header = T, row.names = 1, stringsAsFactors = F))
test.samples <- rownames(read.delim(paste0('files/model_select1/test.dat_zscore_',num,'.txt'),
                                    header = T, row.names = 1, stringsAsFactors = F))
tra.data <- tcga_model_data[tra.samples, ]
dim(tra.data)
test.data <- tcga_model_data[test.samples, ]
dim(test.data)

write.csv(data.frame(cohort=rep(c('Train','Test'),c(length(tra.samples),length(test.samples))),
                     tcga_cli[c(tra.samples,test.samples),
                                    ]),
          '04_model/TCGA_clinical.csv',quote = F,row.names = F)


tra.cox=cox_batch(dat = tcga_exp[corMatrix$gene,tra.samples],
                  time = tcga_cli[tra.samples,]$OS.time,event = tcga_cli[tra.samples,]$OS)
tra.cox=na.omit(tra.cox)
head(tra.cox)

rownames(tra.cox)=gsub('-','__',rownames(tra.cox))
p_cutoff=0.05
table(tra.cox$p.value<p_cutoff)
# FALSE  TRUE 
# 423    21  
tra.cox.fit=tra.cox[tra.cox$p.value<p_cutoff,]
tra.cox.fit
write.csv(tra.cox.fit,file = "04_model/tra.cox.fit.csv")

#########lasso
library(glmnet)
set.seed(num)
fit1=glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
            ,cbind(time=tra.data$OS.time,
                   status=tra.data$OS)
            ,family="cox",nlambda=100, alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
                  ,cbind(time=tra.data$OS.time,
                         status=tra.data$OS)
                  ,family="cox",nlambda=100, alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
names(sig.coef)
length(names(sig.coef))
pdf('04_model/LASSO.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

fmla <- as.formula(paste0("Surv(OS.time, OS) ~",paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tra.data))
cox=step(cox)

module.coxforest=ggforest(cox, data = tra.data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest
ggsave('04_model/gene_forest.pdf',module.coxforest,height = 4,width = 9)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"0.455*KCNK15+0.26*CXCL3+-0.715*PIH1D2+-0.524*DTX1+0.295*LIF+-0.443*TCN2+-0.26*SERPINF2+-0.286*GCNT2"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)
####
###########
lan.dataframe <- as.data.frame(lan)
lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])
# 
lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)
# 
p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FC4100", "Negative" = "#00215E")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
ggsave('04_model/gene_ Coefficients.pdf',p1,height = 4,width = 9)

##4.1########
risktype.col=c('#FFB4C2',"#667BC6")
risk.tra=as.numeric(lan%*%as.matrix(t(tra.data[tra.samples,names(lan)])))
tra.risktype.cli=data.frame(tcga_cli[tra.samples,],Riskscore=risk.tra)

#######
tra.data.point <- surv_cutpoint(tra.risktype.cli, time = "OS.time", event = "OS",
                                variables = 'Riskscore')
tra.cutoff <- as.numeric(summary(tra.data.point)[1])
tra.cutoff
tra.risktype.cli$Risktype=ifelse(tra.risktype.cli$Riskscore>tra.cutoff,'High','Low')
#tra.risktype.cli$Risktype=ifelse(tra.risktype.cli$Riskscore>median(risk.tra),'High','Low')
tra.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                               data = tra.risktype.cli),
                  data=tra.risktype.cli,
                  conf.int = T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                  surv.median.line = 'hv',title='Train cohort',
                  linetype = c("solid", "dashed","strata")[1],
                  palette = risktype.col,ggtheme = custom_theme(),
                  legend = c(0.8,0.75), # 
                  legend.title = "Risktype",legend.labs=c('High','Low'))
tra.km=mg_merge_plot(tra.km$plot,tra.km$table,nrow=2,heights = c(2.5,1),align = 'v')
tra.km


tra.roc=ggplotTimeROC(tra.risktype.cli$OS.time,
                      tra.risktype.cli$OS,
                      tra.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tra.roc


write.table(tcga.risktype.cli,file = "04_model/tcga.risktype.cli.txt",row.names = T,col.names = T,quote = F,sep = "\t")
########4.2########
risk.test=as.numeric(lan%*%as.matrix(t(test.data[test.samples,names(lan)])))
test.risktype.cli=data.frame(tcga_cli[test.samples,],Riskscore=risk.test)
###
test.data.point <- surv_cutpoint(test.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
test.cutoff <- as.numeric(summary(test.data.point)[1])
test.cutoff
test.risktype.cli$Risktype=ifelse(test.risktype.cli$Riskscore>test.cutoff,'High','Low')
#test.risktype.cli$Risktype=ifelse(test.risktype.cli$Riskscore>median(risk.test),'High','Low')
test.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = test.risktype.cli),
                   data=test.risktype.cli,
                   conf.int = T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='Test cohort',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 
                   legend.title = "Risktype",legend.labs=c('High','Low'))
test.km=mg_merge_plot(test.km$plot,test.km$table,nrow=2,heights = c(2.5,1),align = 'v')
test.km


test.roc=ggplotTimeROC(test.risktype.cli$OS.time,
                       test.risktype.cli$OS,
                       test.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
test.roc









##4.3TCGA########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga_cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga_cli,Riskscore=risk.tcga)
########
tcga.data.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
tcga.cutoff <- as.numeric(summary(tcga.data.point)[1])
tcga.cutoff
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
#tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')
tcga.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga.risktype.cli),
                   data=tcga.risktype.cli,
                   conf.int =T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='TCGA',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 
                   legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km=mg_merge_plot(tcga.km$plot,tcga.km$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km

pdf("04_model/all_tcga_roc.pdf",height = 8,width = 8)
mg_surv_pROC(tcga.risktype.cli$OS.time,
             tcga.risktype.cli$OS,
             tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
dev.off()

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tcga.roc
########4.4 GSE44001##########
GSE44001_model_data <- data.frame(GSE44001_cli[, c("OS.time", "OS")],
                                  t(GSE44001_exp[intersect(names(lan),rownames(GSE44001_exp)), GSE44001_cli$Samples]))
colnames(GSE44001_model_data) <- gsub('-', '_', colnames(GSE44001_model_data))

risk.GSE44001=as.numeric(lan%*%as.matrix(t(GSE44001_model_data[GSE44001_cli$Samples,names(lan)])))

GSE44001.risktype.cli=data.frame(GSE44001_cli,Riskscore=risk.GSE44001)
#GSE44001.risktype.cli$Risktype=ifelse(GSE44001.risktype.cli$Riskscore>median(risk.GSE44001),'High','Low')
GSE44001.data.point <- surv_cutpoint(GSE44001.risktype.cli, time = "OS.time", event = "OS",
                                     variables = 'Riskscore')
GSE44001.cutoff <- as.numeric(summary(GSE44001.data.point)[1])
GSE44001.cutoff
GSE44001.risktype.cli$Risktype=ifelse(GSE44001.risktype.cli$Riskscore>GSE44001.cutoff,'High','Low')
GSE44001.roc=ggplotTimeROC(GSE44001.risktype.cli$OS.time,
                           GSE44001.risktype.cli$OS,
                           GSE44001.risktype.cli$Riskscore,mks = c(1,2,3,4))
GSE44001.roc
GSE44001.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = GSE44001.risktype.cli),
                       data=GSE44001.risktype.cli,
                       conf.int = T,pval = T,risk.table = T,
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE44001',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Overall Survival(OS)',
                       legend=c(0.85,0.25),#
                       ggtheme = theme_bw(base_size = 12))
GSE44001.km=mg_merge_plot(GSE44001.km$plot,GSE44001.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE44001.km

####TCGA#####
my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)[7:9]){
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
tcga.risktype.cli$Status <- ifelse(tcga.risktype.cli$OS==0,"alive","dead")
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
##########
p4=mg_merge_plot(mg_merge_plot(tra.km,test.km,tcga.km,common.legend = T,labels = c('C','D','E'),ncol=3),
                 mg_merge_plot(tra.roc,test.roc,tcga.roc,labels = c('F','G','H'),ncol=3),
                 mg_merge_plot(tcga.barplot,GSE44001.km,GSE44001.roc,labels = c('I','J','K'),ncol=3),
                 nrow = 3)

p4
ggsave('04_model/P4c-k.pdf',height = 15,width = 15)



#05.######
dir.create('05_immu')
##estimate####
tcga_estimate <- immu_estimate(exp = tcga_exp)
p5a <-  my_mutiboxplot(tcga_estimate ,group = tcga.risktype.cli$Risktype,group_cols=risktype.col,legend.pos='top',ylab = "Estimate score")
#####TIMER#####
tcga.timer=immu_timer(tcga_exp)
p5b<- mg_PlotMutiBoxplot(data = tcga.timer[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                         legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'Timer score')
p5ab<- mg_merge_plot(p5a,p5b,common.legend =T,widths = c(1,1.3),labels =c ("A","B"),legend = "top")
savePDF(p5ab,file = "05_immu/p5ab.pdf",height = 4,width = 10)
###mcp###
tcga.mcp <- immu_MCPcounter(exp =tcga_exp,isTCGA = T)
saveRDS(tcga.mcp,file ='05_immu/taga.mcp.RDS')
mg_PlotMutiBoxplot(data = tcga.mcp[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                   legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'ssgsea Immune Score')
library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                    y=tcga.mcp[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',MCP_count=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
corrmat.color=colorRampPalette(c('#36BA98', 'white','#F4A261'))(100)

p5c<-quickcor(tcga.mcp[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")
p5c
ggsave('05_immu/p5c.pdf',height = 6,width = 6)
# ICGs####
tcga.icgs=immu_ICGs(tcga_exp)
colnames(tcga.icgs)
icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,tcga_model_data[tcga.risktype.cli$Samples,names(lan)]
                 ,tcga.icgs[tcga.risktype.cli$Samples,c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')])
#c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))

pdf('05_immu/p5d.pdf',height = 6,width = 7,onefile = F)
pheatmap(icg_cor_res$r[-c(1:6),c(names(lan),'Riskcsore')],
             color = circlize::colorRamp2(c(-1, 0, 1), c('#071952', 'white', '#FFB1B1')),
             main="Heatmap", # 
             display_numbers = icg_cor_res.p[-c(1:6),c(names(lan),'Riskcsore')], # 
             cluster_cols = F, # 
             cluster_rows = F,
             show_rownames = T, #
             show_colnames = T,
             fontsize_row = 12, # 
             fontsize_col = 16)
dev.off()

####6####
dir.create('06_risktype.mut')
##########
#tcga.maf=getTCGAMAFByCode('CESC')#
tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.use.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.use.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.use.high,file='06_risktype.mut/tcga.risktype.use.high.txt')
write.table(tcga.risktype.use.low,file='06_risktype.mut/tcga.risktype.use.low.txt')

tcga.maf.high=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.high$Tumor_Sample_Barcode))
tcga.maf.high<-read.maf(tcga.maf.high@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.high.txt')
tcga.maf.high@clinical.data

tcga.maf.low=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.low$Tumor_Sample_Barcode))
tcga.maf.low<-read.maf(tcga.maf.low@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.low.txt')
tcga.maf.low@clinical.data

#######
pdf('06_risktype.mut/Fig6A.pdf',height = 4,width = 7,onefile = F)
oncoplot(maf = tcga.maf.high,top = 10,sortByAnnotation = T)
dev.off()
pdf('06_risktype.mut/Fig6B.pdf',height = 4,width = 7,onefile = F)
oncoplot(maf = tcga.maf.low,top = 10,sortByAnnotation = T)
dev.off()

######### Somatic Interactions
pdf('06_risktype.mut/Fig6C.pdf',height = 5,width = 5,onefile = F)
somaticInteractions(maf = tcga.maf.high,pvalue = c(0.05,0.1),top = 10)
dev.off()

pdf('06_risktype.mut/Fig6D.pdf',height = 5,width = 5,onefile = F)
somaticInteractions(maf = tcga.maf.low,pvalue = c(0.05,0.1),top = 10)
dev.off()
######ssGSEA#####
library(GSEABase)
library(GSVA)
# h.gmt <- getGmt("C:/Users/Administrator/Desktop/data/h.all.v2023.2.Hs.symbols.gmt",
#                 collectionType=BroadCollection(category="h"),
#                 geneIdType=SymbolIdentifier())
# tcga.hall.ssGSEA <- gsva(as.matrix(tcga_exp),h.gmt, method='ssgsea',verbose=TRUE)
# 
# saveRDS(tcga.hall.ssGSEA,file = "06_risktype.mut/tcga.hall.ssGSEA.RDS")
tcga.hall.ssGSEA <- readRDS("06_risktype.mut/tcga.hall.ssGSEA.RDS")
rownames(tcga.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.hall.ssGSEA))
rownames(tcga.hall.ssGSEA)=gsub('_',' ',rownames(tcga.hall.ssGSEA))
tcga.hall.ssGSEA[1:5,1:5]

tcga.pathway.ssgsea=tcga.hall.ssGSEA
tcga.pathwat.dat=cbind.data.frame(Riskscore=tcga.risktype.cli[,'Riskscore'],
                                  t(tcga.pathway.ssgsea[,tcga.risktype.cli$Samples]))
tcga.pathwat.dat[1:5,1:5]
pathway.cor.RS=Hmisc::rcorr(as.matrix(tcga.pathwat.dat),type = 'spearman')
pathway.cor.RS.res=data.frame(Pathway=names(pathway.cor.RS$r['Riskscore',]),
                              cor=as.numeric(pathway.cor.RS$r['Riskscore',]),
                              p.val=as.numeric(pathway.cor.RS$P['Riskscore',]))
head(pathway.cor.RS.res)
pathway.cor.RS.res=pathway.cor.RS.res[-1,]
colnames(pathway.cor.RS.res)=c('Pathway','cor','pvalue')

pathway.cor.RS.res$Pathway=gsub('_',' ',pathway.cor.RS.res$Pathway)
pathway.cor.RS.res$Pathway=factor(pathway.cor.RS.res$Pathway,
                                  levels = pathway.cor.RS.res$Pathway[order(pathway.cor.RS.res$cor,decreasing = T)], ordered=TRUE)

pathway.cor.RS.res$pvalue=ifelse(pathway.cor.RS.res$pvalue==0,1e-16,pathway.cor.RS.res$pvalue)
head(pathway.cor.RS.res)

pathway.cor.RS.res$group=rep('no signif',nrow(pathway.cor.RS.res))
pathway.cor.RS.res$group[which(pathway.cor.RS.res$cor>0 & pathway.cor.RS.res$pvalue<0.05)]='Positive'
pathway.cor.RS.res$group[which(pathway.cor.RS.res$cor<0 & pathway.cor.RS.res$pvalue<0.05)]='Negative'

pdf('06_risktype.mut/Fig6e.pdf',height = 10,width = 10,onefile = F)
ggplot(pathway.cor.RS.res,aes(x=cor,y=Pathway))+
  # geom_bar()+
  scale_fill_gradientn(colors=  c("#219C90",  "#FFF455"))+
  xlab('Cor')+#ylab('HALLMARK Pathway')+
  geom_bar(aes(fill=-log10(pvalue)),stat="identity",position = "dodge") +
  # scale_fill_continuous_sequential(palette = "Viridis",rev = T)+
  #theme_classic()+
  theme(axis.text.y = element_text(size =10),
        legend.position = 'right')
dev.off()




##########GSEA###############
tcga.geneList=getGeneFC(gene.exp=tcga_exp,group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')

gmt2list <- function(annofile){
  
  if (!file.exists(annofile)) {
    stop("There is no such a gmt file!")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
  return(annoList)
}

names(Hall_gmt)=gsub('HALLMARK_','',names(Hall_gmt))
names(Hall_gmt)=gsub('_',' ',names(Hall_gmt))
library(fgsea)
fgseaRes <- fgsea(pathways = Hall_gmt,
                  stats =tcga.geneList ,
                  minSize=10,
                  maxSize=500,
                  nperm=1000)
head(fgseaRes)

nrow(fgseaRes[padj<0.05 & NES > 0])
#5
nrow(fgseaRes[padj<0.05 & NES < 0])
#6

topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=5), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=6), pathway]
topPathways <- c(topPathwaysUp,topPathwaysDown)

library(forcats)
library(ggplot2)
library(ggstance)
fgseaRes$sign<-ifelse(fgseaRes$NES>0,"Activated","Suppressed") ## 

####
pdf('06_risktype.mut/Fig6f_Up.pdf',height = 5,width = 8,onefile = F)
plotGseaTable(Hall_gmt[topPathwaysUp ], tcga.geneList, fgseaRes, 
              gseaParam = 0.5)
dev.off()

pdf('06_risktype.mut/Fig6f_down.pdf',height = 5,width = 8,onefile = F)
plotGseaTable(Hall_gmt[topPathwaysDown ], tcga.geneList, fgseaRes, 
              gseaParam = 0.5)
dev.off()


##########07.###########
dir.create('07_drug')
##7.1tide#####
tcga_tide_dat <- t(scale(t(tcga_exp),scale = F))
dim(tcga_tide_dat)
write.table(tcga_tide_dat,file = '07_drug/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('07_drug/CESC_TIDE.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
tcga_tide_res=cbind(tcga_tide_res[tcga.risktype.cli$Samples,],tcga.risktype.cli)

tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')
tcga_tide_list <- list()
for (fea in c("TIDE","Dysfunction","Exclusion","MDSC","CAF","TAM.M2")) {
  print(fea)
  tmp_plot <- mg_violin_1(data.frame(tcga.risktype.cli$Risktype
                                   ,tcga_tide_res[tcga.risktype.cli$Samples, fea])
                        ,melt = T
                        ,ylab = fea
                        ,jitter=T
                         ,group_col = risktype.col#pal_jco()(9)[1:2]
                        ,test_method = 'wilcox.test'
                        ,cmp_test_method = 'wilcox.test'
                        ,legend.pos = 'tr'
                        ,show_compare = T)
  tcga_tide_list[[fea]] <- tmp_plot
}
fig7a <- cowplot::plot_grid(plotlist = tcga_tide_list,
                                      ncol = 6)
fig7a
# ggsave(plot = fig7a,
#        filename = '07_drug/fig7A_tide_plot.pdf',
#        width = 15, height = 5)


fig7b=my_mutiboxplot(dat = t(tcga_exp[c('CTLA4','HAVCR2','LAG3','TIGIT','PDCD1','CD274'),tcga.risktype.cli$Samples]),
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,
                     ylab = 'Expression',angle = 0,hjust = 0.5)
fig7ab=mg_merge_plot(fig7a,fig7b,nrow = 2,labels = c('A','B'),align = 'h',common.legend = F)
savePDF('07_Immunotherapy/Fig7ab.pdf',fig7ab,height = 4,width = 12)

# 7.3tcga drug ######
library(pRRophetic)
library(ggplot2)


## Cisplatin,
set.seed(12345)
tcga_Cisplatin <- pRRopheticPredict(as.matrix(tcga_exp[, tcga.risktype.cli$Samples])
                                    , "Cisplatin"
                                    , "urogenital_system"
                                    , selection=1
                                    ,dataset = "cgp2016")
tcga_Cisplatin <- data.frame(tcga_Cisplatin)

tcga_durg_ic50_res <- tcga_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
dim(tcga_durg_ic50_res)
for (drug in drugs) {
  print(drug)
  set.seed(12345)
  tmpic50 <- pRRopheticPredict(as.matrix(tcga_exp[, tcga.risktype.cli$Samples])
                               , drug
                               , "urogenital_system"
                               , selection=1
                               , dataset = "cgp2016")
  tmpic50 <- data.frame(tmpic50)
  colnames(tmpic50) <- drug
  tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
}

colnames(tcga_durg_ic50_res)
tcga_durg_ic50_res <- tcga_durg_ic50_res[, -1]
class(tcga_durg_ic50_res)
dim(tcga_durg_ic50_res)
num <- 11:20
mg_PlotMutiBoxplot(tcga_durg_ic50_res[tcga.risktype.cli$Samples, num]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   , add = 'boxplot'
                   , xangle =45
                   , ylab = 'Estimated IC50'
                   # , xlab = 'drug'
                   , group_cols = risktype.col
                   , test_method = 'wilcox.test')
tcga_durg_ic50_res[, c("Cisplatin", "Erlotinib", "Sunitinib", "Paclitaxel", "Sorafenib", "Crizotinib")]


tcga_drug_plot_list <- list()
for (dr in c("Cisplatin", "Erlotinib", "Sunitinib", "Paclitaxel", "Sorafenib", "Crizotinib")) {
  tmp <- mg_violin_1(data.frame(tcga.risktype.cli$Risktype,
                              tcga_durg_ic50_res[tcga.risktype.cli$Samples, dr])
                   ,melt = T
                   ,ylab = 'Estimated IC50'
                   ,xlab = dr
                   # ,ylim = c(0,100)
                   ,test_method = 'wilcox.test'
                   ,cmp_test_method = 'wilcox.test'
                   ,legend.pos = 'tr'
                   ,jitter=T
                   #,#Cluster.col = risktype.col#pal_aaas()(9)[c(2,4,3)]
                   ,show_compare = T,
                   group_col = risktype.col
                   )
  tcga_drug_plot_list[[dr]] <- tmp
}


Fig7c <- cowplot::plot_grid(plotlist = tcga_drug_plot_list,
                                     ncol = 6)
Fig7c 
ggsave(plot = Fig7c ,
       filename = '07_drug/Fig7c .pdf',
       width = 15, height = 5)


pdf('07_drug/icg_tide_drug.pdf', width = 18, height = 15)
cowplot::plot_grid(fig7ab ,
                   Fig7c ,
                   nrow = 2,
                   ncol = 1,
                   rel_heights = c(2,1),
                   labels = c("","C"))
dev.off()




save.image(file = 'CESC_SiaTs2.RData')


