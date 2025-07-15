library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(micropan)
library(RColorBrewer)
library(iNEXT)
library(viridis)
#library(breakaway)
library(dpseg)
library(KneeArrower)
library(ggpmisc)

setwd("C:/Users/Xiaoqian/Desktop/pop_gen/UHGG_plus4/SGBdiffH_500_1pois_1nb_finda/")
lb_list <- Sys.glob("SGB*diffH.txt/SGBdiffH.length_bias_500.txt")

lb_list_names=sapply(strsplit(lb_list,split = 'diffH'),function(x) x[1])
lbfile_list=lapply(lb_list, function(x) read.table(x,sep='\t',header = T))
lbfile_list=mapply(function(x,y) {x$cluster=y;x},lbfile_list,lb_list_names,SIMPLIFY=F)

taxa_list=read.table("../fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice(which.max(n))
lbfile_list_summary=bind_rows(lbfile_list)
lbfile_list_summary=merge(lbfile_list_summary,taxa_list_summary,by='cluster')

#

lbfile_list_summary_n20=lbfile_list_summary %>% group_by(cluster) %>% as.data.frame()
lbfile_list_summary_n20$species=sapply(strsplit(lbfile_list_summary_n20$fgspecies,split = ';'),'[[',3)
lbfile_list_summary_n20$genus=sapply(strsplit(lbfile_list_summary_n20$fgspecies,split = ';'),'[[',2)
lbfile_list_summary_n20=lbfile_list_summary_n20 %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                            TRUE ~ paste0(species)))
lbfile_list_summary_n20$cluster_taxa=paste(lbfile_list_summary_n20$taxa,lbfile_list_summary_n20$cluster,sep = "_")
lbfile_list_summary_n20=lbfile_list_summary_n20 %>% arrange(genus)
lbfile_list_summary_n20[is.na(lbfile_list_summary_n20$Strain.1),]

lbfile_list_summary_n20_list=lbfile_list_summary_n20 %>% group_by(cluster_taxa) %>% group_split()

##Use the clonal divergence after doing the bi-partitioning

lbfile_list_summary_n20_filter=lbfile_list_summary_n20 %>% 
  filter(mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) %>%
  mutate(totalR=sim_fr,div=mu_div,Rdiv=mnb_div) %>%
  mutate(totalR=case_when(sim_fr/hmm_fr>2 ~ hmm_fr, TRUE ~ sim_fr),
         div=case_when(sim_fr/hmm_fr>2 ~ hmm_mu, TRUE ~ mu_div),
         Rdiv=case_when(sim_fr/hmm_fr>2 ~ hmm_mnb, TRUE ~ mnb_div),
         type=case_when(sim_fr/hmm_fr>2 ~ 'C', TRUE ~ 'NC')) %>%
  mutate(totalR=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1  ~ 0, TRUE ~ totalR),
         div=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 10^-5, TRUE ~ div),
         Rdiv=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 0, TRUE ~ Rdiv ))

lbfile_list_summary_n20_filter = lbfile_list_summary_n20_filter %>% group_by(cluster_taxa) %>%
  filter(div>=10^-5) %>%
  distinct(div,totalR, .keep_all = TRUE) %>%
  filter(n()>=5) %>%
 filter(sum(totalR)>0)

#add this step since we are measuring recomb. fraction to mutations
lbfile_list_summary_n20_filter$div=lbfile_list_summary_n20_filter$div*lbfile_list_summary_n20_filter$Alignment.size.iter1*(100-lbfile_list_summary_n20_filter$totalR)/100
lbfile_list_summary_n20_filter$totalR=lbfile_list_summary_n20_filter$totalR/100
#filter list for all taxa with sweeps
sweep_taxa=read.csv('../instrain/confirmed_sweep_from_pairwise.csv')
lbfile_list_summary_n20_filter=lbfile_list_summary_n20_filter %>% 
  filter(cluster %in% sweep_taxa$SGB)
lbfile_list_summary_n20_filter_list=split(lbfile_list_summary_n20_filter,f=lbfile_list_summary_n20_filter$cluster_taxa)

p1=c(0.2,20,40)
p2=c(0.1,10,20)
p3=c(0.1,5,10)

segs_list={}

for (j in 1:length(lbfile_list_summary_n20_filter_list)){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  if (sum(test$x<0.001*2000000)>100){
    ss=seq(1, nrow(test), round(nrow(test)/100))
    test=test[ss,]
  }
  segs={}
  segs[[1]] <- dpseg(x=test$x, y=test$y, jumps=FALSE, P=p1[1],minl=p1[2],maxl=p1[3])
  segs[[2]] <- dpseg(x=test$x, y=test$y, jumps=FALSE, P=p2[1],minl=p2[2],maxl=p2[3])
  segs[[3]] <- dpseg(x=test$x, y=test$y, jumps=FALSE, P=p3[1],minl=p3[2],maxl=p3[3])
  segs[[4]] <- dpseg(x=test$x, y=test$y, jumps=FALSE, P=p1[1],minl=p1[2])
  if (segs[[4]]$segments$r2[1]>0.80 & segs[[4]]$segments$slope[1]>0){
    c=4
  }else {
    c1=0;c2=0;c3=0
    if (segs[[1]]$segments$slope[1]>0){
      c1=segs[[1]]$segments$r2[1]
    }
    if (segs[[2]]$segments$slope[1]>0){
      c2=segs[[2]]$segments$r2[1]
    }
    if (segs[[3]]$segments$slope[1]>0){
      c3=segs[[3]]$segments$r2[1]
    }
    c=which.max(c(c1,c2,c3))
    print (c)
  }
  segs_list[[j]]=segs[[c]]
}

#Use curves that have R2>0.8 directly, otherwise 
# if slope of segment 2 is similar to segment 1, select both 1& 2
#else select the first segment 
#also then use a smaller max length??

redo_list=c()
r1_list=c()
png('new_fits_sweepy_only/segs_fit_r1s.png',res=300,width = 8000,height=5000)
par(mfrow=c(8,8))
for (j in 1:length(segs_list)){
  if (segs_list[[j]]$segments$r2[1]>0.8 & segs_list[[j]]$segments$slope[1]>0){
    r2=round(segs_list[[j]]$segments$r2[1],2)
    plot(segs_list[[j]],main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    r1_list=c(r1_list,j)
  }else {
    redo_list=c(redo_list,j)
  }
}
dev.off()

find_cutoff=function(x){
  ss=x$segments$slope
  cc=1
  if (length(ss)>1){
    while (cc<=length(ss)-1){
      if (ss[cc]*0.75<= ss[cc+1]) {
        cc=cc+1
      } else if (cc<=length(ss)-2){
        if (ss[cc]*0.75<= (ss[cc+1]+ss[cc+2])/2){
          cc=cc+2
        } else {break} 
      } else {break}
    }
  } 
  return (cc)
}

segs_list_2={}
for (j in redo_list ){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  if (sum(test$x<0.001*2000000)>100){
    ss=seq(1, nrow(test), round(nrow(test)/100))
    test=test[ss,]
  }
  if (nrow(segs_list[[j]]$segments)==1){
    test=test
    segs_list_2[[j]]=segs_list[[j]]
  } else {
    n=find_cutoff(segs_list[[j]])
    n1=segs_list[[j]]$segments$start[1]
    n2=segs_list[[j]]$segments$end[n]
    test=test[n1:n2,]
    segs_list_2[[j]]=dpseg(x=test$x, y=test$y, jumps=FALSE, P=0.1,minl=(n2-1),maxl=n2)
  }
}


r2_list=c()
png('new_fits_sweepy_only/segs_fit_rs_p1s.png',res=300,width = 8000,height=5000)
par(mfrow=c(10,8))
for (j in redo_list){
  r2=round(segs_list_2[[j]]$segments$r2[1],2)
  if (r2>0.33 & segs_list_2[[j]]$segments$slope[1]>0){
    plot(segs_list_2[[j]],main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    r2_list=c(r2_list,j)
  }
}
dev.off()

#check ones that are eliminated
r12_list=c(r1_list,r2_list)
new_redo_list=c()
png('new_fits_sweepy_only/r2_omit_summary_2s.png',width=8000,height=6000,res=300)
par(mfrow=c(10,7))
for (j in 1:length(lbfile_list_summary_n20_filter_list)){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  test_o=test
  if (nrow(test)>100){
    if (!(j %in% r12_list)){
      plot(test_o,pch=20,col='grey60',
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],"_",j),xlab="",ylab="",ylim=c(0,1))
      new_redo_list=c(new_redo_list,j)
    }
  }
}
dev.off()


fit_lists=c(r1_list,r2_list)
all_lists=1:length(lbfile_list_summary_n20_filter_list)
unfit_lists=all_lists[!(all_lists %in% fit_lists)]

#we manually check what are the best fits for each figure
cutoff_list=c(1500,900,10000,5000,5000,5000,700,2500,8000,3500,1000,3500,2000,
              900,1000,1000)

#for (j in unfit_lists ){
#  l=lbfile_list_summary_n20_filter_list[[j]]
#  l=l %>% arrange(div)
#  test=data.frame(x=l$div,y=l$totalR)
#  cutoff=cutoff_list[c]
#  test_o=test[test$x<cutoff,]
#  test_o=rbind(c(0,0),test_o)
#  ft=lm(y~0+x,test_o)
#  r2=round(summary(ft)$r.squared,2)
#  plot(test,pch=20,col='grey60',
#       main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
#  points(test_o,col='red',pch=20)
#  abline(ft)
#  c=c+1
#}




seg_slope_list={}
c=1
d=1
png('new_fits_sweepy_only/r2_good_fits_summary.png',width=8000,height=6000,res=300)
par(mfrow=c(10,7))
for (j in 1:length(lbfile_list_summary_n20_filter_list)){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  
  if (j %in% r1_list){
    test_o=test
#    if (sum(test$x<0.001*2000000)>100){
#      ss=seq(1, nrow(test), round(nrow(test)/100))
#      test=test[ss,]
#    }
    
    n1=segs_list[[j]]$segments$start[1]
    n2=segs_list[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      plot(test_o,pch=20,col='grey60',
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
      points(test1,col='red',pch=20)
      abline(ft)
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if (j %in% r2_list) {
    test_o=test
  #  if (sum(test$x<0.001*2000000)>100){
 #     ss=seq(1, nrow(test), round(nrow(test)/100))
 #     test=test[ss,]
 #   }
    n1=segs_list_2[[j]]$segments$start[1]
    n2=segs_list_2[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      plot(test_o,pch=20,col='grey60',
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
      points(test1,col='red',pch=20)
      abline(ft)
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if  (j %in% unfit_lists) {
    test_o=test
    cutoff=cutoff_list[d]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    plot(test_o,pch=20,col='grey60',
         main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
    points(test1,col='red',pch=20)
    abline(ft)
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    d=d+1
  }
}
dev.off()

#check for plots that didn't fit well
badfit_lists=c(22,25,32,47,57,58,59)
cutoff_list_2=c(1200,1200,1000,2000,1000,2500,2500)
#c=1
#for (c in 1:length(badfit_lists)){
#  j=badfit_lists[[c]]
#  l=lbfile_list_summary_n20_filter_list[[j]]
#  l=l %>% arrange(div)
#  test=data.frame(x=l$div,y=l$totalR)
#  cutoff=cutoff_list_2[c]
#  test_o=test[test$x<cutoff,]
#  test_o=rbind(c(0,0),test_o)
#  ft=lm(y~0+x,test_o)
#  r2=round(summary(ft)$r.squared,2)
#  plot(test,pch=20,col='grey60',
#       main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
#  points(test_o,col='red',pch=20)
# abline(ft)
# c=c+1
#}

r1_list=r1_list[!(r1_list %in% badfit_lists)]
r2_list=r2_list[!(r2_list %in% badfit_lists)]


#read in taxa classification
taxa_class=read.csv('../instrain/taxa_class_all.csv')
taxa_class_include=c(taxa_class$commensals)

seg_slope_list={}
c=1
d=1
e=1

png('new_fits_sweepy_only/r2_good_fits_summary_2a_commensals.png',width=5500,height=5500/10*12,res=300)
par(mfrow=c(12,4),mai = c(0.3, 0.3, 0.3, 0.1))
for (j in 1:65){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
 
  
  if (j %in% r1_list){
    test_o=test
    #    if (sum(test$x<0.001*2000000)>100){
    #      ss=seq(1, nrow(test), round(nrow(test)/100))
    #      test=test[ss,]
    #    }
    
    n1=segs_list[[j]]$segments$start[1]
    n2=segs_list[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],' ',r2))
      title_text=gsub('_',' ',mt)
      # Split the title into parts before and after "SGB"
      italic_part <- gsub("(.*?)\\sSGB\\s*(.*)", "\\1", title_text)
      remaining_text <- sub(".*SGB\\s*", "(SGB", title_text)
      remaining_text=gsub(' ',') ',remaining_text)
      if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
      plot(test_o,pch=20,col='grey60',
           main=bquote(italic(.(italic_part)) ~ .(remaining_text)),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.8,cex.axis = 1.8)
      points(test1,col='red',pch=20)
      abline(ft)
      }
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if (j %in% r2_list) {
    test_o=test
    #  if (sum(test$x<0.001*2000000)>100){
    #     ss=seq(1, nrow(test), round(nrow(test)/100))
    #     test=test[ss,]
    #   }
    n1=segs_list_2[[j]]$segments$start[1]
    n2=segs_list_2[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      
      mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],' ',r2))
      title_text=gsub('_',' ',mt)
      # Split the title into parts before and after "SGB"
      italic_part <- gsub("(.*?)\\sSGB\\s*(.*)", "\\1", title_text)
      remaining_text <- sub(".*SGB\\s*", "(SGB", title_text)
      remaining_text=gsub(' ',') ',remaining_text)
      if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
      plot(test_o,pch=20,col='grey60',
           main=bquote(italic(.(italic_part)) ~ .(remaining_text)),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.8,cex.axis = 1.8)
      points(test1,col='red',pch=20)
      abline(ft)
      }
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if  (j %in% unfit_lists) {
    test_o=test
    cutoff=cutoff_list[d]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],' ',r2))
    title_text=gsub('_',' ',mt)
    # Split the title into parts before and after "SGB"
    italic_part <- gsub("(.*?)\\sSGB\\s*(.*)", "\\1", title_text)
    remaining_text <- sub(".*SGB\\s*", "(SGB", title_text)
    remaining_text=gsub(' ',') ',remaining_text)
    if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
    plot(test_o,pch=20,col='grey60',
         main=bquote(italic(.(italic_part)) ~ .(remaining_text)),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.8,cex.axis = 1.8)
    points(test1,col='red',pch=20)
    abline(ft)
    }
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    d=d+1
  } else if  (j %in% badfit_lists) {
    test_o=test
    cutoff=cutoff_list_2[e]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],' ',r2))
    title_text=gsub('_',' ',mt)
    # Split the title into parts before and after "SGB"
    italic_part <- gsub("(.*?)\\sSGB\\s*(.*)", "\\1", title_text)
    remaining_text <- sub(".*SGB\\s*", "(SGB", title_text)
    remaining_text=gsub(' ',') ',remaining_text)
    if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
    plot(test_o,pch=20,col='grey60',
         main=bquote(italic(.(italic_part)) ~ .(remaining_text)),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.8,cex.axis = 1.8)
    points(test1,col='red',pch=20)
    abline(ft)
    }
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    e=e+1
  }
}
dev.off()


png('new_fits_sweepy_only/r2_good_fits_summary_2b_new.png',width=5000,height=5000,res=300)
par(mfrow=c(10,4),mai = c(0.3, 0.3, 0.3, 0.1))
for (j in 32:length(lbfile_list_summary_n20_filter_list)){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  
  if (j %in% r1_list){
    test_o=test
    #    if (sum(test$x<0.001*2000000)>100){
    #      ss=seq(1, nrow(test), round(nrow(test)/100))
    #      test=test[ss,]
    #    }
    
    n1=segs_list[[j]]$segments$start[1]
    n2=segs_list[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
      if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
      plot(test_o,pch=20,col='grey60',
           main=mt,xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.5,cex.axis = 1.5)
      points(test1,col='red',pch=20)
      abline(ft)
      }
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if (j %in% r2_list) {
    test_o=test
    #  if (sum(test$x<0.001*2000000)>100){
    #     ss=seq(1, nrow(test), round(nrow(test)/100))
    #     test=test[ss,]
    #   }
    n1=segs_list_2[[j]]$segments$start[1]
    n2=segs_list_2[[j]]$segments$end[1]
    test1=test[n1:n2,]
    test1=rbind(c(0,0),test1)
    if (max(test$y)/min(test$y)>1.5){
      ft=lm(y ~ 0+x,data=test1)
      r2=round(summary(ft)$r.squared,2)
      mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
      if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
      plot(test_o,pch=20,col='grey60',
           main=mt,xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.5,cex.axis = 1.5)
      points(test1,col='red',pch=20)
      abline(ft)
      }
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if  (j %in% unfit_lists) {
    test_o=test
    cutoff=cutoff_list[d]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
    plot(test_o,pch=20,col='grey60',
         main=mt,xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.5,cex.axis = 1.5)
    points(test1,col='red',pch=20)
    abline(ft)
    }
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    d=d+1
  } else if  (j %in% badfit_lists) {
    test_o=test
    cutoff=cutoff_list_2[e]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    mt=gsub('s__','',paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    if(names(lbfile_list_summary_n20_filter_list)[j] %in% taxa_class_include){
    plot(test_o,pch=20,col='grey60',
         main=mt,xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1),cex.main = 1.5,cex.axis = 1.5)
    points(test1,col='red',pch=20)
    abline(ft)
    }
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    e=e+1
  }
}
dev.off()

all_sweep_info_total_list_SGB_summary=read.csv('../instrain/all_sweep_info_total_list_SGB_summary.csv',header = T)

seg_slope_all=bind_rows(seg_slope_list) 
write.csv(seg_slope_all,'seg_slope_all_nonrm.csv',row.names = F)
seg_slope_all=read.csv('seg_slope_all_nonrm.csv')
seg_slope_all$SGB=sapply(strsplit(seg_slope_all$taxa,split = "_"), function(x) x[length(x)])
seg_slope_all=merge(seg_slope_all,all_sweep_info_total_list_SGB_summary,by='SGB',all.x=T)
seg_slope_all$sweepy=seg_slope_all$n_in_sweep/seg_slope_all$total_n
seg_slope_all$sweepy[is.na(seg_slope_all$sweepy)]=0
SGB_class=read.csv('SGB_classify_list.csv')
seg_slope_all=seg_slope_all %>% mutate(class=case_when(taxa %in% SGB_class$pathogen ~ 'pathogens',
                                                       taxa %in% SGB_class$probiotics ~ 'probiotics',
                                                       TRUE ~ 'commensals'))

seg_slope_all_2=seg_slope_all %>% filter(class!='probiotics')
names(seg_slope_all_2)=c('taxa','recombination_rate','class')

write.csv(seg_slope_all_2,'seg_slope_all_nonrm_2.csv',row.names = F)


taxa_list=read.table("../fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice(which.max(n))
taxa_list_summary$species=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',3)
taxa_list_summary$genus=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',2)
taxa_list_summary$family=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',1)
taxa_list_summary=taxa_list_summary %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                TRUE ~ paste0(species)))
taxa_list_summary$cluster_taxa=paste(taxa_list_summary$taxa,taxa_list_summary$cluster,sep = "_")
seg_slope_all=merge(taxa_list_summary,seg_slope_all,by.x = 'cluster_taxa',by.y='taxa')

num_colors <- length(unique(seg_slope_all$family))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

seg_slope_all$log10_total_n=log10(seg_slope_all$total_n)





seg_slope_all$family=as.factor(seg_slope_all$family)
seg_slope_all$single_ratio=seg_slope_all$MG_n/(seg_slope_all$MG_mixed_n+seg_slope_all$MG_n)
seg_slope_all$single_ratio_2=seg_slope_all$MG_n/(seg_slope_all$MG_exist_n)

seg_slope_all=seg_slope_all %>% 
  filter(sweepy>0)
seg_slope_all_c=seg_slope_all %>% filter(class=='commensals') %>% 
  filter(sweepy>0)
seg_slope_all_c_2=seg_slope_all_c %>% group_by(family) %>% 
  filter(n() > 1)
ggplot(seg_slope_all_c_2,aes(x=family,y=slope))+geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

seg_slope_all_anova=aov(slope ~ family, data = seg_slope_all_c_2)
TukeyHSD(seg_slope_all_anova)

seg_sweepy_all_anova=aov(sweepy ~ family, data = seg_slope_all_c_2)
TukeyHSD(seg_sweepy_all_anova)


cor.test(seg_slope_all_c$sweepy,seg_slope_all_c$single_ratio,method='spearman')
cor.test(seg_slope_all_c$sweepy,seg_slope_all_c$slope,method='spearman')
cor.test(seg_slope_all_c$single_ratio,seg_slope_all_c$slope,method='spearman')


#write.csv(seg_slope_all,'seg_slope_all.csv',row.names = F)

#seg_slope_all_norm=read.csv('seg_slope_all_nonrm.csv')

#slope/100 for non-rm, slope for rm
cor.test(seg_slope_all_c$sweepy,seg_slope_all_c$slope,method = 'spearman')
seg_slope_all$family=gsub('f__','',seg_slope_all$family)
seg_slope_all_2=seg_slope_all %>% filter(class!='probiotics')
p2=ggplot(seg_slope_all_2, aes(x = sweepy, y = slope, color = family,group=class,
                             size=log10(total_n))) +
  scale_color_manual(values = getPalette(num_colors)) +
  facet_wrap(~class)+
  #  scale_y_log10()+
  geom_point() +
  stat_correlation(method='spearman',use_label(c("r", "p")),small.r = TRUE,small.p = TRUE,label.x = "right",size=5) +
  #  stat_poly_line(se=F,lty=2)+
  theme_bw()+
  theme(legend.position = "right",text = element_text(size = 18))+
  guides(color=guide_legend(ncol = 3, byrow=TRUE,title.position = "top",title = "Family",
                            override.aes = list(label = "\u25A0", size = 7)),
         size=guide_legend(ncol=4,title="log10 (sweep size)",override.aes = list(linetype = 0)))+
  xlab('Dominance of sweeps')+
  ylab('Fraction genome recombined per mutation')

png('new_fits_sweepy_only/sweep_slope_correlation_nolog_2_2.png', res=300, width=4000, height =1800)
p2
dev.off()

#for commensals only
seg_slope_all_c$family=gsub('f__','',seg_slope_all_c$family)
p2c=ggplot(seg_slope_all_c, aes(x = sweepy, y = slope, color = family,group=class,
                                size=log10(total_n))) +
  scale_color_manual(values = getPalette(num_colors)) +
  # facet_wrap(~class)+
  #  scale_y_log10()+
  geom_point() +
  stat_correlation(method='spearman',use_label(c("r", "p")),small.r = TRUE,small.p = TRUE,label.x = "right",size=8) +
  #  stat_poly_line(se=F,lty=2)+
  theme_bw()+
  theme(legend.position = "right",text = element_text(size = 22))+
  guides(color=guide_legend(ncol = 2, byrow=TRUE,title.position = "top",title = "Family",
                            override.aes = list(label = "\u25A0", size = 7)),
         size=guide_legend(ncol=4,title="log10 (sweep size)",override.aes = list(linetype = 0)))+
  xlab('Dominance of sweeps')+
  ylab('Fraction genome recombined per mutation')

png('new_fits_sweepy_only/sweep_slope_correlation_nolog_2_commensal_only_2.png', res=300, width=4000, height =2000)
p2c
dev.off()


seg_slope_all_c

 
