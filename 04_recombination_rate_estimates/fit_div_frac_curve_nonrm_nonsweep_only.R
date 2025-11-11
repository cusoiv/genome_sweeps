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
library(patchwork)
library(svglite)

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
#filter list for all taxa with no sweeps
sweep_taxa=read.csv('../instrain/confirmed_sweep_from_pairwise.csv')
taxa_class=read.csv('../SGB_classify_list.csv')

lbfile_list_summary_n20_filter_1=lbfile_list_summary_n20_filter %>% 
  filter(!cluster %in% sweep_taxa$SGB)
lbfile_list_summary_n20_filter=lbfile_list_summary_n20_filter %>% 
  filter(cluster %in% sweep_taxa$SGB)

lbfile_list_summary_n20_filter_full=lbfile_list_summary_n20 %>% 
  filter(!(cluster %in% sweep_taxa$SGB))
lbfile_list_summary_n20_filter_full_1=lbfile_list_summary_n20 %>% 
  filter((cluster %in% sweep_taxa$SGB))
lbfile_list_summary_n20_filter_full_1_n=lbfile_list_summary_n20_filter_full_1%>% 
 group_by(cluster_taxa) %>% summarise(n=n()) %>% filter(!(cluster_taxa %in% c(taxa_class$pathogen,taxa_class$probiotics)))

lbfile_list_summary_n20_filter_full_n=lbfile_list_summary_n20_filter_full%>% 
  group_by(cluster_taxa) %>% summarise(n=n())%>% filter(!(cluster_taxa %in% c(taxa_class$pathogen,taxa_class$probiotics)))


lbfile_list_summary_n20_filter_1_n=lbfile_list_summary_n20_filter_1 %>% group_by(cluster_taxa) %>%
  summarize(n=n()) %>% filter(!cluster_taxa %in% c(taxa_class$pathogen,taxa_class$probiotics))
lbfile_list_summary_n20_filter_n=lbfile_list_summary_n20_filter %>% group_by(cluster_taxa) %>%
  summarize(n=n()) %>% filter(!cluster_taxa %in% c(taxa_class$pathogen,taxa_class$probiotics))
t.test(lbfile_list_summary_n20_filter_1_n$n,lbfile_list_summary_n20_filter_n$n)

lbfile_list_summary_n20_filter_list=split(lbfile_list_summary_n20_filter_1,f=lbfile_list_summary_n20_filter_1$cluster_taxa)

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
#png('new_fits_nonsweepy_only/segs_fit_r1s.png',res=300,width = 8000,height=5000)
#par(mfrow=c(8,8))
for (j in 1:length(segs_list)){
  if (segs_list[[j]]$segments$r2[1]>0.8 & segs_list[[j]]$segments$slope[1]>0){
    r2=round(segs_list[[j]]$segments$r2[1],2)
    plot(segs_list[[j]],main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    r1_list=c(r1_list,j)
  }else {
    redo_list=c(redo_list,j)
  }
}
#dev.off()

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
#png('new_fits_nonsweepy_only/segs_fit_rs_p1s.png',res=300,width = 8000,height=5000)
#par(mfrow=c(10,8))
for (j in redo_list){
  r2=round(segs_list_2[[j]]$segments$r2[1],2)
  if (r2>0.33 & segs_list_2[[j]]$segments$slope[1]>0){
    plot(segs_list_2[[j]],main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2))
    r2_list=c(r2_list,j)
  }
}
#dev.off()

#check ones that are eliminated
r12_list=c(r1_list,r2_list)
new_redo_list=c()
#png('new_fits_nonsweepy_only/r2_omit_summary_2s.png',width=8000,height=6000,res=300)
#par(mfrow=c(10,7))
for (j in 1:length(lbfile_list_summary_n20_filter_list)){
  l=lbfile_list_summary_n20_filter_list[[j]]
  l=l %>% arrange(div)
  test=data.frame(x=l$div,y=l$totalR)
  test_o=test

    if (!(j %in% r12_list)){
      print(j)
      plot(test_o,pch=20,col='grey60',
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],"_",j),xlab="",ylab="",xlim=c(0,5000),ylim=c(0,1))
      new_redo_list=c(new_redo_list,j)
    }
}
#dev.off()

fit_lists=c(r1_list,r2_list)
all_lists=1:length(lbfile_list_summary_n20_filter_list)
unfit_lists=all_lists[!(all_lists %in% fit_lists)]
unfit_lists_select=c(9,10,18,30,35,43,46,67,68,78,85,86,91,94)

cutoff_list=c(2000,1000,1000,800,3000,500,500,500,500,500,500,5000,1000,600)

#c=1
#for (c in 1:length(unfit_lists_select) ){
#  j=unfit_lists_select[[c]]
#  l=lbfile_list_summary_n20_filter_list[[j]]
#  l=l %>% arrange(div)
#  test=data.frame(x=l$div,y=l$totalR)
#  cutoff=cutoff_list[c]
#  test_o=test[test$x<cutoff,]
#  test_o=rbind(c(0,0),test_o)
#  ft=lm(y~0+x,test_o)
#  r2=round(summary(ft)$r.squared,2)
#  plot(test,pch=20,col='grey60',
#       main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',r2),xlab="",ylab="",xlim=c(0,0.001*4000000),ylim=c(0,1))
#  points(test_o,col='red',pch=20)
#  abline(ft)
#  c=c+1
#}

seg_slope_list={}
c=1
d=1
#png('new_fits_nonsweepy_only/r2_good_fits_summary.png',width=8000,height=6000,res=300)
#par(mfrow=c(10,8))
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
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',j),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
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
           main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',j),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
      points(test1,col='red',pch=20)
      abline(ft)
      seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                     slope=ft$coefficients[1])
      c=c+1
    }
  } else if  (j %in% unfit_lists_select) {
    test_o=test
    cutoff=cutoff_list[d]
    test1=test[test$x<cutoff,]
    test1=rbind(c(0,0),test1)
    ft=lm(y~0+x,test1)
    r2=round(summary(ft)$r.squared,2)
    plot(test_o,pch=20,col='grey60',
         main=paste0(names(lbfile_list_summary_n20_filter_list)[[j]],'_',j),xlab="",ylab="",xlim=c(0,0.004*4000000),ylim=c(0,1))
    points(test1,col='red',pch=20)
    abline(ft)
    seg_slope_list[[j]]=data.frame(taxa=names(lbfile_list_summary_n20_filter_list)[[j]],
                                   slope=ft$coefficients[1])
    c=c+1
    d=d+1
  }
}
#dev.off()

badfit_lists=c(13,50,60)
cutoff_list_2=c(1000,1000,2000)
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
#  abline(ft)
#  c=c+1
#}

r1_list=r1_list[!(r1_list %in% badfit_lists)]
r2_list=r2_list[!(r2_list %in% badfit_lists)]

taxa_class_include=names(lbfile_list_summary_n20_filter_list)[!(names(lbfile_list_summary_n20_filter_list) %in% c(taxa_class$pathogen,taxa_class$probiotics))]
taxa_class_total=unique(lbfile_list_summary_n20_filter_full$cluster_taxa)[!(unique(lbfile_list_summary_n20_filter_full$cluster_taxa) %in% c(taxa_class$pathogen,taxa_class$probiotics))]
seg_slope_list={}
c=1
d=1
e=1

#Figure S1a
png('new_fits_nonsweepy_only/r2_good_fits_summary_2a_commensals.png',width=5500,height=5500/10*13,res=300)
par(mfrow=c(13,4),mai = c(0.3, 0.3, 0.3, 0.1))
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
  } else if  (j %in% unfit_lists_select) {
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

seg_slope_all_nonsweepy=bind_rows(seg_slope_list) 
write.csv(seg_slope_all_nonsweepy,'seg_slope_all_nonrm_nonsweepy.csv',row.names = F)

## Compare the recombination rate between sweeps with 

seg_slope_all_nonsweepy=read.csv('seg_slope_all_nonrm_nonsweepy.csv')
seg_slope_all=read.csv('seg_slope_all_nonrm.csv')
seg_slope_all_nonsweepy_commensals=seg_slope_all_nonsweepy %>% 
  filter(!(taxa %in% c(taxa_class$pathogen,taxa_class$probiotics)))
seg_slope_all_sweepy_commensals=seg_slope_all %>% 
  filter(!(taxa %in% c(taxa_class$pathogen,taxa_class$probiotics)))

shapiro.test(log10(seg_slope_all_sweepy_commensals$slope))
shapiro.test(log10(seg_slope_all_nonsweepy_commensals$slope))
t.test(log10(seg_slope_all_sweepy_commensals$slope),
       log10(seg_slope_all_nonsweepy_commensals$slope))

wilcox.test(seg_slope_all_sweepy_commensals$slope,
            seg_slope_all_nonsweepy_commensals$slope)


seg_slope_all_nonsweepy_commensals$Type='SGBs without sweep'
seg_slope_all_sweepy_commensals$Type='SGBs with sweep'
seg_slope_all_commensals=bind_rows(seg_slope_all_nonsweepy_commensals,
                                   seg_slope_all_sweepy_commensals)
seg_slope_all_commensals=seg_slope_all_commensals %>%
  mutate(Genus=case_when(grepl('Blautia',taxa)~'Blautia',TRUE ~ 'Others'))
seg_slope_all_commensals$Type=factor(seg_slope_all_commensals$Type,levels=c('SGBs without sweep','SGBs with sweep'))

seg_slope_all_commensals_median=seg_slope_all_commensals %>% group_by(Type) %>%
  summarise(median_slope=median(slope))

p=ggplot(seg_slope_all_commensals,aes(x=log10(slope),fill=Type))+
  geom_histogram(alpha=0.6, position="identity")+
  facet_wrap(Type ~.,nrow=2)+
  geom_vline(data=seg_slope_all_commensals_median, aes(xintercept=log10(median_slope), color=Type),
             linetype="dashed")+
  scale_fill_manual(values=c("grey80","#cc747c"))+
  scale_color_manual(values=c("grey60","#cc747c"))+
  scale_y_continuous(breaks = seq(0, 10, by = 5))+
  xlab('')+
  ylab('# of SGBs')+
  theme_minimal()+
  theme(text = element_text(size = 16),legend.title = element_blank(),strip.text = element_blank())+
  theme(
    strip.placement = "outside",  # places labels outside axes
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x =element_line(color = "black"),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 20))
    )+
  xlim(c(-5.5,-2))

seg_slope_all_blautia=seg_slope_all_commensals %>% filter(Genus=='Blautia')
seg_slope_all_blautia=seg_slope_all_blautia %>% mutate(class=case_when(taxa=='s__Blautia_massiliensis_SGB4826' ~ 'Blautia massiliensis\n(SGB4826)', TRUE ~ 'Other Blautia'))

pb=ggplot(seg_slope_all_blautia,aes(fill=factor(Type),color=factor(class),x=log10(slope),y='1'))+
  geom_point(size=5,shape = 21,stroke = 2)+
  scale_fill_manual(values=c("grey80","#cc747c"))+
  theme_minimal()+scale_color_manual(values=c("darkred","transparent"),labels = c(bquote(atop(italic("Blautia massiliensis"), "(SGB4826)"))))+
  xlab('log10(Fraction genome recombined per mutation)')+
  xlim(c(-5.5,-2))+
  theme(axis.ticks.y = element_blank(),axis.title.y = element_blank(),axis.text =element_blank(),legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 16),
        legend.text = element_text(size = 16),
        # Make horizontal grid lines (y-axis) black
        panel.grid.major.y = element_line(color = "black"),
        panel.grid.minor.y = element_line(color = "black"))+
  guides(fill = "none")

#Figure 3
png('new_fits_nonsweepy_only/compare_sweep_nonsweep_Blautia.png', res=300, width=2500, height =1200)
wrap_plots(p, pb, ncol = 1) +
  plot_layout(heights = c(3, 1)) 
dev.off()

combined_plots=wrap_plots(p, pb, ncol = 1) +
  plot_layout(heights = c(3, 1)) 
ggsave('new_fits_nonsweepy_only/compare_sweep_nonsweep_Blautia.svg',plot=combined_plots,width=2500, height =1200,units='px')




seg_slope_all_commensals_export=seg_slope_all_commensals[,1:3]
seg_slope_all_commensals_export$Type=gsub('SGBs\n','Commensal ',
                                          seg_slope_all_commensals_export$Type)
seg_slope_all_commensals_export$taxa=gsub('s__','',seg_slope_all_commensals_export$taxa)
seg_slope_all_commensals_export$taxa=gsub('_',' ',seg_slope_all_commensals_export$taxa)  

write.csv(seg_slope_all_commensals_export,'seg_slope_all_commensals_export.csv',row.names = F)

taxa_list_commensals=unique(lbfile_list_summary_n20$cluster_taxa)[!(unique(lbfile_list_summary_n20$cluster_taxa) %in% 
                                                                      c(taxa_class$pathogen,taxa_class$probiotics))]
taxa_list_pathogenes=unique(lbfile_list_summary_n20$cluster_taxa)[(unique(lbfile_list_summary_n20$cluster_taxa) %in% 
                                                                      c(taxa_class$pathogen))]

taxa_list_probiotics=unique(lbfile_list_summary_n20$cluster_taxa)[(unique(lbfile_list_summary_n20$cluster_taxa) %in% 
                                                                     c(taxa_class$probiotics))]


taxa_list_commensals_df=data.frame(taxa=taxa_list_commensals,type='commensals')
taxa_list_pathogens_df=data.frame(taxa=taxa_list_pathogenes,type='pathogen')
taxa_list_probiotics_df=data.frame(taxa=taxa_list_probiotics,type='probiotics')

taxa_list_all_export=bind_rows(taxa_list_commensals_df,taxa_list_pathogens_df,taxa_list_probiotics_df)
taxa_list_all_export$taxa=gsub('s__','',taxa_list_all_export$taxa)
taxa_list_all_export$taxa=gsub('_',' ',taxa_list_all_export$taxa)
write.csv(taxa_list_all_export,'../taxa_list_all_export.csv',row.names = F)
