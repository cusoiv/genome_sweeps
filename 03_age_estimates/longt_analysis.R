library(dplyr)
library(ggplot2)
library(data.table)
library(pegas)

setwd("C:/Users/Xiaoqian/Desktop/pop_gen/UHGG_plus4/instrain/")

all_sweeps=read.table('sweep_db_all.list')
MG_longt_list=Sys.glob('db_*_MQ1/MG_longt/*.IS/output/*_scaffold_info.tsv')
MG_longt=lapply(MG_longt_list,function(x) fread(x,sep='\t'))
MG_longt_names=sapply(MG_longt_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
names(MG_longt)=MG_longt_names
MG_longt=mapply(function(x,y) {x$sample=y;x},MG_longt,MG_longt_names,SIMPLIFY = F)
MG_longt_all=bind_rows(MG_longt)
single_list=Sys.glob('MG_longt/single_SGB*.txt')

taxa_list=read.table("../fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice(which.max(n))
taxa_list_summary$species=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',3)
taxa_list_summary$genus=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',2)
taxa_list_summary=taxa_list_summary %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                TRUE ~ paste0(species)))
taxa_list_summary$cluster_taxa=paste(taxa_list_summary$taxa,taxa_list_summary$cluster,sep = "_")

sweep_info_list={}
sis_sweep_info_list={}
MG_longt_sweep_info_list={}
single_MG_longt_n_list={}
isolate_list={}

for (i in 1:nrow(all_sweeps)){
  
  sweep=all_sweeps$V1[i]
  sweep_1=basename(sweep)
  SGB=strsplit(sweep,split='diffH')[[1]][1]
  
  isolate=read.table(paste0('SGBdiffH/',SGB,'diffH.txt'),header=F)
  isolate_list[[i]]=data.frame(SGB=SGB,isolate_n=nrow(isolate))
  sis_sweep=paste0(all_sweeps$V1[i],'_sister')
  sweep_list=Sys.glob(paste0('db_*_MQ1/',sweep,'/*.IS/output/*_scaffold_info.tsv'))
  sweep_info=lapply(sweep_list,function(x) read.table(x,sep='\t',header=T))
  
  
  genome_name=sapply(sweep_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
  
  sweep_info=mapply(function(x,y) {x$sample=y;x},sweep_info,genome_name,SIMPLIFY = F)
  
  sweep_name=paste0(gsub('/','_',sweep),'_consensus')
  sweep_info=lapply(sweep_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.1)) 
  
  sweep_info_merge=bind_rows(sweep_info)
  
  if (nrow(sweep_info_merge)>0){
    sweep_info_merge$type='sweep'
    sweep_info_merge$SGB=SGB
    sweep_info_merge$sweep=sweep_1
    sweep_info_list[[i]]=sweep_info_merge
  }
  
  
  sis_sweep_list=Sys.glob(paste0('db_*_MQ1/',sis_sweep,'/*.IS/output/*_scaffold_info.tsv'))
  sis_sweep_info=lapply(sis_sweep_list,function(x) read.table(x,sep='\t',header=T))
  genome_name=sapply(sis_sweep_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
  
  sis_sweep_info=mapply(function(x,y) {x$sample=y;x},sis_sweep_info,genome_name,SIMPLIFY = F)
  sis_sweep_name=paste0(gsub('/','_',sis_sweep),'_consensus')
  
  sis_sweep_info=lapply(sis_sweep_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.1)) 
  sis_sweep_info_merge=bind_rows(sis_sweep_info)
  
  if (nrow(sis_sweep_info_merge)>0){
    sis_sweep_info_merge$type='sis_sweep'
    sis_sweep_info_merge$SGB=SGB
    sis_sweep_info_merge$sweep=sweep_1
    sis_sweep_info_list[[i]]=sis_sweep_info_merge
    
  }
  
  single_MG_longt_name=single_list[grep(SGB,single_list)]
  
  if (length(single_MG_longt_name)!=0 ){
    single_MG_longt_list=lapply(single_MG_longt_name,function(x) if (file.size(x)!=0L) read.table(x))
    single_MG_longt=bind_rows(single_MG_longt_list)
    single_MG_longt_n_list[[i]]= data.frame(SGB=SGB,MG_longt_n=nrow(single_MG_longt))
    MG_longt_sweep_info_merge=MG_longt_all %>% filter(scaffold==sweep_name) %>% 
      filter(sample %in% single_MG_longt$V1) %>% filter(breadth>0.5) 
    
    if (nrow(MG_longt_sweep_info_merge)>0){
      MG_longt_sweep_info_merge$type='MG'
      MG_longt_sweep_info_merge$SGB=SGB
      MG_longt_sweep_info_merge$sweep=sweep_1
      MG_longt_sweep_info_list[[i]]=MG_longt_sweep_info_merge
    }
  }
}

sweep_info_total=bind_rows(sweep_info_list)
sis_sweep_info_total=bind_rows(sis_sweep_info_list)
MG_longt_sweep_info_total=bind_rows(MG_longt_sweep_info_list)
single_MG_longt_n_total=bind_rows(single_MG_longt_n_list)
isolate_n_total=bind_rows(isolate_list)


all_sweep_info_total=rbind(sweep_info_total,sis_sweep_info_total)
all_sweep_info_total=rbind(all_sweep_info_total,MG_longt_sweep_info_total)
all_sweep_info_total$type=as.factor(all_sweep_info_total$type)
all_sweep_info_total$dist=1-all_sweep_info_total$popANI_reference
all_sweep_info_total$dist[all_sweep_info_total$dist==0]=10^-6
all_sweep_info_total_list=split(all_sweep_info_total,f=all_sweep_info_total$SGB)
all_sweep_info_total_list_sweep=split(all_sweep_info_total,f=all_sweep_info_total$scaffold)

confirmed_sweeps_commensal=read.csv('confirmed_sweeps_commensal.csv')
all_sweep_info_total_s=all_sweep_info_total %>%
  filter(scaffold %in% confirmed_sweeps_commensal$scaffold)
all_sweep_info_total_s=left_join(all_sweep_info_total_s,confirmed_sweeps_commensal,
                                 by='scaffold')

MG_longt_metadata=read.csv('MG_longt_metadata.csv') %>% select(sample_id,subject_id,
                                                               days_from_first_collection)


all_sweep_info_total_s=inner_join(all_sweep_info_total_s,MG_longt_metadata,
                                 by = c('sample'='sample_id'))

all_sweep_info_total_ss=all_sweep_info_total_s %>% group_by(scaffold) %>%
  filter(dist<max_dist)

all_sweep_info_total_ss_sf=all_sweep_info_total_ss %>% 
  group_by(scaffold,subject_id) %>% summarise(n=n())

all_sweep_info_total_ss_subject_id=inner_join(all_sweep_info_total_s,all_sweep_info_total_ss_sf,
                                          by=c('scaffold','subject_id')) %>%
  filter(dist<0.5*min_sis_dist)

all_sweep_info_total_ss_subject_id_n2=all_sweep_info_total_ss_subject_id %>%
  group_by(scaffold,subject_id) %>% summarise(n=n()) %>% filter(n>=2)

all_sweep_info_total_ss_subject_id_2=inner_join(all_sweep_info_total_s,all_sweep_info_total_ss_subject_id_n2,
                                            by=c('scaffold','subject_id'))


all_sweep_info_total_ss_subject_id_dist=all_sweep_info_total_ss_subject_id_2 %>%
  group_by(scaffold,subject_id) %>% summarise(subject_id_dist=max(dist)-min(dist))
all_sweep_info_total_ss_subject_id_dist=all_sweep_info_total_ss_subject_id_dist %>%
  group_by(scaffold) %>% summarise(mean_subject_id_dist=mean(subject_id_dist),n_subject_id=n())

all_sweep_info_total_ss_dist=merge(all_sweep_info_total_ss_subject_id_dist,confirmed_sweeps_commensal,by='scaffold')
all_sweep_info_total_ss_dist$mean_subject_id_dist[all_sweep_info_total_ss_dist$mean_subject_id_dist==0]=10^-7

p=ggplot(all_sweep_info_total_ss_dist,aes(x=mean_dist,y=mean_subject_id_dist,col=scaffold,
                                          size=n_subject_id))+
  geom_point()+geom_abline()+
  xlim(c(0,0.005))+
  ylim(c(0,0.005))+theme_bw()


png('person_vs_sweep.png',res=300,width=2000,height=1600)
print(p)
dev.off()


all_sweep_info_total_ss_subject_id_2_list=split(all_sweep_info_total_ss_subject_id_2,
                                            f=all_sweep_info_total_ss_subject_id_2$scaffold)

dir.create('longt_MG_consensus_select')
for (i in 1:length(all_sweep_info_total_ss_subject_id_2_list)){
  cs=names(all_sweep_info_total_ss_subject_id_2_list)[i]
  cg=all_sweep_info_total_ss_subject_id_2_list[[i]]$sample
  write.table(cg,paste0('longt_MG_consensus_select/',cs,'.txt'),row.names = F,
              col.names = F,quote = F)
}



