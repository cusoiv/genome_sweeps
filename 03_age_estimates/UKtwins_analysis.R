library(dplyr)
library(ggplot2)
library(data.table)
library(pegas)

setwd("C:/Users/Xiaoqian/Desktop/pop_gen/UHGG_plus4/instrain/")

all_sweeps=read.table('sweep_db_all.list')
MG_UKtwins_list=Sys.glob('db_*_MQ1/MG_UKtwins/*.IS/output/*_scaffold_info.tsv')
MG_UKtwins=lapply(MG_UKtwins_list,function(x) fread(x,sep='\t'))
MG_UKtwins_names=sapply(MG_UKtwins_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
names(MG_UKtwins)=MG_UKtwins_names
MG_UKtwins=mapply(function(x,y) {x$sample=y;x},MG_UKtwins,MG_UKtwins_names,SIMPLIFY = F)
MG_UKtwins_all=bind_rows(MG_UKtwins)
single_list=Sys.glob('MG_UKtwins/single_SGB*.txt')


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
MG_UKtwins_sweep_info_list={}
single_MG_UKtwins_n_list={}
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
  
  single_MG_UKtwins_name=single_list[grep(SGB,single_list)]
  
  if (length(single_MG_UKtwins_name)!=0 ){
    single_MG_UKtwins_list=lapply(single_MG_UKtwins_name,function(x) if (file.size(x)!=0L) read.table(x))
    single_MG_UKtwins=bind_rows(single_MG_UKtwins_list)
    single_MG_UKtwins_n_list[[i]]= data.frame(SGB=SGB,MG_UKtwins_n=nrow(single_MG_UKtwins))
    MG_UKtwins_sweep_info_merge=MG_UKtwins_all %>% filter(scaffold==sweep_name) %>% 
      filter(sample %in% single_MG_UKtwins$V1) %>% filter(breadth>0.5) 
    
    if (nrow(MG_UKtwins_sweep_info_merge)>0){
      MG_UKtwins_sweep_info_merge$type='MG'
      MG_UKtwins_sweep_info_merge$SGB=SGB
      MG_UKtwins_sweep_info_merge$sweep=sweep_1
      MG_UKtwins_sweep_info_list[[i]]=MG_UKtwins_sweep_info_merge
    }
  }
}

sweep_info_total=bind_rows(sweep_info_list)
sis_sweep_info_total=bind_rows(sis_sweep_info_list)
MG_UKtwins_sweep_info_total=bind_rows(MG_UKtwins_sweep_info_list)
single_MG_UKtwins_n_total=bind_rows(single_MG_UKtwins_n_list)
isolate_n_total=bind_rows(isolate_list)

all_sweep_info_total=rbind(sweep_info_total,sis_sweep_info_total)
all_sweep_info_total=rbind(all_sweep_info_total,MG_UKtwins_sweep_info_total)
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

MG_UKtwins_metadata=read.csv('MG_UKtwins_metadata.csv') %>% select(sample_id,age,
                                                                   family,zigosity)

all_sweep_info_total_s=left_join(all_sweep_info_total_s,MG_UKtwins_metadata,
                                  by = c('sample'='sample_id'))

all_sweep_info_total_ss=all_sweep_info_total_s %>% group_by(scaffold) %>%
  filter(dist<max_dist)

all_sweep_info_total_ss$family[is.na(all_sweep_info_total_ss$family)]='isolates'

all_sweep_info_total_ss_sf=all_sweep_info_total_ss %>% 
  group_by(scaffold,family) %>% summarise(n=n()) %>% filter(family!='isolates') 

all_sweep_info_total_ss_family=inner_join(all_sweep_info_total_s,all_sweep_info_total_ss_sf,
                                          by=c('scaffold','family')) %>%
  filter(dist<0.5*min_sis_dist)

all_sweep_info_total_ss_family_n2=all_sweep_info_total_ss_family %>%
  group_by(scaffold,family) %>% summarise(n=n()) %>% filter(n>=2)

all_sweep_info_total_ss_family_2=inner_join(all_sweep_info_total_s,all_sweep_info_total_ss_family_n2,
                                          by=c('scaffold','family'))



all_sweep_info_total_ss_si=all_sweep_info_total_ss %>% 
  group_by(scaffold,family) %>% summarise(n=n()) %>% filter(family=='isolates') 

all_sweep_info_total_ss_isolates=inner_join(all_sweep_info_total_ss,all_sweep_info_total_ss_si,
                                          by=c('scaffold','family'))

all_sweep_info_total_ss_isolates=all_sweep_info_total_ss_isolates %>% 
  filter(scaffold %in% all_sweep_info_total_ss_family_2$scaffold)

all_sweep_info_total_ss_all=bind_rows(all_sweep_info_total_ss_isolates,
                                      all_sweep_info_total_ss_family_2)

all_sweep_info_total_ss_family_dist=all_sweep_info_total_ss_family_2 %>%
  group_by(scaffold,family) %>% summarise(family_dist=max(dist)-min(dist))
all_sweep_info_total_ss_family_dist=all_sweep_info_total_ss_family_dist %>%
  group_by(scaffold) %>% summarise(mean_family_dist=mean(family_dist),n_family=n())

all_sweep_info_total_ss_dist=merge(all_sweep_info_total_ss_family_dist,confirmed_sweeps_commensal,by='scaffold')


p=ggplot(all_sweep_info_total_ss_dist,aes(x=mean_dist,y=mean_family_dist,col=scaffold,
                                        size=n_family))+
  geom_point()+geom_abline()+xlim(c(0,0.005))+
  ylim(c(0,0.005))+theme_bw()
  
all_sweep_info_total_ss_dist$mean_dist/all_sweep_info_total_ss_dist$mean_family_dist

png('family_vs_sweep.png',res=300,width=2000,height=1500)
print(p)
dev.off()


ps=ggplot(all_sweep_info_total_ss_all,aes(x=dist,y=scaffold,col=family))+
  geom_point(alpha=0.4,size=2)+
  theme_bw()+scale_x_log10()+theme(legend.position = "none")



all_sweep_info_total_ss_family_2_list=split(all_sweep_info_total_ss_family_2,
                                            f=all_sweep_info_total_ss_family_2$scaffold)

dir.create('UKtwins_MG_consensus_select')
for (i in 1:length(all_sweep_info_total_ss_family_2_list)){
  cs=names(all_sweep_info_total_ss_family_2_list)[i]
  cg=all_sweep_info_total_ss_family_2_list[[i]]$sample
  write.table(cg,paste0('UKtwins_MG_consensus_select/',cs,'.txt'),row.names = F,
              col.names = F,quote = F)
}


