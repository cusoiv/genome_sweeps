library(dplyr)
library(ggplot2)
#library(micropan)
library(data.table)
library(pegas,lib='~/Rlibs_4')

setwd("/lisc/scratch/dome/xy43/UHGG_plus4/Results/instrain")

all_sweeps=read.table('sweep_db_all.list')
MG_list=Sys.glob('db_*_MQ1/MG*/*.IS/output/*_scaffold_info.tsv')
MG_gene_list=Sys.glob('db_*_MQ1/MG*/*.IS/output/*_gene_info.tsv')
MG_SNP_list=Sys.glob('db_*_MQ1/MG*/*.IS/output/*IS_SNVs.tsv')
MG=lapply(MG_list,function(x) fread(x,sep='\t'))
MG_names=sapply(MG_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
names(MG)=MG_names
MG=mapply(function(x,y) {x$sample=y;x},MG,MG_names,SIMPLIFY = F)
MG_all=bind_rows(MG)
single_list=Sys.glob('MG/single_SGB*.txt')
isolate_gene_list=Sys.glob(paste0('db_*_MQ1/','SGB*/*/*.IS/output/*_gene_info.tsv'))
isolate_MG_gene_list=c(isolate_gene_list,MG_gene_list)
isolate_SNP_list=Sys.glob(paste0('db_*_MQ1/','SGB*/*/*.IS/output/*IS_SNVs.tsv'))

taxa_list=read.table("fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice(which.max(n))
taxa_list_summary$species=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',3)
taxa_list_summary$genus=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',2)
taxa_list_summary=taxa_list_summary %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                            TRUE ~ paste0(species)))
taxa_list_summary$cluster_taxa=paste(taxa_list_summary$taxa,taxa_list_summary$cluster,sep = "_")

sweep_info_list={}
sis_sweep_info_list={}
sweep_gene_info_list={}
sis_sweep_gene_info_list={}
MG_sweep_info_list={}
#MG_sweep_gene_info_list={}
single_MG_n_list={}
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
  sweep_gene_list=Sys.glob(paste0('db_*_MQ1/',sweep,'/*.IS/output/*_gene_info.tsv'))
  sweep_gene_info=lapply(sweep_gene_list,function(x) read.table(x,sep='\t',header=T))
  genome_name=sapply(sweep_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
  
  sweep_info=mapply(function(x,y) {x$sample=y;x},sweep_info,genome_name,SIMPLIFY = F)
  sweep_gene_info=mapply(function(x,y) {x$sample=y;x},sweep_gene_info,genome_name,SIMPLIFY = F)
  sweep_name=paste0(gsub('/','_',sweep),'_consensus')
  sweep_info=lapply(sweep_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.5) %>% filter(coverage>=2)) 
  sweep_gene_info=lapply(sweep_gene_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.5) %>% filter(dNdS_substitutions>2) %>% filter(coverage>=2)) 
  sweep_info_merge=bind_rows(sweep_info)
  sweep_gene_info_merge=bind_rows(sweep_gene_info)
  sweep_gene_info_merge_summary=sweep_gene_info_merge %>% group_by(gene) %>%
    summarise(n=n())
  sweep_gene_info_merge_summary_count=sweep_gene_info_merge_summary %>% filter(n>=length(sweep_gene_info)/2)
  
  if (nrow(sweep_info_merge)>0){
    sweep_info_merge$type='sweep'
    sweep_info_merge$SGB=SGB
    sweep_info_merge$sweep=sweep_1
    sweep_info_list[[i]]=sweep_info_merge
    sweep_gene_info_merge_summary_count$type='sweep'
    sweep_gene_info_merge_summary_count$SGB=SGB
    sweep_gene_info_merge_summary_count$sweep=sweep_1
    sweep_gene_info_list[[i]]=sweep_gene_info_merge_summary_count
  }
  
  
  sis_sweep_list=Sys.glob(paste0('db_*_MQ1/',sis_sweep,'/*.IS/output/*_scaffold_info.tsv'))
  sis_sweep_info=lapply(sis_sweep_list,function(x) read.table(x,sep='\t',header=T))
  genome_name=sapply(sis_sweep_list,function(x) gsub('.IS_scaffold_info.tsv','',basename(x)))
  sis_sweep_gene_list=Sys.glob(paste0('db_*_MQ1/',sis_sweep,'/*.IS/output/*_gene_info.tsv'))
  sis_sweep_gene_info=lapply(sis_sweep_gene_list,function(x) read.table(x,sep='\t',header=T))
  sis_sweep_gene_info=mapply(function(x,y) {x$sample=y;x},sis_sweep_gene_info,genome_name,SIMPLIFY = F)
  sis_sweep_info=mapply(function(x,y) {x$sample=y;x},sis_sweep_info,genome_name,SIMPLIFY = F)
  sis_sweep_name=paste0(gsub('/','_',sis_sweep),'_consensus')
  
  sis_sweep_info=lapply(sis_sweep_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.25) %>% filter(coverage>=2)) #keep low breadth for sister sweep since it is needed as anchor 
  sis_sweep_info_merge=bind_rows(sis_sweep_info)
  sis_sweep_gene_info=lapply(sis_sweep_gene_info,function(x) x%>% filter(scaffold==sweep_name) %>% filter(breadth>0.5) %>% filter(dNdS_substitutions>2) %>% filter(coverage>=2)) 
  
  sis_sweep_gene_info_merge=bind_rows(sis_sweep_gene_info)
  sis_sweep_gene_info_merge_summary=sis_sweep_gene_info_merge %>% group_by(gene) %>%
    summarise(n=n())
  sis_sweep_gene_info_merge_summary_count=sis_sweep_gene_info_merge_summary %>% filter(n>=length(sis_sweep_gene_info)/2)
  
  if (nrow(sis_sweep_info_merge)>0){
    sis_sweep_info_merge$type='sis_sweep'
    sis_sweep_info_merge$SGB=SGB
    sis_sweep_info_merge$sweep=sweep_1
    sis_sweep_info_list[[i]]=sis_sweep_info_merge
    sis_sweep_gene_info_merge_summary_count$type='sis_sweep'
    sis_sweep_gene_info_merge_summary_count$SGB=SGB
    sis_sweep_gene_info_merge_summary_count$sweep=sweep_1
    sis_sweep_gene_info_list[[i]]=sis_sweep_gene_info_merge_summary_count
  }
 
  single_MG_name=single_list[grep(SGB,single_list)]

  if (length(single_MG_name)!=0 ){
    single_MG_list=lapply(single_MG_name,function(x) if (file.size(x)!=0L) read.table(x))
    single_MG=bind_rows(single_MG_list)
    single_MG_n_list[[i]]= data.frame(SGB=SGB,MG_n=nrow(single_MG))
   MG_sweep_info_merge=MG_all %>% filter(scaffold==sweep_name) %>% 
        filter(sample %in% single_MG$V1) %>% filter(breadth>0.25) %>% filter(coverage>=2)
   
      if (nrow(MG_sweep_info_merge)>0){
        MG_sweep_info_merge$type='MG'
        MG_sweep_info_merge$SGB=SGB
        MG_sweep_info_merge$sweep=sweep_1
        MG_sweep_info_list[[i]]=MG_sweep_info_merge
      }
    }
  }

sweep_info_total=bind_rows(sweep_info_list)
sis_sweep_info_total=bind_rows(sis_sweep_info_list)
sweep_gene_info_total=bind_rows(sweep_gene_info_list)
sis_sweep_gene_info_total=bind_rows(sis_sweep_gene_info_list)
sweep_gene_info_total$scaffold=paste0(sweep_gene_info_total$SGB,'diffH_',sweep_gene_info_total$sweep)
sis_sweep_gene_info_total$scaffold=paste0(sis_sweep_gene_info_total$SGB,'diffH_',sis_sweep_gene_info_total$sweep)


MG_sweep_info_total=bind_rows(MG_sweep_info_list)
single_MG_n_total=bind_rows(single_MG_n_list)
isolate_n_total=bind_rows(isolate_list)
sweep_gene_pos_total=sweep_gene_info_total %>% group_by(scaffold) %>% summarise(n=n())
sis_sweep_gene_pos_total=sis_sweep_gene_info_total %>% group_by(scaffold) %>% summarise(n=n())

all_sweep_info_total=rbind(sweep_info_total,sis_sweep_info_total)
all_sweep_info_total=rbind(all_sweep_info_total,MG_sweep_info_total)
all_sweep_info_total$type=as.factor(all_sweep_info_total$type)
all_sweep_info_total$dist=1-all_sweep_info_total$popANI_reference
all_sweep_info_total$dist[all_sweep_info_total$dist==0]=10^-6
all_sweep_info_total_list=split(all_sweep_info_total,f=all_sweep_info_total$SGB)
all_sweep_info_total_list_sweep=split(all_sweep_info_total,f=all_sweep_info_total$scaffold)

ratio_ref=function(s){
  r=rep(1,length(s))
  if (length(s)>=4){
    for (i in 4:length(s)){
      r[i]=s[i]/mean(s[1:(i-1)])
    }
  }
  return(r)
}

all_sweep_info_total_list_sweep_list=list()

for (i in 1:length(all_sweep_info_total_list_sweep)){
       sf=names(all_sweep_info_total_list_sweep)[i]
       all_sweep_info=all_sweep_info_total_list_sweep[[i]]
       all_sweep_info=all_sweep_info %>% distinct(sample,.keep_all = TRUE)
       all_sweep_info=all_sweep_info %>% arrange(dist)
       all_sweep_info$ratio=ratio_ref(all_sweep_info$dist)
       boundary=min(which(all_sweep_info$ratio>=5))
       all_sweep_info$group='group'
       all_sweep_info$bd=FALSE
      if (is.finite(boundary)){
      
        all_sweep_info$group[1:(boundary-1)]='new_sweep'
        all_sweep_info$group[boundary:nrow(all_sweep_info)]='new_sis_sweep'
        if (sum(all_sweep_info$group=='new_sis_sweep')<=2){
          all_sweep_info$group='no_sweep'
        } else {
          all_sweep_info$bd[boundary]=TRUE
          if (all_sweep_info$dist[all_sweep_info$bd==TRUE]>1.5*mean(all_sweep_info$dist[all_sweep_info$type=='sis_sweep']) & sum(all_sweep_info$group=='new_sweep')/nrow(all_sweep_info)>=0.67){
            all_sweep_info$group='no_sweep'
            all_sweep_info$bd[boundary]=FALSE
          } 
	}
      } else {
        all_sweep_info$group='no_sweep'
      }

       all_sweep_info_total_list_sweep_list[[i]]=all_sweep_info
}

names(all_sweep_info_total_list_sweep_list)=names(all_sweep_info_total_list_sweep)

all_sweep_info_total_list_sweep_new=bind_rows(all_sweep_info_total_list_sweep_list)


all_sweep_info_total_list_sweep_new_sweep=all_sweep_info_total_list_sweep_new %>% filter(group=='new_sweep') %>% filter(breadth>0.5 & coverage>2)
all_sweep_info_total_list_sweep_new_no_sweep=all_sweep_info_total_list_sweep_new %>% filter(group!='new_sweep')
all_sweep_info_total_list_sweep_new_sweep=all_sweep_info_total_list_sweep_new_sweep %>% group_by(SGB,sample) %>%
  filter(breadth==max(breadth))
all_sweep_info_total=bind_rows(all_sweep_info_total_list_sweep_new_sweep,all_sweep_info_total_list_sweep_new_no_sweep)
all_sweep_info_total_list=split(all_sweep_info_total,f=all_sweep_info_total$SGB)
all_sweep_info_total_list_sweep=split(all_sweep_info_total,f=all_sweep_info_total$scaffold)

all_sweep_info_total_list_sweep_list=list()
new_sis_sweep_gene_info_list={}
new_sweep_SNV_consensus_tajima_list={}
for (i in 1:length(all_sweep_info_total_list_sweep)){
       sf=names(all_sweep_info_total_list_sweep)[i]
       all_sweep_info=all_sweep_info_total_list_sweep[[i]]
       all_sweep_info=all_sweep_info %>% distinct(sample,.keep_all = TRUE)
       all_sweep_info=all_sweep_info %>% arrange(dist)
       all_sweep_info$ratio=ratio_ref(all_sweep_info$dist)
       boundary=min(which(all_sweep_info$ratio>=5))
       all_sweep_info$group='group'
       all_sweep_info$bd=FALSE
      if (is.finite(boundary)){
        all_sweep_info$group[1:(boundary-1)]='new_sweep'
        all_sweep_info$group[boundary:nrow(all_sweep_info)]='new_sis_sweep'
        if (sum(all_sweep_info$group=='new_sis_sweep')<=2){
          all_sweep_info$group='no_sweep'
        } else {
          all_sweep_info$bd[boundary]=TRUE
          if (all_sweep_info$dist[all_sweep_info$bd==TRUE]>1.5*mean(all_sweep_info$dist[all_sweep_info$type=='sis_sweep']) & sum(all_sweep_info$group=='new_sweep')/nrow(all_sweep_info)>=0.67){
            all_sweep_info$group='no_sweep'
            all_sweep_info$bd[boundary]=FALSE
          }
        }
        new_sis_sweep_list=all_sweep_info$sample[boundary:(boundary+2)]
        new_sis_sweep_list=new_sis_sweep_list[!is.na(new_sis_sweep_list)]
        print (new_sis_sweep_list)
	sn=unique(all_sweep_info$sweep)
        new_sweep_list=all_sweep_info$sample[1:(boundary-1)]
        new_sweep_SNV_list_iso=lapply(new_sweep_list,function(x) isolate_SNP_list[grepl(paste0(".*", sn, ".*", x , ".*"),isolate_SNP_list)])
        new_sweep_SNV_list_MG=lapply(new_sweep_list,function(x) MG_SNP_list[grep(paste0(x,'.IS_SNVs'),MG_SNP_list)])
        new_sweep_SNV_list_loc=mapply(function(x,y) c(x,y),new_sweep_SNV_list_iso,new_sweep_SNV_list_MG)
        new_sweep_SNV_info=lapply(new_sweep_SNV_list_loc,function(x) bind_rows(lapply(x,function(y) fread(y,sep='\t'))))
        new_sweep_SNV_info=mapply(function(x,y) {x$sample=y;x},new_sweep_SNV_info,new_sweep_list,SIMPLIFY = F)
        new_sweep_SNV_info=lapply(new_sweep_SNV_info,function(x) x%>% filter(scaffold==sf))
        new_sweep_SNV_info_merge=bind_rows(new_sweep_SNV_info)
        new_sweep_SNV_info_merge=new_sweep_SNV_info_merge %>% distinct(sample,position,.keep_all = TRUE)
        new_sweep_SNV_info_merge_position=split(new_sweep_SNV_info_merge,f=new_sweep_SNV_info_merge$position)
        
        new_sweep_SNV_consensus=lapply(new_sweep_list, function(x) paste0(unlist(sapply(new_sweep_SNV_info_merge_position, 
                                                                          function(y) {
                                                                            if (x %in% y$sample) {
                                                                              return(y$con_base[which(x==y$sample)])
                                                                            } else {
                                                                              return(y$ref_base[1])
                                                                            }
                                                                          })))
          )
       
        new_sweep_SNV_consensus_bin=as.DNAbin(new_sweep_SNV_consensus)
        new_sweep_SNV_consensus_tajima=as.data.frame(tajima.test(new_sweep_SNV_consensus_bin))
        new_sweep_SNV_consensus_tajima$scaffold=gsub('_consensus','',sf)
        new_sweep_SNV_consensus_tajima_list[[i]]=new_sweep_SNV_consensus_tajima
        
        new_sis_sweep_list_iso=lapply(new_sis_sweep_list,function(x) isolate_gene_list[grepl(paste0(".*", sn, ".*", x , ".*"),isolate_gene_list)])
        new_sis_sweep_list_MG=lapply(new_sis_sweep_list,function(x) MG_gene_list[grep(paste0(x,'.IS_gene_info.tsv'),MG_gene_list)])
        new_sis_sweep_list_loc=mapply(function(x,y) c(x,y),new_sis_sweep_list_MG,new_sis_sweep_list_iso)
        print(new_sis_sweep_list_iso)
	print (new_sis_sweep_list_MG)
        new_sis_sweep_gene_info=lapply(new_sis_sweep_list_loc,function(x) bind_rows(lapply(x,function(y) fread(y,sep='\t'))))
        new_sis_sweep_gene_info=mapply(function(x,y) {x$sample=y;x},new_sis_sweep_gene_info,new_sis_sweep_list,SIMPLIFY = F)
        new_sis_sweep_gene_info=lapply(new_sis_sweep_gene_info,function(x) x%>% filter(scaffold==sf) %>% filter(breadth>0.5) %>% filter(dNdS_substitutions>2) %>% filter(coverage>=2)) 
        
# 	print (new_sis_sweep_gene_info)
	       
        new_sis_sweep_gene_info_merge=bind_rows(new_sis_sweep_gene_info)
        new_sis_sweep_gene_info_merge_summary=new_sis_sweep_gene_info_merge %>% group_by(gene) %>%
          summarise(n=n())
	print (new_sis_sweep_gene_info_merge_summary)
        print (length(new_sis_sweep_gene_info))
        new_sis_sweep_gene_info_merge_summary_count=new_sis_sweep_gene_info_merge_summary #%>% filter(n>=length(new_sis_sweep_gene_info)/2)
        new_sis_sweep_gene_info_merge_summary_count$scaffold=gsub('_consensus','',sf)
        new_sis_sweep_gene_info_list[[i]]=new_sis_sweep_gene_info_merge_summary_count
        
      } else {
        all_sweep_info$group='no_sweep'
      }
       all_sweep_info_total_list_sweep_list[[i]]=all_sweep_info
}
names(all_sweep_info_total_list_sweep_list)=names(all_sweep_info_total_list_sweep)
saveRDS(all_sweep_info_total_list_sweep_list, file="instrain_analysis_rds_redo_06142025/all_sweep_info_total_list_sweep_list.rds")
saveRDS(new_sis_sweep_gene_info_list, file="instrain_analysis_rds_redo_06142025/new_sis_sweep_gene_info_list.rds")
saveRDS(new_sweep_SNV_consensus_tajima_list, file="instrain_analysis_rds_redo_06142025/new_sweep_SNV_consensus_tajima_list.rds")
saveRDS(single_MG_n_total, file="instrain_analysis_rds_redo_06142025/single_MG_n_list.rds")
saveRDS(isolate_n_total,file="instrain_analysis_rds_redo_06142025/isolate_list.rds")
