library(dplyr)
library(ggplot2)
library(ggpomological)
library(tidyr)
library(vegan)
library(RColorBrewer)

setwd("C:/Users/Xiaoqian/Desktop/pop_gen/UHGG_plus4/instrain/")
sweep_list=read.csv('confirmed_sweeps_from_consensus.csv')
UKtwins_metadata=read.csv('MG_UKtwins_metadata.csv') %>% select(sample_id,family,age,age_apart)
longt_metadata=read.csv('MG_longt_metadata.csv') %>% select(sample_id,subject_id,days_from_first_collection)
compare_list=list()
compare_list_summary=list()
compare_list_summary_continent=list()
compare_list_summary_info=list()
compare_UKtwins_list=list()
compare_longt_list=list()
select_total_list=list()

metadata_isolate=read.csv('../UHGG_plus4_metadata.csv') %>% select(Strain,Country)
metadata_MG=read.csv('../metagenome_strainphlan_metadata_s2000.csv') %>% select(sample_id,country)
names(metadata_MG)=c("Strain","Country")
metadata_all=rbind(metadata_isolate,metadata_MG)
country_to_continent=read.csv('../country_to_continent.csv',header = T)


for (i in 1:nrow(sweep_list)){
  sf=paste0(sweep_list$x[i],'_consensus')
  sf_table=paste0('compare_select_lists/',sf,'_sweep.txt')
  if (file.exists(sf_table)){
    
    select_total=read.table(sf_table,header = F)
    select_total=left_join(select_total,metadata_all,by=c("V1"="Strain"))
    select_total=left_join(select_total,country_to_continent,by='Country')
    select_total$Continent_2[is.na(select_total$Continent_2)]='Unknown'
    select_total$scaffold=sf
    select_total_list[[i]]=select_total
    compare_name=gsub('_consensus','_compare/instrainComparer_genomeWide_compare.tsv',gsub('diffH_','diffH/',sf))
    compare_name=Sys.glob(paste0('db*_MQ1/',compare_name))
    UKtwins_name=Sys.glob(gsub('_compare/','_UKtwins_compare/',compare_name))
    longt_name=Sys.glob(gsub('_compare/','_longt_compare/',compare_name))
   
    if (file.size(compare_name) > 0){
      compare=read.table(compare_name,header = T,sep='\t')
      compare$name1=gsub('.sorted.bam','',compare$name1)
      compare$name2=gsub('.sorted.bam','',compare$name2)
      compare=compare %>% mutate(group=case_when(name1 %in% select_total$V1 & name2 %in% select_total$V1 ~ 'new_sweep',
                                                 name1 %in% select_total$V1 & !(name2 %in% select_total$V1) ~ 'new_sis_sweep',
                                                 !(name1 %in% select_total$V1) & name2 %in% select_total$V1 ~ 'new_sis_sweep',
                                                 TRUE ~ 'ns'))
      
      compare=compare%>% filter(group!='ns') %>% filter(percent_compared >0.25)
      compare$dist=1-compare$popANI
      compare_s=compare %>% filter(group=='new_sweep')
      compare_ss=compare %>% filter(group=='new_sis_sweep')
      compare_s=left_join(compare_s,metadata_all,by=c("name1"="Strain"))
      compare_s=left_join(compare_s,metadata_all,by=c("name2"="Strain"))
      compare_s=left_join(compare_s,country_to_continent,by=c("Country.x"="Country"))
      compare_s=left_join(compare_s,country_to_continent,by=c("Country.y"="Country"))
      #compare_sc is where we remove all the samples with unknown geography
      compare_sc=compare_s %>% filter(!(is.na(Continent_2.x)|is.na(Continent_2.y)))
      compare_sc=compare_sc %>% mutate(same_continent=case_when(Continent_2.x==Continent_2.y ~ 'same',
                                                              TRUE ~ 'not same'))
      compare_list_summary_info[[i]]=compare_sc

      compare_sc_2a=compare_sc[,c('name1','name2','dist')]
      compare_sc_2b=compare_sc[,c('name2','name1','dist')]
      names(compare_sc_2b)=c('name1','name2','dist')
      compare_sc_2=bind_rows(compare_sc_2a,compare_sc_2b)
      wide_compare_s <- pivot_wider(
        data = compare_sc_2,
        names_from = name2,
        values_from = dist
      )
      d_same=compare_sc_2$dist[compare_sc$same_continent=='same']
      d_diff=compare_sc$dist[compare_sc$same_continent=='not same']
      if (length(d_same)==0){
        d_same=0
      }
      compare_s_mean_by_continent=data.frame(pvalue=1,delta=mean(d_diff)-mean(d_same))
      
      if (nrow(wide_compare_s)>1){
        wide_compare_s=wide_compare_s[,c('name1',wide_compare_s$name1)]
        max_value <- max(wide_compare_s[2:ncol(wide_compare_s)], na.rm = TRUE)  # Find the maximum value, excluding NA
        wide_compare_s[is.na(wide_compare_s)] <- max_value
        wide_compare_s_m=as.matrix(wide_compare_s[2:ncol(wide_compare_s)])
        diag(wide_compare_s_m) <-NA
        wide_compare_s2=left_join(wide_compare_s,metadata_all,by=c("name1"="Strain"))
        wide_compare_s2=left_join(wide_compare_s2,country_to_continent,by=c("Country"))
        wide_compare_s2$Continent_2[is.na(wide_compare_s2$Continent_2)]='Unknown'
        wide_compare_s2c=wide_compare_s2 %>% group_by(Continent_2) %>% summarise(n=n())
        if (length(unique(wide_compare_s2$Continent_2))>1 & max(wide_compare_s2c$n)>=2){
          compare_s_mean_by_continent$pvalue=anosim(wide_compare_s_m,wide_compare_s2$Continent_2)$signif
        }
      }
  
      compare_s_mean_by_continent$scaffold=sf
      compare_s=compare_s %>% mutate(norm_SNP=population_SNPs/coverage_overlap)
      compare_s_mean=compare_s %>% summarise(mean_dist=mean(dist),mean_SNP=mean(norm_SNP),max_dist=max(dist),max_SNP=max(norm_SNP)) 
      compare_ss_1=data.frame(compare_ss)
      compare_ss_1$name1=compare_ss$name2
      compare_ss_1$name2=compare_ss$name1
      compare_ss=bind_rows(compare_ss,compare_ss_1)
      compare_ss_summary= compare_ss %>% group_by(name1) %>% summarise(min_sis_dist=min(dist)) %>% 
        filter(name1 %in% select_total$V1)
      compare_s_mean$min_sis_dist=mean(compare_ss_summary$min_sis_dist)
      compare_s_mean=compare_s_mean %>% mutate(ratio=min_sis_dist/mean_dist)
      compare_s_mean$scaffold=sf
      compare_list[[i]]=compare
      compare_list_summary[[i]]=compare_s_mean
      compare_list_summary_continent[[i]]=compare_s_mean_by_continent
    }
    
    if (length(UKtwins_name)>0){
      if (file.size(UKtwins_name)>0){
        UKtwins=read.table(UKtwins_name,header = T,sep='\t')
        UKtwins$name1=gsub('.sorted.bam','',UKtwins$name1)
        UKtwins$name2=gsub('.sorted.bam','',UKtwins$name2)
        UKtwins=left_join(UKtwins,UKtwins_metadata,by=c("name1"="sample_id"))
        UKtwins=left_join(UKtwins,UKtwins_metadata,by=c("name2"="sample_id"))
        UKtwins=UKtwins %>% filter(family.x==family.y) %>% filter(percent_compared>0.25)
        UKtwins$dist=1-UKtwins$popANI
        UKtwins$SNP=UKtwins$population_SNPs/UKtwins$percent_compared
        compare_UKtwins_list[[i]]=UKtwins
      }
    }
    if (length(longt_name)>0){
      if (file.size(longt_name)>0){
        longt=read.table(longt_name,header = T,sep='\t')
        longt$name1=gsub('.sorted.bam','',longt$name1)
        longt$name2=gsub('.sorted.bam','',longt$name2)
        longt=left_join(longt,longt_metadata,by=c("name1"="sample_id"))
        longt=left_join(longt,longt_metadata,by=c("name2"="sample_id"))
        longt=longt %>% filter(subject_id.x==subject_id.y) %>% filter(percent_compared>0.25)
        longt$years=abs(longt$days_from_first_collection.y-longt$days_from_first_collection.x)/365
        longt$dist=1-longt$popANI
        longt$SNP=longt$population_SNPs/longt$percent_compared
        compare_longt_list[[i]]=longt
      }
    }
  }
}

taxa_list=read.table("../fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice(which.max(n))
taxa_list_summary$species=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',3)
taxa_list_summary$genus=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',2)
taxa_list_summary$family=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',1)
taxa_list_summary=taxa_list_summary %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                TRUE ~ paste0(species)))
taxa_list_summary$cluster_taxa=paste(taxa_list_summary$taxa,taxa_list_summary$cluster,sep = "_")


compare_table=bind_rows(compare_list)
compare_table$scaffold=gsub('.fasta','',compare_table$genome)
compare_table$scaffold=gsub('cf_','',compare_table$scaffold)

compare_table$SGB=sapply(strsplit(compare_table$scaffold,'diffH'),'[[',1)
compare_table=merge(compare_table,taxa_list_summary,by.x = 'SGB',by.y = 'cluster')
compare_table$taxa=paste0(compare_table$taxa,
                          "_",compare_table$scaffold)

compare_summary=bind_rows(compare_list_summary)
compare_summary$SGB=sapply(strsplit(compare_summary$scaffold,'diffH'),'[[',1)
compare_summary=merge(compare_summary,taxa_list_summary,by.x = 'SGB',by.y = 'cluster')
compare_summary$taxa=paste0(compare_summary$taxa,
                          "_",compare_summary$scaffold)
compare_summary_5=compare_summary %>% filter(ratio>=5 | is.na(ratio) & !(is.na(mean_dist))) 
compare_summary_5s=compare_summary_5 %>% select(SGB,scaffold,mean_dist,mean_SNP,max_dist,max_SNP)
#filter out problematic sweep
compare_summary_5s=compare_summary_5s %>% filter(scaffold!='SGB1855diffH_sweep_2_consensus')
#write.csv(compare_summary_5s,'confirmed_sweep_from_pairwise.csv',row.names = F)

compare_summary_5s_s=compare_summary_5s %>% filter(mean_SNP<500)
compare_summary_info=bind_rows(compare_list_summary_info)
compare_summary_info$genome=gsub('cf_','',compare_summary_info$genome)
compare_summary_info$genome=gsub('.fasta','',compare_summary_info$genome)
compare_summary_info_s=merge(compare_summary_info,compare_summary_5s_s,by.x='genome',by.y='scaffold')

compare_summary_5s=read.csv('confirmed_sweep_from_pairwise.csv')
SGB_class=read.csv('../SGBdiffH_500_1pois_1nb_finda/SGB_classify_list_sweepy.csv')

SGB_class_pathogen=sapply(strsplit(SGB_class$pathogen,split = '_'),'[[',5)
SGB_class_probiotics=sapply(strsplit(SGB_class$probiotics[1:11],split = '_'),'[[',5)
compare_summary_info_s=compare_summary_info_s %>% filter(!(SGB %in% c(SGB_class_probiotics,SGB_class_pathogen)))
compare_summary_info_s= compare_summary_info_s %>% filter(same_continent=='not same')
compare_summary_continent=bind_rows(compare_list_summary_continent)

compare_summary_continent$SGB=sapply(strsplit(compare_summary_continent$scaffold,'diffH'),'[[',1)
compare_summary_continent=merge(compare_summary_continent,taxa_list_summary,by.x = 'SGB',by.y = 'cluster')

compare_summary_continent=compare_summary_continent %>% 
  mutate(class=case_when(cluster_taxa %in% SGB_class$pathogen ~ 'pathogens',
                         cluster_taxa %in% SGB_class$probiotics ~ 'probiotics',
                         TRUE ~ 'commensals'))

#Check for sweeps that are not 
#& pvalue>= 0.1
compare_summary_continent_s=compare_summary_continent %>% 
  filter(class=='commensals'& !is.na(delta) & pvalue>= 0.1) %>% 
  filter(scaffold %in% compare_summary_5s$scaffold)

sweep_df_country_count$taxa=as.character(sweep_df_country_count$taxa)
sweep_df_country_count$scaffold=sapply(strsplit(sweep_df_country_count$taxa,split = 'SGB'),function(x) paste0('SGB',x[2]))
sweep_df_country_count$scaffold=gsub('_sweep','diffH_sweep',sweep_df_country_count$scaffold)
sweep_df_country_count$scaffold=gsub('_cfsweep','diffH_cfsweep',sweep_df_country_count$scaffold)
sweep_df_country_count$scaffold=paste0(sweep_df_country_count$scaffold,'_consensus')
sweep_df_country_count=merge(sweep_df_country_count,compare_summary_continent_s,by='scaffold',all.x=T) %>% filter(distinct_continents>1)
sweep_df_country_count_s2=sweep_df_country_count %>% 
  filter(class=='commensals' & pvalue>= 0.1) 


#SGB_sweep=unique(sapply(strsplit(compare_summary_5s$scaffold,split ='diffH'),'[[',1))

#write.table(SGB_sweep,'SGB_sweep.list',row.names = F,col.names = F,quote=F)

compare_UKtwins=bind_rows(compare_UKtwins_list)
compare_UKtwins$scaffold=gsub('.fasta','',compare_UKtwins$genome)
compare_UKtwins$scaffold=gsub('cf_','',compare_UKtwins$scaffold)
compare_UKtwins$SGB=sapply(strsplit(compare_UKtwins$scaffold,split ='diffH'),'[[',1)

compare_UKtwins_summary=compare_UKtwins %>% group_by(scaffold) %>% summarise(mean_family_SNP=mean(dist),se_family_SNP=sd(SNP),n_family=n())
compare_UKtwins_summary=left_join(compare_UKtwins_summary,compare_summary,by='scaffold')
compare_UKtwins_summary$family_ratio=compare_UKtwins_summary$mean_SNP/compare_UKtwins_summary$mean_family_SNP
compare_UKtwins_summary$family_ratio_se=compare_UKtwins_summary$mean_SNP/compare_UKtwins_summary$mean_family_SNP*compare_UKtwins_summary$se_family_SNP

compare_UKtwins_summary_SGB=compare_UKtwins %>% group_by(SGB) %>% summarise(mean_family_SNP=mean(SNP),n_family=n())
compare_UKtwins_summary_SGB=left_join(compare_UKtwins_summary_SGB,compare_summary,by='SGB')
compare_UKtwins_summary_SGB$family_ratio=compare_UKtwins_summary_SGB$mean_SNP/compare_UKtwins_summary_SGB$mean_family_SNP

compare_longt=bind_rows(compare_longt_list)
compare_longt$scaffold=gsub('.fasta','',compare_longt$genome)
compare_longt$scaffold=gsub('cf_','',compare_longt$scaffold)
compare_longt$SGB=sapply(strsplit(compare_longt$scaffold,split ='diffH'),'[[',1)

compare_longt_summary=compare_longt %>% group_by(scaffold) %>% summarise(mean_person_SNP=mean(SNP),n_person=n())
compare_longt_summary=left_join(compare_longt_summary,compare_summary,by='scaffold')
compare_longt_summary$person_ratio=compare_longt_summary$mean_SNP/compare_longt_summary$mean_person_SNP

compare_longt_summary_SGB=compare_longt %>% group_by(SGB) %>% summarise(mean_person_SNP=mean(SNP),n_person=n())
compare_longt_summary_SGB=left_join(compare_longt_summary_SGB,compare_summary,by='SGB')
compare_longt_summary_SGB$person_ratio=compare_longt_summary_SGB$mean_SNP/compare_longt_summary_SGB$mean_person_SNP

compare_longt_UKtwins_summary=merge(compare_longt_summary,compare_UKtwins_summary,by='scaffold',all=T)
compare_table_s=compare_table %>% filter(scaffold %in% compare_summary_5$scaffold)
compare_table_fs=compare_table %>% filter(!(scaffold %in% compare_summary_5$scaffold))

compare_longt_UKtwins_summary_SGB=merge(compare_longt_summary_SGB,compare_UKtwins_summary_SGB,by='scaffold',all=T)


#Plot for scaffolds that can be found in both longt and twin datasets
compare_longt_UKtwins_summary_2=merge(compare_longt_summary,compare_UKtwins_summary,by='scaffold') %>% 
  filter(ratio.x>=5)
compare_UKtwins$age_apart.x[is.na(compare_UKtwins$age_apart.x)]=18
compare_UKtwins$years=compare_UKtwins$age.x-compare_UKtwins$age_apart.x
compare_UKtwins_select=compare_UKtwins %>% select(scaffold,SNP,years)
compare_UKtwins_select$dataset='UKtwins'
compare_longt_select=compare_longt %>% select(scaffold,SNP,years)
names(compare_longt_select)=c('scaffold','SNP','years')
compare_longt_select$dataset='longt'
compare_UKtwins_longt_select=bind_rows(compare_UKtwins_select,compare_longt_select)

compare_UKtwins_longt_select_2=compare_UKtwins_longt_select %>%
  filter(scaffold %in% compare_longt_UKtwins_summary_2$scaffold)


pg=ggplot()+
  geom_point(data=compare_UKtwins_longt_select_2,aes(x=SNP,y=years,col=dataset))+
  geom_smooth(data=compare_UKtwins_longt_select_2,aes(x=SNP,y=years),method = "lm", se = FALSE, color = "blue")+
  facet_wrap(~scaffold)

#fit to get tmrca for sweep
fit_group1 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1814diffH_cfsweep_1_consensus", ])
fit_group2 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1860diffH_sweep_1_consensus", ])
fit_group3 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB2295diffH_cfsweep_1_consensus", ])

known_x_group1 <- compare_longt_UKtwins_summary_2$max_SNP.x[compare_longt_UKtwins_summary_2$scaffold=='SGB1814diffH_cfsweep_1_consensus']
known_x_group2 <- compare_longt_UKtwins_summary_2$max_SNP.x[compare_longt_UKtwins_summary_2$scaffold=='SGB1860diffH_sweep_1_consensus']
known_x_group3 <- compare_longt_UKtwins_summary_2$max_SNP.x[compare_longt_UKtwins_summary_2$scaffold=='SGB2295diffH_cfsweep_1_consensus']


predicted_values_1 <- predict(fit_group1, newdata = data.frame(SNP = known_x_group1), se.fit = TRUE)
predicted_values_2 <- predict(fit_group2, newdata = data.frame(SNP = known_x_group2), se.fit = TRUE)
predicted_values_3 <- predict(fit_group3, newdata = data.frame(SNP = known_x_group3), se.fit = TRUE)

fit_params_group1 <- coef(summary(lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1814diffH_cfsweep_1_consensus", ])))
fit_params_group2 <- coef(summary(lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1860diffH_sweep_1_consensus", ])))
fit_params_group3 <- coef(summary(lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB2295diffH_cfsweep_1_consensus", ])))



known_y_group1 <- fit_params_group1["(Intercept)", "Estimate"] + fit_params_group1["SNP", "Estimate"] * known_x_group1
known_y_group2 <- fit_params_group2["(Intercept)", "Estimate"] + fit_params_group2["SNP", "Estimate"] * known_x_group2
known_y_group3 <- fit_params_group3["(Intercept)", "Estimate"] + fit_params_group3["SNP", "Estimate"] * known_x_group3

compare_UKtwins_longt_both_fit=data.frame(scaffold=c("SGB1814diffH_cfsweep_1_consensus",
                                                     "SGB1860diffH_sweep_1_consensus","SGB2295diffH_cfsweep_1_consensus"),
                                          mean_estimated_time=c(predicted_values_1$fit,predicted_values_2$fit,predicted_values_3$fit),
                                          se_estimated_time=c(predicted_values_1$se.fit,predicted_values_2$se.fit,predicted_values_3$se.fit),
                                          dataset='both')

pg <- pg +
  geom_point(data = data.frame(SNP = known_x_group1, years = known_y_group1, scaffold = "SGB1814diffH_cfsweep_1_consensus",dataset='average sweep'), 
             aes(x=SNP,y=years,col=dataset)) +
  geom_point(data = data.frame(SNP = known_x_group2, years = known_y_group2, scaffold = "SGB1860diffH_sweep_1_consensus",dataset='average sweep'), 
             aes(x=SNP,y=years,col=dataset))+theme_bw()

#png(paste0('sweep_summary/','tMRCA_longt_UKtwins_both_fit_SNP.png'),width=2000,height=1500,res=300)
#print(pg)
#dev.off()


#plot only for UKtwins datasets
compare_UKtwins_summary_select=compare_UKtwins_summary%>% filter(ratio>=5) %>% 
  filter(!(scaffold %in% compare_UKtwins_longt_select_2$scaffold))
compare_UKtwins_select=compare_UKtwins %>% filter(scaffold %in% compare_UKtwins_summary_select$scaffold)
compare_UKtwins_select=merge(compare_UKtwins_select,compare_summary_5s,by='scaffold')
compare_UKtwins_select=compare_UKtwins_select %>% mutate(estimated_time=max_SNP/SNP*age.x) %>% group_by(scaffold) %>% 
  summarise(mean_estimated_time=mean(estimated_time),se_estimated_time=sd(estimated_time))

compare_longt_summary_select=compare_longt_summary%>% filter(ratio>=5) %>% 
  filter(!(scaffold %in% compare_UKtwins_longt_select_2$scaffold))
compare_longt_select=compare_longt %>% filter(scaffold %in% compare_longt_summary_select$scaffold)
compare_longt_select=merge(compare_longt_select,compare_summary_5s,by='scaffold')
compare_longt_select$SNP[compare_longt_select$SNP==0]=2
compare_longt_select=compare_longt_select %>% mutate(estimated_time=max_SNP/SNP*years) %>% group_by(scaffold) %>% 
  summarise(mean_estimated_time=mean(estimated_time),se_estimated_time=sd(estimated_time))

compare_longt_select$dataset='time series'
compare_UKtwins_select$dataset='UKtwins'
compare_UKtwins_longt_select=bind_rows(compare_longt_select,compare_UKtwins_select)


ps=ggplot(compare_UKtwins_longt_select,aes(x=mean_estimated_time,y=scaffold))+
  geom_bar(stat = 'identity')+
  theme_bw()+xlab('estimated tMRCA of sweep')+ylab('sweep')

#png(paste0('sweep_summary/','tMRCA_longt_UKtwins_only_SNP.png'),width=2000,height=2000,res=300)
#print(ps)
#dev.off()

compare_UKtwins_longt_all=bind_rows(compare_longt_select,compare_UKtwins_select,compare_UKtwins_longt_both_fit)
compare_UKtwins_longt_all$SGB=sapply(strsplit(compare_UKtwins_longt_all$scaffold,'diffH'),'[[',1)
compare_UKtwins_longt_all=left_join(compare_UKtwins_longt_all,taxa_list_summary,by=c("SGB"="cluster"))
compare_UKtwins_longt_all$sweep=sapply(strsplit(compare_UKtwins_longt_all$scaffold,'_'),function(x) paste(x[2], x[3], sep = "_"))
compare_UKtwins_longt_all$sweep=paste(compare_UKtwins_longt_all$cluster_taxa,compare_UKtwins_longt_all$sweep,sep='_')

ps2=ggplot(compare_UKtwins_longt_all,aes(x=mean_estimated_time,y=reorder(sweep,-mean_estimated_time),fill=dataset))+
  geom_bar(stat = 'identity',alpha=0.8)+
  scale_fill_brewer(palette="Dark2")+
  geom_errorbar(aes(xmin = mean_estimated_time - se_estimated_time, xmax = mean_estimated_time + se_estimated_time,color = dataset), width = 0.2)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+xlab('estimated tMRCA of sweep')+ylab('sweep')+
  theme(axis.title = element_text(size = 20),  # Font size for axis titles
axis.text = element_text(size = 20),  # Font size for axis labels
legend.text = element_text(size = 20),  # Font size for legend labels
legend.title = element_text(size = 24))+
  scale_x_log10()

#png(paste0('sweep_summary/','tMRCA_longt_UKtwins_all_SNP.png'),width=3500,height=2000,res=300)
#print(ps2)
#dev.off()


#if we do the same analysis as before but for SGBs
#Plot for scaffolds that can be found in both longt and twin datasets
compare_longt_UKtwins_summary_SGB_2=merge(compare_longt_summary_SGB,compare_UKtwins_summary_SGB,by='scaffold') %>% 
  filter(ratio.x>=5)
#fit to get tmrca for sweep
#The SGBs are the same as the sweep by sweep dataset

fit_group1 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1814diffH_cfsweep_1_consensus", ])
fit_group2 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB1860diffH_sweep_1_consensus", ])
fit_group3 <- lm(years ~ SNP, data = compare_UKtwins_longt_select_2[compare_UKtwins_longt_select_2$scaffold == "SGB2295diffH_cfsweep_1_consensus", ])

known_x_SGB_group1 <- compare_longt_UKtwins_summary_SGB_2$max_SNP.x[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Phocaeicola_vulgatus_SGB1814']
known_x_SGB_group2 <- compare_longt_UKtwins_summary_SGB_2$max_SNP.x[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Bacteroides_faecis_SGB1860']
known_x_SGB_group3 <- compare_longt_UKtwins_summary_SGB_2$max_SNP.x[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Alistipes_shahii_SGB2295']

known_SGB_group1 <- compare_longt_UKtwins_summary_SGB_2$scaffold[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Phocaeicola_vulgatus_SGB1814']
known_SGB_group2 <- compare_longt_UKtwins_summary_SGB_2$scaffold[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Bacteroides_faecis_SGB1860']
known_SGB_group3 <- compare_longt_UKtwins_summary_SGB_2$scaffold[compare_longt_UKtwins_summary_SGB_2$cluster_taxa.x=='s__Alistipes_shahii_SGB2295']

predicted_values_SGB_1 <- predict(fit_group1, newdata = data.frame(SNP = known_x_SGB_group1), se.fit = TRUE)
predicted_values_SGB_2 <- predict(fit_group2, newdata = data.frame(SNP = known_x_SGB_group2), se.fit = TRUE)
predicted_values_SGB_3 <- predict(fit_group3, newdata = data.frame(SNP = known_x_SGB_group3), se.fit = TRUE)


compare_UKtwins_longt_both_SGB_fit=data.frame(scaffold=c(known_SGB_group1,known_SGB_group2,known_SGB_group3),
                                                    
                                          mean_estimated_time=c(predicted_values_SGB_1$fit,predicted_values_SGB_2$fit,predicted_values_SGB_3$fit),
                                          se_estimated_time=c(predicted_values_SGB_1$se.fit,predicted_values_SGB_2$se.fit,predicted_values_SGB_3$se.fit),
                                          dataset='both')


#plot only for UKtwins/longt datasets
compare_summary_5s$SGB=sapply(strsplit(compare_summary_5s$scaffold,'diffH'),'[[',1)
compare_UKtwins_summary_select=compare_UKtwins_summary%>% filter(ratio>=5) %>% 
  filter(!(scaffold %in% compare_UKtwins_longt_select_2$scaffold))
compare_UKtwins_select=compare_UKtwins %>% filter(scaffold %in% compare_UKtwins_summary_select$scaffold)


compare_UKtwins_SGB_select=merge(compare_UKtwins_select,compare_summary_5s,by='SGB')
compare_UKtwins_SGB_select=compare_UKtwins_SGB_select %>% mutate(estimated_time=max_SNP/SNP*age.x) %>% group_by(scaffold.y) %>% 
  summarise(mean_estimated_time=mean(estimated_time),se_estimated_time=sd(estimated_time))
compare_longt_summary_select=compare_longt_summary%>% filter(ratio>=5) %>% 
  filter(!(scaffold %in% compare_UKtwins_longt_select_2$scaffold))
compare_longt_select=compare_longt %>% filter(scaffold %in% compare_longt_summary_select$scaffold)
compare_longt_SGB_select=merge(compare_longt_select,compare_summary_5s,by='SGB')
compare_longt_SGB_select$SNP[compare_longt_SGB_select$SNP==0]=2
compare_longt_SGB_select=compare_longt_SGB_select %>% mutate(estimated_time=max_SNP/SNP*years) %>% group_by(scaffold.y) %>% 
  summarise(mean_estimated_time=mean(estimated_time),se_estimated_time=sd(estimated_time))

compare_longt_SGB_select$dataset='time series'
compare_UKtwins_SGB_select$dataset='UKtwins'
compare_UKtwins_longt_SGB_select=bind_rows(compare_longt_SGB_select,compare_UKtwins_SGB_select)
names(compare_UKtwins_longt_SGB_select)[1]='scaffold'

compare_UKtwins_longt_SGB_all=bind_rows(compare_UKtwins_longt_SGB_select,compare_UKtwins_longt_both_SGB_fit)
compare_UKtwins_longt_SGB_all$SGB=sapply(strsplit(compare_UKtwins_longt_SGB_all$scaffold,'diffH'),'[[',1)
compare_UKtwins_longt_SGB_all=left_join(compare_UKtwins_longt_SGB_all,taxa_list_summary,by=c("SGB"="cluster"))
compare_UKtwins_longt_SGB_all$sweep=sapply(strsplit(compare_UKtwins_longt_SGB_all$scaffold,'_'),function(x) paste(x[2], x[3], sep = "_"))
compare_UKtwins_longt_SGB_all$sweep=paste(compare_UKtwins_longt_SGB_all$cluster_taxa,compare_UKtwins_longt_SGB_all$sweep,sep='_')

ps2=ggplot(compare_UKtwins_longt_SGB_all,aes(x=mean_estimated_time,y=reorder(sweep,-mean_estimated_time),fill=dataset))+
  geom_bar(stat = 'identity',alpha=0.8)+
  scale_fill_brewer(palette="Dark2")+
  geom_errorbar(aes(xmin = mean_estimated_time - se_estimated_time, xmax = mean_estimated_time + se_estimated_time,color = dataset), width = 0.2)+
  scale_color_brewer(palette="Dark2")+
 
  theme_bw()+
  theme(axis.title = element_text(size = 20),  # Font size for axis titles
        axis.text = element_text(size = 20),  # Font size for axis labels
        legend.text = element_text(size = 20),  # Font size for legend labels
        legend.title = element_text(size = 24))+
  xlab('estimated tMRCA of sweep')+ylab('sweep')+
  scale_x_log10()

#png(paste0('sweep_summary/','tMRCA_longt_UKtwins_by_SGB_all_SNP.png'),width=3600,height=2000,res=300)
#print(ps2)
#dev.off()

#compare metagenome data to sweep SNPs
compare_UKtwins_longt_SGB_all=left_join(compare_UKtwins_longt_SGB_all,
                                        compare_summary_5,by='scaffold')
compare_UKtwins_longt_SGB_all=compare_UKtwins_longt_SGB_all %>% mutate(SNP5=max_SNP/5/2,SNP10=max_SNP/10/2)

#compare_UKtwins_longt_SGB_all_2=compare_UKtwins_longt_SGB_all %>% mutate(max_SNP=max_SNP/10)
#compare_UKtwins_longt_SGB_all$mc='1 SNP/year'
#compare_UKtwins_longt_SGB_all_2$mc='10 SNPs/year'
#compare_UKtwins_longt_SGB_all=bind_rows(compare_UKtwins_longt_SGB_all,compare_UKtwins_longt_SGB_all_2)


compare_UKtwins_longt_SGB_all$cluster_taxa.x=gsub('s__','',compare_UKtwins_longt_SGB_all$cluster_taxa.x)
compare_UKtwins_longt_SGB_all$cluster_taxa.x=gsub('_',' ',compare_UKtwins_longt_SGB_all$cluster_taxa.x)
colors <- colorRampPalette(brewer.pal(name="Dark2",n=8))(length(unique(compare_UKtwins_longt_SGB_all$cluster_taxa.x)))


pcompare=ggplot(compare_UKtwins_longt_SGB_all,aes(x=mean_estimated_time,y=SNP5,color=cluster_taxa.x))+
#  geom_line(aes(group=scaffold),linewidth=2,alpha=0.6)+
  geom_errorbar(aes(ymin=SNP10,ymax=max_SNP/2,color = cluster_taxa.x), width = 0,linewidth=1.5,alpha=0.6)+
  geom_errorbar(aes(xmin = mean_estimated_time - se_estimated_time, xmax = mean_estimated_time + se_estimated_time,color = cluster_taxa.x), width = 0,linewidth=1.5,alpha=0.6)+
  geom_point(size=6)+
  theme_bw()+
  scale_color_manual(values = colors, labels = function(x) {
                         sapply(x, function(label) {
                           parts <- strsplit(label, "SGB")[[1]]
                           parts[[2]]=paste0("(SGB",parts[[2]],')')
                           italicized_parts <- lapply(parts, function(part) bquote(italic(.(part))))
                           bquote(.(italicized_parts[[1]])~.(parts[[2]]))
                         })
                       })+
  theme(axis.title = element_text(size = 24),  # Font size for axis titles
        axis.text = element_text(size = 24),  # Font size for axis labels
        legend.text = element_text(size = 24),  # Font size for legend labels
        legend.title = element_text(size = 24))+
  xlim(c(0,25000))+
   ylim(c(0,25000))+
  labs(color='Taxa')+
 scale_x_log10(limits=c(1,25000))+scale_y_log10(limits=c(1,25000))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey",linewidth=2)+
  geom_abline(intercept = -1, slope = 1 , linetype = "dashed", color = "grey50")+
  geom_abline(intercept = 1, slope = 1 , linetype = "dashed", color = "grey50")+
  xlab('Estimated sweep age\n(Time series or twins data)')+
  ylab('Estimated sweep age\n(Constant molecular clock)')+
  guides(color = guide_legend(byrow = TRUE))
    
png('sweep_summary/tMRCA_longt_UKtwins_by_SGB_compare_log_2_e_SNP.png',res=300,width=4500,height=2500)
pcompare
dev.off()

#Add geographical info to the sweep_consensus_5
select_total_all=bind_rows(select_total_list)
select_total_all$SGB=sapply(strsplit(select_total_all$scaffold,'diffH'),'[[',1)
select_total_all=merge(select_total_all,taxa_list_summary,by.x = 'SGB',by.y = 'cluster')
select_total_all$taxa=paste0(select_total_all$taxa,
                            "_",select_total_all$scaffold)

select_total_sweep_5=select_total_all %>% filter(scaffold %in% compare_summary_5$scaffold)
select_total_sweep_5_summary=select_total_sweep_5 %>% group_by(taxa,Continent) %>%
  summarise(nc=n())

#load tree for all clonal frames
library(ape)
library(castor)
library(phangorn)
library(ggtree)
library(RColorBrewer)

ct=read_tree(file='gtdbtk.unrooted.itol.tree')
ct$tip.label=gsub('_cf_consensus','',ct$tip.label)

##all plot panels for commensals
#tree
compare_summary_5$scaffold=gsub("_consensus","",compare_summary_5$scaffold)
cts=keep.tip(ct,compare_summary_5$scaffold)
cts=midpoint(cts)
newtip=inner_join(data.frame(tip.label=cts$tip.label),compare_summary_5,by=c('tip.label'='scaffold'))
cts$tip.label=newtip$taxa
ct_tipseq=get_taxa_name(ggtree(cts))

p=ggtree(cts) + xlim_expand(c(0, 4), 'Tree')
sweep_df_family=compare_summary_5 %>% select(taxa,family)
sweep_df_family$taxa=as.character(sweep_df_family$taxa)
sweep_df_family$new_label=label_pad(sweep_df_family$taxa)

colourCount = length(unique(sweep_df_family$family))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

p=p %<+% sweep_df_family +
  geom_tiplab(aes(label=new_label,color=family),align = T,
              linetype = NULL,size=4,fontface='bold',family='mono')+
  scale_color_manual(values = getPalette(colourCount)) +
  theme(legend.position = "bottom",legend.title=element_text(size=18,face = "bold"),
        legend.text=element_text(size=14),
        legend.background = element_rect(size = 0.5),
        axis.text.x=element_blank(),axis.title.x=element_blank(),
        plot.margin = unit(c(1,-4,1,0), "cm"))+
  guides(color=guide_legend(ncol = 2, byrow=TRUE,title.position = "top",title = "Family",
                            override.aes = list(label = "\u25A0", size = 7)))

p_continent=ggplot(select_total_sweep_5_summary,aes(x=Continent,y=reorder(taxa,desc(taxa)),
                                          col=Continent,size=nc))+
  geom_point()+theme_classic()+
  theme(legend.position = "bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=14),
        axis.title.x=element_blank(),
        legend.text=element_text(size=18),legend.title=element_text(size=18,face = "bold"),
        legend.background = element_rect(size = 0.5))+
  guides(color=guide_legend(ncol = 2, byrow=TRUE,title.position = "top",
                            override.aes = list(size = 5),title='Geography'),
         size=guide_legend(ncol = 1, byrow=TRUE,title.position = "top",title='No. of samples'))+
  scale_color_pomological()



sweep_df_list=split(sweep_df,f=sweep_df$taxa)
sweep_df_p=list()
for (i in 1:length(sweep_df_list)){
  sweep_tab=sweep_df_list[[i]] %>% count(Continent,group) %>% spread(group,n,fill = 0)
  # sweep_tab = sweep_tab %>% mutate(Country = case_when(new_sweep==0 ~ 'others', TRUE ~ Country))
  sweep_tab= sweep_tab %>% group_by(Continent) %>% summarise(new_sis_sweep=sum(new_sis_sweep),
                                                             new_sweep=sum(new_sweep)) %>% 
    select(-Continent)
  if (nrow(sweep_tab)>1){
    sweep_df_p[[i]]=fisher.test(sweep_tab,workspace = 4000000)$p.value
  } else {
    sweep_df_p[[i]]=1
  }
}


names(sweep_df_p)=names(sweep_df_list)
sweep_df_p=as.data.frame(unlist(sweep_df_p))
sweep_df_p$taxa=as.factor(rownames(sweep_df_p))
names(sweep_df_p)=c('pvalue','taxa')
sweep_df_p=sweep_df_p %>% mutate(local=case_when(pvalue<=0.05 ~ TRUE,
                                                 pvalue>0.05 ~ FALSE))
sweep_df_p$taxa=factor(sweep_df_p$taxa,levels=ct_tipseq)
p_local=ggplot(sweep_df_p,aes(x='local',y=reorder(taxa,desc(taxa)),fill=local))+
  geom_tile()+theme_classic()+
  scale_fill_manual(values = c('grey80','burlywood1'))+
  theme(legend.position = "bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.line.y = element_blank(),axis.text.x=element_text(size=18,angle = 90),axis.title.x=element_blank(),
        legend.text=element_text(size=18),legend.title=element_text(size=18,face = "bold"),
        legend.direction='vertical',
        legend.background = element_rect(size = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))+
  labs(fill="Local")+
  coord_cartesian(expand = FALSE)
p_local_leg=get_legend(p_local)
p_continent_leg=get_legend(p_continent)

ggplot(compare_table_s,aes(x=SNP,y=scaffold,col=group))+geom_point()+
  scale_x_log10()


