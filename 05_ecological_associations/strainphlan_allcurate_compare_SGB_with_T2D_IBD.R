library(phangorn)
library(ggtree)
library(ggplot2)
library(castor)
library(RColorBrewer)
library(dplyr)
library(ggnewscale)
library(reshape2)
library(purrr)
library(ComplexHeatmap)
library(tidyr)
library(SignifReg)
library(tibble)
library(ggtext)
library(forcats)

setwd("C:/Users/Xiaoqian/Desktop/pop_gen/UHGG_plus4/strainphlan/allcurated_with_IBD_T2D_meta4_beta1_vJan21/")

taxa_list=read.table("../../fastANI_clust.average.diffH.94s.fgspecies.txt",
                     sep='\t',header=T)
taxa_list_summary=taxa_list %>% group_by(cluster) %>% count(fgspecies) %>% slice_max(order_by=n,with_ties = F)
taxa_list_summary$species=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',3)
taxa_list_summary$genus=sapply(strsplit(taxa_list_summary$fgspecies,split = ';'),'[[',2)
taxa_list_summary=taxa_list_summary %>% mutate(taxa = case_when(species=='s__' ~ paste0(genus,'_',species),
                                                                TRUE ~ paste0(species)))
taxa_list_summary$cluster_taxa=paste(taxa_list_summary$taxa,
                                     taxa_list_summary$cluster,sep = "_")
taxa_list_summary=taxa_list_summary %>% arrange(genus)


SGB_exist_list=Sys.glob("allcurated_sams/*")

SGBstrain_list=list()

for (i in 1: length(SGB_exist_list)){
  
  SGB=gsub('allcurated_sams/','',SGB_exist_list[[i]])
  SGB=gsub('_exist.txt','',SGB)
  SGB=gsub('_group','',SGB)
  SGB=gsub('10068','10068s',SGB)
  taxa=taxa_list_summary$cluster_taxa[taxa_list_summary$cluster==SGB]
  if (file.size(SGB_exist_list[[i]]) > 0) {
    SGBstrain=read.table(SGB_exist_list[[i]])
    SGBstrain$V1=gsub('_marker','', SGBstrain$V1)
    SGBstrain_list[[i]]=SGBstrain
    names(SGBstrain_list)[i]=taxa
  }
}

MG_metadata=read.csv("../../metagenome_strainphlan_metadata.csv")
SGBstrain_list=Filter(Negate(is.null), SGBstrain_list)
SGBstrain_list_merge=bind_rows(SGBstrain_list)
SGBstrain_list=Filter(function(x) nrow(x) > 1, SGBstrain_list)

baseline_data_age=read.csv('baseline_data_age.csv')
baseline_data_CRC=read.csv('baseline_data_CRC.csv')
baseline_data_IBD=read.csv('baseline_data_IBD.csv')
baseline_data_UC=read.csv('baseline_data_UC.csv')
baseline_data_CD=read.csv('baseline_data_CD.csv')
baseline_data_T2D=read.csv('baseline_data_T2D.csv')

#for age dataset
SGBstrain_merged_list_1=lapply(SGBstrain_list,function(x)
  merge(baseline_data_age,x,by.x='sample_id',by.y='V1',all.x = T) %>% select(sample_id,subject_id,disease,age,age_category,country,family,non_westernized)
  %>% mutate(isSGB=case_when(sample_id %in% x$V1 ~'SGB', TRUE ~ 'notSGB')))

#this step makes family either family or subject id, so then if we slice for family,
#we end up slicing for both subject and family
SGBstrain_merged_list_1=lapply(SGBstrain_merged_list_1,function(x)
  x %>% 
    group_by(family) %>% slice(1) %>% ungroup())

#for CRC

SGBstrain_merged_list_2=lapply(SGBstrain_list,function(x)
  merge(baseline_data_CRC,x,by.x='sample_id',by.y='V1',all.x = T) %>% select(sample_id,subject_id,disease,age,age_category,country,family,non_westernized)
  %>% mutate(isSGB=case_when(sample_id %in% x$V1 ~'SGB', TRUE ~ 'notSGB')))

#this step makes family either family or subject id, so then if we slice for family,
#we end up slicing for both subject and family
SGBstrain_merged_list_2=lapply(SGBstrain_merged_list_2,function(x)
  x %>% 
    group_by(family) %>% slice(1) %>% ungroup())

#for CD

SGBstrain_merged_list_3a=lapply(SGBstrain_list,function(x)
  merge(baseline_data_CD,x,by.x='sample_id',by.y='V1',all.x = T) %>% select(sample_id,subject_id,disease,age,age_category,country,family,non_westernized)
  %>% mutate(isSGB=case_when(sample_id %in% x$V1 ~'SGB', TRUE ~ 'notSGB')))

#this step makes family either family or subject id, so then if we slice for family,
#we end up slicing for both subject and family
SGBstrain_merged_list_3a=lapply(SGBstrain_merged_list_3a,function(x)
  x %>% 
    group_by(family) %>% slice(1) %>% ungroup())

#for UC

SGBstrain_merged_list_3b=lapply(SGBstrain_list,function(x)
  merge(baseline_data_UC,x,by.x='sample_id',by.y='V1',all.x = T) %>% select(sample_id,subject_id,disease,age,age_category,country,family,non_westernized)
  %>% mutate(isSGB=case_when(sample_id %in% x$V1 ~'SGB', TRUE ~ 'notSGB')))

#this step makes family either family or subject id, so then if we slice for family,
#we end up slicing for both subject and family
SGBstrain_merged_list_3b=lapply(SGBstrain_merged_list_3b,function(x)
  x %>% 
    group_by(family) %>% slice(1) %>% ungroup())

#for T2D

SGBstrain_merged_list_4=lapply(SGBstrain_list,function(x)
  merge(baseline_data_T2D,x,by.x='sample_id',by.y='V1',all.x = T) %>% select(sample_id,subject_id,disease,age,age_category,country,family,non_westernized)
  %>% mutate(isSGB=case_when(sample_id %in% x$V1 ~'SGB', TRUE ~ 'notSGB')))

#this step makes family either family or subject id, so then if we slice for family,
#we end up slicing for both subject and family
SGBstrain_merged_list_4=lapply(SGBstrain_merged_list_4,function(x)
  x %>% 
    group_by(family) %>% slice(1) %>% ungroup())


#only select for commensal SGBs
SGB_class=read.csv('../../SGBdiffH_500_1pois_1nb_finda/SGB_classify_list_sweepy.csv')
commensal_select=names(SGBstrain_merged_list_1)[!(names(SGBstrain_merged_list_1) %in% c(SGB_class$pathogen,SGB_class$probiotics))]

#filter for commensals
SGBstrain_merged_list_1=SGBstrain_merged_list_1[commensal_select]

new_isSGB_list_unmerged_1s=lapply(SGBstrain_merged_list_1, function(x) x%>% select(age_category, isSGB))
new_isSGB_list_unmerged_1s=do.call(cbind,new_isSGB_list_unmerged_1s)
isSGB_names=grepl('isSGB',names(new_isSGB_list_unmerged_1s))
new_isSGB_list_unmerged_1s=cbind(new_isSGB_list_unmerged_1s[,1],
                                     new_isSGB_list_unmerged_1s[isSGB_names])
names(new_isSGB_list_unmerged_1s)[1]='age_category'
new_isSGB_list_unmerged_1s=new_isSGB_list_unmerged_1s %>% select(-where(~ all(. == "notSGB")))
new_isSGB_list_unmerged_1s=new_isSGB_list_unmerged_1s %>% select(-where(~ all(. == "SGB")))
new_isSGB_list_unmerged_1s=mutate_if(new_isSGB_list_unmerged_1s, is.character, as.factor)
if (ncol(new_isSGB_list_unmerged_1s)>1){
    new_isSGB_list_unmerged_1s_glm=glm(age_category~.,
                                         data=new_isSGB_list_unmerged_1s,family='binomial')
    
    new_isSGB_list_unmerged_1s_glm_null = glm(age_category~1,
                                                data=new_isSGB_list_unmerged_1s,family='binomial')
    
    scope = list(lower=formula(new_isSGB_list_unmerged_1s_glm_null),
                 upper=formula(new_isSGB_list_unmerged_1s_glm))
    
    select.fit = SignifReg(new_isSGB_list_unmerged_1s_glm_null, scope = scope, direction = "forward", trace = TRUE,adjust.method = "fdr")
    R2=1-(logLik(select.fit) / logLik(new_isSGB_list_unmerged_1s_glm_null))
    new_isSGB_R2_unmerged_1=R2
    
    new_isSGB_tab_unmerged_1=lapply(SGBstrain_merged_list_1, function(x) x%>% count(isSGB,age_category) %>% spread(age_category,n,fill = 0))
    new_isSGB_tab_unmerged_1=lapply(new_isSGB_tab_unmerged_1, function(x) x%>% 
                                        select(-isSGB))
    
    df_pSGB=as.data.frame(select.fit$coefficients)
    df_pSGB$SGB=row.names(df_pSGB)
    rownames(df_pSGB)=NULL
    names(df_pSGB)[1]='coef'
    df_pSGB$p.value=select.fit$steps.info$max_pvalue
    df_pSGB$SGB=gsub('.isSGBSGB','',df_pSGB$SGB)
    df_pSGB=df_pSGB %>% filter(!is.na(p.value))
    if (nrow(df_pSGB)>=1){
      for (s in 1:nrow(df_pSGB)){
        ss=df_pSGB$SGB[s]
        df_pSGB$fp.value[s]=fisher.test(new_isSGB_tab_unmerged_1[ss][[1]])$p.value
      }
      df_pSGB$p.adjust=p.adjust(df_pSGB$fp.value,method='fdr')
    }
  }else{
    df_pSGB=data.frame(coef=0,sweep=NA,p.value=NA,fp.value=NA,p.adjust=NA)
  }
  
new_isSGB_p_df_unmerged_1_age=df_pSGB
names(new_isSGB_p_df_unmerged_1_age)[2]='SGB'



new_isSGB_list_unmerged_1_age=list()
new_isSGB_geo_p_unmerged_1_age=list()

for (s in 1:nrow(new_isSGB_p_df_unmerged_1_age)){
  taxa=new_isSGB_p_df_unmerged_1_age$SGB[s]
 
  
  stt=SGBstrain_merged_list_1[taxa]
  SGB_list=stt[[1]]%>% filter(isSGB=='SGB')
 # new_isSGB_list_unmerged_1_age[[s]]=sweep_list
  
#  names(new_isSGB_list_unmerged_1_age)[s]=taxa
  # is sweep dominated by a country, and also check sweep direction (enrich for age or no)
  observed_freq <- table(SGB_list$country)
  if (length(observed_freq)>1){
    new_isSGB_geo_p_unmerged_1_age[[s]]=c(chisq.test(observed_freq)$p.value,new_isSGB_p_df_unmerged_1_age$coef[s]>0)
  }else{
    new_isSGB_geo_p_unmerged_1_age[[s]]=c(0,new_isSGB_p_df_unmerged_1_age$coef[s]>0)
  }
  names(new_isSGB_geo_p_unmerged_1_age)[s]=taxa
}

new_isSGB_geo_p_unmerged_1_age_df=as.data.frame(t(bind_rows(new_isSGB_geo_p_unmerged_1_age,.id='id')))
names(new_isSGB_geo_p_unmerged_1_age_df)=c('geo_enriched','age_enriched')
new_isSGB_geo_p_unmerged_1_age_df=new_isSGB_geo_p_unmerged_1_age_df %>%
  mutate(geo_enriched=case_when(geo_enriched<=0.05 ~ TRUE, TRUE ~ FALSE))
new_isSGB_geo_p_unmerged_1_age_df=new_isSGB_geo_p_unmerged_1_age_df %>%
  mutate(age_enriched=case_when(age_enriched=='1' ~ 'enriched', TRUE ~ 'depleted'))


new_isSGB_geo_p_unmerged_1_age_df_nogeo=new_isSGB_geo_p_unmerged_1_age_df %>%
#  filter(geo_enriched==FALSE) %>% 
select(age_enriched)
names(new_isSGB_geo_p_unmerged_1_age_df_nogeo)='Senior'

colors=c('#4393C3','#D6604D')
p1=Heatmap(new_isSGB_geo_p_unmerged_1_age_df_nogeo,col=colors,row_names_max_width = unit(12, "cm"))
png('SGB_reg_age.png',width=2000,height=2500,res=300)
print(p1)
dev.off()

#filter for commensals
SGBstrain_merged_list_2=SGBstrain_merged_list_2[commensal_select]

new_isSGB_list_unmerged_2s=lapply(SGBstrain_merged_list_2, function(x) x%>% select(disease, isSGB))
new_isSGB_list_unmerged_2s=do.call(cbind,new_isSGB_list_unmerged_2s)
isSGB_names=grepl('isSGB',names(new_isSGB_list_unmerged_2s))
new_isSGB_list_unmerged_2s=cbind(new_isSGB_list_unmerged_2s[,1],
                                 new_isSGB_list_unmerged_2s[isSGB_names])
names(new_isSGB_list_unmerged_2s)[1]='disease'
new_isSGB_list_unmerged_2s=new_isSGB_list_unmerged_2s %>% select(-where(~ all(. == "notSGB")))
new_isSGB_list_unmerged_2s=new_isSGB_list_unmerged_2s %>% select(-where(~ all(. == "SGB")))
new_isSGB_list_unmerged_2s=mutate_if(new_isSGB_list_unmerged_2s, is.character, as.factor)
if (ncol(new_isSGB_list_unmerged_2s)>1){
  new_isSGB_list_unmerged_2s_glm=glm(disease~.,
                                     data=new_isSGB_list_unmerged_2s,family='binomial')
  
  new_isSGB_list_unmerged_2s_glm_null = glm(disease~1,
                                            data=new_isSGB_list_unmerged_2s,family='binomial')
  
  scope = list(lower=formula(new_isSGB_list_unmerged_2s_glm_null),
               upper=formula(new_isSGB_list_unmerged_2s_glm))
  
  select.fit = SignifReg(new_isSGB_list_unmerged_2s_glm_null, scope = scope, direction = "forward", trace = TRUE,adjust.method = "fdr")
  R2=1-(logLik(select.fit) / logLik(new_isSGB_list_unmerged_2s_glm_null))
  new_isSGB_R2_unmerged_2=R2
  
  new_isSGB_tab_unmerged_2=lapply(SGBstrain_merged_list_2, function(x) x%>% count(isSGB,disease) %>% spread(disease,n,fill = 0))
  new_isSGB_tab_unmerged_2=lapply(new_isSGB_tab_unmerged_2, function(x) x%>% 
                                    select(-isSGB))
  
  df_pSGB=as.data.frame(select.fit$coefficients)
  df_pSGB$SGB=row.names(df_pSGB)
  rownames(df_pSGB)=NULL
  names(df_pSGB)[1]='coef'
  df_pSGB$p.value=select.fit$steps.info$max_pvalue
  df_pSGB$SGB=gsub('.isSGBSGB','',df_pSGB$SGB)
  df_pSGB=df_pSGB %>% filter(!is.na(p.value))
  if (nrow(df_pSGB)>=1){
    for (s in 1:nrow(df_pSGB)){
      ss=df_pSGB$SGB[s]
      df_pSGB$fp.value[s]=fisher.test(new_isSGB_tab_unmerged_2[ss][[1]])$p.value
    }
    df_pSGB$p.adjust=p.adjust(df_pSGB$fp.value,method='fdr')
  }
}else{
  df_pSGB=data.frame(coef=0,sweep=NA,p.value=NA,fp.value=NA,p.adjust=NA)
}

new_isSGB_p_df_unmerged_2_CRC=df_pSGB
names(new_isSGB_p_df_unmerged_2_CRC)[2]='SGB'

new_isSGB_list_unmerged_2_CRC=list()
new_isSGB_geo_p_unmerged_2_CRC=list()

for (s in 1:nrow(new_isSGB_p_df_unmerged_2_CRC)){
  taxa=new_isSGB_p_df_unmerged_2_CRC$SGB[s]
  
  
  stt=SGBstrain_merged_list_2[taxa]
  SGB_list=stt[[1]]%>% filter(isSGB=='SGB')
#  new_isSGB_list_unmerged_2_CRC[[s]]=sweep_list
  
 # names(new_isSGB_list_unmerged_2_CRC)[s]=taxa
  # is sweep dominated by a country, and also check sweep direction (enrich for disease or no)
  observed_freq <- table(SGB_list$country)
  if (length(observed_freq)>1){
    new_isSGB_geo_p_unmerged_2_CRC[[s]]=c(chisq.test(observed_freq)$p.value,new_isSGB_p_df_unmerged_2_CRC$coef[s]<0) #here needs to be <0 because the correlation is for non-CRC
  }else{
    new_isSGB_geo_p_unmerged_2_CRC[[s]]=c(0,new_isSGB_p_df_unmerged_2_CRC$coef[s]<0)  #here needs to be <0 because the correlation is for non-CRC
  }
  names(new_isSGB_geo_p_unmerged_2_CRC)[s]=taxa
}

new_isSGB_geo_p_unmerged_2_CRC_df=as.data.frame(t(bind_rows(new_isSGB_geo_p_unmerged_2_CRC,.id='id')))
names(new_isSGB_geo_p_unmerged_2_CRC_df)=c('geo_enriched','CRC_enriched')
new_isSGB_geo_p_unmerged_2_CRC_df=new_isSGB_geo_p_unmerged_2_CRC_df %>%
  mutate(geo_enriched=case_when(geo_enriched<=0.05 ~ TRUE, TRUE ~ FALSE))
new_isSGB_geo_p_unmerged_2_CRC_df=new_isSGB_geo_p_unmerged_2_CRC_df %>%
  mutate(CRC_enriched=case_when(CRC_enriched=='1' ~ 'enriched', TRUE ~ 'depleted'))

#currently not filtering for geo specificity
new_isSGB_geo_p_unmerged_2_CRC_df_nogeo=new_isSGB_geo_p_unmerged_2_CRC_df %>%
#  filter(geo_enriched==FALSE) %>% 
  select(CRC_enriched)
names(new_isSGB_geo_p_unmerged_2_CRC_df_nogeo)='CRC'

colors=c('#4393C3','#D6604D')
p1=Heatmap(new_isSGB_geo_p_unmerged_2_CRC_df_nogeo,col=colors,row_names_max_width = unit(12, "cm"))
png('SGB_reg_CRC.png',width=1800,height=2000,res=300)
print(p1)
dev.off()

# #merge age and CRC datasets and plot with one heatmap
# age_CRC_merged_df=merge(new_isSGB_geo_p_unmerged_1_age_df_nogeo,
#                         new_isSGB_geo_p_unmerged_2_CRC_df_nogeo,by=0,all=T)
# rownames(age_CRC_merged_df)=age_CRC_merged_df$Row.names
# age_CRC_merged_df$Row.names=NULL
# 
# rownames(age_CRC_merged_df)=gsub('s__','',rownames(age_CRC_merged_df))
# 
# write.csv(age_CRC_merged_df,'age_CRC_SGB_df.csv')
# 
# colors2=c('#4393C3','#D6604D','grey60')
# p2=Heatmap(age_CRC_merged_df,col=colors2,na_col = "white",
#            rect_gp = gpar(col = "grey60", lwd = 2),
#            column_gap = unit(5, "mm"),
#            row_names_max_width = unit(12, "cm"),name='Direction')
# 
# png('SGB_reg_age_CRC.png',width=1800,height=2000,res=300)
# print(p2)
# dev.off()

#filter for commensals, CD
SGBstrain_merged_list_3a=SGBstrain_merged_list_3a[commensal_select]

new_isSGB_list_unmerged_3as=lapply(SGBstrain_merged_list_3a, function(x) x%>% select(disease, isSGB))
new_isSGB_list_unmerged_3as=do.call(cbind,new_isSGB_list_unmerged_3as)
isSGB_names=grepl('isSGB',names(new_isSGB_list_unmerged_3as))
new_isSGB_list_unmerged_3as=cbind(new_isSGB_list_unmerged_3as[,1],
                                 new_isSGB_list_unmerged_3as[isSGB_names])
names(new_isSGB_list_unmerged_3as)[1]='disease'
new_isSGB_list_unmerged_3as=new_isSGB_list_unmerged_3as %>% select(-where(~ all(. == "notSGB")))
new_isSGB_list_unmerged_3as=new_isSGB_list_unmerged_3as %>% select(-where(~ all(. == "SGB")))
new_isSGB_list_unmerged_3as=mutate_if(new_isSGB_list_unmerged_3as, is.character, as.factor)
if (ncol(new_isSGB_list_unmerged_3as)>1){
  new_isSGB_list_unmerged_3as_glm=glm(disease~.,
                                     data=new_isSGB_list_unmerged_3as,family='binomial')
  
  new_isSGB_list_unmerged_3as_glm_null = glm(disease~1,
                                            data=new_isSGB_list_unmerged_3as,family='binomial')
  
  scope = list(lower=formula(new_isSGB_list_unmerged_3as_glm_null),
               upper=formula(new_isSGB_list_unmerged_3as_glm))
  
  select.fit = SignifReg(new_isSGB_list_unmerged_3as_glm_null, scope = scope, direction = "forward", trace = TRUE,adjust.method = "fdr")
  R2=1-(logLik(select.fit) / logLik(new_isSGB_list_unmerged_3as_glm_null))
  new_isSGB_R2_unmerged_3a=R2
  
  new_isSGB_tab_unmerged_3a=lapply(SGBstrain_merged_list_3a, function(x) x%>% count(isSGB,disease) %>% spread(disease,n,fill = 0))
  new_isSGB_tab_unmerged_3a=lapply(new_isSGB_tab_unmerged_3a, function(x) x%>% 
                                    select(-isSGB))
  
  df_pSGB=as.data.frame(select.fit$coefficients)
  df_pSGB$SGB=row.names(df_pSGB)
  rownames(df_pSGB)=NULL
  names(df_pSGB)[1]='coef'
  df_pSGB$p.value=select.fit$steps.info$max_pvalue
  df_pSGB$SGB=gsub('.isSGBSGB','',df_pSGB$SGB)
  df_pSGB=df_pSGB %>% filter(!is.na(p.value))
  if (nrow(df_pSGB)>=1){
    for (s in 1:nrow(df_pSGB)){
      ss=df_pSGB$SGB[s]
      df_pSGB$fp.value[s]=fisher.test(new_isSGB_tab_unmerged_3a[ss][[1]])$p.value
    }
    df_pSGB$p.adjust=p.adjust(df_pSGB$fp.value,method='fdr')
  }
}else{
  df_pSGB=data.frame(coef=0,sweep=NA,p.value=NA,fp.value=NA,p.adjust=NA)
}

new_isSGB_p_df_unmerged_3a_CD=df_pSGB
names(new_isSGB_p_df_unmerged_3a_CD)[2]='SGB'

new_isSGB_list_unmerged_3a_CD=list()
new_isSGB_geo_p_unmerged_3a_CD=list()

for (s in 1:nrow(new_isSGB_p_df_unmerged_3a_CD)){
  taxa=new_isSGB_p_df_unmerged_3a_CD$SGB[s]
  
  
  stt=SGBstrain_merged_list_3a[taxa]
  SGB_list=stt[[1]]%>% filter(isSGB=='SGB')
  #  new_isSGB_list_unmerged_3a_CD[[s]]=sweep_list
  
  # names(new_isSGB_list_unmerged_3a_CD)[s]=taxa
  # is sweep dominated by a country, and also check sweep direction (enrich for disease or no)
  observed_freq <- table(SGB_list$country)
  if (length(observed_freq)>1){
    new_isSGB_geo_p_unmerged_3a_CD[[s]]=c(chisq.test(observed_freq)$p.value,new_isSGB_p_df_unmerged_3a_CD$coef[s]<0) #here needs to be <0 because the correlation is for non-CD
  }else{
    new_isSGB_geo_p_unmerged_3a_CD[[s]]=c(0,new_isSGB_p_df_unmerged_3a_CD$coef[s]<0)  #here needs to be <0 because the correlation is for non-CD
  }
  names(new_isSGB_geo_p_unmerged_3a_CD)[s]=taxa
}

new_isSGB_geo_p_unmerged_3a_CD_df=as.data.frame(t(bind_rows(new_isSGB_geo_p_unmerged_3a_CD,.id='id')))
names(new_isSGB_geo_p_unmerged_3a_CD_df)=c('geo_enriched','CD_enriched')
new_isSGB_geo_p_unmerged_3a_CD_df=new_isSGB_geo_p_unmerged_3a_CD_df %>%
  mutate(geo_enriched=case_when(geo_enriched<=0.05 ~ TRUE, TRUE ~ FALSE))
new_isSGB_geo_p_unmerged_3a_CD_df=new_isSGB_geo_p_unmerged_3a_CD_df %>%
  mutate(CD_enriched=case_when(CD_enriched=='1' ~ 'enriched', TRUE ~ 'depleted'))

#currently not filtering for geo specificity
new_isSGB_geo_p_unmerged_3a_CD_df_nogeo=new_isSGB_geo_p_unmerged_3a_CD_df %>%
  #  filter(geo_enriched==FALSE) %>% 
  select(CD_enriched)
names(new_isSGB_geo_p_unmerged_3a_CD_df_nogeo)='CD'

colors=c('#4393C3','#D6604D')
p1=Heatmap(new_isSGB_geo_p_unmerged_3a_CD_df_nogeo,col=colors,row_names_max_width = unit(12, "cm"))
png('SGB_reg_CD.png',width=1800,height=2000,res=300)
print(p1)
dev.off()

#filter for commensals, UC
SGBstrain_merged_list_3b=SGBstrain_merged_list_3b[commensal_select]

new_isSGB_list_unmerged_3bs=lapply(SGBstrain_merged_list_3b, function(x) x%>% select(disease, isSGB))
new_isSGB_list_unmerged_3bs=do.call(cbind,new_isSGB_list_unmerged_3bs)
isSGB_names=grepl('isSGB',names(new_isSGB_list_unmerged_3bs))
new_isSGB_list_unmerged_3bs=cbind(new_isSGB_list_unmerged_3bs[,1],
                                 new_isSGB_list_unmerged_3bs[isSGB_names])
names(new_isSGB_list_unmerged_3bs)[1]='disease'
new_isSGB_list_unmerged_3bs=new_isSGB_list_unmerged_3bs %>% select(-where(~ all(. == "notSGB")))
new_isSGB_list_unmerged_3bs=new_isSGB_list_unmerged_3bs %>% select(-where(~ all(. == "SGB")))
new_isSGB_list_unmerged_3bs=mutate_if(new_isSGB_list_unmerged_3bs, is.character, as.factor)
if (ncol(new_isSGB_list_unmerged_3bs)>1){
  new_isSGB_list_unmerged_3bs_glm=glm(disease~.,
                                     data=new_isSGB_list_unmerged_3bs,family='binomial')
  
  new_isSGB_list_unmerged_3bs_glm_null = glm(disease~1,
                                            data=new_isSGB_list_unmerged_3bs,family='binomial')
  
  scope = list(lower=formula(new_isSGB_list_unmerged_3bs_glm_null),
               upper=formula(new_isSGB_list_unmerged_3bs_glm))
  
  select.fit = SignifReg(new_isSGB_list_unmerged_3bs_glm_null, scope = scope, direction = "forward", trace = TRUE,adjust.method = "fdr")
  R2=1-(logLik(select.fit) / logLik(new_isSGB_list_unmerged_3bs_glm_null))
  new_isSGB_R2_unmerged_3b=R2
  
  new_isSGB_tab_unmerged_3b=lapply(SGBstrain_merged_list_3b, function(x) x%>% count(isSGB,disease) %>% spread(disease,n,fill = 0))
  new_isSGB_tab_unmerged_3b=lapply(new_isSGB_tab_unmerged_3b, function(x) x%>% 
                                    select(-isSGB))
  
  df_pSGB=as.data.frame(select.fit$coefficients)
  df_pSGB$SGB=row.names(df_pSGB)
  rownames(df_pSGB)=NULL
  names(df_pSGB)[1]='coef'
  df_pSGB$p.value=select.fit$steps.info$max_pvalue
  df_pSGB$SGB=gsub('.isSGBSGB','',df_pSGB$SGB)
  df_pSGB=df_pSGB %>% filter(!is.na(p.value))
  if (nrow(df_pSGB)>=1){
    for (s in 1:nrow(df_pSGB)){
      ss=df_pSGB$SGB[s]
      df_pSGB$fp.value[s]=fisher.test(new_isSGB_tab_unmerged_3b[ss][[1]])$p.value
    }
    df_pSGB$p.adjust=p.adjust(df_pSGB$fp.value,method='fdr')
  }
}else{
  df_pSGB=data.frame(coef=0,sweep=NA,p.value=NA,fp.value=NA,p.adjust=NA)
}

new_isSGB_p_df_unmerged_3b_UC=df_pSGB
names(new_isSGB_p_df_unmerged_3b_UC)[2]='SGB'

new_isSGB_list_unmerged_3b_UC=list()
new_isSGB_geo_p_unmerged_3b_UC=list()

for (s in 1:nrow(new_isSGB_p_df_unmerged_3b_UC)){
  taxa=new_isSGB_p_df_unmerged_3b_UC$SGB[s]
  
  
  stt=SGBstrain_merged_list_3b[taxa]
  SGB_list=stt[[1]]%>% filter(isSGB=='SGB')
  #  new_isSGB_list_unmerged_3b_UC[[s]]=sweep_list
  
  # names(new_isSGB_list_unmerged_3b_UC)[s]=taxa
  # is sweep dominated by a country, and also check sweep direction (enrich for disease or no)
  observed_freq <- table(SGB_list$country)
  if (length(observed_freq)>1){
    new_isSGB_geo_p_unmerged_3b_UC[[s]]=c(chisq.test(observed_freq)$p.value,new_isSGB_p_df_unmerged_3b_UC$coef[s]<0) #here needs to be <0 because the correlation is for non-UC
  }else{
    new_isSGB_geo_p_unmerged_3b_UC[[s]]=c(0,new_isSGB_p_df_unmerged_3b_UC$coef[s]<0)  #here needs to be <0 because the correlation is for non-UC
  }
  names(new_isSGB_geo_p_unmerged_3b_UC)[s]=taxa
}

new_isSGB_geo_p_unmerged_3b_UC_df=as.data.frame(t(bind_rows(new_isSGB_geo_p_unmerged_3b_UC,.id='id')))
names(new_isSGB_geo_p_unmerged_3b_UC_df)=c('geo_enriched','UC_enriched')
new_isSGB_geo_p_unmerged_3b_UC_df=new_isSGB_geo_p_unmerged_3b_UC_df %>%
  mutate(geo_enriched=case_when(geo_enriched<=0.05 ~ TRUE, TRUE ~ FALSE))
new_isSGB_geo_p_unmerged_3b_UC_df=new_isSGB_geo_p_unmerged_3b_UC_df %>%
  mutate(UC_enriched=case_when(UC_enriched=='1' ~ 'enriched', TRUE ~ 'depleted'))

#currently not filtering for geo specificity
new_isSGB_geo_p_unmerged_3b_UC_df_nogeo=new_isSGB_geo_p_unmerged_3b_UC_df %>%
  #  filter(geo_enriched==FALSE) %>% 
  select(UC_enriched)
names(new_isSGB_geo_p_unmerged_3b_UC_df_nogeo)='UC'

colors=c('#4393C3','#D6604D')
p1=Heatmap(new_isSGB_geo_p_unmerged_3b_UC_df_nogeo,col=colors,row_names_max_width = unit(12, "cm"))
png('SGB_reg_UC.png',width=1800,height=2000,res=300)
print(p1)
dev.off()

#filter for commensals

SGBstrain_merged_list_4=SGBstrain_merged_list_4[commensal_select]

new_isSGB_list_unmerged_4s=lapply(SGBstrain_merged_list_4, function(x) x%>% select(disease, isSGB))
new_isSGB_list_unmerged_4s=do.call(cbind,new_isSGB_list_unmerged_4s)
isSGB_names=grepl('isSGB',names(new_isSGB_list_unmerged_4s))
new_isSGB_list_unmerged_4s=cbind(new_isSGB_list_unmerged_4s[,1],
                                 new_isSGB_list_unmerged_4s[isSGB_names])
names(new_isSGB_list_unmerged_4s)[1]='disease'
new_isSGB_list_unmerged_4s=new_isSGB_list_unmerged_4s %>% select(-where(~ all(. == "notSGB")))
new_isSGB_list_unmerged_4s=new_isSGB_list_unmerged_4s %>% select(-where(~ all(. == "SGB")))
new_isSGB_list_unmerged_4s=mutate_if(new_isSGB_list_unmerged_4s, is.character, as.factor)
if (ncol(new_isSGB_list_unmerged_4s)>1){
  new_isSGB_list_unmerged_4s_glm=glm(disease~.,
                                     data=new_isSGB_list_unmerged_4s,family='binomial')
  
  new_isSGB_list_unmerged_4s_glm_null = glm(disease~1,
                                            data=new_isSGB_list_unmerged_4s,family='binomial')
  
  scope = list(lower=formula(new_isSGB_list_unmerged_4s_glm_null),
               upper=formula(new_isSGB_list_unmerged_4s_glm))
  
  select.fit = SignifReg(new_isSGB_list_unmerged_4s_glm_null, scope = scope, direction = "forward", trace = TRUE,adjust.method = "fdr")
  R2=1-(logLik(select.fit) / logLik(new_isSGB_list_unmerged_4s_glm_null))
  new_isSGB_R2_unmerged_4=R2
  
  new_isSGB_tab_unmerged_4=lapply(SGBstrain_merged_list_4, function(x) x%>% count(isSGB,disease) %>% spread(disease,n,fill = 0))
  new_isSGB_tab_unmerged_4=lapply(new_isSGB_tab_unmerged_4, function(x) x%>% 
                                    select(-isSGB))
  
  df_pSGB=as.data.frame(select.fit$coefficients)
  df_pSGB$SGB=row.names(df_pSGB)
  rownames(df_pSGB)=NULL
  names(df_pSGB)[1]='coef'
  df_pSGB$p.value=select.fit$steps.info$max_pvalue
  df_pSGB$SGB=gsub('.isSGBSGB','',df_pSGB$SGB)
  df_pSGB=df_pSGB %>% filter(!is.na(p.value))
  if (nrow(df_pSGB)>=1){
    for (s in 1:nrow(df_pSGB)){
      ss=df_pSGB$SGB[s]
      df_pSGB$fp.value[s]=fisher.test(new_isSGB_tab_unmerged_4[ss][[1]])$p.value
    }
    df_pSGB$p.adjust=p.adjust(df_pSGB$fp.value,method='fdr')
  }
}else{
  df_pSGB=data.frame(coef=0,sweep=NA,p.value=NA,fp.value=NA,p.adjust=NA)
}

new_isSGB_p_df_unmerged_4_T2D=df_pSGB
names(new_isSGB_p_df_unmerged_4_T2D)[2]='SGB'

new_isSGB_list_unmerged_4_T2D=list()
new_isSGB_geo_p_unmerged_4_T2D=list()

for (s in 1:nrow(new_isSGB_p_df_unmerged_4_T2D)){
  taxa=new_isSGB_p_df_unmerged_4_T2D$SGB[s]
  
  
  stt=SGBstrain_merged_list_4[taxa]
  SGB_list=stt[[1]]%>% filter(isSGB=='SGB')
  #  new_isSGB_list_unmerged_4_T2D[[s]]=sweep_list
  
  # names(new_isSGB_list_unmerged_4_T2D)[s]=taxa
  # is sweep dominated by a country, and also check sweep direction (enrich for disease or no)
  observed_freq <- table(SGB_list$country)
  if (length(observed_freq)>1){
    new_isSGB_geo_p_unmerged_4_T2D[[s]]=c(chisq.test(observed_freq)$p.value,new_isSGB_p_df_unmerged_4_T2D$coef[s]>0) #here needs to be >0 because the correlation is for non-T2D
  }else{
    new_isSGB_geo_p_unmerged_4_T2D[[s]]=c(0,new_isSGB_p_df_unmerged_4_T2D$coef[s]>0)  #here needs to be <0 because the correlation is for non-T2D
  }
  names(new_isSGB_geo_p_unmerged_4_T2D)[s]=taxa
}

new_isSGB_geo_p_unmerged_4_T2D_df=as.data.frame(t(bind_rows(new_isSGB_geo_p_unmerged_4_T2D,.id='id')))
names(new_isSGB_geo_p_unmerged_4_T2D_df)=c('geo_enriched','T2D_enriched')
new_isSGB_geo_p_unmerged_4_T2D_df=new_isSGB_geo_p_unmerged_4_T2D_df %>%
  mutate(geo_enriched=case_when(geo_enriched<=0.05 ~ TRUE, TRUE ~ FALSE))
new_isSGB_geo_p_unmerged_4_T2D_df=new_isSGB_geo_p_unmerged_4_T2D_df %>%
  mutate(T2D_enriched=case_when(T2D_enriched=='1' ~ 'enriched', TRUE ~ 'depleted'))

#currently not filtering for geo specificity
new_isSGB_geo_p_unmerged_4_T2D_df_nogeo=new_isSGB_geo_p_unmerged_4_T2D_df %>%
  #  filter(geo_enriched==FALSE) %>% 
  select(T2D_enriched)
names(new_isSGB_geo_p_unmerged_4_T2D_df_nogeo)='T2D'

colors=c('#4393C3','#D6604D')
p1=Heatmap(new_isSGB_geo_p_unmerged_4_T2D_df_nogeo,col=colors,row_names_max_width = unit(12, "cm"))
png('SGB_reg_T2D.png',width=1800,height=2000,res=300)
print(p1)
dev.off()


#merge all datasets and plot with one heatmap
age_CRC_CD_UC_T2D_merged_df=list(
  new_isSGB_geo_p_unmerged_1_age_df_nogeo,
  new_isSGB_geo_p_unmerged_2_CRC_df_nogeo,
  new_isSGB_geo_p_unmerged_3a_CD_df_nogeo,
  new_isSGB_geo_p_unmerged_3b_UC_df_nogeo,
  new_isSGB_geo_p_unmerged_4_T2D_df_nogeo
) %>%
  lapply(rownames_to_column, var = "rowname") %>%
  reduce(full_join, by = "rowname") %>%
  column_to_rownames(var = "rowname")


write.csv(age_CRC_CD_UC_T2D_merged_df,'age_CRC_CD_UC_T2D_merged_df.csv')

rownames(age_CRC_CD_UC_T2D_merged_df)=gsub('s__', '*', rownames(age_CRC_CD_UC_T2D_merged_df))
rownames(age_CRC_CD_UC_T2D_merged_df)=gsub('_', ' ', rownames(age_CRC_CD_UC_T2D_merged_df))
rownames(age_CRC_CD_UC_T2D_merged_df)=gsub(' SGB','* (SGB',rownames(age_CRC_CD_UC_T2D_merged_df))
rownames(age_CRC_CD_UC_T2D_merged_df)=paste0(rownames(age_CRC_CD_UC_T2D_merged_df),')')

age_CRC_CD_UC_T2D_merged_df$Species=rownames(age_CRC_CD_UC_T2D_merged_df)
age_CRC_CD_UC_T2D_merged_df_long <- age_CRC_CD_UC_T2D_merged_df %>%
  pivot_longer(-Species, names_to = "Condition", values_to = "Association") %>%
  mutate(Association = factor(Association, levels = c("depleted", "enriched"))) %>%
  mutate(Association = fct_explicit_na(Association, na_level = "no association"))
age_CRC_CD_UC_T2D_merged_df_long$Condition <- factor(age_CRC_CD_UC_T2D_merged_df_long$Condition, levels = c("Senior", "CRC", "CD", "UC","T2D"))
p_all=ggplot(age_CRC_CD_UC_T2D_merged_df_long, aes(x = Condition, y = Species, fill = Association)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("enriched" = "#D6604D", "depleted" = "#4393C3", "no association" = "grey90")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = ggtext::element_markdown(),  # allows italics/bold
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

png('SGB_associations_all.png',res=300,width=2000,height=2000)
p_all
dev.off()
