compare <- function(x, y) {            # Create custom function in R
  sum(x$V1 %in% y$V1)==nrow(x)
}

merge_sweeps=function(sf,cf){
  if (length(sf)>1){
    sfdf=expand.grid(names(sf),names(sf),stringsAsFactors = F)
    sfdf=sfdf %>% filter(Var1!=Var2)
    rm_list=apply(sfdf,1,function(x) compare(sf[[x[1]]],sf[[x[2]]]))
    rm_names=sfdf$Var1[rm_list]
    keep_names=names(sf)[!names(sf) %in% rm_names]
    sf=sf[keep_names]
  }
  if (length(cf)>0 & length(sf)>0 ){
    sfcfdf=expand.grid(names(sf),names(cf),stringsAsFactors = F)
    rm_list=apply(sfcfdf,1,function(x) compare(sf[[x[1]]],cf[[x[2]]]))
    rm_names=sfcfdf$Var1[rm_list]
    keep_names=names(sf)[!names(sf) %in% rm_names]
    sf=sf[keep_names]
  }
  tsf=c(cf,sf)
  return(tsf)
}


sweep_finder=function(t,cutoff){
  clade_list=get_clade_list(t)$clades
  ntips=length(t$tip.label)
  ratio_df=data.frame(node=ntips+1,ratio=NA,mean_outdist=NA)
  for (n in 2:t$Nnode){
    st=get_subtree_at_node(t, n)$subtree
    nt=length(st$tip.label)
    md=sum(get_all_pairwise_distances(st,only_clades=st$tip.label))/((nt-1)*nt)
    if (md<cutoff){
      
      # find last ancestral node
      an=get_ancestral_nodes(t,ntips+n, Nsplits=1)
      
      # find sister node
      sisters=clade_list[ntips+an,c(2,3)]
      sister_n=sisters[sisters!=(ntips+n)]
      
      # if sister node is a node:
      if (sister_n>ntips){
        
        # find all tip to node distance for sister node
        sister_st=get_subtree_at_node(t, sister_n-ntips)$subtree
        all_distances_sister_st = get_all_distances_to_root(sister_st)
        tip_distances_sister_st=all_distances_sister_st[1:Ntip(sister_st)]
        tips_sister_st=sister_st$tip.label
        
        # find all tip to node distance for node
        all_distances_st = get_all_distances_to_root(st)
        tip_distances_st=all_distances_st[1:Ntip(st)]
        
        # find all pairwise distances between tips for node
        all_tip_pairs_st_sq=get_all_pairwise_distances(st,only_clades = st$tip.label)
        all_tip_pairs_st=all_tip_pairs_st_sq[lower.tri(all_tip_pairs_st_sq)]
        tips_st=st$tip.label
        
        # all combinations of tip distances 
        all_diff_clade_dist=expand.grid(tip_distances_sister_st,tip_distances_st)
        all_diff_clade_tips=expand.grid(tips_sister_st,tips_st)
        names(all_diff_clade_tips)=c('Strain.1','Strain.2')
        ad=get_pairwise_distances(t,sister_n,ntips+n)  #add in distance between nodes
        all_diff_clade_dist$dist=all_diff_clade_dist$Var1+all_diff_clade_dist$Var2+ad
        ratio=mean(all_diff_clade_dist$dist)/mean(all_tip_pairs_st)
    
      }else{
        # find sister tip
        tips_sister_st=t$tip.label[sister_n]
        
        # find all tip to node distance for node
   
        all_distances_st = get_all_distances_to_root(st)
        tip_distances_st=all_distances_st[1:Ntip(st)]
        
        # find all pairwise distances between tips for node
        all_tip_pairs_st_sq=get_all_pairwise_distances(st,only_clades = st$tip.label)
        all_tip_pairs_st=all_tip_pairs_st_sq[lower.tri(all_tip_pairs_st_sq)]
        tips_st=st$tip.label
        
        # all combinations of tip distances 
        ad=get_pairwise_distances(t,sister_n,ntips+n) #sister tip to node distance
        all_diff_clade_dist=expand.grid(ad,tip_distances_st)
        all_diff_clade_tips=expand.grid(tips_sister_st,tips_st)
        all_diff_clade_dist$dist=all_diff_clade_dist$Var1+all_diff_clade_dist$Var2
        ratio=mean(all_diff_clade_dist$dist)/mean(all_tip_pairs_st)
      }
      ratio_df[n,]=c(ntips+n,ratio,mean(all_diff_clade_dist$dist))
    }
  }
  return (ratio_df)
}


compare_sweeps=function(sweep1,sweep2){
  sweep1=lapply(sweep1,function(x) x=filter(x, !grepl("marker", V1)))
  sweep2=lapply(sweep2,function(x) x=filter(x, !grepl("marker", V1)))
  sweep1=Filter(function(x) nrow(x) > 0, sweep1)
  sweep2=Filter(function(x) nrow(x) > 0, sweep2)
  
  if (length(sweep1)>0 & length(sweep2)>0 ){
    sfcfdf=expand.grid(names(sweep1),names(sweep2),stringsAsFactors = F)
    keep_list=apply(sfcfdf,1,function(x) compare(sweep1[[x[1]]],sweep2[[x[2]]]))
    keep_names=unique(sfcfdf$Var1[keep_list])
  }else(
    keep_names=c()
  )
  return(keep_names)
}

compare_sweeps_2=function(sweep1,sweep2){
  sweep1=lapply(sweep1,function(x) x=filter(x, !grepl("marker", V1)))
  sweep2=lapply(sweep2,function(x) x=filter(x, !grepl("marker", V1)))
  sweep1=Filter(function(x) nrow(x) > 0, sweep1)
  sweep2=Filter(function(x) nrow(x) > 0, sweep2)
  sweep1_total=bind_rows(sweep1)
  sweep2_total=bind_rows(sweep2)
  isweep=intersect(sweep1_total$V1,sweep2_total$V1)
  if (nrow(sweep1_total)>0 & nrow(sweep2_total)>0){
   rf=c(length(isweep)/nrow(sweep1_total),length(isweep)/nrow(sweep2_total))
  }else{
    rf=c(0,0)
  }
  return (rf)
}


compare_sweeps_withMG=function(sweep1,sweep2){
  if (length(sweep1)>0 & length(sweep2)>0 ){
    sfcfdf=expand.grid(names(sweep1),names(sweep2),stringsAsFactors = F)
    keep_list=apply(sfcfdf,1,function(x) compare(sweep1[[x[1]]],sweep2[[x[2]]]))
    keep_names=unique(sfcfdf$Var1[keep_list])
  }else(
    keep_names=c()
  )
  return(keep_names)
}

compare_sweeps_2_withMG=function(sweep1,sweep2){
  sweep1_total=bind_rows(sweep1)
  sweep2_total=bind_rows(sweep2)
  isweep=intersect(sweep1_total$V1,sweep2_total$V1)
  if (nrow(sweep1_total)>0 & nrow(sweep2_total)>0){
    rf=c(length(isweep)/nrow(sweep1_total),length(isweep)/nrow(sweep2_total))
  }else{
    rf=c(0,0)
  }
  return (rf)
}

get_new_mean= function(x) {y=get_subtree_at_node(t, x-Ntip(t))$subtree
new_mean=NA
if (Ntip(y)>2){
  ratio_df_1=sweep_finder(y,cutoff)
  ratio_df_1_5=ratio_df_1 %>% filter(ratio>=5)
  
  if (nrow(ratio_df_1_5)>0){
    tips_to_drop=lapply(ratio_df_1_5$node, function(x1) {z=as.data.frame(get_subtree_at_node(y, x1-Ntip(y))$subtree$tip.label);names(z)="V1";z})
    names(tips_to_drop)=paste0('tips',seq_along(tips_to_drop))
    tips_to_drop_m=merge_sweeps(tips_to_drop,list())
    tips_to_drop_m_f=lapply(tips_to_drop_m,function(t) t[1:(nrow(t)-1),])
    y1=drop.tip(y, unlist(tips_to_drop_m_f))
    all_tip_pairs_st_sq=get_all_pairwise_distances(y1,only_clades = y1$tip.label)
    all_tip_pairs_st=all_tip_pairs_st_sq[lower.tri(all_tip_pairs_st_sq)]
    new_mean=mean(all_tip_pairs_st)}};new_mean}
