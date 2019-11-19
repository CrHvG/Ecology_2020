# 0 functions ------

  # nested_cvs_plots - finds cvs plots that are 1000m, based on five modules (4 intensive and 1 R)
    cvs_1000_400_100m_plots<-function(x){
      require(moments)
      require(ineq)
      plots<-sort(unique(x[x$plot_stemObservationArea==1000,"OBSERVATION_ID"])) # restruct plots to those at 100, 400 and 1000m2
      stems2<-x[x$OBSERVATION_ID%in%plots,]
      mod_mat<-as.matrix(table(stems2$OBSERVATION_ID,stems2$Stratum.Index))
      mod_mat[mod_mat>0]<-1
      mod_num<-as.matrix(cbind(as.numeric(row.names(mod_mat)),apply(mod_mat,1,sum)))
      final_plots<-as.vector(mod_num[mod_num[,2]==5,1])
      final_stems<-x[x$OBSERVATION_ID%in%final_plots,c("OBSERVATION_ID","currentTaxonNameWithSp","Stratum.Index","stemDiameter","stem_count" ,"moduleArea")]
      final_stems<-final_stems[complete.cases(final_stems),]
      if (!all( unique(final_stems$moduleArea) %in% c(100,600))){stop ("Weird module configuration. Check it out.")}
      final_stems2<-final_stems[,-6]
      names(final_stems2)<-c("plot","spp","mod","dbh","count")
      final_stems2$mod<-as.factor(gsub("[mod ]","",final_stems2$mod))
      final_stems2[(final_stems2$mod=="r" | final_stems2$mod=="R"),"mod"]<-"R"
      final_stems2$ba<-round(((final_stems2$dbh/2)^2)*(3.1419),2)
      return(final_stems2)
    }

  # structure at 100m2 (now each of 100m)
    strct_div_fun<-function(x){
      rows<-length(unique(x$plot))*8 # i.e. all plots get seven rows (1000m, 600m, 400m, and 4 100m modules)
      misc_vars<-c("scale","mod","plot")
      div_vars<-c("SR","H_ct","D_ct","iD_ct","J_ct","H_ba","D_ba","iD_ba","J_ba","rar_SR","total_stems")
      strct_vars<-c("qmd","max_ba","dense_st","dense_ba","cv_ba","gc_ba","la_ba","skew_ba","kurt_ba")
      tmp_mat<-data.frame(matrix(0,nrow=rows, ncol=(length(misc_vars)+length(div_vars)+length(strct_vars))))
      colnames(tmp_mat)<-c(misc_vars,div_vars,strct_vars)
      progress.bar <- create_progress_bar("text")
      progress.bar$init(length(unique(x$plot)))

      m<-1
      for (i in unique(x$plot)){
        # 1000m
        tmp_mat[m,3:14]<-calc_div(tmp_mat,x[x$plot==i,])
        tmp_mat[m,c(1,15:23)]<-c(1000,calc_strct(tmp_mat,x[x$plot==i,])) 
        tmp_mat[m,c("dense_st","dense_ba")]<-tmp_mat[m,c("dense_st","dense_ba")]*(1/10)# adjust from 100m density
        # 600m
        tmp_mat[m+1,3:14]<-calc_div(tmp_mat,x[x$plot==i & x$mod=="R",])
        tmp_mat[m+1,c(1,15:23)]<-c(600,calc_strct(tmp_mat,x[x$plot==i & x$mod=="R",])) 
        tmp_mat[m+1,c("dense_st","dense_ba")]<-tmp_mat[m+1,c("dense_st","dense_ba")]*(1/6)# adjust from 100m density
        # 400m
        tmp_mat[m+2,3:14]<-calc_div(tmp_mat,x[x$plot==i & x$mod!="R",])
        tmp_mat[m+2,c(1,15:23)]<-c(400,calc_strct(tmp_mat,x[x$plot==i & x$mod!="R",])) 
        tmp_mat[m+2,c("dense_st","dense_ba")]<-tmp_mat[m+2,c("dense_st","dense_ba")]*(1/4)# adjust from 100m density
        # 100m
        tmp<-x[x$plot==i & x$mod!="R",] 
        tmp$mod<-as.numeric(as.character(tmp$mod))
        mod<-matrix(NA,4,20)
        n<-1
        for (j in unique(tmp$mod)){
          #j<-unique(tmp$mod)[4]
          tmp2<-tmp[tmp$mod==j,]
          mod[n,1:11]<-calc_div(mod,tmp2)
          mod[n,12:20]<-calc_strct(mod,tmp2)
          n<-n+1
        }
        tmp_mat[(m+3):(m+6),3:23]<-mod
        tmp_mat[(m+3):(m+7),1]<-100
        tmp_mat[(m+3):(m+6),2]<-unique(tmp$mod)
        #100m average
        tmp_mat[(m+7),2]<-0
        tmp_mat[(m+7),3:22]<-apply(mod,2,mean)
        m<-m+8
        progress.bar$step()
      }
      return(tmp_mat)
    }
    
    calc_strct<-function(strct_x,tmp_x){
      BAs<-unlist(mapply(rep, tmp_x$ba, tmp_x$count,SIMPLIFY=F))
      if (sd(BAs)==0 | is.na(sd(BAs))) BAs<-rnorm(length(BAs)*2,BAs,.001) # jitter if uniform
      strct_tmp<-matrix(NA,1,9)
      strct_tmp[,1]<-sqrt(sum((tmp_x$dbh*tmp_x$count)^2)/sum(tmp_x$count)) #qmd
      strct_tmp[,2]<-max(BAs) #max_ba
      strct_tmp[,3]<-sum(tmp_x$count)/100 #dense_st
      strct_tmp[,4]<-sum(BAs)/100 #dense_ba
      strct_tmp[,5]<-sd(BAs)/mean(BAs) # cv ba
      strct_tmp[,6]<-Gini(BAs) #gc_ba 
      strct_tmp[,7]<-Lasym(BAs)
      strct_tmp[,8]<-skewness(BAs)
      strct_tmp[,9]<-kurtosis(BAs) 
      return(strct_tmp)
    }
    
    calc_div<-function(strct_x,tmp_x){
      tmp_ct<-dcast(tmp_x, plot ~ spp,value.var='count',fill=0,fun.aggregate=sum)[,-1]
      tmp_ba<-dcast(tmp_x, plot ~ spp,value.var='ba',fill=0,fun.aggregate=sum)[,-1]
      div_tmp<-matrix(NA,1,12)
      div_tmp[,1]<-mean(tmp_x$plot) #plot
      if(length(unique(tmp_x$spp))==1){div_tmp[,2]<-1}else{div_tmp[,2]<-apply(tmp_ct, 1, function(c)sum(c!=0))}
      div_tmp[,3]<-diversity(tmp_ct, "shannon")  
      div_tmp[,4]<-diversity(tmp_ct, "simpson") 
      div_tmp[,5]<-diversity(tmp_ct, "invsimpson")  
      div_tmp[,6]<-ifelse(specnumber(tmp_ct)==1,0,diversity(tmp_ct)/log(specnumber(tmp_ct)))
      div_tmp[,7]<-diversity(tmp_ba, "shannon")  
      div_tmp[,8]<-diversity(tmp_ba, "simpson") 
      div_tmp[,9]<-diversity(tmp_ba, "invsimpson")  
      div_tmp[,10]<-ifelse(specnumber(tmp_ba)==1,0,diversity(tmp_ba)/log(specnumber(tmp_ba)))
      inext_output<-iNEXT(t(tmp_ct), q=0, datatype="abundance",endpoint=rarefied_sample_size,knots=20)  
      div_tmp[,11]<-inext_output$iNextEst[[1]][inext_output$iNextEst[[1]]$m==rarefied_sample_size,"qD"][1]
      div_tmp[,12]<-sum(tmp_ct)
      return(div_tmp)
    }
    
# 1 - input stems data
  setwd("...0_stems")
  stems_files<-list.files() # list files
  read.csv(stems_files[[1]])
  stems_raw<-lapply(stems_files,read.csv) # reads file list
  stems<-lapply(stems_raw,cvs_1000_400_100m_plots) # list of all plots w 100,400,600,1000m nested sampling

  # strct and div (with median at each scale for target rarefied abundance)
    g<-1;i<-2693
    stems_strct_list<-list()
    for (g in 1:5){
      x<-stems[[g]]
      rows<-length(unique(x$plot))*8 
      misc_vars<-c("scale","mod","plot")
      div_vars<-c("SR","H_ct","D_ct","iD_ct","J_ct","H_ba","D_ba","iD_ba","J_ba","rar_SR","total_stems")
      strct_vars<-c("qmd","max_ba","dense_st","dense_ba","cv_ba","gc_ba","la_ba","skew_ba","kurt_ba")
      tmp_mat<-data.frame(matrix(0,nrow=rows, ncol=(length(misc_vars)+length(div_vars)+length(strct_vars))))
      colnames(tmp_mat)<-c(misc_vars,div_vars,strct_vars)
      progress.bar <- create_progress_bar("text")
      progress.bar$init(length(unique(x$plot)))
      
      #plot loop
      m<-1
      for (i in unique(x$plot)){
        #i<-2693
        # 1000m
        tmp_mat[m,3:14]<-calc_div(tmp_mat,x[x$plot==i,])
        tmp_mat[m,c(1,15:23)]<-c(1000,calc_strct(tmp_mat,x[x$plot==i,])) 
        tmp_mat[m,c("dense_st","dense_ba")]<-tmp_mat[m,c("dense_st","dense_ba")]*(1/10)# adjust from 100m density
        # 600m
        tmp_mat[m+1,3:14]<-calc_div(tmp_mat,x[x$plot==i & x$mod=="R",])
        tmp_mat[m+1,c(1,15:23)]<-c(600,calc_strct(tmp_mat,x[x$plot==i & x$mod=="R",])) 
        tmp_mat[m+1,c("dense_st","dense_ba")]<-tmp_mat[m+1,c("dense_st","dense_ba")]*(1/6)# adjust from 100m density
        # 400m
        tmp_mat[m+2,3:14]<-calc_div(tmp_mat,x[x$plot==i & x$mod!="R",])
        tmp_mat[m+2,c(1,15:23)]<-c(400,calc_strct(tmp_mat,x[x$plot==i & x$mod!="R",])) 
        tmp_mat[m+2,c("dense_st","dense_ba")]<-tmp_mat[m+2,c("dense_st","dense_ba")]*(1/4)# adjust from 100m density
        # 100m
        tmp<-x[x$plot==i & x$mod!="R",] 
        tmp$mod<-as.numeric(as.character(tmp$mod))
        mod<-matrix(NA,4,21)
        n<-1
        for (j in unique(tmp$mod)){
          tmp2<-tmp[tmp$mod==j,]
          mod[n,1:12]<-calc_div(mod,tmp2)
          mod[n,13:21]<-calc_strct(mod,tmp2)
          n<-n+1
        }
        tmp_mat[(m+3):(m+6),3:23]<-mod
        tmp_mat[(m+3):(m+7),1]<-100
        tmp_mat[(m+3):(m+6),2]<-unique(tmp$mod)
        tmp_mat[(m+7),2]<-0
        tmp_mat[(m+7),3:23]<-apply(mod,2,mean)
        m<-m+8
        progress.bar$step()
      }
      stems_strct_list[[g]]<-tmp_mat
      print(g)
    }
    
    stems_strct<-do.call(rbind, stems_strct_list) # puts list into a data frame
    write.csv(stems_strct,"/Users/crhvg/Dropbox/chris/UNC/abd/4_CVS/data/1_woody_full.csv",row.names = F)


# 2. get cover values of only stems above 1cm (400m indicative of plot) sums all species for each of 4 100m subplots, and then takes mean -------

  rm(list=ls())
  library(plyr)
  library(dplyr)
  
  load('...all_woody.r')
  setwd("...0_multiscale_presence/")
  mp_files<-list.files()
  
  cvs_full_data<-read.csv("...s/18_CVS/data/cvs_full_data.csv")
  
  head(cvs_full_data)
  cvs_full_data[cvs_full_data$plot==1205,]
  plot_list<-unique(cvs_full_data$plot)
  
  setwd("...0_stems")
  stems_files<-list.files()
  
  cc<-as.data.frame(matrix(NA,length(plot_list),5))
  colnames(cc)<-c("plot","cc","tree_sr_1000","tree_sr_400","tree_sr_100")
  
  gm_mean = function(a){prod(a)^(1/length(a))}
  
  m<-1;depth<-1;i<-3;j<-1140;stem_mods<-"mod 1"
  for (i in 1:length(stems_files)){
    mp_tmp_state<-read.csv(paste("...0_multiscale_presence/",mp_files[i],sep=""))
    stems_tmp_state<-read.csv(paste("...0_stems/",stems_files[i],sep=""))
    state_plots<-unique(stems_tmp_state$OBSERVATION_ID)
    state_plots<-state_plots[state_plots%in%plot_list]
    
    for (j in 1:length(state_plots)){
      stems_tmp_plot<-stems_tmp_state[stems_tmp_state$OBSERVATION_ID==state_plots[j],]
      stems<-unique(stems_tmp_plot[stems_tmp_plot$stemDiameter>=5,4])
      my_summary<-stems_tmp_plot[stems_tmp_plot$Stratum.Index!="mod R",] %>% group_by(Stratum.Index) %>%  dplyr::summarize(x=length(unique(currentTaxonNameWithSp)))
      mp_tmp_plot<-mp_tmp_state[mp_tmp_state$OBSERVATION_ID==state_plots[j],]
      tmp_plot2<-mp_tmp_plot[is.na(mp_tmp_plot$cornerNumber) & mp_tmp_plot$currentTaxonNameWithSp%in%stems,]
      tmp_plot3<-ddply(tmp_plot2,~modNumber,summarise,sum=sum(coverPercent))
      cc[m,1]<-state_plots[j]
      cc[m,2]<-round(gm_mean(tmp_plot3[,2]))
      cc[m,3]<-length(stems)
      cc[m,4]<-length(unique(stems_tmp_plot[stems_tmp_plot$Stratum.Index!="mod R" & stems_tmp_plot$stemDiameter>=1,4]))
      cc[m,5]<-sum(my_summary[,2])/4
      print(m<-m+1)
    }
  }
  
  write.csv(cc,"...stems_cc_sr_5cm.csv",row.names=F)

  
# 3 prepare diversity, structure, EPA, and soils data ------------
    rm(list=ls())
    library(ecodist)
    home<-"...data"
    setwd(home)
  
  # diversity
    SR_woody<-read.csv("multiscale_sr_woody.csv")
    colnames(SR_woody)<-c("plot","w_sr_1000","w_sr_400","w_sr_100","w_sr_10","w_sr_1","w_sr_01","w_sr_001")
    SR_nonwoody<-read.csv("multiscale_sr_nonwoody.csv")
    colnames(SR_nonwoody)<-c("plot","nw_sr_1000","nw_sr_400","nw_sr_100","nw_sr_10","nw_sr_1","nw_sr_01","nw_sr_001")
  
  # strct
    stems_strct<-read.csv("1_woody_full.csv")
    strct<-stems_strct[stems_strct$mod==0,c(1:3,15:23)]
    strct_full<-strct[order(strct$plot),]
    
    m1000_strc<-strct_full[strct_full$scale==1000,-c(1,2)]
    colnames(m1000_strc)[-1]<-paste(colnames(m1000_strc)[-1],"_1000m",sep="")
    m400_strc<-strct_full[strct_full$scale==400,-c(1,2)]
    colnames(m400_strc)[-1]<-paste(colnames(m400_strc)[-1],"_400m",sep="")
    m100_strc<-strct_full[strct_full$scale==100,-c(1,2)]
    colnames(m100_strc)[-1]<-paste(colnames(m100_strc)[-1],"_100m",sep="")
  
  # cc and LAI
    ccp<-read.csv("stems_cc_sr_5cm.csv")
    names(ccp)[2]<-"cc"
    ccp$ccp<-ccp$cc
    ccp$cc<-ccp$cc/100
    ccp$ccp[ccp$ccp>100]<-100  # everything above 100 is 100%
  
  # load compositional and epa ecoregion data
    setwd(home)
    comp_classes<-read.csv("2_vege_classes.csv")
    test<-comp_classes[comp_classes$formation=="F008",]
    
    epa<-read.csv("0_epa.csv")
    table(epa$US_L3CODE)
    epa$US_L2CODE<-epa$US_L3CODE
    epa[epa$US_L3CODE==66,8]<-84
    epa[epa$US_L3CODE==45,8]<-83
    epa[epa$US_L3CODE==65,8]<-83
    epa[epa$US_L3CODE==63,8]<-85
    epa[epa$US_L3CODE==75,8]<-85
    epa2<-epa[,c(1,6,8)]
  
  # SR
    tmp<-merge(SR_nonwoody,SR_woody,by="plot")
    tmp2<-cbind(tmp[,1],tmp[,2:8]+tmp[,9:15],tmp[,-1])
    colnames(tmp2)[1:8]<-c("plot","all_sr_1000","all_sr_400","all_sr_100","all_sr_10","all_sr_1","all_sr_01","all_sr_001")
  
  # structure
    tmp4<-merge(tmp2,m1000_strc,by="plot")
    tmp5<-merge(tmp4,m400_strc,by="plot")
    tmp6<-merge(tmp5,m100_strc,by="plot")
    tmp8<-merge(tmp6,ccp[,c(1,2,6)],by="plot")
  
  # EPA ecoregion physiognomy
    tmp9<-merge(tmp8,comp_classes[,c(2,5:7)],by="plot")
    tmp10<-tmp9[tmp9$division %in% c("D006","D008","D011","D062"), ]
    cvs_full_data<-merge(tmp10,epa2,by="plot")
    cvs_full_data$comp<-NA
    cvs_full_data[cvs_full_data$formation=="F008",57]<-"Cool Temperate"
    cvs_full_data[cvs_full_data$formation=="F018",57]<-"Warm Temperate"
    cvs_full_data[cvs_full_data$formation=="F026",57]<-"Swamp Floodplains"
    cvs_full_data<-cvs_full_data[,-54]
    names(cvs_full_data)[56]<-"formation"
    
    cvs_full_data$ecogeo3<-NA
    cvs_full_data[cvs_full_data$US_L3CODE%in%c(63,65,75),57]<-"Coastal Forests"
    cvs_full_data[cvs_full_data$US_L3CODE%in%c(66,67),57]<-"App Forests"
    cvs_full_data[cvs_full_data$US_L3CODE==45,57]<-"Piedmont Forests"
    
    cvs_full_data$ecogeo4<-cvs_full_data$ecogeo3
    cvs_full_data[cvs_full_data$formation=="Swamp Floodplains",58]<-"Swamp Floodplains"
    
    cvs_full_data$ecogeo5<-cvs_full_data$ecogeo3
    cvs_full_data[cvs_full_data$macrogroup=="M007",59]<-"Longleaf Woodlands"
    
    cvs_full_data$ecogeo6<-cvs_full_data$ecogeo4
    cvs_full_data[cvs_full_data$macrogroup=="M007",60]<-"Longleaf Woodlands"
  
  # soils ---
    plot_list<-unique(cvs_full_data$plot)
    setwd(paste(home,"/0_soils",sep=""))
    soils_files<-list.files() # list files
    soils_raw<-lapply(soils_files,read.csv) 
    soils_full<-do.call(rbind, soils_raw)
    soils<-soils_full[soils_full$OBSERVATION.ID %in% plot_list,]
    names(soils)[30]<-"plot"
    soils2<-soils[,-c(5,31)]
    
    names(soils2)[c(1,3,4,6,9,10,11,12,13,14,16,17,18,21,22,23,24,25,28)]
    i<-1
    for (i in 1:28){
      hist(soils2[,i],main=names(soils2)[i])
    }
          # log transform: soilOrganic_avg      soilSilt_avg         soilClay_avg         exchangeCapacity_avg
          # S_avg                P_avg                Ca_ppm_avg            Mg_ppm_avg            K_ppm_avg             Na_ppm_avg
          # percent_Mg_avg       percent_K_avg        percent_Na_avg       B_ppm_avg            Fe_ppm_avg            Mn_ppm_avg            Cu_ppm_avg
          # Zn_ppm_avg           Ca_over_Mg_ppm_avg
    
    nonnormal<-c(1,3,4,6,9,10,11,12,13,14,16,17,18,21,22,23,24,25,28)
    
    for (i in nonnormal){
      soils2[,i]<-log(soils2[,i]+0.01)
      print(shapiro.test(soils2[,i]+0.01))
    }
    
    soils3<-soils2[complete.cases(soils2),]
    names(soils3) # RKP says drop %Ca, %Mg, %K, %Na & %other
    my_pca <- princomp(formula = ~., data = soils3[,-c(15:19,29)],na.action=na.exclude,cor=T,scores=T)          
    
    cor_mat<-mycor2m(soils3[,-c(15:19,29)],my_pca$scores[,1:3])
    
    cor_mat[order(cor_mat[,1],decreasing = T),] #soil fertility (base cation availablity)
  
  # write the table
    best_vars<-c("Ca_ppm_avg","Mg_ppm_avg","K_ppm_avg","exchangeCapacity_avg","soilPH_avg","Cu_ppm_avg","Zn_ppm_avg","baseSaturation_avg","soilSand_avg","soilSilt_avg","soilOrganic_avg","N_avg","percent_H_avg","Mn_ppm_avg")
    cor_mat<-mycor2m(soils3,my_pca$scores[,1:3])
    load <- with(my_pca, unclass(loadings))
    aload <- abs(load)
    percent_contrib<-sweep(aload, 2, colSums(aload), "/")
    
    percent_contrib<-round(percent_contrib[row.names(percent_contrib)%in%best_vars,1:3],3)
    cor_mat<-round(cor_mat[row.names(cor_mat)%in%best_vars,],3)
    
    pca_mat<-as.data.frame(cbind(cor_mat[,1],percent_contrib[,1],cor_mat[,2],percent_contrib[,2],cor_mat[,3],percent_contrib[,3]))
    
    cor_mat_pvals<-mycor2m_pval(soils3,my_pca$scores[,1:3])
    cor_mat_pvals<-round(cor_mat_pvals[row.names(cor_mat_pvals)%in%best_vars,],3)
    
    row.names(cor_mat_pvals)<-row.names(pca_mat)<-c("Organic matter (%)","sand (%)","silt (%)",
                                                    "pH","CEC","Base Saturation","N (%)","Ca (ppm)","Mg (ppm)",
                                                    "K (ppm)","H (%)","Mn (ppm)","Cu (ppm)", "Zn (ppm)")
    names(cor_mat_pvals)<-names(cor_mat)<-c("PCA 1","PCA 2","PCA 3")
    cor_mat_pvals<-cbind(cor_mat_pvals[,1],rep(0,14),cor_mat_pvals[,2],rep(0,14),cor_mat_pvals[,3],rep(0,14))
    pca_mat[which(cor_mat_pvals>0.01,arr.ind =T)]<-0.00
    pca_mat<-pca_mat[order(pca_mat[,1],decreasing = T),]
    names(pca_mat)<-c("vl","pc","vl","pc","vl","pc")
    write.csv(pca_mat,"...Table_1_cor_mat.csv")
  
  
  # merge and write
    soils4<-as.data.frame(cbind(soils3$plot,my_pca$scores[,1],my_pca$scores[,2]))
    soils5<-soils4[complete.cases(soils4),]
    colnames(soils5)<-c("plot","soil_fert1","soil_fert2")
    cvs_full_data2<-merge(cvs_full_data,soils5,by="plot")
  
  # add tree SRs
    cvs_full_data3<-merge(cvs_full_data2,ccp[,c(1,3:5)],by="plot")
    cvs_full_data3$ground_sr_1000<-cvs_full_data3$all_sr_1000-cvs_full_data3$tree_sr_1000
    cvs_full_data3$ground_sr_400<-cvs_full_data3$all_sr_400-cvs_full_data3$tree_sr_400
    cvs_full_data3$ground_sr_100<-cvs_full_data3$all_sr_100-cvs_full_data3$tree_sr_100
    write.csv(cvs_full_data3,"...cvs_full_data_5cm.csv",row.names =F)  
    