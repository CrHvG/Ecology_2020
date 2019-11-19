# Fig. 1 maps ---------

  # ArcGIS

# Fig 2. schematic ----------

  dev.off()
  rm(list=ls())
  library(ggplot2)
  library(cowplot)
  library(stringr)
  library(reshape2)
  library(ggthemes)
  library(RColorBrewer)
  
  setwd("...6_soil_maxes")
  cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  
  names(cvs_full_data)
  
  mydata<-cvs_full_data[,colnames(cvs_full_data)%in%c("all_sr_1000","all_sr_100","all_sr_10","cc","soil_fert1")]
  mydata2<-melt(mydata, id.vars=c("cc", "soil_fert1"))
  mydata2$variable<-str_split_fixed(mydata2$variable, "_", 3)[,3]
  mydata2$variable <- factor(mydata2$variable, levels = c("1000","100","10"))
  
  sz<-10
  mycols<-brewer.pal(n = 11, name = "RdYlBu")[c(2,4,10)]
  
  sf<-ggplot(mydata2, aes(soil_fert1,value,group=variable,colour = variable)) + 
    geom_point(size=.1,alpha=0.5) +
    geom_smooth(aes(fill = variable),method="loess", se=T, fullrange=T, level=0.95,span=.9,size=.5) +
    scale_color_manual(values=mycols,name = bquote('scale ('*m^2*")")) + 
    scale_fill_manual(values=mycols,name = bquote('scale ('*m^2*")")) +
    theme_bw() + 
    ggtitle("(a)")+
    theme(legend.position = c(0.2, 0.82),
          plot.title=element_text(hjust=-0.07,face="bold",vjust=-.5),
          text=element_text(size=sz,family = "Times"),
          axis.text=element_text(size=sz,family = "Times"),
          legend.background = element_blank(),
          legend.title=element_text(face="bold"),
          legend.key = element_blank(), 
          plot.background = element_blank(), 
          panel.grid = element_blank()) +
    ylab("species richness (SR)") + xlab("soil fertility (SF)") 
  
  
  cc<-ggplot(mydata2, aes(cc,value,group=variable,colour = variable)) + 
    geom_point(size=.1, alpha=0.5) +
    geom_smooth(aes(fill = variable),method="loess", se=TRUE, fullrange=FALSE, level=0.95,span=1,size=.5) +
    scale_color_manual(values=mycols) + 
    scale_fill_manual(values=mycols) +
    theme_bw() + 
    ggtitle("(b)")+
    theme(plot.title=element_text(hjust=-0.07,face="bold",vjust=-.5),
          text=element_text(size=sz,family = "Times"),
          axis.text.x=element_text(size=sz,family = "Times"),
          legend.background = element_blank(),
          legend.key = element_blank(), 
          axis.text.y=element_blank(),
          plot.background = element_blank(), 
          panel.grid = element_blank()) +
    ylab("")+xlab("canopy cover (CC)")
  
  legend_b <- get_legend(sf + theme(plot.margin=margin(l=0.1,r=0.5,unit="cm"),legend.text=element_text(family="Times"),legend.position="right"))
  
  
  prow3<-plot_grid(sf+theme(plot.margin=margin(c(.2,0,.2,0.2),unit="cm")),
                   cc+theme(legend.position="none",plot.margin=margin(c(.2,0.2,.2,0),unit="cm")),
                   ncol = 2, 
                   rel_widths = c(1,.9))
  prow3
  
  pdf("...Fig_2_loess.pdf",height=4,width=5)
  prow3
  dev.off()


# Fig. 3. Five-way polynomial - SF1 - all plots across scales -------
  
  rm(list=ls())
  library(ggplot2)
  library(cowplot)
  library(ecodist)
  library(quantreg)
  library(Hmisc)
  library(MASS)
  library(dplyr)
  library(plyr)
  library(rstan)
  library(rstanarm)
  
  setwd("...5_interaction_models")
  cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  # or the "all" one where SR is also scaled
  names(cvs_full_data)
  tmp1<-cvs_full_data[,c(2:22,63:68,61,50)]
  tmp1<-tmp1[complete.cases(tmp1),]
  tmp1[,1:27]<-apply(tmp1[,1:27],2,function(x){(x-min(x)) /(max(x)-min(x))})
  tmp1[,28:29]<-scale(tmp1[,28:29])
  
  poly_model_results<-as.data.frame(matrix(NA,27,28))
  names(poly_model_results)<-c("plot_subset","scale","species_subset","soil_2.5","soil_25","soil_50","soil_75","soil_97.5",
                               "cc_2.5","cc_25","cc_50","cc_75","cc_97.5",
                               "soil_quad_2.5","soil_quad_25","soil_quad_50","soil_quad_75","soil_quad_97.5",
                               "cc_quad_2.5","cc_quad_25","cc_quad_50","cc_quad_75","cc_quad_97.5",
                               "soil_cc_2.5","soil_cc_25","soil_cc_50","soil_cc_75","soil_cc_97.5")
  
  g<-1;i<-22
  for (i in 1:27){
    tmp<-tmp1[,c(i,28,29)]
    names(tmp)<-c("SR","SF","CC")
    poly_model<-stan_glm(SR~SF*CC+I(SF^2)+I(CC^2), data = tmp, QR = TRUE,chains = 4, iter = 10000) 
    beta_hat_CIs<-as.data.frame(summary(poly_model))[2:6,4:8]
    poly_model_results[g,1]<-"all"
    poly_model_results[g,2]<-strsplit(names(tmp1)[i], "_")[[1]][3]
    poly_model_results[g,3]<-strsplit(names(tmp1)[i], "_")[[1]][1]
    poly_model_results[g,4:28]<-cbind(beta_hat_CIs[1,],beta_hat_CIs[2,],beta_hat_CIs[3,],beta_hat_CIs[4,],beta_hat_CIs[5,])
    print(i);g<-g+1
  }
  poly_model_results$scale <- factor(poly_model_results$scale, levels = rev(c("1000","400","100","10","1","01","001")))
  poly_model_results$scale<-revalue(poly_model_results$scale, c("01"="0.1", "001"="0.01"))
  poly_model_results$species_subset <- factor(poly_model_results$species_subset, levels = c( "ground","nw","all","w","tree"))
  poly_model_results$species_subset<-revalue(poly_model_results$species_subset, c("nw"="herbaceous", "w"="woody","all"="total"))
  write.csv(poly_model_results,"...poly_model_results.csv")
  poly_model_results<-read.csv("...poly_model_results.csv")
  
  library(RColorBrewer)
  mycols<-c("chartreuse4", brewer.pal(n = 11, name = "PiYG")[8],"deepskyblue3",
            brewer.pal(n = 11, name = "Spectral")[c(4,2)]) #11 before the 4
  
  mysizes<-c(.25,.5,1,1.25,1.5)
  
  poly_model_results$soil_cc_sig<-poly_model_results$cc_quad_sig<-poly_model_results$soil_quad_sig<-poly_model_results$cc_sig<-poly_model_results$soil_sig<-as.character(poly_model_results$species_subset)
  i<-1
  
  for (i in 1:dim(poly_model_results)[1]){
    if (poly_model_results[i,4]<0 & poly_model_results[i,8]>0) {poly_model_results[i,29]<-"non"}
    if (poly_model_results[i,9]<0 & poly_model_results[i,13]>0) {poly_model_results[i,30]<-"non"}
    if (poly_model_results[i,14]<0 & poly_model_results[i,18]>0) {poly_model_results[i,31]<-"non"}
    if (poly_model_results[i,19]<0 & poly_model_results[i,23]>0) {poly_model_results[i,32]<-"non"}
    if (poly_model_results[i,24]<0 & poly_model_results[i,28]>0) {poly_model_results[i,33]<-"non"}
  }
  size_tmp<-14
  soil_fert1_plot<-ggplot(poly_model_results, aes(x = as.factor(scale),y=soil_50,colour=factor(species_subset),shape=factor(species_subset),fill=soil_sig)) +
    theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
    geom_errorbar(aes(ymin=soil_2.5, ymax=soil_97.5), width=0,position=position_dodge(width=-0.6)) +
    geom_point(position=position_dodge(width=-0.6),size=1.5) +
    guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +    
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = c("non"="white", "ground"=mycols[1], "herbaceous"=mycols[2], "total"=mycols[3], "woody"=mycols[4], "tree"=mycols[5])) +
    scale_colour_manual(values=mycols) + 
    ggtitle("(a)")+
    theme(plot.title=element_text(hjust=-0.13,face="bold",vjust=-1.5,family="Times",size=size_tmp),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(family="Times",size=(size_tmp+2)),
          legend.title=element_blank(),
          axis.title=element_text(family="Times",size=size_tmp),
          axis.text.y=element_text(family="Times"),
          axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
    xlab(bquote('scale ('*m^2*")")) + 
    ylab(expression(paste("SF (", beta[1],")")))
  
  
  soil2_plot<-ggplot(poly_model_results, aes(x = as.factor(scale),y=soil_quad_50,colour=factor(species_subset),shape=factor(species_subset),fill=soil_quad_sig)) +
    theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
    geom_errorbar(aes(ymin=soil_quad_2.5, ymax=soil_quad_97.5), width=0,position=position_dodge(width =-0.6)) +
    geom_point(position=position_dodge(width=-0.6),size=1.5) +
    guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +    
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
    scale_size_manual(values=mysizes) + scale_colour_manual(values=mycols) + 
    ggtitle("(b)")+
    theme(plot.title=element_text(hjust=-0.13,face="bold",vjust=-1.5,family="Times",size=size_tmp),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(family="Times",size=size_tmp),
          axis.text.y=element_text(family="Times"),
          axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
    xlab(" ") + ylab(expression(paste(SF^2," (", beta[2],")")))+labs(colour="Species")
  
  cc_plot<-ggplot(poly_model_results, aes(x = as.factor(scale),y=cc_50,colour=factor(species_subset),shape=factor(species_subset),fill=cc_sig)) +
    theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
    geom_errorbar(aes(ymin=cc_2.5, ymax=cc_97.5), width=0,position=position_dodge(width =-0.6)) +
    geom_point(position=position_dodge(width=-0.6),size=1.5) +
    guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +    
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
    scale_size_manual(values=mysizes) + scale_colour_manual(values=mycols) + 
    ggtitle("(c)")+
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=size_tmp),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(family="Times",size=size_tmp),
          axis.text.y=element_text(family="Times"),
          axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
    xlab(" ") + ylab(expression(paste("CC (", beta[3],")")))+labs(colour="Species")
  
  cc2_plot<-ggplot(poly_model_results, aes(x = as.factor(scale),y=cc_quad_50,colour=factor(species_subset),shape=factor(species_subset),fill=cc_quad_sig)) +
    theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
    geom_errorbar(aes(ymin=cc_quad_2.5, ymax=cc_quad_97.5), width=0,position=position_dodge(width =-0.6)) +
    geom_point(position=position_dodge(width=-0.6),size=1.5) +
    guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +    
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
    scale_size_manual(values=mysizes) + scale_colour_manual(values=mycols) + 
    ggtitle("(d)")+
    theme(plot.title=element_text(hjust=-0.13,face="bold",vjust=-1.5,family="Times",size=size_tmp),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(family="Times",size=size_tmp),
          axis.text.y=element_text(family="Times"),
          axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
    xlab(" ") + ylab(expression(paste(CC^2," (", beta[4],")")))+labs(colour="Species")
  
  soil_cc_plot<-ggplot(poly_model_results, aes(x = as.factor(scale),y=soil_cc_50,colour=factor(species_subset),shape=factor(species_subset),fill=soil_cc_sig)) +
    theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
    geom_errorbar(aes(ymin=soil_cc_2.5, ymax=soil_cc_97.5), width=0,position=position_dodge(width =-0.6)) +
    geom_point(position=position_dodge(width=-0.6),size=1.5) +
    guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +    
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
    scale_size_manual(values=mysizes) + scale_colour_manual(values=mycols) + 
    ggtitle("(e)")+
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=size_tmp),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(family="Times",size=size_tmp),
          axis.text.y=element_text(family="Times"),
          axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
    xlab(" ") + ylab(expression(paste("SF:CC (", beta[5],")")))+labs(colour="Species")
  
  setwd("...figs_tabs")
  
  prow <- plot_grid(soil_fert1_plot + theme(legend.position="none",plot.margin=margin(l=0.1,r=0,unit="cm")),
                    soil2_plot + theme(legend.position="none",axis.text.y=element_blank(),plot.margin=margin(l=0,unit="cm")),
                    cc_plot + theme(legend.position="none",axis.text.y=element_blank(),plot.margin=margin(l=0,unit="cm")),
                    cc2_plot + theme(legend.position="none",axis.text.y=element_blank(),plot.margin=margin(l=0,unit="cm")),
                    soil_cc_plot + theme(legend.position="none",axis.text.y=element_blank(),plot.margin=margin(l=0,r=0.1,unit="cm")),
                    nrow = 1,rel_widths = c(1,.85,.85,.85,.85),align="h")
  
  legend_b <- get_legend(soil_fert1_plot + theme(legend.title=element_blank(),legend.text=element_text(family="Times"),legend.position="bottom"))
  
  pdf("Fig_3_five_way_poly_scale_class4.pdf",height=4.5,width=11)
  plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
  dev.off()

# Fig. 4 SR~SF+SF2 across cc and scale --------

  dev.off()
  rm(list=ls())
  library(ecodist)
  library(quantreg)
  library(forcats)
  library(Hmisc)
  library(dplyr)
  library(RColorBrewer)
  library(cowplot)
  library(plyr)
  library(rstan)
  library(rstanarm)
  
    setwd("...6_soil_maxes")
    cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  
    group_num<-3
    cvs.dat <- cvs_full_data %>%  mutate(cc_grp = cut2(cc, g=group_num))
    my_groups<-unique(cvs.dat$cc_grp)
    
    cvs_full_data<-cvs_full_data[,c(2:22,61,50)]
    cvs_full_data<-cvs_full_data[complete.cases(cvs_full_data),]
    cvs_full_data[,1:21]<-apply(cvs_full_data[,1:21],2,function(x){(x-min(x)) /(max(x)-min(x))})
    cvs_full_data[,22:23]<-scale(cvs_full_data[,22:23])
    
    
    sz=12
    names(cvs_full_data)
    soil_fert_plots<-list();models_range<-list();poly_models<-list()
    scale<-1
    for(scale in 1:7){
      poly_model_results<-as.data.frame(matrix(NA,group_num*3,6))
      names(poly_model_results)<-c("scale","species_subset","cc_level","soil_2.5","soil_50","soil_97.5")
      tmp<-cvs_full_data[,c(scale,scale+7,scale+14,22,23)]
      
      groups<-split(tmp, cut2(tmp$cc, g=group_num))
      
      k<-3;i<-2;g<-1;gf<-2
      for (k in 1:group_num){
        for (gf in 1:3){
          tmp_subset<-groups[[k]][,c(gf,4,5)]
          names(tmp_subset)<-c("SR","SF","CC")
          poly_model<-stan_glm(SR~SF+I(SF^2), data = tmp_subset, QR = TRUE,chains = 4, iter = 10000) 
          print(plot(poly_model, plotfun = "combo", prob = 0.5)) 
          poly_model_results[g,1]<-strsplit(names(groups[[k]][,c(gf,4,5)])[1], "_")[[1]][3]
          poly_model_results[g,2]<-strsplit(names(groups[[k]][,c(gf,4,5)])[1], "_")[[1]][1]
          poly_model_results[g,3]<-as.character(sort(my_groups)[k])
          poly_model_results[g,4:6]<-as.data.frame(summary(poly_model))[2,c(4,6,8)]
          g<-g+1
        }
      }
      poly_model_results$species_subset[poly_model_results$species_subset=="w"]<-"woody"
      poly_model_results$species_subset[poly_model_results$species_subset=="nw"]<-"nonwoody"
      poly_models[[scale]]<-poly_model_results
      models_range[[scale]]<-c(min(poly_model_results$soil_2.5),max(poly_model_results$soil_97.5))
    }
    
    soil_betas<-do.call(rbind,poly_models)
    names(soil_betas)<-c("scale","species_subset","group","min","mean","max")
    write.csv(soil_betas,"...soil_betas.csv")
    soil_betas<-read.csv("...oil_betas.csv")
    soil_betas<-soil_betas[,-1]
    soil_betas$scale[37:63]<-rep(c(1,0.1,0.01),each=9)
    soil_betas$scale<-as.factor(soil_betas$scale)
    
    soil_betas$scale <- factor(soil_betas$scale, levels = c(1000,400,100,10,1,0.1,0.01))
    soil_betas_all<-subset(soil_betas,species_subset=="all")
    
    mycolors<- brewer.pal(n = 11, name = "RdYlBu")[c(1:4,9:11)]
    soil_betas_all_tmp<-soil_betas_all
    
    SF<- ggplot(soil_betas_all, aes(x = as.factor(group),y=mean,col=factor(scale),shape=factor(scale))) +
      theme_bw() +  coord_flip() + 
      geom_point(size=1.5,position=position_dodge(width = -0.6)) +
      guides(shape=FALSE,colour = guide_legend(override.aes = list(shape =c(15:19,15,18),fill=mycolors))) +    
      
      scale_shape_manual(values = c(15:19,15,18)) +
      geom_errorbar(aes(ymin=min, ymax=max), size=0.5,position=position_dodge(width = -0.6), width=0) +
      theme(plot.title=element_text(hjust=-0.1,vjust=-1.5,size=11,face="bold"),
            text=element_text(size=sz,family = "Times"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.height = unit(.2, "cm"))+
      geom_hline(yintercept = 0, linetype = "longdash") +
      ggtitle("(a)") + 
      xlab("canopy cover levels")  + labs(col=bquote('scale ('*m^2*")"),y=expression(SF,"~",SR (beta)))+
      ylab(expression(paste("SF (", beta[1],")"))) + scale_color_manual(values = mycolors) +
      geom_point(aes(x=1.255, y=soil_betas_all_tmp[1,5]),  fill="white",colour=mycolors[1],size=1.5,shape=22) +
      geom_point(aes(x=1.17, y=soil_betas_all_tmp[4,5]),  fill="white",colour=mycolors[2],size=1.5,shape=21) 

  # SR~SF2 across cc and scale

    for(scale in 1:7){
      poly_model_results<-as.data.frame(matrix(NA,group_num*3,6))
      names(poly_model_results)<-c("scale","species_subset","cc_level","soil_2.5","soil_50","soil_97.5")
      tmp<-cvs_full_data[,c(scale,scale+7,scale+14,22,23)]
      
      groups<-split(tmp, cut2(tmp$cc, g=group_num))
      
      k<-3;i<-2;g<-1;gf<-2
      for (k in 1:group_num){
        for (gf in 1:3){
          tmp_subset<-groups[[k]][,c(gf,4,5)]
          names(tmp_subset)<-c("SR","SF","CC")
          poly_model<-stan_glm(SR~SF+I(SF^2), data = tmp_subset, QR = TRUE,chains = 4, iter = 10000) 
          print(plot(poly_model, plotfun = "combo", prob = 0.5)) 
          poly_model_results[g,1]<-strsplit(names(groups[[k]][,c(gf,4,5)])[1], "_")[[1]][3]
          poly_model_results[g,2]<-strsplit(names(groups[[k]][,c(gf,4,5)])[1], "_")[[1]][1]
          poly_model_results[g,3]<-as.character(sort(my_groups)[k])
          poly_model_results[g,4:6]<-as.data.frame(summary(poly_model))[3,c(4,6,8)]
          g<-g+1
        }
      }
      poly_model_results$species_subset[poly_model_results$species_subset=="w"]<-"woody"
      poly_model_results$species_subset[poly_model_results$species_subset=="nw"]<-"nonwoody"
      poly_models[[scale]]<-poly_model_results
      models_range[[scale]]<-c(min(poly_model_results$soil_2.5),max(poly_model_results$soil_97.5))
    }
    
    soil_betas2<-do.call(rbind,poly_models)
    names(soil_betas2)<-c("scale","species_subset","group","min","mean","max")
    write.csv(soil_betas2,"..soil_betas2.csv")
    soil_betas2<-read.csv("...soil_betas2.csv")
    soil_betas2<-soil_betas2[,-1]
    soil_betas2$scale[37:63]<-rep(c(1,0.1,0.01),each=9)
    soil_betas2$scale<-as.factor(soil_betas2$scale)
    
    soil_betas2$scale <- factor(soil_betas2$scale, levels = c(1000,400,100,10,1,0.1,0.01))
    soil_betas_all<-subset(soil_betas2,species_subset=="all")
    
    mycolors<- brewer.pal(n = 11, name = "RdYlBu")[c(1:4,9:11)]
    soil_betas_all_tmp2<-soil_betas_all
    
    SF2<- ggplot(soil_betas_all, aes(x = as.factor(group),y=mean,col=factor(scale),shape=factor(scale))) +
      theme_bw() +  coord_flip() + 
      geom_point(size=1.5,position=position_dodge(width = -0.6)) +
      guides(shape=FALSE,colour = guide_legend(override.aes = list(shape =c(15:19,15,18),fill=mycolors))) +    
      
      scale_shape_manual(values = c(15:19,15,18)) +
      geom_errorbar(aes(ymin=min, ymax=max), size=0.5,position=position_dodge(width = -0.6), width=0) +
      theme(plot.title=element_text(hjust=-0.1,vjust=-1.5,size=11,face="bold"),
            text=element_text(size=sz,family = "Times"),
            axis.text.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.height = unit(.2, "cm"))+
      geom_hline(yintercept = 0, linetype = "longdash") +
      ggtitle("(b)") + 
      xlab("")  +
      ylab(expression(paste(SF^2," (", beta[2],")"))) + scale_color_manual(values = mycolors)+
      
      geom_point(aes(x=2.83, y=soil_betas_all_tmp2[18,5]),  fill="white",colour=mycolors[6],size=1.5,shape=22) +
      geom_point(aes(x=2.745, y=soil_betas_all_tmp2[21,5]),  fill="white",colour=mycolors[7],size=1,shape=23) +
      geom_point(aes(x=1.745, y=soil_betas_all_tmp2[20,5]),  fill="white",colour=mycolors[7],size=1,shape=23) 

    setwd("...figs_tabs")

    prow <- plot_grid(SF + theme(legend.position="none",plot.margin=margin(l=0.1,r=.1,unit="cm")),
                      SF2 + theme(legend.position="none",plot.margin=margin(l=0.1,r=.1,unit="cm")),rel_widths = c(1,0.8))
    legend_b <- get_legend(SF + theme(plot.margin=margin(l=0,r=0,unit="cm"),legend.text=element_text(family="Times"),legend.position="right"))
    
    pdf("Fig_4_five_way_poly_scale_class3.pdf",height=4,width=6)
    plot_grid( prow, legend_b, ncol = 2, rel_widths = c(1, .2))
    dev.off()


# Fig. 5  1000m2 comp across levels of cc (quartile) -----------------

  dev.off()
  rm(list=ls())
  library(ecodist)
  library(ggplot2)
  library(quantreg)
  library(Hmisc)
  library(reshape2)
  library(grid)
  library(cowplot)
  library(MASS)
  library(dplyr)
  library(plyr)
  
  setwd("...6_soil_maxes")
  cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  
  group_num<-3
  cvs.dat <- cvs_full_data %>%  mutate(cc_grp = cut2(cc, g=group_num))
  my_groups<-unique(cvs.dat$cc_grp)

  cvs_full_data[,c(61,50)]<-scale(cvs_full_data[,c(61,50)])
  tmp<-cvs_full_data[,c(2,61,50)]
  
  cvs.dat <- tmp %>%  mutate(cc_grp = cut2(cc, g=group_num))
  groups <- cvs.dat %>% split(.$cc_grp) 
  plot_list<-list();my_max1<-list();my_max2<-list()
  mycols<- brewer.pal(n = 11, name = "RdYlBu")[c(1,9,10)]
  sz<-10
  g.num<-2
  plots<-list()

    for (g.num in 1:group_num){
      tmp_data<-groups[[g.num]]
      
      quart_model<-rq(tmp_data[,1]~tmp_data[,2]+I(tmp_data[,2]^2), tau=.9)
      quad_model<-lm(all_sr_1000~soil_fert1+I(soil_fert1^2), data=tmp_data)      
      
      pframe <- with(tmp_data,expand.grid(soil_fert1=unique(soil_fert1)))
      tmp_data$quadratic <- round(predict(quad_model,newdata=pframe,type="response"))
      
      pframe <- with(tmp_data,expand.grid(soil_fert1=unique(soil_fert1)))
      tmp_data$quartile <- round(predict(quart_model,newdata=pframe,type="response"))
      
      tmp_data2<-tmp_data[,-c(3,4)]
      head(tmp_data)
      tmp_data3 <- melt(tmp_data2, id=c("all_sr_1000","soil_fert1"))
      names(tmp_data3)[3]<-"Model"
      
      tmp_data3$Model<-revalue(tmp_data3$Model, c("quadratic"="mean", "quartile"="90%"))
      
      my_max1[[g.num]]<-round(c(sort(tmp_data[which.max(fitted(quart_model)),2]),fitted(quart_model)[which(fitted(quart_model)==max(fitted(quart_model)))]),1)
      my_max2[[g.num]]<-round(c(sort(tmp_data[which.max(fitted(quad_model)),2]),fitted(quad_model)[which(fitted(quad_model)==max(fitted(quad_model)))]),1)
      
      tmp_data3$Model <- factor(tmp_data3$Model, levels = c("90%", "mean"))
      
      plots[[g.num]]<- ggplot(tmp_data3, aes(soil_fert1,all_sr_1000,group=Model)) + 
        geom_point(size=.1,color="grey80") +
        geom_smooth(aes(y=value,group = Model,color=Model),size=0.5) +
        scale_color_manual(values=mycols) +
        ylim(min(tmp$all_sr_1000),max(tmp$all_sr_1000)) +
        ylab("") + xlab("") +
        coord_cartesian(clip = 'off') +
        theme(text=element_text(size=sz,family = "Times"),
              axis.text=element_text(size=sz,family = "Times"),
              legend.position=c(.1,.8),legend.title=element_blank(),
              legend.text=element_text(size=sz,family = "Times"),
              legend.key.height = unit(.1, "cm"),
              plot.title = element_text(size=sz,hjust = 0.5,face="plain",margin = margin(t = 0, b = -10))) +
        scale_x_continuous(expand = c(0,0),limits = c(min(tmp$soil_fert1),max(tmp$soil_fert1))) +
        guides(color=guide_legend(override.aes=list(fill=NA)))
    }
    
    dev.off()
    prow <- plot_grid( plots[[3]]+ylab("")+theme(axis.text.x = element_blank(),plot.margin = unit(c(.2, .8, 0, .3), "cm")) + ggtitle(paste("CC: ",sort(my_groups)[3],sep="")) +
                         geom_segment(x = my_max1[[3]][1], y = my_max1[[3]][2], xend = my_max1[[3]][1], yend = 0,color=mycols[1], linetype = 2) +
                         geom_segment(x = my_max2[[3]][1], y = my_max2[[3]][2], xend = my_max2[[3]][1], yend = 0,color=mycols[2], linetype = 2) +
                         geom_point(aes(x=my_max1[[3]][1], y=my_max1[[3]][2]), colour=mycols[1],size=1) +
                         geom_point(aes(x=my_max2[[3]][1], y=my_max2[[3]][2]), colour=mycols[2],size=1) +
                         annotate("text", x = my_max1[[3]][1], y = my_max1[[3]][2]+20, label = round(my_max1[[3]][1],1),size=sz/3,color=mycols[1]) +
                         annotate("text", x = my_max2[[3]][1], y = my_max2[[3]][2]+20, label = round(my_max2[[3]][1],1),size=sz/3,color=mycols[3]),
                       
                       plots[[2]]+  ylab(bquote('SR (1000'~m^2*")")) + theme(axis.text.x = element_blank(),legend.position = "none",plot.margin = unit(c(-.3, .8, 0, .15), "cm"))+ ggtitle(paste("CC: ",sort(my_groups)[2],sep="")) +
                         geom_segment(x = my_max1[[2]][1], y = my_max1[[2]][2], xend = my_max1[[2]][1], yend = 0,color=mycols[1], linetype = 2) +
                         geom_segment(x = my_max2[[2]][1], y = my_max2[[2]][2], xend = my_max2[[2]][1], yend = 0,color=mycols[2], linetype = 2) +
                         geom_point(aes(x=my_max1[[2]][1], y=my_max1[[2]][2]), colour=mycols[1],size=1) +
                         geom_point(aes(x=my_max2[[2]][1], y=my_max2[[2]][2]), colour=mycols[2],size=1) +
                         annotate("text", x = my_max1[[2]][1], y = my_max1[[2]][2]+20, label = round(my_max1[[2]][1],1),size=sz/3,color=mycols[1]) +
                         annotate("text", x = my_max2[[2]][1], y = my_max2[[2]][2]+20, label = round(my_max2[[2]][1],1),size=sz/3,color=mycols[3]),
                       
                       plots[[1]]+ ylab("")+xlab("SF")+ theme(legend.position = "none",plot.margin = unit(c(-.3, .8, 0, .3), "cm"))+ ggtitle(paste("CC: ",sort(my_groups)[1],sep="")) +
                         geom_segment(x = my_max1[[1]][1], y = my_max1[[1]][2], xend = my_max1[[1]][1], yend = 0,color=mycols[1], linetype = 2) +
                         geom_segment(x = my_max2[[1]][1], y = my_max2[[1]][2], xend = my_max2[[1]][1], yend = 0,color=mycols[2], linetype = 2) +
                         geom_point(aes(x=my_max1[[1]][1], y=my_max1[[1]][2]), colour=mycols[1],size=1) +
                         geom_point(aes(x=my_max2[[1]][1], y=my_max2[[1]][2]), colour=mycols[2],size=1) +
                         annotate("text", x = my_max1[[1]][1], y = my_max1[[1]][2]+20, label = round(my_max1[[1]][1],1),size=sz/3,color=mycols[1]) +
                         annotate("text", x = my_max2[[1]][1], y = my_max2[[1]][2]+20, label = round(my_max2[[1]][1],1),size=sz/3,color=mycols[3]),
                       
                       rel_heights=c(1,.8,.8),ncol = 1)
    prow
    pdf("...Fig_5_opt_3_levels_quart.pdf",height=4,width=5)
    prow
    dev.off()

#Fig. 6 SF_opt across scales -----------

  rm(list = ls(all = TRUE))
  library(tidyverse)
  library(rstan)
  library(rstanarm)
  library(Hmisc)
  
  cvs_full_data1<-read_csv("...cvs_full_data_5cm.csv")  # or the "all" one where SR is also scaled
  cvs_full_data1<-cvs_full_data1[,c(1:8,61,50)]
  cvs_full_data<-cvs_full_data1[complete.cases(cvs_full_data1),]
  cvs_full_data[,2:8]<-apply(cvs_full_data[,2:8],2,function(x){(x-min(x)) /(max(x)-min(x))})
  cvs_full_data[,9:10]<-scale(cvs_full_data[,9:10])
  
  n.grp<-3;scale<-6
  my_boxes<-list();results<-list()
  for (scale in c(1:7)){
    
    cvs.dat<-cvs_full_data[,c(1,(scale+1),9,10)]
    names(cvs.dat)[2]<-"SR"
    
    cvs.dat <- cvs.dat %>%  mutate(cc_grp = cut2(cc, g=n.grp))
    cvs.dist <- cvs.dat %>% split(.$cc_grp) 
    
    fit_stan_lm <- function(dat) {stan_glm(SR ~ soil_fert1 + I(soil_fert1^2),  data = dat)}

    get_opt_fert_NB <- function(stan) {
      stan.mcmc <- as.data.frame(stan) %>%
        as_tibble()
      colnames(stan.mcmc) <- c('int', 'soil_fert1', 'soil_fert12','reciprocal_dispersion')
      stan.mcmc %>%
        mutate(opt_fert = -0.5*soil_fert1/soil_fert12)
    }
    
    cvs.post <- cvs.dat %>% 
      split(.$cc_grp) %>% 
      map(. %>% fit_stan_lm())
    
    
    cvs.fert_opt<- cvs.post %>% 
      map(. %>% get_opt_fert_NB())
    
    nonsig_index<-lapply(cvs.post,function(x){
      all(as.data.frame(summary(x))[3,4]<0,as.data.frame(summary(x))[3,8]>0)
    })
    
    
    pe<-cvs.fert_opt %>%
      lapply(., function(x) {t(apply(x , 2 , quantile , probs = c(0.025, .5, 0.975) , na.rm = TRUE ))} )  %>%
      lapply(., `[`,5,)
    pe2<-as.data.frame(do.call(rbind,pe))
    pe3<-cbind(row.names(pe2),pe2)
    pe3[which(nonsig_index==TRUE),2:4]<-c(NA,NA,NA)
    names(pe3)<-c("group","min","mean","max")
    
    my_boxes[[scale]]<-pe3
    results[[scale]]<-lapply(cvs.post,as.data.frame(summary))
  }

    pe4<-do.call(rbind,my_boxes)
    pe5<-cbind(pe4,rep(c(1000,400,100,10,1,0.1,0.01), each=3))
    unique(pe5$group)
    names(pe5)[5]<-"scale"
    pe5$scale<-as.factor(pe5$scale)
    pe5$scale<-fct_rev(pe5$scale)
    
    cvs.dat <- cvs_full_data1 %>%  mutate(cc_grp = cut2(cc, g=n.grp))
    
    
    levels(pe5$group)<-rev(unique(cvs.dat$cc_grp))
    pe5$group<- levels(pe5$group)[c(2,3,1)]
    write.csv(pe5,"...figs_tabs/Fig_6.csv")
    
    library(RColorBrewer)
    mycolors<- brewer.pal(n = 11, name = "RdYlBu")[c(1:4,9:11)]
    
    sz=10
    fert_opt_plot<-ggplot(pe5, aes(x = as.factor(group),y=mean,col=as.factor(scale),shape=factor(scale))) +
      theme_bw() +  coord_flip() + 
      geom_point(size=1,position=position_dodge(width = -0.6)) +
      guides(shape=FALSE,colour = guide_legend(override.aes = list(shape =c(15:19,15,18),fill=mycolors))) +    
      geom_errorbar(aes(ymin=min, ymax=max), size=0.3,position=position_dodge(width = -0.6), width=0) +
      scale_shape_manual(values = c(15:19,15,18)) +
      theme(plot.title=element_text(hjust=0.5,size=11),
            text=element_text(size=sz,family = "Times"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.height = unit(.2, "cm"))+
      scale_color_manual(values = mycolors) +
      xlab("canopy cover levels") + labs(col=bquote('scale ('*m^2*")"),y=expression(SF[maxSR]))
    
    pdf("...Fig_6_opt.pdf",height=3,width=3)
    fert_opt_plot
    dev.off()


# Fig 7 Soil and cc by plot segmentation ########

    dev.off()
    rm(list=ls())
    library(ecodist)
    library(quantreg)
    library(Hmisc)
    library(MASS)
    library(stringr)
    library(reshape)
    library(cowplot)
    
    cvs_full_data<-read.csv("...cvs_full_data_5cm.csv") 
    setwd("...figs_tabs")
    names(cvs_full_data)
    test<-cvs_full_data[,c(61,57)]   #23,51
    
    test2<-test %>% melt(id.vars=c("ecogeo3")) 
    test2
    
    library(RColorBrewer)
    mycols<- brewer.pal(n = 11, name = "Spectral")[c(1,3,5)]
    sz1<-11
    sz2<-12
    
    means<-aggregate(test2$value, by=list(Category=test2$ecogeo3), FUN=median)
    means[,2]<-round(means[,2],1)
    v4<-
      
      ggplot(aes(y = value, x = factor(ecogeo3), fill = ecogeo3,colour=ecogeo3), data = subset(test2,variable=="soil_fert1")) + 
      geom_violin() +
      xlab("")+ ylab("soil fertility (index)") + 
      theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.2, 0, 0.2), "cm"),
            axis.title=element_text(size = sz2,family="Times"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols) + 
      scale_fill_manual(values=mycols) + 
      ggtitle("(a)") +
      annotate("text",x=1,y=means[1,2],label=means[1,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=2,y=means[2,2],label=means[2,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=3,y=means[3,2],label=means[3,2],colour="black",hjust=0.5,family="Times") 
    
    
    test<-cvs_full_data[,c(50,57)]   #23,51
    test2<-test %>% melt(id.vars=c("ecogeo3")) 
    test2$ecogeo3<-revalue(test2$ecogeo3, c("App Forests"="Appalachian Forests","Piedmont Forests"="Piedmont Forests"))
    test2$ecogeo3<-factor(test2$ecogeo3, levels = c("Coastal Forests", "Piedmont Forests", "Appalachian Forests"))
    
    means2<-aggregate(test2$value, by=list(Category=test2$ecogeo3), FUN=median)
    means2[,2]<-round(means2[,2],1)
    
    v5<- ggplot(aes(y = value, x = factor(ecogeo3), fill = ecogeo3,colour=ecogeo3), data = subset(test2,variable=="cc")) + 
      geom_violin() +    
      xlab("")+ ylab("canopy closure") + 
      theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
            legend.position = "none",
            plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"),
            axis.title=element_text(size = sz2,family="Times"),
            axis.text.x = element_text(angle = 15,hjust = 1,size = sz2,family="Times"),
            axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols) + 
      scale_fill_manual(values=mycols) +
      ggtitle("(b)") +
      annotate("text",x=1,y=means2[1,2],label=means2[1,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=2,y=means2[2,2],label=means2[2,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=3,y=means2[3,2],label=means2[3,2],colour="black",hjust=0.5,family="Times") 
    
    test<-cvs_full_data[,c(61,60)]   #23,51
    test2<-test %>% melt(id.vars=c("ecogeo6")) 
    test2$ecogeo6<-as.character(test2$ecogeo6)
    test2[test2$ecogeo6!="Longleaf Woodlands",1]<-"Not Longleaf"
    mycols<- brewer.pal(n = 11, name = "Spectral")[c(8,11)]
    
    means3<-aggregate(test2$value, by=list(Category=test2$ecogeo6), FUN=median)
    means3[,2]<-round(means3[,2],1)
    
    
    v6<- ggplot(aes(y = value, x = factor(ecogeo6), fill = ecogeo6,colour=ecogeo6), data = subset(test2,variable=="soil_fert1")) + 
      geom_violin(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9)) +
      xlab("")+ ylab("") + 
      theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.2, 0, 0.2), "cm"),
            axis.title=element_text(size = sz2,family="Times"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols) + 
      scale_fill_manual(values=mycols) +
      ggtitle("(c)") + 
      annotate("text",x=1,y=means3[1,2],label=means3[1,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=2,y=means3[2,2],label=means3[2,2],colour="black",hjust=0.5,family="Times")
    
    
    test<-cvs_full_data[,c(50,60)]   #23,51
    test2<-test %>% melt(id.vars=c("ecogeo6")) 
    test2$ecogeo6<-as.character(test2$ecogeo6)
    test2[test2$ecogeo6!="Longleaf Woodlands",1]<-"Not Longleaf"
    test2$ecogeo6<-as.factor(test2$ecogeo6)
    
    means4<-aggregate(test2$value, by=list(Category=test2$ecogeo6), FUN=median)
    means4[,2]<-round(means4[,2],1)
    
    
    v7<-ggplot(aes(y = value, x = ecogeo6, fill = ecogeo6,colour=ecogeo6), data = subset(test2,variable=="cc")) + 
      geom_violin(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9)) +
      xlab("")+ ylab("") + 
      theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
            legend.position = "none",
            plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"),
            axis.title=element_text(size = sz2,family="Times"),
            axis.text.x = element_text(angle = 15,hjust = 1,size = sz2,family="Times"),
            axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols) + 
      scale_fill_manual(values=mycols) +
      ggtitle("(d)") + 
      annotate("text",x=1,y=means4[1,2],label=means4[1,2],colour="black",hjust=0.5,family="Times") +
      annotate("text",x=2,y=means4[2,2],label=means4[2,2],colour="black",hjust=0.5,family="Times")
    
    
    prow<-plot_grid(v4,v6,v5,v7,align="v",hjust = -.2, ncol = 2,rel_heights=c(.83,1))
    
    pdf("fig_7_plots_2_box_plots_ecogeo3_2way.pdf",height=5,width=6)
    prow
    dev.off()

# Fig. 8  boxplots and violins of ecoregions ecogeo3 ------

  dev.off()
  rm(list=ls())
  library(ecodist)
  library(quantreg)
  library(Hmisc)
  library(MASS)
  library(stringr)
  library(reshape)
  library(cowplot)
  library(plyr)
  
  cvs_full_data<-read.csv("...cvs_full_data_5cm.csv") 
  setwd("...figs_tabs")
  names(cvs_full_data)
  cvs_full_data$cc<-cvs_full_data$cc/100
  table(cvs_full_data$ecogeo6)
  annotate_size<-2.8
  
  levels(cvs_full_data$ecogeo6)[5]<-"Riparian"
  
  op <- par(mfrow = c(4,2), mar=c(1, 4, 1, 1) + 0.1)
  
  test<-cvs_full_data[,c(57,2,4,8,9,11,15,16,18,22,63,65,66,68)]   #23,51
  test2<-test %>% melt(id.vars=c("ecogeo3")) 
  test3<-cbind(test2[,1],as.data.frame(str_split_fixed(test2$variable, "_", 3))[,-2],test2[,3])
  head(test3)
  names(test3)<-c("ecoregion","class","scale","value")
  
  test3$class <- factor(test3$class, levels = c("all","ground","nw","w","tree"))
  test3$class<-revalue(test3$class, c("nw"="herbaceous", "w"="woody"))
  
  test3$ecoregion <- factor(test3$ecoregion, levels = c("Longleaf Woodlands","Coastal Forests","Riparian","Piedmont Forests","App Forests"))
  test3$ecoregion<-revalue(test3$ecoregion, c("Riparian"="Swamp Floodplains", "App Forests"="Appalachian Forests"))
  
  
  library(RColorBrewer)
  mycols<-c("deepskyblue3","chartreuse4", brewer.pal(n = 11, name = "PiYG")[8],
            brewer.pal(n = 11, name = "Spectral")[c(4,2)]) #11 before the 4
  
  mycols<-brewer.pal(n = 11, name = "Spectral")[7:11]
  
  ghj<-subset(test3,scale=="1000")
  max(ghj$value)
  sz1<-11
  sz2<-11
  test4<-test3[test3$scale=="1000",]
  means<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
  maxes<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
  
  v1<-ggplot(aes(y = value, x = ecoregion, fill = class,colour=class), data = subset(test3,scale=="1000")) + 
    geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
    xlab("")+ ylab(bquote('SR (1000'~m^2*")")) + 
    #coord_cartesian(ylim = c(0, max(maxes[,3])))+
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          plot.margin = unit(c(0, 0.2, 0, -.1), "cm"),
          axis.text.x = element_blank(),
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols)+ 
    scale_fill_manual(values=mycols) +
    annotate("text",x=.64,y=maxes[1,3],label=maxes[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.64,y=maxes[2,3],label=maxes[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.64,y=maxes[3,3],label=maxes[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=.82,y=maxes[4,3],label=maxes[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.82,y=maxes[5,3],label=maxes[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.82,y=maxes[6,3],label=maxes[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1,y=maxes[7,3],label=maxes[7,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes[8,3],label=maxes[8,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3,y=maxes[9,3],label=maxes[9,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.18,y=maxes[10,3],label=maxes[10,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.18,y=maxes[11,3],label=maxes[11,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3.18,y=maxes[12,3],label=maxes[12,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.36,y=maxes[13,3],label=maxes[13,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.36,y=maxes[14,3],label=maxes[14,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3.36,y=maxes[15,3],label=maxes[15,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    
    annotate("segment", x = .565, xend = .715, y = means[1,3], yend =  means[1,3],colour = "black",size=.3) +
    annotate("segment", x = 1.565, xend = 1.715, y = means[2,3], yend =  means[2,3],colour = "black",size=.3) +
    annotate("segment", x = 2.565, xend = 2.715, y = means[3,3], yend =  means[3,3],colour = "black",size=.3) +
    
    annotate("segment", x = .745, xend = .895, y = means[4,3], yend =  means[4,3],colour = "black",size=.3) +
    annotate("segment", x = 1.745, xend = 1.895, y = means[5,3], yend =  means[5,3],colour = "black",size=.3) +
    annotate("segment", x = 2.745, xend = 2.895, y = means[6,3], yend =  means[6,3],colour = "black",size=.3) +
    
    annotate("segment", x = .925, xend = 1.075, y = means[7,3], yend =  means[7,3],colour = "black",size=.3) +
    annotate("segment", x = 1.925, xend = 2.075, y = means[8,3], yend =  means[8,3],colour = "black",size=.3) +
    annotate("segment", x = 2.925, xend = 3.075, y = means[9,3], yend =  means[9,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.105, xend = 1.255, y = means[10,3], yend =  means[10,3],colour = "black",size=.3) +
    annotate("segment", x = 2.105, xend = 2.255, y = means[11,3], yend =  means[11,3],colour = "black",size=.3) +
    annotate("segment", x = 3.105, xend = 3.255, y = means[12,3], yend =  means[12,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.285, xend = 1.435, y = means[13,3], yend =  means[13,3],colour = "black",size=.3) +
    annotate("segment", x = 2.285, xend = 2.435, y = means[14,3], yend =  means[14,3],colour = "black",size=.3) +
    annotate("segment", x = 3.285, xend = 3.435, y = means[15,3], yend =  means[15,3],colour = "black",size=.3) +   ggtitle("(a)")
  
  test4<-test3[test3$scale=="100",]
  means2<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
  maxes2<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
  maxes2[,3]<-round(maxes2[,3])
  
  v2<-ggplot(aes(y = value, x = ecoregion, fill = class,colour=class), data = subset(test3,scale=="100")) + 
    geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
    #coord_cartesian(ylim = c(0, 100))+
    xlab("")+ ylab(bquote('SR (100'~m^2*")")) +
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.x = element_blank(),
          plot.margin = unit(c(-.3, 0.2, 0, 0.1), "cm"),
          axis.text.y = element_text(size = sz2,family="Times"))  +
    scale_colour_manual(values=mycols)+ scale_fill_manual(values=mycols) +
    
    annotate("text",x=.64,y=maxes2[1,3],label=maxes2[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.64,y=maxes2[2,3],label=maxes2[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.64,y=maxes2[3,3],label=maxes2[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=.82,y=maxes2[4,3],label=maxes2[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.82,y=maxes2[5,3],label=maxes2[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.82,y=maxes2[6,3],label=maxes2[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1,y=maxes2[7,3],label=maxes2[7,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes2[8,3],label=maxes2[8,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3,y=maxes2[9,3],label=maxes2[9,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.18,y=maxes2[10,3],label=maxes2[10,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.18,y=maxes2[11,3],label=maxes2[11,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3.18,y=maxes2[12,3],label=maxes2[12,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.36,y=maxes2[13,3],label=maxes2[13,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.36,y=maxes2[14,3],label=maxes2[14,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3.36,y=maxes2[15,3],label=maxes2[15,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    annotate("segment", x = .565, xend = .715, y = means2[1,3], yend =  means2[1,3],colour = "black",size=.3) +
    annotate("segment", x = 1.565, xend = 1.715, y = means2[2,3], yend =  means2[2,3],colour = "black",size=.3) +
    annotate("segment", x = 2.565, xend = 2.715, y = means2[3,3], yend =  means2[3,3],colour = "black",size=.3) +
    
    annotate("segment", x = .745, xend = .895, y = means2[4,3], yend =  means2[4,3],colour = "black",size=.3) +
    annotate("segment", x = 1.745, xend = 1.895, y = means2[5,3], yend =  means2[5,3],colour = "black",size=.3) +
    annotate("segment", x = 2.745, xend = 2.895, y = means2[6,3], yend =  means2[6,3],colour = "black",size=.3) +
    
    annotate("segment", x = .925, xend = 1.075, y = means2[7,3], yend =  means2[7,3],colour = "black",size=.3) +
    annotate("segment", x = 1.925, xend = 2.075, y = means2[8,3], yend =  means2[8,3],colour = "black",size=.3) +
    annotate("segment", x = 2.925, xend = 3.075, y = means2[9,3], yend =  means2[9,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.105, xend = 1.255, y = means2[10,3], yend =  means2[10,3],colour = "black",size=.3) +
    annotate("segment", x = 2.105, xend = 2.255, y = means2[11,3], yend =  means2[11,3],colour = "black",size=.3) +
    annotate("segment", x = 3.105, xend = 3.255, y = means2[12,3], yend =  means2[12,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.285, xend = 1.435, y = means2[13,3], yend =  means2[13,3],colour = "black",size=.3) +
    annotate("segment", x = 2.285, xend = 2.435, y = means2[14,3], yend =  means2[14,3],colour = "black",size=.3) +
    annotate("segment", x = 3.285, xend = 3.435, y = means2[15,3], yend =  means2[15,3],colour = "black",size=.3) +   ggtitle("(b)")
  
  test4<-test3[test3$scale=="001",]
  means3<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
  maxes3<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
  maxes3[,3]<-round(maxes3[,3],1)
  
  v3<-ggplot(aes(y = value, x = ecoregion, fill = class,colour=class), data = subset(test3,scale=="001")) + 
    geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
    #coord_cartesian(ylim = c(0, 3))+
    xlab("")+ ylab(bquote('SR (0.01'~m^2*")")) +
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          plot.margin = unit(c(-.3, 0.2, 0, 0.2), "cm"),
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.x = element_text(angle = 10,hjust = 1,size = 10,family="Times"),
          axis.text.y = element_text(size = sz2,family="Times")) + scale_colour_manual(values=mycols[c(1,3,4)]) + 
    scale_fill_manual(values=mycols[c(1,3,4)]) +
    
    annotate("text",x=.7,y=maxes3[1,3],label=maxes3[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.7,y=maxes3[2,3],label=maxes3[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.7,y=maxes3[3,3],label=maxes3[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1,y=maxes3[4,3],label=maxes3[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes3[5,3],label=maxes3[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3,y=maxes3[6,3],label=maxes3[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.3,y=maxes3[7,3],label=maxes3[7,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.3,y=maxes3[8,3],label=maxes3[8,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=3.3,y=maxes3[9,3],label=maxes3[9,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    
    annotate("segment", x = .5725, xend = .8275, y = means3[1,3], yend =  means3[1,3],colour = "black",size=.3) +
    annotate("segment", x = 1.5725, xend = 1.8275, y = means3[2,3], yend =  means3[2,3],colour = "black",size=.3) +
    annotate("segment", x = 2.5725, xend = 2.8275, y = means3[3,3], yend =  means3[3,3],colour = "black",size=.3) +
    
    annotate("segment", x = .8725, xend = 1.1275, y = means3[4,3], yend =  means3[4,3],colour = "black",size=.3) +
    annotate("segment", x = 1.875, xend = 2.1275, y = means3[5,3], yend =  means3[5,3],colour = "black",size=.3) +
    annotate("segment", x = 2.875, xend = 3.1275, y = means3[6,3], yend =  means3[6,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.1725, xend = 1.4275, y = means3[7,3], yend =  means3[7,3],colour = "black",size=.3) +
    annotate("segment", x = 2.1725, xend = 2.4275, y = means3[8,3], yend =  means3[8,3],colour = "black",size=.3) +
    annotate("segment", x = 3.1725, xend = 3.4275, y = means3[9,3], yend =  means3[9,3],colour = "black",size=.3) +   ggtitle("(c)")
  
  
    prow1 <- plot_grid(v1,v2,v3,hjust = 0,vjust=1,rel_heights=c(.8,.75,.88), nrow= 3,align="v")
    
    test<-cvs_full_data[,c(60,2,4,8,9,11,15,16,18,22,63,65,66,68)]   #23,51
    test2<-test %>% melt(id.vars=c("ecogeo6")) 
    test3<-cbind(test2[,1],as.data.frame(str_split_fixed(test2$variable, "_", 3))[,-2],test2[,3])
    head(test3)
    names(test3)<-c("ecoregion","class","scale","value")
    
    test3$class <- factor(test3$class, levels = c("all","ground","nw","w","tree"))
    test3$class<-revalue(test3$class, c("nw"="herbaceous", "w"="woody"))
    
    test3$ecoregion <- factor(test3$ecoregion, levels = c("Longleaf Woodlands","Coastal Forests","Riparian","Piedmont Forests","App Forests"))
    test3$ecoregion<-revalue(test3$ecoregion, c("Riparian"="Swamp Floodplains", "App Forests"="Appalachian Forests"))
    
    head(test3)
    library(RColorBrewer)
  
    test3$ecoregion<-as.character(test3$ecoregion)
    test3[test3$ecoregion!="Longleaf Woodlands",1]<-"Not Longleaf"
    
    test4<-test3[test3$scale=="1000",]
    means4<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
    maxes4<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
    maxes4[,3]<-round(maxes4[,3])
    
    v4<-ggplot(aes(y = value, x = factor(ecoregion), fill = class,colour=class), data = subset(test3,scale=="1000")) + 
      geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
      xlab("")+ ylab("")+ 
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          legend.text=element_text(size=14),
          plot.margin = unit(c(0, 0.2, 0, -.1), "cm"),
          axis.text.x = element_blank(),
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.y = element_text(size = sz2,family="Times")) + 
    scale_colour_manual(values=mycols)+ scale_fill_manual(values=mycols) +
    
    annotate("text",x=.64,y=maxes4[1,3],label=maxes4[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.64,y=maxes4[2,3],label=maxes4[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=.82,y=maxes4[3,3],label=maxes4[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.82,y=maxes4[4,3],label=maxes4[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1,y=maxes4[5,3],label=maxes4[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes4[6,3],label=maxes4[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.18,y=maxes4[7,3],label=maxes4[7,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.18,y=maxes4[8,3],label=maxes4[8,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.36,y=maxes4[9,3],label=maxes4[9,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.36,y=maxes4[10,3],label=maxes4[10,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    annotate("segment", x = .565, xend = .715, y = means4[1,3], yend =  means4[1,3],colour = "black",size=.3) +
    annotate("segment", x = 1.565, xend = 1.715, y = means4[2,3], yend =  means4[2,3],colour = "black",size=.3) +
    annotate("segment", x = .745, xend = .895, y = means4[3,3], yend =  means4[3,3],colour = "black",size=.3) +
    annotate("segment", x = 1.745, xend = 1.895, y = means4[4,3], yend =  means4[4,3],colour = "black",size=.3) +
    annotate("segment", x = .925, xend = 1.075, y = means4[5,3], yend =  means4[5,3],colour = "black",size=.3) +
    annotate("segment", x = 1.925, xend = 2.075, y = means4[6,3], yend =  means4[6,3],colour = "black",size=.3) +
    annotate("segment", x = 1.105, xend = 1.255, y = means4[7,3], yend =  means4[7,3],colour = "black",size=.3) +
    annotate("segment", x = 2.105, xend = 2.255, y = means4[8,3], yend =  means4[8,3],colour = "black",size=.3) +
    annotate("segment", x = 1.285, xend = 1.435, y = means4[9,3], yend =  means4[9,3],colour = "black",size=.3) +
    annotate("segment", x = 2.285, xend = 2.435, y = means4[10,3], yend =  means4[10,3],colour = "black",size=.3) +   ggtitle("(d)")

    test4<-test3[test3$scale=="100",]
    means5<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
    maxes5<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
    maxes5[,3]<-round(maxes5[,3])

  
  v5<-ggplot(aes(y = value, x = ecoregion, fill = class,colour=class), data = subset(test3,scale=="100")) + 
    geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
    xlab("")+ ylab("")+ 
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.x = element_blank(),
          plot.margin = unit(c(-.3, 0.2, 0, 0.1), "cm"),
          axis.text.y = element_text(size = sz2,family="Times"))  +
    scale_colour_manual(values=mycols)+ scale_fill_manual(values=mycols)+
    
    annotate("text",x=.64,y=maxes5[1,3],label=maxes5[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.64,y=maxes5[2,3],label=maxes5[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=.82,y=maxes5[3,3],label=maxes5[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.82,y=maxes5[4,3],label=maxes5[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1,y=maxes5[5,3],label=maxes5[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes5[6,3],label=maxes5[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.18,y=maxes5[7,3],label=maxes5[7,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.18,y=maxes5[8,3],label=maxes5[8,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.36,y=maxes5[9,3],label=maxes5[9,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.36,y=maxes5[10,3],label=maxes5[10,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    annotate("segment", x = .565, xend = .715, y = means5[1,3], yend =  means5[1,3],colour = "black",size=.3) +
    annotate("segment", x = 1.565, xend = 1.715, y = means5[2,3], yend =  means5[2,3],colour = "black",size=.3) +
    annotate("segment", x = .745, xend = .895, y = means5[3,3], yend =  means5[3,3],colour = "black",size=.3) +
    annotate("segment", x = 1.745, xend = 1.895, y = means5[4,3], yend =  means5[4,3],colour = "black",size=.3) +
    annotate("segment", x = .925, xend = 1.075, y = means5[5,3], yend =  means5[5,3],colour = "black",size=.3) +
    annotate("segment", x = 1.925, xend = 2.075, y = means5[6,3], yend =  means5[6,3],colour = "black",size=.3) +
    annotate("segment", x = 1.105, xend = 1.255, y = means5[7,3], yend =  means5[7,3],colour = "black",size=.3) +
    annotate("segment", x = 2.105, xend = 2.255, y = means5[8,3], yend =  means5[8,3],colour = "black",size=.3) +
    annotate("segment", x = 1.285, xend = 1.435, y = means5[9,3], yend =  means5[9,3],colour = "black",size=.3) +
    annotate("segment", x = 2.285, xend = 2.435, y = means5[10,3], yend =  means5[10,3],colour = "black",size=.3) +   ggtitle("(e)")
  
  test4<-test3[test3$scale=="001",]
  means6<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=median)
  maxes6<-aggregate(test4$value, by=list(test4$ecoregion,test4$class), FUN=max)
  maxes6[,3]<-round(maxes6[,3],1)
  
  v6<-ggplot(aes(y = value, x = ecoregion, fill = class,colour=class), data = subset(test3,scale=="001")) + 
    geom_boxplot(outlier.colour="grey",outlier.size=.1,position=position_dodge(.9))+
    xlab("")+ylab("")+ 
    theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz1),
          legend.position = "none",
          plot.margin = unit(c(-.3, 0.2, 0, 0.2), "cm"),
          axis.title=element_text(size = sz1,family="Times"),
          axis.text.x = element_text(angle = 10,hjust = 1,size = 10,family="Times"),
          axis.text.y = element_text(size = sz2,family="Times")) + 
    scale_colour_manual(values=mycols[c(1,3,6)]) + 
    scale_fill_manual(values=mycols[c(1,3,4)]) +
    
    annotate("text",x=.7,y=maxes6[1,3],label=maxes6[1,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.,y=maxes6[3,3],label=maxes6[3,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.3,y=maxes6[5,3],label=maxes6[5,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=1.7,y=maxes6[2,3],label=maxes6[2,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2,y=maxes6[4,3],label=maxes6[4,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    annotate("text",x=2.3,y=maxes6[6,3],label=maxes6[6,3],colour="grey20",hjust=0.5,size=annotate_size,family="Times") +
    
    annotate("segment", x = .5725, xend = .8275, y = means6[1,3], yend =  means6[1,3],colour = "black",size=.3) +
    annotate("segment", x = .8725, xend = 1.1275, y = means6[3,3], yend =  means6[3,3],colour = "black",size=.3) +
    annotate("segment", x = 1.1725, xend = 1.4275, y = means6[5,3], yend =  means6[5,3],colour = "black",size=.3) +
    
    annotate("segment", x = 1.5725, xend = 1.8275, y = means6[2,3], yend =  means6[2,3],colour = "black",size=.3) +
    annotate("segment", x = 1.8725, xend = 2.1275, y = means6[4,3], yend =  means6[4,3],colour = "black",size=.3) +
    annotate("segment", x = 2.1725, xend = 2.4275, y = means6[6,3], yend =  means6[6,3],colour = "black",size=.3) +   ggtitle("(f)")
  
  
  prow2 <- plot_grid(v4,v5,v6,hjust = .2,vjust=1,rel_heights=c(.8,.75,.88), nrow= 3,align="v")
  
  legend_b <- get_legend(v4 + theme(legend.title=element_blank(),
                                    legend.text=element_text(family="Times",size=11),
                                    legend.justification="center",
                                    plot.margin = unit(c(0, 0, -.5, 0), "cm"),
                                    legend.position="bottom"))
  
  
  prow3<-plot_grid(prow1,prow2, hjust=1,ncol = 2)
  
  pdf("Fig_8_28.pdf",height=6,width=8)
  print(plot_grid(prow3,legend_b, hjust=1,ncol = 1, rel_heights = c(1, .05)))
  dev.off()


# Fig 9 - all five variables by ecoregion ------
    
    rm(list=ls())
    library(ggplot2)
    library(cowplot)
    library(ecodist)
    library(quantreg)
    library(Hmisc)
    library(MASS)
    library(dplyr)
    library(plyr)
    library(grid)
    library(gridExtra)
    library(RColorBrewer)
    library(rstan)
    library(rstanarm)
    
    cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  # or the "all" one where SR is also scaled
    
    test<-cvs_full_data[cvs_full_data$ecogeo6=="Longleaf Woodlands",]
    table(test$ecogeo3)
    
    head(cvs_full_data[,c(2:4,61,50,57,60)])
    cvs_full_data[,c(61,50)]<-scale(cvs_full_data[,c(61,50)])
    
    poly_model_results<-as.data.frame(matrix(NA,28,29))
    names(poly_model_results)<-c("subset_type","plot_subset","scale","species_subset","soil_2.5","soil_25","soil_50","soil_75","soil_97.5",
                                 "cc_2.5","cc_25","cc_50","cc_75","cc_97.5",
                                 "soil_quad_2.5","soil_quad_25","soil_quad_50","soil_quad_75","soil_quad_97.5",
                                 "cc_quad_2.5","cc_quad_25","cc_quad_50","cc_quad_75","cc_quad_97.5",
                                 "soil_cc_2.5","soil_cc_25","soil_cc_50","soil_cc_75","soil_cc_97.5")
    
    
    g<-1;i<-7;eco3<-"Piedmont Forests"
    for (i in c(2:8)){
      tmp<-cvs_full_data[,c(i,61,50,57,60)]
      tmp<-tmp[complete.cases(tmp),]
      names(tmp)[1:3]<-c("SR","SF","CC")
      
      # eco3  
      for (eco3 in (unique(cvs_full_data$ecogeo3))){
        tmp2<-tmp[tmp$ecogeo3==eco3,]
        tmp2<-tmp2[complete.cases(tmp2),]
        tmp2[,1]<-(tmp2[,1]-min(tmp2[,1])) /(max(tmp2[,1])-min(tmp2[,1]))
        poly_model<-stan_glm(SR~SF*CC+I(SF^2)+I(CC^2), data = tmp2, QR = TRUE,chains = 4, iter = 10000) 
        beta_hat_CIs<-as.data.frame(summary(poly_model))[2:6,4:8]
        poly_model_results[g,1]<-"ecogeo3"
        poly_model_results[g,2]<-eco3
        poly_model_results[g,3]<-strsplit(names(cvs_full_data)[i], "_")[[1]][3]
        poly_model_results[g,4]<-strsplit(names(cvs_full_data)[i], "_")[[1]][1]
        poly_model_results[g,5:29]<-cbind(beta_hat_CIs[1,],beta_hat_CIs[2,],beta_hat_CIs[3,],beta_hat_CIs[4,],beta_hat_CIs[5,])
        print(i);g<-g+1
      }
    }
    
    cvs_full_data[,c(61,50)]<-scale(cvs_full_data[,c(61,50)])
    
    for (i in c(2:8)){
      tmp<-cvs_full_data[,c(i,61,50,57,60)]
      tmp<-tmp[complete.cases(tmp),]
      names(tmp)[1:3]<-c("SR","SF","CC")
      tmp2<-tmp[tmp$ecogeo6=="Longleaf Woodlands",]
      tmp2[,1]<-(tmp2[,1]-min(tmp2[,1])) /(max(tmp2[,1])-min(tmp2[,1]))
      poly_model<-stan_glm(SR~SF*CC+I(SF^2)+I(CC^2), data = tmp2, QR = TRUE,chains = 4, iter = 10000) 
      beta_hat_CIs<-as.data.frame(summary(poly_model))[2:6,4:8]
      poly_model_results[g,1]<-"ecogeo8"
      poly_model_results[g,2]<-"Longleaf Woodlands"
      poly_model_results[g,3]<-strsplit(names(cvs_full_data)[i], "_")[[1]][3]
      poly_model_results[g,4]<-strsplit(names(cvs_full_data)[i], "_")[[1]][1]
      poly_model_results[g,4]<-strsplit(names(cvs_full_data)[i], "_")[[1]][1]
      poly_model_results[g,5:29]<-cbind(beta_hat_CIs[1,],beta_hat_CIs[2,],beta_hat_CIs[3,],beta_hat_CIs[4,],beta_hat_CIs[5,])
      print(i);g<-g+1
    }
    
    poly_model_results$plot_subset <- factor(poly_model_results$plot_subset, levels = c("App Forests","Piedmont Forests","Coastal Forests","Longleaf Woodlands"))
    poly_model_results$plot_subset<-revalue(poly_model_results$plot_subset, c("App Forests"="Appalachian Forests", "Coastal Forests"="Coastal Forests"))
    poly_model_results$scale <- factor(poly_model_results$scale, levels = rev(c("1000","400","100","10","1","01","001")))
    poly_model_results$scale<-revalue(poly_model_results$scale, c("01"="0.1", "001"="0.01"))
    write.csv(poly_model_results,"...poly_model_results_ecoregion.csv")
    poly_model_results<-read.csv("...poly_model_results_ecoregion.csv")
    poly_model_results<-poly_model_results[,-1]
    
    poly_model_results$soil_cc_sig<-poly_model_results$cc_quad_sig<-poly_model_results$soil_quad_sig<-poly_model_results$cc_sig<-poly_model_results$soil_sig<-as.character(poly_model_results$species_subset)
    i<-1
    
    for (i in 1:dim(poly_model_results)[1]){
      if (poly_model_results[i,5]<0 & poly_model_results[i,9]>0) {poly_model_results[i,30]<-"non"}
      if (poly_model_results[i,10]<0 & poly_model_results[i,14]>0) {poly_model_results[i,31]<-"non"}
      if (poly_model_results[i,15]<0 & poly_model_results[i,19]>0) {poly_model_results[i,32]<-"non"}
      if (poly_model_results[i,20]<0 & poly_model_results[i,24]>0) {poly_model_results[i,33]<-"non"}
      if (poly_model_results[i,25]<0 & poly_model_results[i,29]>0) {poly_model_results[i,34]<-"non"}
    }
    
    
    test<-cbind(poly_model_results[,2:4],poly_model_results[,c(5,7,9)],poly_model_results[30])
    test[,3]<-"SF"
    colnames(test)[4:7]<-c("my025","my50","my975","sig")
    
    test2<-cbind(poly_model_results[,2:4],poly_model_results[,c(10,12,14)],poly_model_results[31])
    test2[,3]<-"CC"
    colnames(test2)[4:7]<-c("my025","my50","my975","sig")
    
    test3<-cbind(poly_model_results[,2:4],poly_model_results[,c(15,17,19)],poly_model_results[32])
    test3[,3]<-"SF2"
    colnames(test3)[4:7]<-c("my025","my50","my975","sig")
    
    test4<-cbind(poly_model_results[,2:4],poly_model_results[,c(20,22,24)],poly_model_results[33])
    test4[,3]<-"CC2"
    colnames(test4)[4:7]<-c("my025","my50","my975","sig")
    
    test5<-cbind(poly_model_results[,2:4],poly_model_results[,c(25,27,29)],poly_model_results[34])
    test5[,3]<-"SF:CC"
    colnames(test5)[4:7]<-c("my025","my50","my975","sig")
    
    
    nonsigs<-rbind(test,test2,test3,test4,test5)
    
    nonsigs[which(nonsigs$sig=="all"),7]<-nonsigs[which(nonsigs$sig=="all"),3]
    names(nonsigs)[3]<-"variable"
    
    nonsigs$sig <- factor(nonsigs$sig , levels = c("SF","SF2","CC","CC2","SF:CC","non"))
    nonsigs$variable <- factor(nonsigs$variable , levels = c("SF","SF2","CC","CC2","SF:CC","non"))
    
    
    library(RColorBrewer)
    
    mycols<-c(brewer.pal(n = 9, name = "YlOrRd")[7], #red
              brewer.pal(n = 9, name = "YlOrRd")[4],# orange
              brewer.pal(n = 9, name = "Blues")[4], #blue
              brewer.pal(n = 11, name = "PiYG")[10], #green
              brewer.pal(n = 11, name = "Spectral")[11]) #purple
    
    ecogeo<-"Piedmont Forests"
    myplot<-list()
    mynum<-1
    sz<-13
    
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(nonsigs,plot_subset==ecogeo), aes(x = as.factor(scale),y=my50,colour=factor(variable),shape=factor(variable),fill=sig)) + 
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=my025, ymax=my975), width=0,position=position_dodge(width = -0.6)) +
        geom_point(size=1.5,position=position_dodge(width = -0.6)) +
        scale_fill_manual(values = c("non"="white","SF"=mycols[1],"SF2"=mycols[2],"CC"=mycols[3],"CC2"=mycols[4],"SF:CC"=mycols[5])) +
        guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:25,fill=mycols))) +   
        scale_shape_manual(values = c(21:25)) +
        scale_colour_manual(values=mycols,labels=c("SF",expression(paste(SF^2)),"CC",expression(paste(CC^2)),"SF:CC")) + 
        theme(plot.title=element_text(hjust=-0.12,face="bold",vjust=-1.5,family="Times",size=sz),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text=element_text(family="Times",size=sz),
              axis.title=element_text(family="Times",size=sz),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        xlab(bquote('scale ('~m^2~ ")")) +
        ylab(ecogeo)+labs(colour="")
      
      mynum<-mynum+1
    }
    
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               legend.text=element_text(family="Times",size=sz),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    plot_grid( myplot[[3]], legend_b, ncol = 1, rel_heights = c(1, .07))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none") +      ggtitle("(a)"),
                      myplot[[2]] + theme(legend.position="none")+xlab("")+      ggtitle("(b)"),
                      myplot[[1]] + theme(legend.position="none")+xlab("")+      ggtitle("(c)"),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=4),
                      myplot[[4]] + theme(legend.position="none")+      ggtitle("(d)"),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    
    pdf(paste("...Fig_9.pdf",sep=""),height=5,width=11)
    prow2
    dev.off()


# Appendix 1 plot design---------

  # manually drawn in pptx to scale

# Appendix 2 Species richness distribution across all plots, split by growth form and strata--------
  dev.off()
  rm(list=ls())
  library(ecodist)
  library(quantreg)
  library(Hmisc)
  library(MASS)
  library(cowplot)
  library(RColorBrewer)
  library(reshape)
  
  cvs_full_data<-read.csv("...cvs_full_data_5cm.csv") 
  setwd("...figs_tabs")
  table(cvs_full_data$ecogeo6)
  
  nw_model<-lm(all_sr_1000~nw_sr_1000,data=cvs_full_data)
  summary(nw_model)
  
  ground_model<-lm(all_sr_1000~ground_sr_1000,data=cvs_full_data)
  summary(ground_model)
  
  w_model<-lm(all_sr_1000~w_sr_1000,data=cvs_full_data)
  summary(w_model)
  
  tree_model<-lm(all_sr_1000~tree_sr_1000,data=cvs_full_data)
  summary(tree_model)
  
  test<-cvs_full_data[,c(2,9,16,63,66)]
  sz<-10
  mycols<-c("chartreuse4", brewer.pal(n = 11, name = "PiYG")[8],brewer.pal(n = 11, name = "Spectral")[c(4,2)]) #11 before the 4
  
  # wnw
    all_ws<- aggregate(test$nw_sr_1000, by=list(Category=test$all_sr_1000), FUN=sum)
    all_nws<-aggregate(test$w_sr_1000, by=list(Category=test$all_sr_1000), FUN=sum)
    mat<-merge(all_ws,all_nws,by="Category",all = TRUE)
    mat$sum<-mat[,2]+mat[,3]
    all_all<-as.data.frame(table(cvs_full_data$all_sr_1000))
    names(all_all)<-c("Category","Freq")
    mat2<-merge(mat,all_all,by="Category",all = TRUE)
    
    mat2$nonwoody<-mat2[,2]/mat2$sum*mat2$Freq
    mat2$woody<-mat2[,3]/mat2$sum*mat2$Freq
    mat3<-mat2[,c(1,6,7)]
    mat4 <- melt(mat3, id.var="Category")
    mat4$value<-mat4$value/sum(mat4$value)
    head(mat4)
    mat4$variable<-revalue(mat4$variable, c("nonwoody"="herbaceous", "w"="woody"))
    
    wnw<-ggplot(mat4, aes(x = Category, y = value, fill = variable)) + 
      theme(text=element_text(size=sz,family = "Times"),
            legend.position = c(0.5, .75),
            axis.text=element_text(size=sz,family = "Times")) +
      geom_bar(width=1,stat = "identity") + xlab("SR") + ylab("Frequency") +
      guides(fill=guide_legend(title="Growth Form"))+
      scale_fill_manual(values=mycols[c(2,3)]) + scale_y_continuous(labels = scales::percent_format(accuracy = .5))

  # tree-ground
    all_ws<- aggregate(test$ground_sr_1000, by=list(Category=test$all_sr_1000), FUN=sum)
    all_nws<-aggregate(test$tree_sr_1000, by=list(Category=test$all_sr_1000), FUN=sum)
    mat<-merge(all_ws,all_nws,by="Category",all = TRUE)
    mat$sum<-mat[,2]+mat[,3]
    all_all<-as.data.frame(table(cvs_full_data$all_sr_1000))
    names(all_all)<-c("Category","Freq")
    mat2<-merge(mat,all_all,by="Category",all = TRUE)
    
    mat2$ground<-mat2[,2]/mat2$sum*mat2$Freq
    mat2$tree<-mat2[,3]/mat2$sum*mat2$Freq
    mat3<-mat2[,c(1,6,7)]
    mat4 <- melt(mat3, id.var="Category")
    mat4$value<-mat4$value/sum(mat4$value)
    
    head(mat4)
    tg<-ggplot(mat4, aes(x = Category, y = value, fill = variable)) + 
      theme(text=element_text(size=sz,family = "Times"),
            legend.position = c(0.5, .75),
            legend.background = element_rect(colour = 'black', fill = 'white', size = 1),
            axis.text=element_text(size=sz,family = "Times")) +
      geom_bar(width=1,stat = "identity") + xlab("SR") + ylab("Frequency") +
      guides(fill=guide_legend(title="Strata"))+
      scale_fill_manual(values=mycols[c(1,4)]) + scale_y_continuous(labels = scales::percent_format(accuracy = .5))

    prow <- plot_grid(wnw,tg,labels = c("   (a)", "   (b)"),hjust = 0, nrow= 1,align="v")
    
    pdf("append_4_histograms.pdf",height=3,width=6)
    print(prow)
    dev.off()

# Appendix 3 compare AICs---------
  dev.off()
  rm(list=ls())
  library(ecodist)
  library(quantreg)
  library(Hmisc)
  library(MASS)
  library(AER)
  
  setwd("....5_interaction_models")
  cvs_full_data<-read.csv("...vs_full_data_5cm.csv")  # or the "all" one where SR is also scaled
  names(cvs_full_data)
  
  library(MuMIn)
  glm1<-glm(all_sr_1000~soil_fert1*cc+I(soil_fert1^2)+I(cc^2),data=cvs_full_data,na.action = "na.fail")
  
  test<-dredge(glm1)
  test
  
  my_table<-as.data.frame(matrix(NA,27,5))
  names(my_table)<-c("response","scale","model","adj R2","sig vars")
  g<-1;i<-21;vari<-3
  for (vari in 1:27){
    tmp<-tmp1[,c(vari,28,29)]
    
    
    a<-lm(tmp[,1]~tmp[,2])
    b<-lm(tmp[,1]~I(tmp[,2]^2))
    c<-lm(tmp[,1]~tmp[,2]+I(tmp[,2]^2))
    d<-lm(tmp[,1]~tmp[,3])
    e<-lm(tmp[,1]~I(tmp[,3]^2))
    f<-lm(tmp[,1]~tmp[,3]+I(tmp[,3]^2))
    g<-lm(tmp[,1]~tmp[,2]*tmp[,3])
    h<-lm(tmp[,1]~tmp[,2]*tmp[,3]+I(tmp[,2]^2))
    i<-lm(tmp[,1]~tmp[,2]*tmp[,3]+I(tmp[,3]^2))
    j<-lm(tmp[,1]~tmp[,2]*tmp[,3]+I(tmp[,2]^2)+I(tmp[,3]^2))
    k<-lm(tmp[,1]~I(tmp[,2]^2)+I(tmp[,3]^2))
    l<-lm(tmp[,1]~I(tmp[,2]^2)*I(tmp[,3]^2))
    
    
    AIC_results<-AIC(a,b,c,d,e,f,g,h,i,j,k,l)
    my_table[vari,1]<-strsplit(names(tmp)[1], "_")[[1]][1]
    my_table[vari,2]<-strsplit(names(tmp)[1], "_")[[1]][3]
    my_table[vari,3]<-print(row.names(AIC_results)[which.min(AIC_results$AIC)])
    if (my_table[vari,3]=="h"){
      my_table[vari,4]<-round(summary(h)$adj.r.squared,2)
      my_table[vari,5]<-paste(row.names(summary(h)$coeff)[-1][summary(h)$coeff[-1,4] < 0.05],collapse=", ")
    }
    if (my_table[vari,3]=="i"){ 
      my_table[vari,4]<-round(summary(i)$adj.r.squared,2)
      my_table[vari,5]<-paste(row.names(summary(i)$coeff)[-1][summary(i)$coeff[-1,4] < 0.05],collapse=", ")
    }
    if (my_table[vari,3]=="j"){ 
      my_table[vari,4]<-round(summary(j)$adj.r.squared,2)
      my_table[vari,5]<-paste(row.names(summary(j)$coeff)[-1][summary(j)$coeff[-1,4] < 0.05],collapse=", ")
    }
    
  }
  
  write.csv(my_table,"...AIC_R2.csv",row.names=F)

# Appendix 4. Each of the five betas by ecogeo ------

    rm(list=ls())
    library(ggplot2)
    library(cowplot)
    library(ecodist)
    library(quantreg)
    library(Hmisc)
    library(MASS)
    library(dplyr)
    library(plyr)
    library(grid)
    library(gridExtra)
    library(RColorBrewer)
    
    cvs_full_data<-read.csv("...cvs_full_data_5cm.csv")  # or the "all" one where SR is also scaled
    
    test<-cvs_full_data[cvs_full_data$ecogeo6=="Longleaf Woodlands",]
    table(test$ecogeo3)
    
    head(cvs_full_data[,c(2:4,61,50,57,60)])
    cvs_full_data[,c(61,50)]<-scale(cvs_full_data[,c(61,50)])
    
    poly_model_results<-as.data.frame(matrix(NA,(27*4),19))
    names(poly_model_results)<-c("subset_type","plot_subset","scale","species_subset","soil_2.5","soil_50","soil_97.5",
                                 "cc_2.5","cc_50","cc_97.5",
                                 "soil_quad_2.5","soil_quad_50","soil_quad_97.5",
                                 "cc_quad_2.5","cc_quad_50","cc_quad_97.5",
                                 "soil_cc_2.5","soil_cc_50","soil_cc_97.5")
    
    g<-1;i<-7;eco3<-"Longleaf Woodlands"
    for (i in c(2:22,63:68)){
      tmp<-cvs_full_data[,c(i,61,50,57,60)]
      tmp<-tmp[complete.cases(tmp),]
      
    # eco3  
      for (eco3 in (unique(cvs_full_data$ecogeo3))){
        tmp2<-tmp[tmp$ecogeo3==eco3,]
        tmp2<-tmp2[complete.cases(tmp2),]
        tmp2[,1]<-(tmp2[,1]-min(tmp2[,1])) /(max(tmp2[,1])-min(tmp2[,1]))
        poly_model<-lm(tmp2[,1]~tmp2[,2]*tmp2[,3]+I(tmp2[,2]^2)+I(tmp2[,3]^2))
        beta_hats<-summary(poly_model)$coefficients[,1]
        se<-summary(poly_model)$coefficients[,2]
        beta_hat_CIs<- cbind(beta_hats-1.96*se,beta_hats,beta_hats+1.96*se)
        poly_model_results[g,1]<-"ecogeo3"
        poly_model_results[g,2]<-eco3
        poly_model_results[g,3]<-strsplit(names(cvs_full_data)[i], "_")[[1]][3]
        poly_model_results[g,4]<-strsplit(names(cvs_full_data)[i], "_")[[1]][1]
        poly_model_results[g,5:19]<-cbind(beta_hat_CIs[2,],beta_hat_CIs[3,],beta_hat_CIs[4,],beta_hat_CIs[5,],beta_hat_CIs[6,])
        g<-g+1
      }
    }
    
    cvs_full_data[,c(61,50)]<-scale(cvs_full_data[,c(61,50)])
    
    for (i in c(2:22,63:68)){
      tmp<-cvs_full_data[,c(i,61,50,57,60)]
      tmp<-tmp[complete.cases(tmp),]
      tmp2<-tmp[tmp$ecogeo6=="Longleaf Woodlands",]
      tmp2[,1]<-(tmp2[,1]-min(tmp2[,1])) /(max(tmp2[,1])-min(tmp2[,1]))
      poly_model<-lm(tmp2[,1]~tmp2[,2]*tmp2[,3]+I(tmp2[,2]^2)+I(tmp2[,3]^2))
      beta_hats<-summary(poly_model)$coefficients[,1]
      se<-summary(poly_model)$coefficients[,2]
      beta_hat_CIs<- cbind(beta_hats-1.96*se,beta_hats,beta_hats+1.96*se)
      poly_model_results[g,1]<-"ecogeo8"
      poly_model_results[g,2]<-"Longleaf Woodlands"
      poly_model_results[g,3]<-strsplit(names(cvs_full_data)[i], "_")[[1]][3]
      poly_model_results[g,4]<-strsplit(names(cvs_full_data)[i], "_")[[1]][1]
      poly_model_results[g,5:19]<-cbind(beta_hat_CIs[2,],beta_hat_CIs[3,],beta_hat_CIs[4,],beta_hat_CIs[5,],beta_hat_CIs[6,])
      g<-g+1
    }
    
    poly_model_results$species_subset <- factor(poly_model_results$species_subset, levels = c( "ground","nw","all","w","tree"))
    poly_model_results$species_subset<-revalue(poly_model_results$species_subset, c("nw"="herbaceous", "w"="woody","all"="total"))
    poly_model_results$plot_subset <- factor(poly_model_results$plot_subset, levels = c("App Forests","Piedmont Forests","Coastal Forests","Longleaf Woodlands"))
    poly_model_results$plot_subset<-revalue(poly_model_results$plot_subset, c("App Forests"="Appalachian Forests", "Coastal Forests"="Coastal Forests"))
    poly_model_results$scale <- factor(poly_model_results$scale, levels = rev(c("1000","400","100","10","1","01","001")))
    poly_model_results$scale<-revalue(poly_model_results$scale, c("01"="0.1", "001"="0.01"))
    
    
    poly_model_results$soil_cc_sig<-poly_model_results$cc_quad_sig<-poly_model_results$soil_quad_sig<-poly_model_results$cc_sig<-poly_model_results$soil_sig<-as.character(poly_model_results$species_subset)
    i<-1
    
    for (i in 1:dim(poly_model_results)[1]){
      if (poly_model_results[i,5]<0 & poly_model_results[i,7]>0) {poly_model_results[i,20]<-"non"}
      if (poly_model_results[i,8]<0 & poly_model_results[i,10]>0) {poly_model_results[i,21]<-"non"}
      if (poly_model_results[i,11]<0 & poly_model_results[i,13]>0) {poly_model_results[i,22]<-"non"}
      if (poly_model_results[i,14]<0 & poly_model_results[i,16]>0) {poly_model_results[i,23]<-"non"}
      if (poly_model_results[i,17]<0 & poly_model_results[i,19]>0) {poly_model_results[i,24]<-"non"}
    }
    
    
    sz=10
    mycols<-c("chartreuse4", brewer.pal(n = 11, name = "PiYG")[8],"deepskyblue3",brewer.pal(n = 11, name = "Spectral")[c(4,2)])
    ecogeo<-"Piedmont Forests"
    myplot<-list()
    mynum<-1
    
  # SF
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(poly_model_results,plot_subset==ecogeo), aes(x = as.factor(scale),y=soil_50,colour=factor(species_subset),fill=soil_sig)) +
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=soil_2.5, ymax=soil_97.5), width=0,position=position_dodge(width = -0.5)) +
        geom_point(position=position_dodge(width = -0.5),shape=21) +
        scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
        guides(fill=FALSE,colour = guide_legend(override.aes = list(shape =19))) +    
        scale_colour_manual(values=mycols) + ggtitle(ecogeo) +
        theme(plot.title=element_text(hjust=0.5,family="Times",size=11),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(family="Times"),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        xlab(bquote('scale ('~m^2~ ")")) +
        ylab(expression(paste("slope (", beta[1],")")))+labs(colour="Species")  
      mynum<-mynum+1
    }
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none"),
                      myplot[[2]] + theme(legend.position="none")+xlab(""),
                      myplot[[1]] + theme(legend.position="none")+xlab(""),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=2),
                      myplot[[4]] + theme(legend.position="none"),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),labels = c("(a)", "(b)", "(c)","","(d)"),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    y.grob <- textGrob("SR ~ SF", gp=gpar(family="Times", fontsize=15), rot=90)
    
    pdf(paste("...App_4_SF.pdf",sep=""),height=4,width=11)
    print(grid.arrange(arrangeGrob(prow2, left = y.grob)))
    dev.off()
    
    
  # SF2
    myplot<-list()  ;mynum<-1
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(poly_model_results,plot_subset==ecogeo), aes(x = as.factor(scale),y=soil_quad_50,colour=factor(species_subset),fill=soil_quad_sig)) +
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=soil_quad_2.5, ymax=soil_quad_97.5), width=0,position=position_dodge(width = -0.5)) +
        geom_point(position=position_dodge(width = -0.5),shape=21) +
        scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
        guides(fill=FALSE,colour = guide_legend(override.aes = list(shape =19))) +    
        scale_colour_manual(values=mycols) + 
        ggtitle(ecogeo)+xlab(bquote('scale ('~m^2~ ")")) +
        theme(plot.title=element_text(hjust=0.5,family="Times",size=11),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(family="Times"),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        ylab(expression(paste("slope (", beta[2],")")))+labs(colour="Species")  
      mynum<-mynum+1
    }
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none"),
                      myplot[[2]] + theme(legend.position="none")+xlab(""),
                      myplot[[1]] + theme(legend.position="none")+xlab(""),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=2),
                      myplot[[4]] + theme(legend.position="none")+xlab(""),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),labels = c("(a)", "(b)", "(c)","","(d)"),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    unit <- quote(SF^2)
    y.grob <- textGrob( bquote('SR ~ ' ~ .(unit)), gp=gpar(family="Times", fontsize=15), rot=90)
    
    pdf(paste("...App_4_SF2.pdf",sep=""),height=4,width=11)
    print(grid.arrange(arrangeGrob(prow2, left = y.grob)))
    dev.off()
    
    
  # cc 
    myplot<-list();mynum<-1
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(poly_model_results,plot_subset==ecogeo), aes(x = as.factor(scale),y=cc_50,colour=factor(species_subset),fill=cc_sig)) +
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=cc_2.5, ymax=cc_97.5), width=0,position=position_dodge(width = -0.5)) +
        geom_point(position=position_dodge(width = -0.5),shape=21) +
        scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
        guides(fill=FALSE,colour = guide_legend(override.aes = list(shape =19))) +    
        scale_colour_manual(values=mycols) + 
        ggtitle(ecogeo) + xlab(bquote('scale ('~m^2~ ")")) +
        theme(plot.title=element_text(hjust=0.5,family="Times",size=11),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(family="Times"),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        ylab(expression(paste("slope (", beta[3],")")))+labs(colour="Species")  
      mynum<-mynum+1
    }
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none"),
                      myplot[[2]] + theme(legend.position="none")+xlab(""),
                      myplot[[1]] + theme(legend.position="none")+xlab(""),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=2),
                      myplot[[4]] + theme(legend.position="none")+xlab(""),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),labels = c("(a)", "(b)", "(c)","","(d)"),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    y.grob <- textGrob("SR ~ CC", gp=gpar(family="Times", fontsize=15), rot=90)
    
    pdf(paste("...App_4_CC.pdf",sep=""),height=4,width=11)
    print(grid.arrange(arrangeGrob(prow2, left = y.grob)))
    dev.off()
    
  # CC2
    myplot<-list();mynum<-1
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(poly_model_results,plot_subset==ecogeo), aes(x = as.factor(scale),y=cc_quad_50,colour=factor(species_subset),fill=cc_quad_sig)) +
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=cc_quad_2.5, ymax=cc_quad_97.5), width=0,position=position_dodge(width = -0.5)) +
        geom_point(position=position_dodge(width = -0.5),shape=21) +
        scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
        guides(fill=FALSE,colour = guide_legend(override.aes = list(shape =19))) +    
        scale_colour_manual(values=mycols) + 
        ggtitle(ecogeo)+ xlab(bquote('scale ('~m^2~ ")")) +
        theme(plot.title=element_text(hjust=0.5,family="Times",size=11),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(family="Times"),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        ylab(expression(paste("slope (", beta[4],")")))+labs(colour="Species") 
      mynum<-mynum+1
    }
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none"),
                      myplot[[2]] + theme(legend.position="none")+xlab(""),
                      myplot[[1]] + theme(legend.position="none")+xlab(""),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=2),
                      myplot[[4]] + theme(legend.position="none")+xlab(""),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),labels = c("(a)", "(b)", "(c)","","(d)"),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    unit <- quote(CC^2)
    y.grob <- textGrob( bquote('SR ~ ' ~ .(unit)), gp=gpar(family="Times", fontsize=15), rot=90)
    
    pdf(paste("...App_4_CC2.pdf",sep=""),height=4,width=11)
    print(grid.arrange(arrangeGrob(prow2, left = y.grob)))
    dev.off()
    
    
  # SF:CC
    myplot<-list();mynum<-1
    for (ecogeo in unique(poly_model_results$plot_subset)){
      
      myplot[[mynum]]<-ggplot(subset(poly_model_results,plot_subset==ecogeo), aes(x = as.factor(scale),y=soil_cc_50,colour=factor(species_subset),fill=soil_cc_sig)) +
        theme_bw()+geom_hline(yintercept = 0, linetype = "longdash")  +  coord_flip()+ 
        geom_errorbar(aes(ymin=soil_cc_2.5, ymax=soil_cc_97.5), width=0,position=position_dodge(width = -0.5)) +
        geom_point(position=position_dodge(width = -0.5),shape=21) +
        scale_fill_manual(values = c("non"="white","ground"=mycols[1],"herbaceous"=mycols[2],"total"=mycols[3],"woody"=mycols[4],"tree"=mycols[5])) +
        guides(fill=FALSE,colour = guide_legend(override.aes = list(shape =19))) +    
        scale_colour_manual(values=mycols) + 
        ggtitle(ecogeo)+ xlab(bquote('scale ('~m^2~ ")")) +
        theme(plot.title=element_text(hjust=0.5,family="Times",size=11),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(family="Times"),
              axis.text.y=element_text(family="Times"),
              axis.text.x=element_text(family="Times",angle=25, hjust = 1))+
        ylab(expression(paste("slope (", beta[5],")")))+labs(colour="Species")  
      mynum<-mynum+1
    }
    
    legend_b <- get_legend(myplot[[1]] + theme(legend.title=element_blank(),
                                               plot.margin=unit(c(0,0,.5,0), "cm"),
                                               legend.position="bottom"))
    
    prow <- plot_grid(myplot[[3]] + theme(legend.position="none"),
                      myplot[[2]] + theme(legend.position="none")+xlab(""),
                      myplot[[1]] + theme(legend.position="none")+xlab(""),
                      ggdraw() + draw_line(x=c(0.5,0.5),y=c(0,1),color="black",size=2),
                      myplot[[4]] + theme(legend.position="none")+xlab(""),
                      align = 'vh',rel_widths=c(1,1,1,.3,1),labels = c("(a)", "(b)", "(c)","","(d)"),hjust = -1, nrow = 1)
    
    prow2<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .07))
    y.grob <- textGrob("SR ~ SF:CC", gp=gpar(family="Times", fontsize=15), rot=90)
    
    pdf(paste("...App_4_SF_CC.pdf",sep=""),height=4,width=11)
    print(grid.arrange(arrangeGrob(prow2, left = y.grob)))
    dev.off()

# Appendix 5  ---------
    
    #list of plots