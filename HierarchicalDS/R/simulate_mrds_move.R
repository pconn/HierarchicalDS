#' simulate MRDS double observer data subject to movement and measurement error; assumes equal area bins
#' @param n_species Number of species to simulate data for
#' @param n_bins Number of bins for complete data (note that measurement error can result in some of the animals outside of strip width being detected)
#' @param Obs_bins A vector giving which of the n_bins are within the strip width
#' @param measure_par Gives measurement error precision (if gaussian = TRUE) or measurement error rate (if gaussian=FALSE)
#' @param move_par Gives movement precision (if gaussian = TRUE) or rate (if gaussian=FALSE)
#' @param gaussian If TRUE (default), uses a Gaussian distribution for measurement error and a half-normal for movement; if FALSE (default), uses exponential / half exponential (Laplace dist)
#' @param one_side If TRUE (default), assumes movement can only be away from the center line
#' @param first_obs_correct If TRUE, measurement error for observer 1 is assumed to be zero (default is FALSE)
#' @return a community mrds dataset
#' @export
#' @keywords simulation, mrds
#' @author Paul B. Conn
simulate_mrds_move <- function(n_species=4,n_bins=10,Obs_bins=c(2,3,4,5,6),measure_par=2.0,move_par=0.75,gaussian=TRUE,one_side=TRUE,first_obs_correct=FALSE){
  N<-round(runif(n_species,19.5,100.5))
  n_obs_bins = length(Obs_bins)
  
  Beta_det=c(1,-0.2,-0.2*c(1:(n_obs_bins-1)),0.02,0.5,rnorm(n_species-1,0,0.25))
   #obs 1 (bin 1), obs 2 offset, offset for bin 2, ..., offset for bin n.bins, grp size,fly,species
  cor_par=1.0 #correlation in max distance bin (linear from zero)
  
  
  ### group model - poisson-normal mixture
  Mu_group=rnorm(n_species,log(1),0.1)
  sigma_grp=0.000001
  #sigma_grp = 1.2
  
  ### measurement error model
  #measure_par = 1.5  #par parameter for exponential distribution

  ### proportion flying
  Prop_fly = runif(n_species,0.4,0.8)
  
  Complete_data = matrix(0,2000,9)
  colnames(Complete_data) = c("species","d1_true","d2_true","d1_obs","d2_obs","g_size","fly","det1","det2")
  counter=1
  
  Prob_measure = Prob_move = matrix(0,n_bins,n_bins) 
  if(gaussian==TRUE){
    for(ibin in 1:n_bins){
      Bin_distance = rep(1:n_bins)-ibin
      Prob_measure[ibin,] = dnorm(Bin_distance,0,sqrt(1/measure_par))
      Prob_move[ibin,]=dnorm(Bin_distance,0,sqrt(1/move_par))
      if(one_side==TRUE & ibin>1)Prob_move[ibin,1:(ibin-1)]=0
    }   
  }
  else{  #double exponential
    for(ibin in 1:n_bins){
      Bin_distance = rep(1:n_bins)-ibin
      Prob_measure[ibin,] = dexp(abs(Bin_distance),measure_par) #double exponential random variates
      Prob_move[ibin,] = dexp(abs(Bin_distance),move_par)
      if(one_side==TRUE & ibin>1)Prob_move[ibin,1:(ibin-1)]=0
    }
  }
  Prob_measure = Prob_measure/apply(Prob_measure,1,'sum')
  Prob_move = Prob_move/apply(Prob_move,1,'sum')
  
  Sigma = diag(2)
  for(isp in 1:n_species){
    Cur_rows=counter:(counter+N[isp]-1)
    Complete_data[Cur_rows,"species"]=isp
    Complete_data[Cur_rows,"fly"] = rbinom(N[isp],1,Prop_fly[isp])
    Complete_data[Cur_rows,"g_size"] = rpois(N[isp],exp(rnorm(N[isp],Mu_group[isp],sigma_grp)))+1  
    Complete_data[Cur_rows,"d1_true"] = round(runif(N[isp],0.5,n_bins+0.5))
    for(irow in 1:N[isp]){
      Complete_data[Cur_rows[irow],"d1_obs"]=which(rmultinom(1,1,Prob_measure[Complete_data[Cur_rows[irow],"d1_true"],])==1)
      Complete_data[Cur_rows[irow],"d2_true"]=which(rmultinom(1,1,Prob_move[Complete_data[Cur_rows[irow],"d1_true"],])==1)
      Complete_data[Cur_rows[irow],"d2_obs"]=which(rmultinom(1,1,Prob_measure[Complete_data[Cur_rows[irow],"d2_true"],])==1)
    }
    if(first_obs_correct==TRUE)Complete_data[Cur_rows,"d1_obs"] =  Complete_data[Cur_rows,"d1_true"]
    
    
    #detection model
    X_obs1 = matrix(0,N[isp],length(Beta_det))
    X_obs1[,1]=1
    if(n_obs_bins>1){
      for(iind in 1:N[isp]){
        if(Complete_data[Cur_rows[iind],"d1_obs"] %in% Obs_bins[-1])X_obs1[iind,2+Complete_data[Cur_rows[iind],"d1_obs"]-Obs_bins[1]]=1
      }
    }
    X_obs1[,2+n_obs_bins]=Complete_data[Cur_rows,"g_size"]
    X_obs1[,3+n_obs_bins]=Complete_data[Cur_rows,"fly"]
    if(isp>1)X_obs1[,2+n_obs_bins+isp]=1
    X_obs2=X_obs1
    if(n_obs_bins>1){
      X_obs2[,3:(1+n_obs_bins)]=0
      for(iind in 1:N[isp]){
        if(Complete_data[Cur_rows[iind],"d2_obs"] %in% Obs_bins[-1])X_obs2[iind,2+Complete_data[Cur_rows[iind],"d2_obs"]-Obs_bins[1]]=1
      }
    }
    X_obs2[,2]=1
    
    Mu_1 = X_obs1%*%Beta_det
    Mu_2 = X_obs2%*%Beta_det
    Probit_det = matrix(0,N[isp],2)
    for(iind in 1:N[isp]){
      Sigma[1,2]=Sigma[2,1]=abs(Complete_data[Cur_rows[iind],"d1_true"]-Obs_bins[1])/(n_bins-Obs_bins[1])*cor_par
      if(Complete_data[Cur_rows[iind],"d1_obs"])
      Probit_det[iind,]=rmvnorm(1,c(Mu_1[iind],Mu_2[iind]),Sigma)
    }
    #ensure animals with perceived distances outside of strip have zero detection 
    Probit_det[which(!(Complete_data[Cur_rows,"d1_obs"]%in%Obs_bins)),1]=-1
    Probit_det[which(!(Complete_data[Cur_rows,"d2_obs"]%in%Obs_bins)),2]=-1
    ### determine detections
    Complete_data[Cur_rows,c("det1","det2")]=1.0*(Probit_det>0)
    counter=counter+N[isp]
    
  }
  Complete_data=Complete_data[1:(counter-1),]
  Obs_data=Complete_data[,-c(2,3)]
  Obs_data[which(!(Obs_data[,"d1_obs"]%in% Obs_bins)),"det1"]=0
  Obs_data[which(!(Obs_data[,"d2_obs"]%in% Obs_bins)),"det2"]=0
  Obs_data=Obs_data[-which(Obs_data[,"det1"]==0 & Obs_data[,"det2"]==0),]
  
  #put into format for hierarchical_DS_move.R  
  HDS_data = matrix(0,nrow(Obs_data)*2,9)
  colnames(HDS_data)=c("Transect","Match","Observer","Obs","Species","Seat","Distance","Group","Fly")
  for(irec in 1:nrow(Obs_data)){
    HDS_data[(irec-1)*2+1,]=c(1,irec,1,Obs_data[irec,"det1"],Obs_data[irec,"species"],1,Obs_data[irec,"d1_obs"],Obs_data[irec,"g_size"],Obs_data[irec,"fly"])
    HDS_data[(irec-1)*2+2,]=c(1,irec,2,Obs_data[irec,"det2"],Obs_data[irec,"species"],2,Obs_data[irec,"d2_obs"],Obs_data[irec,"g_size"],Obs_data[irec,"fly"])
  }
  HDS_data[which(HDS_data[,"Obs"]==0),"Distance"]=NA
  
  
  Data_in_strip = Complete_data[which(Complete_data[,"d1_true"]%in%Obs_bins),]
  G_true=N_true=rep(0,n_species)
  for(isp in 1:n_species){
    Which_sp = which(Data_in_strip[,"species"]==isp)
    G_true[isp] = length(Which_sp)
    N_true[isp]=sum(Data_in_strip[Which_sp,"g_size"])
  }

  return(list(Dat=HDS_data,Complete_data=Complete_data,G_true=G_true,N_true=N_true,Mu_group=Mu_group,Beta=Beta_det))
}