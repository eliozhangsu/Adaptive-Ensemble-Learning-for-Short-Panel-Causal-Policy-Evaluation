rm(list=ls())

 

print(Sys.time())

# Load the packages (comes with base R)
library(parallel)
library(doParallel)
library(foreach)




 

# Print the result
print(detectCores())



##Main References: https://github.com/susanathey/MCPanel/blob/master/tests/examples_from_paper/california/california_smoking.R






# Set the number of cores to use
num_cores <- 10 # You can replace this with the desired number of cores



print(num_cores)


# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)
  
  

# Define your R function to be run in parallel
par_function <- function(x) {
  
######################################################################################################################

  
  unlink(".RData")
  options(scipen=999)
  
  # Load necessary libraries
  library(dplyr)
  library(plm)  # For Fixed Effects models
  library(reshape2)  # For reshaping data
  library(softImpute)  # For matrix completion
  library(tidyr)
  library(foreign)
  library(Synth)
  library(did)
  #library(Rcpp)
  library(MCPanel)
  library(glmnet)
  library(ggplot2)
  library(latex2exp)
  
  #install.packages("devtools")
  #install.packages("latex2exp")
  library(devtools) 
  #install_github("susanathey/MCPanel")
  
  #d <- read.dta("repgermany.dta")
  
  #install.packages("devtools")
  #install.packages("latex2exp")
  #library(devtools) 
  #install_github("susanathey/MCPanel")
  #library(MCPanel)
  
  # Function to simulate panel data
  #There are totally "units" units and "periods" periods.
  #unit 1 is the only treated unit and the treatment is applied after (periods-T0) (>).
  
  
  simulate_panel_data <- function(units, periods, T0, num_covariates,sigma,SNR) {
    
    #units<-2; periods<-3;num_covariates<-4
    data <- expand.grid(Unit = 1:units, Period = 1:periods)
    for (i in 1:num_covariates) {
      data[[paste0("Covariate_", i)]] <- rep(0,nrow(data) )
    }
    for (j in 1:units){
      for (k in 1:num_covariates){
        data[data$Unit==j,paste("Covariate_",k,sep="")]<-arima.sim(n = periods, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
              sd = sqrt(0.1796))
      }
    }
  
      
    for (j in 1:units){
      for (k in 1:periods){
        data[data$Unit==j&data$Period==k,paste0("Covariate_",1:num_covariates,sep="")]<-data[data$Unit==j&data$Period==k,paste0("Covariate_",1:num_covariates,sep="")]+rnorm(1)
        }
    }
    
    
    #### DGP Settings: Synthetic  ####
    data<-data[order(data$Unit,data$Period),]
    data$Eta <- rep(rnorm(units),each=periods )

    data$Mean <- 0.5*(data$Covariate_2) + 1.5*(data$Covariate_4) + data$Eta   ##Signal without treatment
    
    data[data$Unit==1,"Mean"] <- data[data$Unit==3,"Mean"]*0.5 +data[data$Unit==5,"Mean"]*0.5
    
    data$Mean <- data$Mean  * SNR 
    #### DGP Settings End ####


   
    data$Nonlinear_Effect <- sin (data$Covariate_1^2) + 0.5 * (data$Covariate_3)^2
    data$Treatment <- ifelse( (data$Unit == 1) & (data$Period >= (periods - T0 + 1)), 1, 0)
    data$Treatment_Effect <- data$Treatment * data$Nonlinear_Effect
  
    data$Random_error <- rnorm(nrow(data),0, sd = sigma)
    data$Outcome <- data$Mean+data$Treatment_Effect + data$Random_error
    
    data$Outcome_wi_treatment_mean <- data$Mean + data$Nonlinear_Effect 
    data$Outcome_wo_treatment_mean <- data$Mean 
    data$Outcome_wo_treatment <- data$Mean + data$Random_error
    return(data)
  }
  
  

  
  
  

  
  
  
  
  ##############################################################################################
  
  # We pull the data first
  
  #Panel <- read.dta("http://dss.princeton.edu/training/Panel101.dta")
  
  #fixed <- plm(y ~ x1, data=Panel, index=c("country", "year"), model="within")  #fixed model
  #random <- plm(y ~ x1, data=Panel, index=c("country", "year"), model="random")  #random model
  #phtest(fixed,random) #Hausman test
  
  
    effect_control <- function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
    #train_length = 50
    train_data <- panel_data[(panel_data$Period<=train_length)&(panel_data$Unit %in%units_select), ]
    fixed <-  plm(as.formula(paste("Outcome ~", paste0("Covariate_",1:num_covariates,collapse="+"),sep="")), data=train_data, index=c("Unit", "Period"), model="within")  #fixed model
    random <- plm(as.formula(paste("Outcome ~", paste0("Covariate_",1:num_covariates,collapse="+"),sep="")), data=train_data, index=c("Unit", "Period"), model="random")  #random model
    
    
    phtest(fixed,random) #Hausman test
    
     
    test_data<- panel_data[panel_data$Period>train_length&panel_data$Period<=(train_length+test_length)&panel_data$Unit==test_unit, ] #
    newdata.p <- pdata.frame(test_data, index=c("Unit", "Period") )
    
    prediction_fixed  <- predict(fixed,newdata.p)
    
    prediction_random <- predict(random,newdata.p)
    
    prediction <-data.frame(prediction_fixed,prediction_random)
    return(prediction)
  
  }
  
  
  
  ####################################################################################
  
  # Synthetic Control
  # test_length=5; units_select<-c(2:10);train_length=190;test_unit<-2;panel_data[panel_data$Unit==test_unit,"Outcome"]<-panel_data[panel_data$Unit==3,"Outcome"] 
  synthetic_control <- function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
  
  #T1<-400; test_length=5; n<-units;train_length=400  
  dataprep_out <-
      dataprep(
        foo = panel_data,
        predictors    = paste0("Covariate_",1:num_covariates,sep=""),
        dependent     = "Outcome",
        unit.variable = "Unit",
        time.variable = "Period",
        #special.predictors = list(
        #  list("industry" ,1981:1990, c("mean")),
        #  list("schooling",c(1980,1985), c("mean")),
        #  list("invest80" ,1980, c("mean"))
        #),
        treatment.identifier = test_unit,
        controls.identifier = setdiff(units_select, test_unit),
        time.predictors.prior = 1:train_length,
        time.optimize.ssr = 1:train_length,
        #unit.names.variable = "country",
        time.plot = (train_length+1):(train_length+test_length)
      )
    synth_out <- synth(dataprep_out)
    synth_out$solution.w
    
    units_train<-setdiff(units_select,test_unit)
    test_data<-panel_data[( panel_data$Unit%in%units_train) & panel_data$Period>train_length & panel_data$Period<=(train_length+test_length),c("Unit","Period","Outcome")]
    test_data_mat<-matrix(test_data[,"Outcome"],length(units_train),test_length)
    
    prediction_sc <- as.vector(t(test_data_mat)%*%(synth_out$solution.w))
    prediction_sc
    
    #panel_data[panel_data$Unit==test_unit& panel_data$Period>train_length,c("Unit","Period","Outcome")]
    #panel_data[panel_data$Unit==3& panel_data$Period>train_length,c("Unit","Period","Outcome")]  
    #dataprep_out$Y1plot
    #pred
    #synth_out$solution.w
    prediction<-data.frame(prediction_sc)
    
    return(prediction)
  }
  
  
  
  
  ################################################################################################
  
  
  # Matrix Completion (MC)
  matrix_completion <- function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
    
    #train_length <-70; test_length<-5;test_unit<-2;units_select<-(2:10)
    
    
    panel_data<-panel_data[order(panel_data$Unit, panel_data$Period),]
  
    outcome_select<-panel_data[panel_data$Unit%in%units_select,"Outcome"]
    outcome_matrix <-matrix(outcome_select, length(units_select), length(outcome_select)/length(units_select) )
    train_data <-outcome_matrix[, 1:(train_length+test_length)]
    
    
    treat_mat <- matrix(1, length(units_select), train_length+test_length )
    treat_mat[test_unit==units_select,(train_length+1):(train_length+test_length)]<-0
    Y_obs <- train_data * treat_mat
  
    ## ------
    ## MC-NNM
    ## ------  
    
    
    est_model_MCPanel <- mcnnm_cv(Y_obs, treat_mat, to_estimate_u = 1, to_estimate_v = 1)
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(train_length+test_length,est_model_MCPanel$u) + t(replicate(length(units_select),est_model_MCPanel$v))
    prediction_MC <- (est_model_MCPanel$Mhat)[treat_mat ==0]
    prediction_MC
    
    
    
    ## -----
    ## Horizontal Regression
    ## EN : It does Not cross validate on alpha (only on lambda) and keep alpha = 1 (LASSO).
    ##      Change num_alpha to a larger number, if you are willing to wait a little longer.
    ## -----
    
    
    
    est_model_EN <- en_mp_rows(Y_obs, treat_mat, num_alpha = 50)
    
    prediction_EN <- (est_model_EN)[treat_mat ==0]
    prediction_EN
    
    
    
    
    ## -----
    ## Vertical Regression
    ## EN_T : It does Not cross validate on alpha (only on lambda) and keep alpha = 1 (LASSO).
    ##        Change num_alpha to a larger number, if you are willing to wait a little longer.
    ## -----
    
    
    
    est_model_ENT <- t(en_mp_rows(t(Y_obs), t(treat_mat), num_alpha = 50))
    
    prediction_ENT <- (est_model_ENT)[treat_mat ==0]
    prediction_ENT
    
    
    ## -----
    ## DID
    ## -----
    
    est_model_DID <- DID(Y_obs, treat_mat)
    prediction_DID <- est_model_DID[treat_mat ==0]
    prediction_DID
    
    ## -----
    ## SC-ADH
    ## -----
    
    est_model_ADH <- adh_mp_rows(Y_obs, treat_mat)
    prediction_ADH <- est_model_ADH[treat_mat ==0]
    prediction_ADH
    
    prediction<-data.frame(prediction_MC,prediction_EN,prediction_ENT,prediction_DID,prediction_ADH)
    return(prediction)
    }
  
  
  
  
  
  
  
  
  counterfactual_pred<-function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
  
  pred1<-effect_control(panel_data, units_select, num_covariates, train_length,test_length,test_unit)
  pred1
  
  
  pred2<-synthetic_control(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
  pred2
  
  pred3<-matrix_completion(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
  pred3
  
  test_outcome<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Outcome"]
  test_outcome
  
  test_outcome_wo_treatment<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Outcome_wo_treatment"]
  test_outcome_wo_treatment
  
  test_mean_wi_treatment<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Outcome_wi_treatment_mean"]
  test_mean_wi_treatment
  
  test_mean_wo_treatment<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Outcome_wo_treatment_mean"]
  test_mean_wo_treatment
                                                                                                                                                       1
  Treatment_Effect<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Treatment_Effect"]
  Treatment_Effect
  
  
  
  pred_matrix<-cbind(pred1, pred2, pred3, test_outcome, test_outcome_wo_treatment, test_mean_wo_treatment, test_mean_wi_treatment,Treatment_Effect)
  return(pred_matrix)
  
  }
  
  ###############################################
  ###############################################
  ###############################################
sim<-function(){
  
 
  
  for (units in c(40)) { 
  for (sigma in c(0.5, 1, 2)) {
  for (periods in c(15)) {
      

      T0 <- 7  # Treated time periods
      num_covariates <- 5  # Number of covariates
      bs_n<- 100  #
      snr <- 1 #signal to noise ratio
      
      
      
      panel_data <- simulate_panel_data(units = units, periods = periods, T0 = T0, num_covariates = num_covariates,sigma = sigma, SNR = snr)
      
      SNRATIO<- sqrt(mean((panel_data$Mean)^2))/sigma
      
      
      
      
      
      
      ######################################################
      ######################################################
      ###################Prediction#########################
      ######################################################
      ######################################################
      
      
      train_length<- (periods-T0) ; test_length<-T0  #The parameters need to be adjusted for cross-validation 
      
      units_select <- c(1:units);  test_unit<-1
      
      output<-counterfactual_pred(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
      
      pred_error<-abs(output[,1:8]-output[,"test_outcome_wo_treatment"])
      pred_error_mean<-apply(pred_error,2,mean)
      er<-pred_error_mean
      
      
      treatment_est<- (output[,"test_outcome"]-output[,1:8] )
      
      treatment_error<-abs(output[,"Treatment_Effect"]-treatment_est)
      apply(treatment_error,2,mean)
      
      pred_mat<-output[,1:9]
      pred_mat      
      
      for (cv_test_length in c(2)) {
          

  
  
  ######################################################
  ######################################################
  ###################Cross-validation###################
  ######################################################
  ######################################################
  test_length <- cv_test_length    ##Cross-validation parameters.
  train_length<- (periods-T0-test_length)  ##Cross-validation parameters.
  
  units_select <- c(1:units);  test_unit<-1
  
  output_cv<-counterfactual_pred(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
  
  cv_error<-abs(output_cv[,1:8]-output_cv[,"test_outcome"])
  cv_error_mean<-apply(cv_error,2,mean)

  cv<-cv_error_mean
  
  
  
  cv_er<- NULL
  cv_er<-(er[cv==min(cv)])[1]  
  
  
  ######################################################
  ######################################################
  ###############Perturbed Cross-Validation ############
  ######################################################
  ######################################################   
  

  
  
  devi_seq <- c(0.0625, 0.125,0.25, 0.5,1,2,4,8,16,32)
  
  ##########################################
  ##########################################
  ##########################################
  ########################################## 
  

  
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  
  
  all_normal_mc_er_vec <- NULL
  for (devi in devi_seq) {
            bs_error_mat<-NULL
            
              
            
            for (k in 1:bs_n){
                    try({ 
                    test_length <- cv_test_length   ##Cross-validation parameters.
                    train_length<- (periods-T0-cv_test_length)  ##Cross-validation parameters.
                    
                    units_select <- c(1:units);  test_unit<-1
                    
                    panel_data0<-NULL #perturbed panel dataframe
                    panel_data0<-panel_data
                    
                    
                    panel_data0<-panel_data0[order(panel_data0$Unit,panel_data0$Period),]
                    

                    dfk<- 1/(1+devi^2*sigma^2)
                    random_pert<- ((rnorm(nrow(panel_data0))*(devi*sigma))^2+1)*dfk
                    panel_data0$Outcome<-panel_data0$Outcome*random_pert
                    
                     
                    output_bs<-counterfactual_pred(panel_data0, units_select,num_covariates,train_length,test_length,test_unit)
                    
                    bs_error0<- abs(output_bs[,"test_outcome"]-output_bs[,1:8] )
                    
                    bs_error_mat<-rbind(bs_error_mat, apply(bs_error0,2,mean))
                    panel_data0<-NULL #perturbed panel dataframe
                    },silent=TRUE)
              
                    } 
            print(bs_error_mat)
            bs_error_mat == apply(bs_error_mat,1,min)
            model_weight<-apply(bs_error_mat == apply(bs_error_mat,1,min),2,mean)
            model_weight
            
            bs_error_mean<-apply(bs_error_mat,2,mean)
            
            bs1<-bs_error_mean;
            bs_er <- NULL
            bs_er<-(er[bs1==min(bs1)])[1]
            
            mc_er <- NULL
            mc_er<-mean(abs(as.matrix(pred_mat[,1:8])%*%as.vector(model_weight)-output[,"test_outcome_wo_treatment"]))
            all_normal_mc_er_vec <- c(all_normal_mc_er_vec, bs_er, mc_er)
  }
  all_normal_mc_er_vec_names <- paste(c("all_normal_bs_er","all_normal_mc_er"),rep(devi_seq,each=2),sep="_")

  
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  
  
  all_chi_mc_er_vec <- NULL
  for (devi in devi_seq) {
    bs_error_mat<-NULL
    
    
    
    for (k in 1:bs_n){
      try({
      test_length <- cv_test_length   ##Cross-validation parameters.
      train_length<- (periods-T0-cv_test_length)  ##Cross-validation parameters.
      
      units_select <- c(1:units);  test_unit<-1
      
      panel_data0<-NULL #perturbed panel dataframe
      panel_data0<-panel_data
      
      
      panel_data0<-panel_data0[order(panel_data0$Unit,panel_data0$Period),]
      
      dfk<- devi*(sigma)
      random_pert <- rchisq(nrow(panel_data0),dfk)*1/(dfk)
      panel_data0$Outcome<-panel_data0$Outcome*random_pert
      
      
      output_bs<-counterfactual_pred(panel_data0, units_select,num_covariates,train_length,test_length,test_unit)
      
      bs_error0<- abs(output_bs[,"test_outcome"]-output_bs[,1:8] )
      
      bs_error_mat<-rbind(bs_error_mat, apply(bs_error0,2,mean))
      panel_data0<-NULL #perturbed panel dataframe
      },silent=TRUE)
      
    } 
    print(bs_error_mat)
    bs_error_mat == apply(bs_error_mat,1,min)
    model_weight<-apply(bs_error_mat == apply(bs_error_mat,1,min),2,mean)
    model_weight
    
    bs_error_mean<-apply(bs_error_mat,2,mean)
    
    bs1<-bs_error_mean;
    
    bs_er <-NULL
    bs_er<-(er[bs1==min(bs1)])[1]
    
    mc_er <- NULL
    mc_er<-mean(abs(as.matrix(pred_mat[,1:8])%*%as.vector(model_weight)-output[,"test_outcome_wo_treatment"]))
    all_chi_mc_er_vec <- c(all_chi_mc_er_vec, bs_er, mc_er)
  }
  all_chi_mc_er_vec_names <- paste(c("all_chi_bs_er","all_chi_mc_er"),rep(devi_seq,each=2),sep="_")
  
  
  
  #########################################################
  #########################################################
  #########################################################
  #########################################################
  #########################################################
  #########################################################  
  

  
  
  result <-c(snr, units, periods, T0, cv_test_length, num_covariates, sigma, SNRATIO, bs_n,              all_chi_mc_er_vec, all_normal_mc_er_vec, cv_er, er)
  names(result)<-c("snr","units","periods","T0", "cv_test_length","num_cov", "sigma","sn_ratio","bs_n",  all_chi_mc_er_vec_names, all_normal_mc_er_vec_names, "cv_er",paste(names(er),"_er",sep=""))
  print(result)

  
  ##path at laptop
  ##path <-"C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/R code/sim_synthetic16/"
  
  ##path at desktop
  path <-"C:/Users/yongl/OneDrive - Seattle University/2024 CSE Summer UG Research Award/R code/sim_synthetic16/"
  
  
  vvv<-round(abs(rnorm(1)*10000000000000000000000000))  
  
  
  
  filename<-paste(path,vvv,".","txt",sep="")
  write.table(t(result), file=filename, row.names=FALSE,col.names=T)
  gc()  
  
  
  }
  }
  }  
  }
  }


  for (k in 1:11)    {  sim()  }

                                                                  z
  
  
############################################################################################################
}


# Generate a sequence of values to be passed to the function
input_values <- 1:num_cores

################################################################################
################################################################################
################################################################################
################################################################################ 

# Run the function in parallel
result_list <- foreach(i = input_values, .combine = c) %dopar% {
  par_function(i)
}

# Stop the parallel backend
stopCluster(cl)

# 'result_list' now contains a list of matrices
#print(result_list)



print(Sys.time())
print('simulation is done')



