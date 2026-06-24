rm(list=ls())



unlink(".RData")
options(scipen=999)

#install.packages("devtools")
#install.packages("latex2exp")
#install_github("susanathey/MCPanel")
#install.packages("devtools")
#install.packages("latex2exp")
#install_github("susanathey/MCPanel")


# Load necessary libraries
library(dplyr)
library(plm)  # For Fixed Effects models
library(reshape2)  # For reshaping data
library(softImpute)  # For matrix completion
library(tidyr)
library(foreign)
library(Synth)
library(did)
library(MCPanel)
library(glmnet)
library(ggplot2)
library(latex2exp)
library(devtools) 







library(doParallel)
library(foreach)





##Data Cleaning##


##Main References: https://github.com/susanathey/MCPanel/blob/master/tests/examples_from_paper/california/california_smoking.R

sodatax_ct<-c("San Francisco") 

weight_matrix<-read.csv("C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/sodatax/adult_population_ratio.csv",sep=",",header=T)
weight_matrix<-weight_matrix[(weight_matrix$County!="Alameda")& (weight_matrix$County!="Tulare"),]

#unique(weight_matrix$County)
#unique(weight_matrix$Year)

panel_data00<-read.csv("C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/sodatax/California_obesity_demographics_COMPLETE_UPDATED_split_obesity.csv",sep=",",header=T)


panel_data00<-panel_data00[(panel_data00$County!="Alameda")& (panel_data00$County!="All"),]

ct_list<-unique(panel_data00$County)




names(panel_data00)[4]<-"Percentage.OverweightandObese"
names(panel_data00)[6]<-"Percentage.Overweight"
names(panel_data00)[8]<-"Percentage.Obese"
names(panel_data00)[10]<-"Pecentange.Latino"
names(panel_data00)[12]<-"Pecentange.White"
names(panel_data00)[36]<-"Percentage.College.Completed"
names(panel_data00)[52]<-"Percentage.Full.Time.Employed"
names(panel_data00)[62]<-"Percentage.Urban"
names(panel_data00)[66]<-"Percentage.Uninsured"
names(panel_data00)[74]<-"Percentage.Medicaid"

panel_data00[panel_data00$County=="San Francisco","Percentage.Obese"]




periods<-13 # Total number of periods.
T1<-6  # Post-treatment time periods
T0 <- periods - T1   # Pre-treatment time periods

##Parameters in Model and Estimation
notax_units<-length(ct_list)-1 # number of untreated units
units<-(notax_units+1); #total number of units
units_select <- (1:units) 

test_unit <- 1
nx <- 4 #number of covariates
bs_n<- 100
pcv_test_length  <- 3
a<-0.5


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
panel_data0<-panel_data00[,c("County","Year","Percentage.Obese",
                          "Pecentange.Latino","Percentage.College.Completed",
                          "Percentage.Full.Time.Employed","Percentage.Urban")]
 
#panel_data0<-panel_data00[,c("County","Year","Percentage.OverweightandObese",
#                             "Pecentange.Latino","Pecentange.White","Percentage.College.Completed",
#                             "Percentage.Full.Time.Employed","Percentage.Urban")]




names(panel_data0)[3]<-"Outcome"
names(panel_data0)[4:ncol(panel_data0)]<-paste("Covariate_", 1:nx, sep="")
panel_data0$Period<- (panel_data0$Year-2010)





notax_ct<-setdiff(ct_list,sodatax_ct)
 
Unit<-c(1:length(ct_list))
County<-c(sodatax_ct,notax_ct)
ct_units<-data.frame(County,Unit)


##Cleaned Panel Data and Output
df <-merge(panel_data0,ct_units, by.x="County",by.y="County")
df0 <- df[order(df$Unit,df$Period),]
write.csv(df0, "C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/sodatax/panel_data_clean.csv", row.names = F)




df0<-read.csv("C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/sodatax/panel_data_clean.csv",sep=",",header=T)
df0$Outcome <- log(  (df0$Outcome/100)/(1- df0$Outcome/100) )

# Estimate sigma from pre-treatment fixed-effect residuals
sigma_data <- df0[df0$Period <= T0, ]

fe_sigma_model <- lm(
  Outcome ~ Covariate_1 + Covariate_2 + Covariate_3 +
    Covariate_4 + factor(Unit) + factor(Year),
  data = sigma_data
)
summary(fe_sigma_model)
sigma <- summary(fe_sigma_model)$sigma
print(sigma)


################################################################################
################################################################################
################################################################################
################################################################################
##Definition of Functions####################################################### 
  
  
  effect_control <- function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
    
    #train_length = 5; units_select = c(2:33); test_unit = 8
    
    train_data <- panel_data[(panel_data$Period<=train_length)&(panel_data$Unit %in%units_select), ]
    fixed <-  plm(as.formula(paste("Outcome ~", paste0("Covariate_",1:num_covariates,collapse="+"),sep="")), data=train_data, index=c("Unit", "Period"), model="within")  #fixed model
    random <- plm(as.formula(paste("Outcome ~", paste0("Covariate_",1:num_covariates,collapse="+"),sep="")), data=train_data, index=c("Unit", "Period"), model="random")  #random model

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
  
  #train_length = 5; units_select = c(2:33); test_unit = 8; T1<-6 
  dataprep_out <-
      dataprep(
        foo = panel_data,
        predictors    = paste0("Covariate_",1:num_covariates,sep=""),
        dependent     = "Outcome",
        unit.variable = "Unit",
        time.variable = "Period",

        treatment.identifier = test_unit,
        controls.identifier = setdiff(units_select, test_unit),
        time.predictors.prior = 1:train_length,
        time.optimize.ssr = 1:train_length,

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
  
  
  
  
 
  
  # Matrix Completion (MC)
  matrix_completion <- function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
    
    #train_length = 5; units_select = c(2:33); test_unit = 8; T1<-6
    
    
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
  
  
  
  
  #train_length = 5; units_select = c(2:33); test_unit = 8; T1<-6
   
  counterfactual_pred<-function(panel_data, units_select, num_covariates, train_length,test_length,test_unit) {
  
  pred1<-effect_control(panel_data, units_select, num_covariates, train_length,test_length,test_unit)
  pred1
  
  
  pred2<-synthetic_control(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
  pred2
  
  pred3<-matrix_completion(panel_data, units_select,num_covariates,train_length,test_length,test_unit)
  pred3
  

  test_outcome<-panel_data[panel_data$Unit==test_unit&panel_data$Period>train_length&panel_data$Period<=(train_length +test_length),"Outcome"]
  test_outcome
  
  
  
  pred_matrix<-cbind(pred1, pred2, pred3, test_outcome)
  return(pred_matrix)
  }


  
##End of Definition of Functions################################################
################################################################################
################################################################################  

  


################################################################################
################################################################################
###################PCV for San Francisco########################################
bs_error_mat<-NULL
model_weight_mat <- NULL

pred_mat<-counterfactual_pred(panel_data=df0, units_select = unique(df0$Unit), num_covariates=nx, train_length=(periods-T1),test_length=T1,test_unit=1)
pred_mat



#bs_n <- 5  
for (k in 1:bs_n){
          try({ 
            
          panel_data0<- NULL  
          panel_data0<-df0
          
          
          panel_data0<-panel_data0[order(panel_data0$Unit,panel_data0$Period),]
          
          devi <- a * sigma
          dfk<- 1/(1+devi^2)
          random_pert<- ((rnorm(nrow(panel_data0))*devi)^2+1)*dfk
          panel_data0$Outcome<-panel_data0$Outcome*random_pert

           
           
          output_bs<-counterfactual_pred(panel_data0, units_select,num_covariates=nx,train_length=( periods- pcv_test_length -T1),test_length=pcv_test_length,test_unit=1)
          
          bs_error0<- abs(output_bs[,"test_outcome"]-output_bs[,1:8] )
          
          bs_error_mat<-rbind(bs_error_mat, apply(bs_error0,2,mean))
          },silent=TRUE)
          } 


model_weight<-apply(bs_error_mat == apply(bs_error_mat,1,min),2,mean)
model_weight
cty_weight = cty = "San Francisco"
model_weight_mat <- rbind(model_weight_mat, data.frame(cty_weight,t(model_weight))  )
 


bs_error_mean<-apply(bs_error_mat,2,mean)
bs_error_mean == min( bs_error_mean )
pred_mat[,1:8][,bs_error_mean == min( bs_error_mean )]
selection_pcv <- pred_mat[,1:8][,bs_error_mean == min( bs_error_mean )]

  
combination_pcv<-as.matrix(pred_mat[,1:8])%*%as.vector(model_weight)
my_pred <- cbind(pred_mat,selection_pcv, combination_pcv)
my_pred <- exp(my_pred)/(1+exp(my_pred))*100
my_pred

my_pred$selection_gap  <- my_pred$selection_pcv - my_pred$test_outcome
my_pred$combination_gap <- my_pred$combination_pcv - my_pred$test_outcome


gap_mat <- NULL
selection_gap <- my_pred$selection_gap
combination_gap <- my_pred$combination_gap
selection_pred <- my_pred$selection_pcv
combination_pred <- my_pred$combination_pcv
actual_rate <- my_pred$test_outcome
county <- df0[df0$Unit==1&df0$Period>=T0+1,"County"]
years <- df0[df0$Unit==1&df0$Period>=T0+1,"Year"]
  
  
gap_mat <-rbind(gap_mat, data.frame(county,years, actual_rate,selection_pred,selection_gap,combination_pred,combination_gap))
gap_mat
 




################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################
################Placebo Analysis ###############################################

panel_data_plb<-df0[df0$County!= "San Francisco",]

units_select_plb <- unique(panel_data_plb$Unit);

plcb_gap <-NULL

for (kk in units_select_plb) { 

        #kk <- 28
        print(kk)
      
        pred_mat_plb<-counterfactual_pred(panel_data_plb, units_select_plb, num_covariates=nx, train_length=(periods-T1),test_length=T1,test_unit=kk)
        print(pred_mat_plb)
      
      ###################PCV################################
      
        
        
        
        #bs_n <- 5
        bs_error_mat<-NULL
        for (k in 1:bs_n){
          try({ 
          
          panel_data0<-  NULL
          panel_data0<-  panel_data_plb
          
          devi <- a * sigma
          dfk<- 1/(1+devi^2)
          random_pert<- ((rnorm(nrow(panel_data0))*devi)^2+1) * dfk
          panel_data0$Outcome<-panel_data0$Outcome*random_pert
        
          
        
          output_bs<-counterfactual_pred(panel_data0, units_select_plb,num_covariates=nx,train_length=(periods-T1-pcv_test_length),test_length=pcv_test_length,test_unit=kk)
          
          bs_error0<- abs(output_bs[,"test_outcome"]-output_bs[,1:8] )
          
          bs_error_mat<-rbind(bs_error_mat, apply(bs_error0,2,mean))
          },silent=TRUE)
          } 

      
      
      
        model_weight_plb<-apply(bs_error_mat == (apply(bs_error_mat,1,min)), 2, mean)
        print(model_weight_plb)
        
        ###
        bs_error_mean<-apply(bs_error_mat,2,mean)
        bs_error_mean == min( bs_error_mean )

        selection_pcv_plb <- pred_mat_plb[,1:8][,bs_error_mean == min( bs_error_mean )]        
        
        ###

        my_pred_plcb<-NULL 
        combination_pcv_plb<-as.matrix(pred_mat_plb[,1:8])%*%as.vector(model_weight_plb)
        
        
        ##
        my_pred_plcb<-cbind(pred_mat_plb, selection_pcv_plb, combination_pcv_plb)
        
        my_pred_plcb <- exp(my_pred_plcb)/(1+exp(my_pred_plcb))*100
         
        
        my_pred_plcb$selection_gap <- my_pred_plcb$selection_pcv - my_pred_plcb$test_outcome
        my_pred_plcb$combination_gap <- my_pred_plcb$combination_pcv - my_pred_plcb$test_outcome 
        print(my_pred_plcb)
        
        

        ##
        actual_rate <- my_pred_plcb$test_outcome
        selection_pred <- my_pred_plcb$selection_pcv
        selection_gap <- my_pred_plcb$selection_gap
        combination_pred <- my_pred_plcb$combination_pcv
        combination_gap <- my_pred_plcb$combination_gap
        
        
        county <- df0[df0$Unit==kk&df0$Period>=T0+1,"County"]
        years <- df0[df0$Unit==kk&df0$Period>=T0+1,"Year"]
        
        
        gap_mat <-rbind(gap_mat,data.frame(county,years, actual_rate,selection_pred,selection_gap,combination_pred,combination_gap))
        cty_weight <- unique(county)
        model_weight_mat <- rbind(model_weight_mat, data.frame(cty_weight,t(model_weight_plb))  )        
        
        
        
        print(aggregate(cbind(selection_gap, combination_gap)~county,data=gap_mat,mean)  )
        print(model_weight_mat)
         
        
        }


print(aggregate(gap~county,gap_mat,mean))


model_weight_mat
gap_mat

path = "C:/Users/ezhang2/OneDrive - Seattle University/2024 CSE Summer UG Research Award/sodatax/"

write.csv(gap_mat, paste(path, "placebo_gap_mat_all_",a,"a_",pcv_test_length,"cv" ,".csv",sep=""), row.names = F)

write.csv(model_weight_mat, paste(path, "model_weight_all_",a,"a_",pcv_test_length,"cv" ,".csv",sep=""), row.names = F)

############## End #############################################################
############## End #############################################################


sodatax_weight<-weight_matrix[! (weight_matrix$County%in% sodatax_ct),]
sodatax_weight_sum<-aggregate(Total.Adult.Population~Year,sodatax_weight,sum)
names(sodatax_weight_sum)[2]<- "Total.Adult.Population.Sum"

sodatax_weight<-merge(sodatax_weight, sodatax_weight_sum,by.x="Year",by.y="Year",all.x=TRUE )
sodatax_weight$Weight<-sodatax_weight$Total.Adult.Population/sodatax_weight$Total.Adult.Population.Sum
#sodatax_weight


sodatax_df<-df[!df$County%in%sodatax_ct,]
sodatax_df<-merge(sodatax_df,sodatax_weight[,c("Year", "County","Weight")], by.x=c("County","Year"), by.y=c("County","Year") )

sodatax_df$weighted.rate<-sodatax_df$Weight* exp(sodatax_df$Outcome)/(1+ exp(sodatax_df$Outcome) )
head(sodatax_df)
#dim(sodatax_df)
nosodatax_rate<-(aggregate(weighted.rate~Year,sodatax_df,sum)[,2])
aggregate(Weight~Year,sodatax_df,sum)

sf_actural_rate<-exp(df0[df0$Unit==1,"Outcome"])/( 1 + exp(df0[df0$Unit==1,"Outcome"]))*100

sf_counterfactual_rate<-c((sf_actural_rate[1:7]),my_pred[,"prediction_pcv"])

output<-(data.frame(nosodatax_rate,sf_actural_rate,sf_counterfactual_rate))
output$gap <- output$sf_actural_rate-output$sf_counterfactual_rate
output 



