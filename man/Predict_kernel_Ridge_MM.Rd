
\name{Predict_kernel_Ridge_MM}
\alias{Predict_kernel_Ridge_MM}
%- Also NEED an '\alias' for EACH other topic documented here.
\encoding{utf8}

\title{
Predict function for Kernel_Ridge_MM object
%%  ~~function to do ... ~~
}

\description{
Predict the value(s) for a vector or a design matrix of covariates (i.e. features)
}

\usage{
		Predict_kernel_Ridge_MM( Model_kernel_Ridge_MM, Matrix_covariates_target,
		
		X_target=as.vector(rep(1,dim(Matrix_covariates_target)[1])),
		
		Z_target=diag(1,dim(Matrix_covariates_target)[1]) ) 
	}
	
%- maybe also 'usage' for other objects documented here.

\arguments{
  
  \item{Model_kernel_Ridge_MM}{a Kernel_Ridge_MM object}
  
  \item{Matrix_covariates_target}{numeric matrix; design matrix of covariates for target data}
  
  \item{X_target}{numeric matrix; design matrix of predictors with fixed effects for target data (default is a vector of ones)}
  
  \item{Z_target}{numeric matrix; design matrix of predictors with random effects for target data (default is identity matrix)}


 }

 \details{ 
The matrix Matrix_covariates_target is mandatory to build the kernel matrix (with Matrix_covariates_train from Model_kernel_Ridge_MM) for prediction.
}

\value{
  \item{f_hat }{ Predicted value for target data, i.e. f_hat = X_target*Beta_hat + Z_target*U_target where U_target=K_target_train*alpha_train
  and alpha_train is the BLUP of alpha for the model, i.e. alpha_train=Cov(alpha,Y_train)*Var(Y_train)^{-1}*(Y_train - E[Y_train]) }
}
 
\author{Laval Jacquin <jacquin.julien@gmail.com>}

\examples{

\dontrun{

library(KRMM)

### SIMULATE DATA 
set.seed(123)
p=200
N=100

beta=rnorm(p, mean=0, sd=1.0)
X=matrix(runif(p*N, min=0, max=1), ncol=p, byrow=TRUE)  #X: covariates (i.e. predictors)

f=X\%*\%beta                    #f: data generating process (i.e. DGP)
E=rnorm(N, mean=0, sd=0.5)

Y=f+E                           #Y: observed response data

hist(f)
hist(beta)
Nb_train=floor((2/3)*N)

###======================================================================###
### CREATE TRAINING AND TARGET SETS FOR RESPONSE AND PREDICTOR VARIABLES ###
###======================================================================###

Index_train=sample(1:N, size=Nb_train, replace=FALSE)

### Covariates (i.e. predictors) for training and target sets

Predictors_train=X[Index_train, ]
Response_train=Y[Index_train]

Predictors_target=X[-Index_train, ]
True_value_target=f[-Index_train]    #True value (generated by DGP) we want to predict


###=================================================================================###
### PREDICTION WITH KERNEL RIDGE REGRESSION SOLVED WITHIN THE MIXED MODEL FRAMEWORK ### 
###=================================================================================###

Gaussian_KRR_model_train = Kernel_Ridge_MM( Y_train=Response_train,
 Matrix_covariates_train=Predictors_train, method="RKHS", rate_decay_kernel=5.0)

### Predict new entries for target set and measure prediction error

f_hat_target_Gaussian_KRR = Predict_kernel_Ridge_MM( Gaussian_KRR_model_train, 
 Matrix_covariates_target=Predictors_target )

plot(f_hat_target_Gaussian_KRR, True_value_target)

}


}
