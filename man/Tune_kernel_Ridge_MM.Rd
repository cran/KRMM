
\name{Tune_kernel_Ridge_MM}
\alias{Tune_kernel_Ridge_MM}
\encoding{utf8}

\title{Tune kernel ridge regression in the mixed model framework}

\description{
Tune_kernel_Ridge_MM tunes the rate of decay parameter of kernels, by K-folds cross-validation, for kernel ridge regression
}

\usage{
  		Tune_kernel_Ridge_MM( Y_train, X_train=as.vector(rep(1,length(Y_train))),
		
		Z_train=diag(1,length(Y_train)), Matrix_covariates_train, 
		
		method="RKHS", kernel="Gaussian", rate_decay_kernel=0.1, 
		
		degree_poly=2, scale_poly=1, offset_poly=1,
		
		degree_anova=3, init_sigma2K=2, init_sigma2E=3,

		convergence_precision=1e-8, nb_iter=1000, display="FALSE",

		rate_decay_grid=seq(0.1,1.0,length.out=10),

		nb_folds=5, loss="mse")		
	}

\arguments{

	\item{rate_decay_grid}{
	Grid over which the rate of decay is tuned by K-folds cross-validation
	}
  
	\item{nb_folds}{
	Number of folds, i.e. K=nb_folds (default is 5)
	}
  
	\item{loss}{ mse (mean square error) or cor (correlation) (default is mse)}
  
	\item{Y_train}{numeric vector; response vector for training data}

	\item{X_train}{numeric matrix; design matrix of predictors with fixed  effects for training data (default is a vector of ones)}
  
	\item{Z_train}{numeric matrix; design matrix of predictors with random effects for training data (default is identity matrix)}

	\item{Matrix_covariates_train}{numeric matrix of entries used to build the kernel matrix}
 
	\item{method}{character string; RKHS, GBLUP or RR-BLUP}

	\item{kernel}{character string; Gaussian, Laplacian or ANOVA (kernels for RKHS regression ONLY, the linear kernel is automatically built for GBLUP and RR-BLUP and hence no kernel is supplied for these methods) }

	\item{rate_decay_kernel}{	numeric scalar;
		hyperparameter of the Gaussian, Laplacian or ANOVA kernel (default is 0.1)
	}

	\item{degree_poly, scale_poly, offset_poly}{numeric scalars; parameters for polynomial kernel (defaults are 2, 1 and 1 respectively)}
  
	\item{degree_anova}{
		numeric scalar; parameter for ANOVA kernel (defaults is 3)
	}
  
	\item{init_sigma2K, init_sigma2E}{
		numeric scalars;
		initial guess values, associated to the mixed model variance parameters, for the EM-REML algorithm
		(defaults are 2 and 3 respectively)
		}
	
	\item{convergence_precision, nb_iter}{		
		numeric scalars; convergence precision (i.e. tolerance) associated to the mixed model variance parameters,
		for the EM-REML algorithm, and number of maximum iterations allowed if 
		convergence is not reached (defaults are 1e-8 and 1000 respectively)
	}	
		
	\item{display}{
		boolean (TRUE or FALSE character string); 
		should estimated components be displayed at each   
		iteration 
	}

}


\value{

\item{tuned_model}{the tuned model (a Kernel_Ridge_MM object)}

\item{expected_loss_grid}{the average loss for each rate of decay tested over the grid}

\item{optimal_h}{the rate of decay minimizing the average loss}

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

f=X\%*\%beta        #f: data generating process (i.e. DGP)
E=rnorm(N, mean=0, sd=0.5)

Y=f+E                                                  #Y: response data

hist(f)
hist(beta)
Nb_train=floor((2/3)*N)

###======================================================================###
###	CREATE TRAINING AND TARGET SETS FOR RESPONSE AND PREDICTOR VARIABLES ###
###======================================================================###

Index_train=sample(1:N, size=Nb_train, replace=FALSE)

### Covariates (i.e. predictors) for training and target sets

Predictors_train=X[Index_train, ]
Response_train=Y[Index_train]

Predictors_target=X[-Index_train, ]
True_value_target=f[-Index_train]    #True value (generated by DGP) we want to predict

###=======================###
### Tuned Gaussian Kernel ###
###=======================###

Tuned_Gaussian_KRR_train = Tune_kernel_Ridge_MM( Y_train=Response_train, Matrix_covariates_train
=Predictors_train, method='RKHS', rate_decay_grid=seq(1,10,length.out=10), nb_folds=5, loss='mse' )

Tuned_Gaussian_KRR_model_train = Tuned_Gaussian_KRR_train$tuned_model
Tuned_Gaussian_KRR_train$optimal_h
Tuned_Gaussian_KRR_train$rate_decay_grid
Tuned_Gaussian_KRR_train$expected_loss_grid

dev.new()
plot(Tuned_Gaussian_KRR_train$rate_decay_grid, Tuned_Gaussian_KRR_train$expected_loss_grid,
 type="l", main="Tuning the rate of decay (for Gaussian kernel) with K-folds cross-validation")

### Predict with tuned model
 
f_hat_target_tuned_Gaussian_KRR = Predict_kernel_Ridge_MM( Tuned_Gaussian_KRR_model_train, 
Matrix_covariates_target=Predictors_target )

mean((f_hat_target_tuned_Gaussian_KRR-True_value_target)^2)
cor(f_hat_target_tuned_Gaussian_KRR,True_value_target)


}


}
