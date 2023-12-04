library(leaps)
library(glmnet)
library(plotmo)
library(crayon)
library(corrplot) # Load the "corrplot" package
library(psych)    # Load the "psych" package
library(car)      # Load the "car" package for Variance Inflation Factor (VIF) calculation



#### Function definitions ####

## print min of coefficients metrics function
min_coefficients <- function(summary) {
  rss_min <- which.min(summary$rss)
  adjr2_min <- which.max(summary$adjr2)
  cp_min <- which.min(summary$cp)
  bic_min <- which.min(summary$bic)
  cat("\nLocation of RSS min: ", rss_min, "\n")
  cat("Location of adj-RSq max:",adjr2_min,"\n")
  cat("Location of Cp min:",cp_min,"\n")
  cat("Location of BIC min:",bic_min,"\n")
  values <- c(rss_min, adjr2_min, cp_min, bic_min)
  return (values)
  
}

## plot metrics function
plot_metrics <- function(summary){
  dev.new()
  par(mfrow=c(2,2))
  
  plot(round(summary$rss, 3), xlab="Number of Variables ", ylab="RSS", type="l")
  points(which.min(round(summary$rss, 3)),min(round(summary$rss,3)), col="red",cex=2,pch=20)
  
  plot(round(summary$adjr2, 3) ,xlab="Number of Variables ", ylab="Adjusted RSq",type="l")
  points(which.max(round(summary$adjr2, 3)),max(round(summary$adjr2,3)), col="red",cex=2,pch=20)
  
  plot(round(summary$cp,3) ,xlab="Number of Variables ",ylab="Cp", type="l")
  points(which.min(round(summary$cp,3)),min(round(summary$cp,3)),col="red",cex=2,pch=20)
  
  plot(round(summary$bic,3) ,xlab="Number of Variables ",ylab="BIC",type="l")
  points(which.min(round(summary$bic,3)),min(round(summary$bic, 3)),col="red",cex=2,pch=20)
  
  cat("Plotted!")
}

## print test error metrics functions
test_errors <- function(test, regfit, summary){
  test.mat=model.matrix(Y~.,data=test)
  name=c('adjr2', 'cp', 'bic')
  val.errors=rep(NA,3)
  val.n_beta=rep(NA,3)
  
  for (i in 1:length(name)){
    if(name[i]=='adjr2'){
      coefi = coef(regfit,id=which.max(get(name[i],summary)))
    }else{
      coefi = coef(regfit,id=which.min(get(name[i],summary)))
    }
    pred = test.mat[,names(coefi)]%*%coefi
    val.n_beta[i] = length(coefi)
    val.errors[i] = mean((test$Y-pred)^2)
  }
  cat("\t\t    ", name, "\nErrors BSS: ", val.errors)
  cat("\nbest BSS: ", name[which.min(val.errors)], " with number of regressors: ", val.n_beta[which.min(val.errors)]-1, "\n")
  
  values <- c(min(val.errors), val.n_beta[which.min(val.errors)]-1, name[which.min(val.errors)])
  return(values)
}




######### READ DATASET ######### 
#setwd('D:/duino/DriveUnisa/2023/Data Analysis/finalProject_DA')
dataset <- read.csv(file = 'D:/duino/DriveUnisa/2023/Data Analysis/finalProject_DA/Gruppo11/RegressionDataset_DA_group11.csv')
head(dataset)
dim(dataset)
n = nrow(dataset); n
data <- na.omit(dataset) # Remove rows with missing values


######### TRAIN AND TEST SPLIT ########
set.seed(1)
#set.seed(2023)        # Use a fixed seed for reproducibility of the experiment
# Save "Y" values in "y" and all values except for "Y" in "x"
train_percentage = 0.8

x <- model.matrix(Y~., data)[, -1] # Save all values except for "Y"
y <- data$Y # Save "Y" values
train <- sample(1:nrow(x), train_percentage * nrow(x)) # Sample 80% of the data for training set
test <- (-train)  # The remaining data will be the test set
y.test <- y[test] # Save test set "Y" values
train_data <- data[train, ]
test_data <- data[test, ]


# Parameters
n <- nrow(data[train, ]) # Number of observations in the training set
p <- ncol(x[train, ])    # Number of regressors


#Dict which contains best MSEs and associated technique
final_mses_dict <- list()


######### Preliminary analysis ######### 
######### Correlation Matrix ######### 
dev.new()
corData <- round(cor(x), digits = 2) # Calculate the correlation matrix of x and round the values to 2 decimal places
# Uncomment the next line to save the plot as a PDF file
# pdf(file="corr.pdf", width=26, height=15)
corPlot(corData, cex = 0.22, show.legend = TRUE, main = "Correlation Matrix") # Plot the correlation matrix with specified font size and display legend



######### Linear Model ######### 
cat(bold('Linear Model\n'))
lm.mod <- lm(Y ~ ., data = data[train, ]) # Fit a linear regression model to the data, using all variables in training set
summary(lm.mod) # Get the summary of the model
lm.predict <- predict(lm.mod, newdata = (data[test, ])) # Make predictions using the model on the test data
lm.coeff <- coef(lm.mod) # Get the coefficients of the model
lm.coeff

lm.mse <- mean((lm.predict - y[test])^2) # Calculate the mean squared error between the predicted values and the actual values in the test data
lm.mse # Display the mean squared error



######### BSS ######### 
cat(bold('Best Subset Selection\n'))
regfit.full=regsubsets(Y~., train_data,  really.big=TRUE, nvmax = 50)
reg.summary=summary(regfit.full)
reg.summary

## Print min number of regressors for all metrics
min_coef_full = min_coefficients(reg.summary)

## Plot RSS, adjusted R2, Cp, and BIC for all of the models at once
plot_metrics(reg.summary)

## Test errors BSS
errors_full = test_errors(test_data, regfit.full, reg.summary)

## Taking best MSE for BSS 
best_bss_mse = errors_full[[1]]
best_bss_name = errors_full[[3]]
key <- paste("bss ", best_bss_name)
final_mses_dict[[key]] <- best_bss_mse



######### FORWARD ######### 
cat(bold('Forward Stepwise Selection\n'))
regfit.fwd=regsubsets(Y~., train_data, nvmax=50, method ="forward")
reg.summaryForward=summary(regfit.fwd)

## Print min number of regressors for all metrics
min_coef_fwd = min_coefficients(reg.summaryForward)

## Plot RSS, adjusted R2, Cp, and BIC for all of the models at once
plot_metrics(reg.summaryForward)

## Test errors Forward
errors_fwd = test_errors(test_data, regfit.fwd,reg.summaryForward)
print(errors_fwd)

## Taking best MSE for Forward BSS 
best_fwd_mse = as.numeric(errors_fwd[[1]])
best_fwd_name = errors_fwd[[3]]
key <- paste("fwd ", best_fwd_name)
final_mses_dict[[key]] <- best_fwd_mse

print(final_mses_dict)



######### BACKWARD ######### 
cat(bold('Backward Stepwise Selection\n'))
regfit.bwd=regsubsets(Y~., train_data, nvmax=50, method ="backward")
reg.summaryBackward=summary(regfit.bwd)

## Print min number of regressors for all metrics
min_coef_bwd = min_coefficients(reg.summaryBackward)

## Plot RSS, adjusted R2, Cp, and BIC for all of the models at once
plot_metrics(reg.summaryBackward)

## Test errors Backward
errors_bwd = test_errors(test_data, regfit.bwd,reg.summaryBackward)

## Taking best MSE for Backward BSS 
best_bwd_mse = as.numeric(errors_bwd[[1]])
best_bwd_name = errors_bwd[[3]]
key <- paste("bwd ", best_bwd_name)
final_mses_dict[[key]] <- best_bwd_mse



######### HYBRID STEPWISE ######### 
cat(bold('Hybrid Stepwise Selection\n'))
regfit.seq=regsubsets(Y~., train_data, nvmax=50, method ="seqrep")
reg.summarySeq=summary(regfit.seq)

## Print min number of regressors for all metrics
min_coef_seq = min_coefficients(reg.summarySeq)

## Plot RSS, adjusted R2, Cp, and BIC for all of the models at once
plot_metrics(reg.summarySeq)

## Test errors hybrid
errors_seq = test_errors(test_data, regfit.seq,reg.summarySeq)

## Taking best MSE for Hybrid BSS 
best_seq_mse = as.numeric(errors_seq[[1]])
best_seq_name = errors_seq[[3]]
key <- paste("seq ", best_seq_name)
final_mses_dict[[key]] <- best_seq_mse



######### RIDGE REGRESSION ######### 
cat(bold('Ridge Regression\n'))
grid=10^seq(6,-2,length=1000)
ridge.mod = glmnet(x[train,], y[train], alpha=0, lambda=grid, thresh=1e-12)
dev.new()
plot_glmnet(ridge.mod, xvar = "lambda")

## best lambda selection
cv.out=cv.glmnet(x[train,], y[train],alpha=0,lambda=grid)
#cv.out=cv.glmnet(x[train,], y[train],alpha=0)
dev.new()
plot(cv.out)
bestlam=cv.out$lambda.min; bestlam; log(bestlam) #best lambda
cv.out$lambda.1se # one standard error rule
log(cv.out$lambda.1se) 
ridge.pred=predict(ridge.mod,s=bestlam ,newx=x[test,])
mse_ridge <- mean((ridge.pred-y[test])^2); mse_ridge

# Refit the ridge regression model to the entire data set
ridge.out <- glmnet(x, y, alpha = 0, lambda = grid) 
ridge.coef <- predict(ridge.out, type = "coefficients", s = bestlam)[1:p + 1, ] # Get the coefficients for the best lambda value
ridge.coef # Print the coefficients

## Taking MSE for Ridge
final_mses_dict[["ridge"]] <- mse_ridge



######### LASSO REGRESSION ######### 
cat(bold('Lasso Regression\n'))
grid=10^seq(6,-2,length=1000)
lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=grid, thresh=1e-12)
dev.new()
plot_glmnet(lasso.mod, xvar = "lambda")

## best lambda selection
cv.out=cv.glmnet(x[train,], y[train],alpha=1,lambda=grid)
dev.new()
plot(cv.out)
bestlam=cv.out$lambda.min; print(bestlam);print(log(bestlam))
print(cv.out$lambda.1se)
print(log(cv.out$lambda.1se))
lasso.pred=predict(lasso.mod,s=cv.out$lambda.1se ,newx=x[test,])
mse_lasso <- mean((lasso.pred-y[test])^2); mse_lasso

# Refit the ridge regression model to the entire data set
#out=glmnet(x[train,],y[train],alpha=1,lambda=grid)
lasso_out <- glmnet(x, y, alpha=1, lambda=grid)
lasso.coef <- predict(lasso_out, type="coefficients", s=cv.out$lambda.1se)[1:p+1, ] # Get the coefficients for the best 1se lambda value
lasso.coef # Print the coefficients
lasso.coef[lasso.coef!=0]
cat("Number of coefficients equal to 0:",sum(lasso.coef==0),"\n")
cat("Number of coefficients:",p-sum(lasso.coef==0),"\n")

## Taking MSE for Lasso
final_mses_dict[["lasso"]] <- mse_lasso



######### ELASTIC NET ######### 
cat(bold('Elastic Net\n'))
grid=10^seq(6,-2,length=1000)
best_mse = NULL
best_lam = NULL
best_enet = NULL
best_alpha = NULL
for (alp in seq(0.05,0.95,by=0.05)) {
  cat(bold('Current Alpha value: ', alp,"\n"))
  enet.mod = glmnet(x[train,], y[train], alpha=alp, lambda=grid)
  #dev.new()
  #plot(enet.mod,label = T)
  #dev.new()
  #plot(enet.mod,label = T, xvar = "lambda")
  #plot_glmnet(enet.mod, xvar = "lambda")
  # best lambda selection
  cv.out=cv.glmnet(x[train,],y[train],alpha=alp,lambda=grid,nfolds = 5)
  #dev.new()
  #plot(cv.out)
  min_lam=cv.out$lambda.min
  lam_1se=cv.out$lambda.1se
  #print(cv.out$lambda.min)
  #print(log(cv.out$lambda.min))
  #print(cv.out$lambda.1se)
  #print(log(cv.out$lambda.1se))
  enet.pred_with_min_lam=predict(enet.mod,s=min_lam ,newx=x[test,])
  enet.pred_with_1selam=predict(enet.mod,s=lam_1se ,newx=x[test,])
  mse_enet_with_min_lam <- mean((enet.pred_with_min_lam-y[test])^2)
  mse_enet_with_1selam <- mean((enet.pred_with_1selam-y[test])^2)
  print(mse_enet_with_min_lam)
  print(mse_enet_with_1selam)
  min_mse = min()
  if (mse_enet_with_min_lam < mse_enet_with_1selam){
    if (is.null(best_mse) || mse_enet_with_min_lam < best_mse){
      best_mse = mse_enet_with_min_lam
      best_lam = min_lam
      best_enet = enet.mod
      best_alpha = alp
    }
  }else{
    if (is.null(best_mse)|| mse_enet_with_1selam < best_mse){
      best_mse = mse_enet_with_min_lam
      best_lam = lam_1se
      best_enet = enet.mod
      best_alpha = alp
    }
  }
}
print(best_lam)
print(best_mse)
print(best_alpha)
enet.coef = predict(best_enet, type="coefficients", s=best_lam)[1:p+1,]
## Taking MSE for Elastic Net
final_mses_dict[["enet"]] <- best_mse


######### BEST TECHNIQUE ######### 
## Taking best MSE between all adopted techniques 
cat(bold('Best MSE between all adopted techniques\n'))
print(final_mses_dict)
values <- unlist(final_mses_dict)
final_mse <- min(values)
best_technique <- names(final_mses_dict)[which.min(final_mses_dict)]

print(final_mse)
print(best_technique)



######### PHRASE FINDER ######### 
cat(bold('Phrase finder\n'))
beta=coef(regfit.fwd ,which.min(reg.summaryForward$cp));
intToUtf8(round(beta/100, digits = 0))


