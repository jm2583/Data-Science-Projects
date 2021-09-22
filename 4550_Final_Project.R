#LoadStock Data
load("C:/Users/Jonathan Mac/Desktop/4550_Notes/Project_export/ARIMA_data/Stocks.RData")
write.table(Stocks, file = "4550ProjectStockData.csv", sep=",", row.names = FALSE,
            col.names = TRUE)
Stocks #name of data


#Libraries
library(forecast) # For Arima function
library(astsa) # For acf2
#install.packages('rugarch')
library(rugarch) # For model fitting
library(MASS) 
install.packages('EnvStats')
library(EnvStats)

############################ PRE-STEPS ############################

#Original time plot
apple = Stocks[,1]
plot(apple, main = "Apple Stock", ylab = "Price") 


#Transform
log_apple = log(apple) 
d_log_apple = diff(log_apple) #to obtain ARIMA(p,1,q) order
plot(d_log_apple)

sq_d_log_apple = d_log_apple^2 #to obtain GARCH(p,q) order
plot(sq_d_log_apple)


#Plot ACFs to determine order
#acf2 function outputs error "lags exceed number of observations"
acf2(d_log_apple,max.lag=30)

acf2(sq_d_log_apple,max.lag=30)
#lags at h = 1,5,6 exceed blue dotted cut-offs at both ACF and PACF
#thus all are candidates for GARCH(p,q) like GARCH(1,1), GARCH(5,5), GARCH(6,6)

#################################################################

############################ Final Project ############################

### GARCH Model ###



# model specification
spec <- ugarchspec( 
  mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
  # "sGARCH" is the usual GARCH model we defined:
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  # We assume a normal model for e_t: 
  distribution.model="norm"
)

# fit
fit <- ugarchfit(spec=spec, data = d_log_apple, out.sample = 59)
fit

# plot
plot(fit) #1, 8, 9, 10, 11: not normal




# fit student t distribution
ehat <- as.numeric(residuals(fit, standardize=TRUE))
fitdistr(x=ehat, densfun='t') 
qqPlot(x=ehat, distribution='t', 
       param.list=list(df=4), add.line=TRUE)

# refit
spec_t <- ugarchspec(
  mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
  variance.model = list(model = "sGARCH", garchOrder = c(5, 5)),
  distribution.model="std"
)

fit_t <- ugarchfit(spec=spec_t, data = d_log_apple, out.sample = 59)
fit_t





#ARMA(2,2) + GARCH(5,5) without student's t distribution
new_spec_t <- ugarchspec(
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  variance.model = list(model = "sGARCH", garchOrder = c(5, 5)),
  distribution.model="std"
)

#AR(1) + GARCH(1,1) with student's t distribution
new_fit_t <- ugarchfit(spec=new_spec_t, data = d_log_apple, out.sample = 89)
new_fit_t





### OPTIONAL: State Space Model ###

plot(Stocks)

y <- t(rbind(Stocks[-510:-568,], matrix(rep(NA, 59*6), nrow=59))) #510 bc inclusive


# B in MARSS notation is equivalent to Phi in our notation. 
Bvec <- c('phi11', 'phi21', 'phi31','phi41','phi51','phi61', 
          'phi12', 'phi22', 'phi32','phi42','phi52','phi62',  
          'phi13', 'phi23', 'phi33','phi43','phi53','phi63', 
          'phi14', 'phi24', 'phi34','phi44','phi54','phi64', 
          'phi15', 'phi25', 'phi35','phi45','phi55','phi65', 
          'phi16', 'phi26', 'phi36','phi46','phi56','phi66')
B <- matrix(data = Bvec, nrow = 6, ncol= 6)

# Z in MARSS notation is equivalent to A in our notation. 

Zvec <- c(1,0,0,0,0,0, 
          0,1,0,0,0,0, 
          0,0,1,0,0,0,
          0,0,0,1,0,0,
          0,0,0,0,1,0,
          0,0,0,0,0,1)
Z <- matrix(Zvec, 6, 6)

# Q in MARSS notation is equivalent to Q in our notation. 
# Remember that Q is symmetric. 
Qvec <- c('q11', 'q12', 'q13', 'q14', 'q15', 'q16',
          'q12', 'q22', 'q23','q24', 'q25', 'q26', 
          'q13', 'q23', 'q33','q34', 'q35', 'q36',
          'q14', 'q24', 'q34','q44', 'q45', 'q46',
          'q15', 'q25', 'q35','q45', 'q55', 'q56',
          'q16', 'q26', 'q36','q46', 'q56', 'q66')

Q <- matrix(Qvec, 6, 6)

# R in MARSS notation is equivalent to R in our notation. 
# Remember that Q is symmetric. 
Rvec <- c('r11', 'r12', 'r13', 'r14', 'r15', 'r16',
          'r12', 'r22', 'r23','r24', 'r25', 'r26', 
          'r13', 'r23', 'r33','r34', 'r35', 'r36',
          'r14', 'r24', 'r34','r44', 'r45', 'r46',
          'r15', 'r25', 'r35','r45', 'r55', 'r56',
          'r16', 'r26', 'r36','r46', 'r56', 'r66')
R <- matrix(Rvec, 6, 6)


# In MARSS notation, u is a drift or intercept parameter. 
# In our notation, this would be included via the set of 
# exogenous variables.

Uvec <- c(0, 0, 0,0,0,0)
U = matrix(Uvec, 6, 1)

# In MARSS notation, a is an intercept parameter. 
# In our notation, this would be included via the set of 
# exogenous variables.

Avec <- c(0, 0, 0,0,0,0)
A = matrix(Avec, 6, 1)

# We put all these together in a list: 

model.list <- list(B=B,U=U,A=A,Q=Q,Z=Z,R=R,tinitx=0)

# tinitx=0 tells MARSS how we're doing the indexing.

# We fit the models. 
# First we use the EM algorithm, but this doesn't converge. 
# If the EM algorithm doesn't converge, we pass the non-converged 
# estimates into the BFGS algorithm as initial values. 

fit_kem <- MARSS(y=y, model=model.list, method='kem')

fit_bfgs_with_kem <- MARSS(y=y, model=model.list, method = "BFGS", 
                           inits = fit_kem)

MARSSparamCIs(fit_bfgs_with_kem)

plot(fit_bfgs_with_kem)
