#' Semiparametric Model-Assisted Estimation under a Bernoulli Sampling Design
#'
#' \code{sreg_ber} is used to estimate the total parameter of a finite population generated from a semi-parametric generalized gamma population under a Bernoulli sampling design.
#' @param location_formula a symbolic description of the systematic component of the location model to be fitted.
#' @param scale_formula a symbolic description of the systematic component of the scale model to be fitted.
#' @param data a data frame, list containing the variables in the model.
#' @param pi numeric, represents the first order probability. Default value is 0.5.
#' @return \code{sampling_design} is the name of the sampling design used in the estimation process.
#' @return \code{N} is the population size.
#' @return \code{n} is the random sample size used in the estimation process.
#' @return \code{first_order_probabilities} vector of the first order probabilities used in the estimation process.
#' @return \code{sample} is the random sample used in the estimation process.
#' @return \code{total_y_sreg} is the SREG estimate of the total parameter of the finite population.
#' @references Sarndal C.E.,  Swensson B., and Wretman J. (2003). Model Assisted Survey Sampling. Springer.
#' @references Cardozo C.A.,  Paula G., and Vanegas L. (2021). Generalized log-gamma semiparametric models with P-spline smoothing. Submitted.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' #This example use the data set 'apipop' of the survey package.
#' library(survey)
#' library(dplyr)
#' data(api)
#' attach(apipop)
#' Apipop <- filter(apipop,full!= 'NA')
#' Apipop <- filter(Apipop, stype == 'H')
#' Apipop <- select(Apipop,c(api00,grad.sch,full))
#' sreg_ber(api00 ~  pb(grad.sch), scale_formula = ~ full - 1, data= Apipop, pi=0.25)
#' sum(Apipop$api00)
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist GG
#' @importFrom dplyr select
#' @importFrom TeachingSampling S.BE
#' @importFrom caret knnreg
#' @importFrom stats model.matrix predict
#' @importFrom methods missingArg
#' @export sreg_ber
#'
sreg_ber = function(location_formula,scale_formula,data,pi){
  if (missingArg(location_formula))
      'The location_formula is missing!'
  if (missingArg(scale_formula))
    'The scale_formula is missing!'
  if (missingArg(pi))
      pi=0.1
  t_y_sreg <- numeric()
  N <- dim(select(data,as.character(location_formula)[2]))[1]
  index <- S.BE(N,pi)
  sample <- data[index,]
  n <- dim(sample)[1]
  factors <- rep(1/pi, n)
  sample <- data.frame(sample,factors)
  fit_sreg <- try(gamlss(location_formula, sigma.formula = scale_formula, family=GG(sigma.link="identity"), data=sample, weights = factors, control=gamlss.control(n.cyc = 100)),silent=TRUE)
  cond <- is.list(fit_sreg)
  if(cond==TRUE){
    sigma_est <- fit_sreg$sigma.coefficients
    lambda_est <-  fit_sreg$nu.coefficients*sigma_est
    X_train <- data.frame(model.matrix(fit_sreg,what='mu'),model.matrix(fit_sreg,what='sigma'))
    y_train <-  as.matrix(select(sample,as.character(location_formula)[2]))
    model_knn <- knnreg(X_train,y_train)
    X_pop <- select(data,-as.character(location_formula)[2])
    if(length(names(X_train))!=length(names(select(data,-as.character(location_formula)[2])))){
       X_pop  <- data.frame(1,select(data,-as.character(location_formula)[2]))
       names(X_pop) <- names(X_train)
    }
    eta_est <- predict(model_knn,X_pop)
    X_s <- select(data,colnames(model.matrix(fit_sreg, what='sigma')))
    sigma_est <- sigma_est*X_s
    if(abs(lambda_est) < 0.1)
       lambda_est <- 0.1
    mu_est <- (abs(lambda_est)^(2*sigma_est/lambda_est))*eta_est*(gamma((1/lambda_est^2) + (sigma_est/lambda_est)))/gamma(1/(lambda_est^2))
    cond2 <- fit_sreg$nu.coefficients*fit_sreg$sigma.coefficients
    t_y_sreg  <- sum(mu_est) + sum(factors*(select(sample,as.character(location_formula)[2]) - mu_est[index]))
    output <- list(sampling_design="Bernoulli", N =N, n=n, first_order_probabilities=pi,  sample=sample, total_y_sreg=t_y_sreg)
    class(output) <- "sregsurvey"
    return(output)
  }
  else{
    return('The regression fitting process was not sucessful. Try again!')
  }
}
