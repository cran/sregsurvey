#' Semiparametric Model-Assisted Estimation under a Proportional to Size Sampling Design
#'
#' \code{sreg_pips} is used to estimate the total parameter of a finite population generated from a semi-parametric generalized gamma population under a proportional to size without-replacement sampling design.
#' @param location_formula a symbolic description of the systematic component of the location model to be fitted.
#' @param scale_formula a symbolic description of the systematic component of the scale model to be fitted.
#' @param data a data frame, list containing the variables in the model.
#' @param x vector, an auxiliary variable to calculate the inclusion probabilities of each unit.
#' @param n numeric, sample size.
#' @param ... further parameters accepted by caret and survey functions.
#' @references Cardozo C.A, Alonso C. (2021) Semi-parametric model assisted estimation in finite populations. In preparation.
#' @references Cardozo C.A.,  Paula G., and Vanegas L. (2021). Generalized log-gamma semiparametric models with P-spline smoothing. Submitted.
#' @references Sarndal C.E.,  Swensson B., and Wretman J. (2003). Model Assisted Survey Sampling. Springer-Verlag.
#' @return \code{sampling_design} is the name of the sampling design used in the estimation process.
#' @return \code{N} is the population size.
#' @return \code{n} is the sample size used in the estimation process.
#' @return \code{first_order_probabilities} vector of the first order probabilities used in the estimation process.
#' @return \code{sample} is the random sample used in the estimation process.
#' @return \code{estimated_total_y_sreg} is the SREG estimate of the total parameter of the finite population.

#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' library(sregsurvey)
#' library(survey)
#' library(dplyr)
#' library(gamlss)
#' data(api)
#' attach(apipop)
#' Apipop <- filter(apipop,full!= 'NA')
#' Apipop <- filter(Apipop, stype == 'H')
#' Apipop <- Apipop %>% dplyr::select(api00,grad.sch,full,api99)
#' n=ceiling(0.25*dim(Apipop)[1])
#' aux_var <- Apipop %>% dplyr::select(api99)
#' sreg_pips(api00 ~  pb(grad.sch), scale_formula = ~ full - 1, data= Apipop, x= aux_var, n=n)
#' # The total population value is
#' sum(Apipop$api00)
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist GG
#' @importFrom dplyr select filter
#' @importFrom TeachingSampling S.piPS
#' @importFrom caret knnreg
#' @importFrom stats model.matrix predict
#' @importFrom methods missingArg
#' @export sreg_pips
sreg_pips = function(location_formula,scale_formula,data,x,n,...){
  if (missingArg(location_formula))
    'The location_formula is missing!'
  if (missingArg(scale_formula))
    'The scale_formula is missing!'
  if (missingArg(x))
    'The auxiliary variable is missing!'
  if (missingArg(n))
    'The sample size is missing!'
  t_y_sreg <- numeric()
  N <- dim(select(data,as.character(location_formula)[2]))[1]
  index <- S.piPS(n,as.matrix(x))
  new_data <- select(data,-names(x))
  sample <- new_data[index[,1],]
  factors <- 1/index[,2]
  sample <- data.frame(sample,factors)
  fit_sreg <- try(gamlss(location_formula, sigma.formula = scale_formula, family=GG(sigma.link="identity"), data=sample, weights = factors, control=gamlss.control(n.cyc = 200,trace = FALSE)),silent=TRUE)
  cond <- is.list(fit_sreg)
  if(cond==TRUE){
    sigma_est <- fit_sreg$sigma.coefficients
    lambda_est <-  fit_sreg$nu.coefficients*sigma_est
    X_train <- data.frame(model.matrix(fit_sreg,what='mu'),model.matrix(fit_sreg,what='sigma'))
    y_train <-  as.matrix(select(sample,as.character(location_formula)[2]))
    model_knn <- knnreg(X_train,y_train)
    X_pop <- select(new_data,-as.character(location_formula)[2])
    if(length(names(X_train))!=length(names(select(new_data,-as.character(location_formula)[2])))){
      X_pop  <- data.frame(1,select(new_data,-as.character(location_formula)[2]))
      names(X_pop) <- names(X_train)
    }
    eta_est <- predict(model_knn,X_pop)
    X_s <- select(new_data,colnames(model.matrix(fit_sreg, what='sigma')))
    sigma_est <- sigma_est*X_s
    if(abs(lambda_est) < 0.1)
      lambda_est <- 0.1
    mu_est <- (abs(lambda_est)^(2*sigma_est/lambda_est))*eta_est*(gamma((1/lambda_est^2) + (sigma_est/lambda_est)))/gamma(1/(lambda_est^2))
    cond2 <- fit_sreg$nu.coefficients*fit_sreg$sigma.coefficients
    t_y_sreg  <- sum(mu_est) + sum(factors*(select(sample,as.character(location_formula)[2]) - mu_est[index[,1],]))
    output <- list(sampling_design="Proportional to Size", N =N, n=n, first_order_probabilities=pi,  sample=sample, estimated_total_y_sreg=t_y_sreg)
    class(output) <- "sregsurvey"
    return(output)
  }
  else{
    return('The regression fitting process was not sucessful!')
  }
}
