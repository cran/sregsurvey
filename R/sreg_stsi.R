#' Semiparametric Model-Assisted Estimation under a Stratified Sampling
#' with Simple Random Sampling Without Replace in each stratum.
#'
#' \code{sreg_stsi} is used to estimate the total parameter of a finite population generated from a semi-parametric generalized gamma population
#' under a stratified sampling with simple random sampling without-replacement in each stratum.
#' @param location_formula a symbolic description of the systematic component of the location model to be fitted.
#' @param scale_formula a symbolic description of the systematic component of the scale model to be fitted.
#' @param data a data frame, list containing the variables in the model.
#' @param stratum vector, represents the strata of each unit in the population
#' @param n integer, represents a fixed sample size.
#' @param ss_sizes vector, represents a vector with the sample size in each stratum. By default is a proportional allocation, Sarndal et. al. (2003).
#' @param ... further parameters accepted by caret and survey functions.
#' @return \code{sampling_design} is the name of the sampling design used in the estimation process.
#' @return \code{N} is the population size.
#' @return \code{n} is the fixed sample size used in the estimation process.
#' @return \code{first_order_probabilities} vector of the first order probabilities used in the estimation process.
#' @return \code{sample} is the random sample used in the estimation process.
#' @return \code{total_y_sreg} is the SREG estimate of the total parameter of the finite population.
#' @references Cardozo C.A, Alonso C. (2021) Semi-parametric model assisted estimation in finite populations. In preparation.
#' @references Cardozo C.A.,  Paula G., and Vanegas L. (2021). Generalized log-gamma semiparametric models with P-spline smoothing. Submitted.
#' @references Sarndal C.E.,  Swensson B., and Wretman J. (2003). Model Assisted Survey Sampling. Springer-Verlag.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' library(sregsurvey)
#' library(survey)
#' library(dplyr)
#' library(magrittr)
#' library(gamlss)
#' data(api)
#' attach(apipop)
#' Apipop <- filter(apipop,full!= 'NA')
#' Apipop <- Apipop %>% dplyr::select(api00,grad.sch,full,stype)
#' dim(Apipop)
#' sreg_stsi(api00 ~  pb(grad.sch),scale_formula = ~ full - 1, n=400, stratum = 'stype',data = Apipop)
#' # The total population value is
#' sum(Apipop$api00)
#'
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist GG
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
#' @importFrom TeachingSampling S.SI
#' @importFrom caret knnreg
#' @importFrom stats model.matrix predict
#' @importFrom methods missingArg
#' @export sreg_stsi
sreg_stsi = function(location_formula, scale_formula, stratum, data, n, ss_sizes,...){
  if (missingArg(location_formula))
    'The formula is missing!'
  if (missingArg(scale_formula))
    'The sigma formula is missing!'
  if (missingArg(stratum))
    'The stratum variable is missing!'
  t_y_sreg <- numeric()
  N <- dim(data)[1]

  var_stratum <- data %>% dplyr::select(stratum)
  t_var_stratum <- table(var_stratum)
  val_stratum <- names(t_var_stratum)
  H <- length(val_stratum)
  Ns <- as.numeric(t_var_stratum)
  if (missingArg(ss_sizes)){
    ss_sizes <- round(n*(Ns/N))
      for(i in 1:H){
         if(Ns[i] >= 60){
           if(ss_sizes[i] < 60 ){
             print('The minimum sample size, 60, to guarantee an adecuate fitting process was allocated!')
             ss_sizes[i] <- 60}
         }
         else
            stop('There are strata with very few units, less than 60!')
      }
  }
  fractions <- round(ss_sizes/Ns,2)
  index_acum <- vector()
  totals <- vector()
  samples <- list()
  for(j in 1:H){
    sample <- data %>% filter(var_stratum==val_stratum[j])
    sample <- sample %>% dplyr::select(-stratum)
    subtotal <- sreg_srswr(api00 ~  pb(grad.sch), scale_formula = ~ full - 1, data= sample, fraction=fractions[j])
    if (is.list(subtotal)){
        print(paste('Stratum',j, 'done!'))
        samples <- append(samples,subtotal$sample)
        totals[j] <- subtotal$total_y_sreg
        }
    else{
      return(paste('The regression fitting process was not sucessful in ',j, ' stratum, Try again!'))
    }
  }
  output <- list(sampling_design="Stratified sampling with simple random sampling without replacement in each stratum.", N = N, H = H, Ns = Ns, ss_sizes = ss_sizes, fracs=fractions, total_strata=totals, total_y_sreg=sum(totals))
  class(output) <- "sregsurvey"
  return(output)
}