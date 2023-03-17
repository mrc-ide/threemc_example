#define TMB_LIB_INIT R_init_threemc
#include <TMB.hpp>

/***************************************************/
/* Function to inverse logit transform for vectors */
/***************************************************/
template <class vector>
vector invlogit_vec(vector x){
  vector y = 1.0 / (1.0 + exp(-x));
  return y;
}
/*******************************************************************/
/* Function to inverse logit transform on a general interval [a,b] */
/*******************************************************************/
template <class Type>
Type geninvlogit(Type x, Type a, Type b){
  Type y; 
  y = 1.0 / (1.0 + exp(-x));
  y = y * (b - a) + a;
  return y;
}
/************************************************************************/
/* Objective function to specify model and to optimize model parameters */
/************************************************************************/
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  ////////////////////////
  /// Data definitions ///
  ////////////////////////
  // Survival analysis matrices
  DATA_SPARSE_MATRIX(A_mmc); // Matrix selecting instantaneous hazard for medically circumcised pop
  DATA_SPARSE_MATRIX(A_tmc); // Matrix selecting instantaneous hazard for traditionally circumcised pop
  DATA_SPARSE_MATRIX(A_mc); // Matrix selecting instantaneous hazard for unknown circumcised pop
  DATA_SPARSE_MATRIX(B); // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
  DATA_SPARSE_MATRIX(C); // Matrix selecting relevant cumulative hazard entry for interval censored pop
  DATA_SPARSE_MATRIX(IntMat1); // Integration matrix for cumulative hazard 
  DATA_SPARSE_MATRIX(IntMat2); // Integration matrix for lagged cumulative hazard 
  
  // Design matrices 
  DATA_SPARSE_MATRIX(X_fixed_mmc); // Design matrix for the fixed effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_time_mmc); // Design matrix for the temporal random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_age_mmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_space_mmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agetime_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agespace_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_spacetime_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_fixed_tmc); // Design matrix for the fixed effects in the traditional circumcision hazard rate
  DATA_SPARSE_MATRIX(X_age_tmc); // Design matrix for the stratification random effects in the traditional circumcision hazard rate
  DATA_SPARSE_MATRIX(X_space_tmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agespace_tmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  
  // Precision matrices 
  DATA_SPARSE_MATRIX(Q_space); // Aggregation matrix for number of circumcisions performed
  
  //////////////////
  /// Parameters ///
  //////////////////
  // Fixed Effects
  PARAMETER_VECTOR(u_fixed_mmc);
  PARAMETER_VECTOR(u_fixed_tmc);
  
  // Age random effect
  PARAMETER_VECTOR(u_age_mmc); 
  PARAMETER_VECTOR(u_age_tmc); 
  
  // Temporal random effects 
  PARAMETER_VECTOR(u_time_mmc);
  
  // Spatial random effects
  PARAMETER_VECTOR(u_space_mmc);
  PARAMETER_VECTOR(u_space_tmc);
  
  // Interactions
  PARAMETER_ARRAY(u_agetime_mmc);
  PARAMETER_ARRAY(u_agespace_mmc);
  PARAMETER_ARRAY(u_spacetime_mmc);
  PARAMETER_ARRAY(u_agespace_tmc);
  
  // Standard deviations 
  PARAMETER(logsigma_age_mmc);       Type sigma_age_mmc       = exp(logsigma_age_mmc);
  PARAMETER(logsigma_time_mmc);      Type sigma_time_mmc      = exp(logsigma_time_mmc);
  PARAMETER(logsigma_space_mmc);     Type sigma_space_mmc     = exp(logsigma_space_mmc);
  PARAMETER(logsigma_agetime_mmc);   Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  PARAMETER(logsigma_agespace_mmc);  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  PARAMETER(logsigma_spacetime_mmc); Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  PARAMETER(logsigma_age_tmc);       Type sigma_age_tmc       = exp(logsigma_age_tmc);
  PARAMETER(logsigma_space_tmc);     Type sigma_space_tmc     = exp(logsigma_space_tmc);
  PARAMETER(logsigma_agespace_tmc);  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);
  
  // Autocorrelation parameters 
  PARAMETER(logitrho_mmc_time1);  Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time2);  Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time3);  Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age1);   Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age2);   Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age3);   Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age1);   Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age2);   Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));
  
  //////////////////////////////////
  /// Prior on the fixed effects ///
  //////////////////////////////////
  // Negative log likelihood definition
  Type nll = Type(0);
  
  // Fixed effects for the medical circumcision rate
  nll -= dnorm(u_fixed_mmc,  Type(0), Type(5), TRUE).sum();
  
  // Fixed effects for the traditional circumcision rate
  nll -= dnorm(u_fixed_tmc, Type(0), Type(5), TRUE).sum();
  
  ////////////////////////////////////////////
  /// Prior on the temporal random effects ///
  ////////////////////////////////////////////
  // AR1 Process 
  nll += AR1(rho_mmc_time1)(u_time_mmc);
  
  // Sum to zero constraint 
  nll -= dnorm(u_time_mmc.sum(), Type(0), Type(0.001) * u_time_mmc.size(), TRUE);
  
  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time_mmc, Type(1), TRUE) + logsigma_time_mmc;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time1, Type(3), Type(3), TRUE); 
  
  ///////////////////////////////////////
  /// Prior on the age random effects ///
  ///////////////////////////////////////
  // AR1 processes
  nll += AR1(rho_mmc_age1)(u_age_mmc);
  nll += AR1(rho_tmc_age1)(u_age_tmc);
  
  // Sum to zero constraint 
  nll -= dnorm(u_age_mmc.sum(), Type(0), Type(0.001) * u_age_mmc.size(), TRUE);
  nll -= dnorm(u_age_tmc.sum(), Type(0), Type(0.001) * u_age_tmc.size(), TRUE);
  
  // Prior on the standard deviation for the aeg random effects
  nll -= dexp(sigma_age_mmc, Type(1), TRUE) + logsigma_age_mmc;
  nll -= dexp(sigma_age_tmc, Type(1), TRUE) + logsigma_age_tmc;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_age1, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age1, Type(3), Type(2), TRUE);
  
  ///////////////////////////////////////////
  /// Prior on the spatial random effects ///
  ///////////////////////////////////////////
  // Gaussian markov random field with prespecified precision matrix 
  nll += GMRF(Q_space)(u_space_mmc);
  nll += GMRF(Q_space)(u_space_tmc);
  
  // Sum to zero constraints
  nll -= dnorm(u_space_mmc.sum(), Type(0), Type(0.001) * u_space_mmc.size(), TRUE);
  nll -= dnorm(u_space_tmc.sum(), Type(0), Type(0.001) * u_space_tmc.size(), TRUE);
  
  // Prior on the standard deviation for the spatial random effects
  nll -= dexp(sigma_space_mmc, Type(1), TRUE) + logsigma_space_mmc;
  nll -= dexp(sigma_space_tmc, Type(1), TRUE) + logsigma_space_tmc;
  
  ///////////////////////////////////////////////
  /// Prior on the interaction random effects ///
  ///////////////////////////////////////////////
  // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
  nll += SEPARABLE(AR1(rho_mmc_time2), AR1(rho_mmc_age2))(u_agetime_mmc);
  nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_age3))(u_agespace_mmc);
  nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_time3))(u_spacetime_mmc);
  nll += SEPARABLE(GMRF(Q_space), AR1(rho_tmc_age2))(u_agespace_tmc);
  
  // Sum-to-zero constraints
  for (int i = 0; i < u_agespace_mmc.cols(); i++) {
    nll -= dnorm(u_agespace_mmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_agetime_mmc.cols(); i++) {
    nll -= dnorm(u_agetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_agetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_spacetime_mmc.cols(); i++) {
    nll -= dnorm(u_spacetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_spacetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_agespace_tmc.cols(); i++) {
    nll -= dnorm(u_agespace_tmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_tmc.col(i).size(), true);
  } 
  
  // Vectorising the interaction
  vector<Type> u_agespace_mmc_v(u_agespace_mmc);
  vector<Type> u_agetime_mmc_v(u_agetime_mmc);
  vector<Type> u_spacetime_mmc_v(u_spacetime_mmc);
  vector<Type> u_agespace_tmc_v(u_agespace_tmc);
  
  // Prior on the standard deviation for the interaction random effects
  nll -= dexp(sigma_agespace_mmc,  Type(1), TRUE) + logsigma_agespace_mmc;
  nll -= dexp(sigma_agetime_mmc,   Type(1), TRUE) + logsigma_agetime_mmc;
  nll -= dexp(sigma_spacetime_mmc, Type(1), TRUE) + logsigma_spacetime_mmc;
  nll -= dexp(sigma_agespace_tmc,  Type(1), TRUE) + logsigma_agespace_tmc;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time2, Type(3), Type(2), TRUE); 
  nll -= dnorm(logitrho_mmc_age2,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_time3, Type(3), Type(2), TRUE); 
  nll -= dnorm(logitrho_mmc_age3,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age2,  Type(3), Type(2), TRUE);
  
  //////////////////////////////
  /// Estimating hazard rate ///
  //////////////////////////////
  // Medical hazard rate
  vector<Type> haz_mmc = X_fixed_mmc * u_fixed_mmc + 
    X_time_mmc * u_time_mmc * sigma_time_mmc + 
    X_space_mmc * u_space_mmc * sigma_space_mmc + 
    X_age_mmc * u_age_mmc * sigma_age_mmc + 
    X_agetime_mmc * u_agetime_mmc_v * sigma_agetime_mmc + 
    X_agespace_mmc * u_agespace_mmc_v * sigma_agespace_mmc + 
    X_spacetime_mmc * u_spacetime_mmc_v * sigma_spacetime_mmc;
  
  // Traditional hazard rate	
  vector<Type> haz_tmc = X_fixed_tmc * u_fixed_tmc +
    X_space_tmc * u_space_tmc * sigma_space_tmc + 
    X_age_tmc * u_age_tmc * sigma_age_tmc + 
    X_agespace_tmc * u_agespace_tmc_v * sigma_agespace_tmc;
  
  // Rates on [0,1] scale 
  haz_tmc = invlogit_vec(haz_tmc);
  haz_mmc = invlogit_vec(haz_mmc);
  
  // Adjustment such that \lambda_mmc + \lambda_tmc \in [0,1]  
  // Medical rate to only take from the remaining proportion 
  // not taken through traditional circumcision (1 - \lambda_tmc)
  haz_mmc = haz_mmc * (1 - haz_tmc);
  
  // Total hazard rate
  vector<Type> haz = haz_mmc + haz_tmc;
  
  // Survival probabilities
  vector<Type> logprob  = log(Type(1.0) - haz);
  vector<Type> surv     = exp(IntMat1 * logprob);
  vector<Type> surv_lag = exp(IntMat2 * logprob);
  vector<Type> leftcens = Type(1.0) - surv;
  
  // Incidence 
  vector<Type> inc_tmc = haz_tmc * surv_lag;
  vector<Type> inc_mmc = haz_mmc * surv_lag;
  vector<Type> inc = haz * surv_lag;
  
  // Cumulative incidence 
  vector<Type> cum_inc_tmc = IntMat1 * inc_tmc;
  vector<Type> cum_inc_mmc = IntMat1 * inc_mmc;
  vector<Type> cum_inc = cum_inc_tmc + cum_inc_mmc;
  
  //////////////////
  /// Likelihood ///
  //////////////////
  // Getting likelihood for those medically circumcised
  nll -= (A_mmc * log(inc_mmc)).sum();
  
  // Getting likelihood for those traditionally circumcised
  nll -= (A_tmc * log(inc_tmc)).sum();
  
  // Getting likelihood for those circumcised of unknown type
  nll -= (A_mc * log(inc)).sum();
  
  // Getting likelihood for those right censored
  nll -= (B * log(surv)).sum();
  
  // Getting likelihood for those left censored
  nll -= (C * log(leftcens)).sum();
  
  ///////////////////////////
  /// Reporting variables ///
  ///////////////////////////
  REPORT(haz_mmc);     // Medical hazard rate
  REPORT(haz_tmc);     // Traditional hazard rate
  REPORT(haz);         // Total hazard rate
  REPORT(inc_tmc);     // Traditional circumcision incidence rate
  REPORT(inc_mmc);     // Medical circumcision incidence rate
  REPORT(inc);         // Total circumcision incidence rate
  REPORT(cum_inc_tmc); // Traditional circumcision cumulative incidence rate
  REPORT(cum_inc_mmc); // Medical circumcision cumulative incidence rate
  REPORT(cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(surv);        // Survival probabilities
  
  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return nll;
}
