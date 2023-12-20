###############################################################################
###                                                                         ###
###                               Appendix A1                               ###
###                                                                         ###
###       Excerpt of code written in Program R to perform statistical       ###
###       population reconstruction of moose (Alces alces) throughout       ###
###       Minnesota using harvest, effort, aerial, and telemetry data       ###
###                                                                         ###
###############################################################################

###############################################################################
############################## Set Up Work Space ##############################
###############################################################################

# Clear global environment and import required packages
{
  rm(list=ls())
  
  require(BB)
  require(pso)
  require(readxl)
  require(ggplot2)
  require(numDeriv)
  require(truncnorm)
}

# Define function to calculate the log of a binomial coefficient
binomial_coeff_log <- function(n, k) {
  log_n <- lgamma(n+1)
  log_k <- lgamma(k+1)
  log_nk <- lgamma(n-k+1)
  return(log_n - (log_k + log_nk))
}

# Define function to calculate the log of a multinomial coefficient
multinomial_coeff_log <- function(n, k) {
  log_n <- lgamma(n + 1)
  log_k <- sum(lgamma(k + 1))
  return(log_n - log_k)
}

###############################################################################
############################### Import Raw Data ###############################
###############################################################################

# Import age-at-harvest matrix for state harvest with each year as a row
{
  stateHarvest <- matrix(c(2004,	0,	24,	127,	245,
                           2005,	0,	27,	136,	284,
                           2006,	0,	28,	133,	269,
                           2007,	0,	0,	115,	229,
                           2008,	0,	0,	110,	245,
                           2009,	0,	0,	103,	223,
                           2010,	0,	0,	109,	212,
                           2011,	0,	0,	53,	103,
                           2012,	0,	0,	46,	87,
                           2013,	0,	0,	0,	0,
                           2014,	0,	0,	0,	0,
                           2015,	0,	0,	0,	0,
                           2016,	0,	0,	0,	0,
                           2017,	0,	0,	0,	0,
                           2018,	0,	0,	0,	0,
                           2019,	0,	0,	0,	0,
                           2020,	0,	0,	0,	0,
                           2021,	0,	0,	0,	0,
                           2022,	0,	0,	0,	0), byrow = T, ncol = 5)

    colnames(stateHarvest) <- c("Year", "Calves", "Cows", "Bulls", "Licenses")
}

# Import age-at-harvest matrix for tribal harvest with each year as a row
{
  tribeHarvest <- matrix(c(2004,	1,	11,	49,	51,
                           2005,	2,	16,	46,	137,
                           2006,	1,	11,	36,	138,
                           2007,	0,	8,	39,	125,
                           2008,	0,	4,	26,	126,
                           2009,	1,	9,	35,	123,
                           2010,	0,	7,	29,	128,
                           2011,	1,	5,	25,	128,
                           2012,	0,	5,	31,	112,
                           2013,	0,	0,	0,	0,
                           2014,	0,	0,	0,	0,
                           2015,	0,	0,	0,	0,
                           2016,	0,	0,	28,	60,
                           2017,	0,	0,	25,	68,
                           2018,	0,	0,	26,	80,
                           2019,	0,	0,	25,	80,
                           2020,	0,	0,	38,	80,
                           2021,	0,	0,	29,	93,
                           2022,	0,	0,	14,	42), byrow = T, ncol = 5)
  
  colnames(tribeHarvest) <- c("Year", "Calves", "Cows", "Bulls", "Licenses")
}

# Import aerial survey data with each year as a row
{
  aerialSurvey <- matrix(c(2005, 8160, 6090, 11410, 0.52, 19, 9, 1.04,
                           2006, 8840, 6790, 11910, 0.34, 13, 5, 1.09,
                           2007, 6860, 5320, 9150, 0.29, 13, 3, 0.89,
                           2008, 7890, 6080, 10600, 0.36, 16, 2, 0.77,
                           2009, 7840, 6270, 10040, 0.32, 14, 2, 0.94,
                           2010, 5700, 4540, 7350, 0.28, 13, 3, 0.83,
                           2011, 4900, 3870, 6380, 0.24, 13, 1, 0.64,
                           2012, 4230, 3250, 5710, 0.36, 15, 6, 1.08,
                           2013, 2760, 2160, 3650, 0.33, 12, 3, 1.23,
                           2014, 4350, 3220, 6210, 0.44, 17, 3, 1.24,
                           2015, 3450, 2610, 4770, 0.29, 13, 3, 0.99,
                           2016, 4020, 3230, 5180, 0.42, 17, 5, 1.03,
                           2017, 3710, 3010, 4710, 0.36, 15, 4, 0.91,
                           2018, 3030, 2320, 4140, 0.37, 15, 4, 1.25,
                           2019, 4180, 3250, 5580, 0.32, 13, 3, 1.24,
                           2020, 3150, 2400, 4320, 0.36, 18, 2, 0.90,
                           2021, NA, NA, NA, NA, NA, NA, NA,
                           2022, 4700, 3440, 6780, 0.45, 19, 3, 0.94,
                           2023, 3290, 2480, 4560, 0.38, 16, 6, 1.26), 
                         byrow = T, ncol = 8)
  
  colnames(aerialSurvey) <- c("Year", "Estimate", "Lo90CI", "Hi90CI", 
                             "Calf:Cow", "%Calves", "%CowsW/Twins", "Bull:Cow")
}

# Import radio telemetry data with each year as a row
{
  telemetryCounts <- matrix(c(2005,  0,  0,  0,  0,  0,  0,
                              2006,  0,  0,  0,  0,  0,  0,
                              2007,  0,  0,  0,  0,  0,  0,
                              2008,  0,  0,  0,  0,  0,  0,
                              2009,  0,  0,  0,  0,  0,  0,
                              2010,  0,  0,  0,  0,  9,  2,
                              2011,  0,  2,  0,  0, 13,  6,
                              2012,  0,  2,  1,  0, 12,  7,
                              2013,  0,  0,  0,  0, 10,  4,
                              2014,  0,  0,  1,  0, 10,  4,
                              2015,  0,  0,  1,  0,  8,  2,
                              2016,  0,  1,  0,  0,  5,  0,
                              2017,  0,  1,  0,  0,  4,  0,
                              2018,  0,  0,  0,  0,  0,  0,
                              2019,  0,  0,  0,  0,  0,  0,
                              2020,  0,  0,  0,  0,  0,  0,
                              2021,  0,  0,  0,  0,  0,  0,
                              2022,  0,  0,  0,  0,  0,  0,
                              2023,  0,  0,  0,  0,  0,  0),
                            byrow = T, ncol = 7)
  
  colnames(telemetryCounts) <- c("Year", "CalvesDied", "CowsDied", "BullsDied",
                                 "CalvesAtRisk","CowsAtRisk","BullsAtRisk")
}

###############################################################################
############################### Synthesize Data ###############################
###############################################################################

# Extract three-tier age-and-sex-at-harvest for both harvest seasons
{
  yearRange <- as.vector(stateHarvest[,"Year"])
  
  h1 <- stateHarvest[,c("Calves", "Cows", "Bulls")]
  h2 <- tribeHarvest[,c("Calves", "Cows", "Bulls")]
  h3 <- h1 + h2
  
  rownames(h1) <- c(yearRange)
  rownames(h2) <- c(yearRange)
  rownames(h3) <- c(yearRange)
}

# Extract estimates of yearly catch-effort and scale to a mean of one
{
  f1 <- as.vector(stateHarvest[,"Licenses"])
  f2 <- as.vector(tribeHarvest[,"Licenses"])
  f3 <- f1 + f2
  
  f1 <- f1 / mean(f1[f1 != 0])
  f2 <- f2 / mean(f2[f2 != 0])
  f3 <- f3 / mean(f3[f3 != 0])
}

# Extract counts and standard errors from raw aerial survey data
{
  a <- cbind((aerialSurvey[,"Calf:Cow"] / 
                (aerialSurvey[,"Calf:Cow"] + aerialSurvey[,"Bull:Cow"] + 1) *
                 aerialSurvey[,"Estimate"]),
             (1 / 
                (aerialSurvey[,"Calf:Cow"] + aerialSurvey[,"Bull:Cow"] + 1) *
                 aerialSurvey[,"Estimate"]),
             (aerialSurvey[,"Bull:Cow"] / 
                (aerialSurvey[,"Calf:Cow"] + aerialSurvey[,"Bull:Cow"] + 1) *
                 aerialSurvey[,"Estimate"]))
  
  rownames(a) <- c(yearRange)
  colnames(a) <- c("Calves", "Cows", "Bulls")
  
  a_se <- (aerialSurvey[,"Hi90CI"] - aerialSurvey[,"Estimate"] + 
           aerialSurvey[,"Estimate"] - aerialSurvey[,"Lo90CI"]) / 2 / 1.645 *
          (a / rowSums(a)) 
}

# Extract counts for number of deaths and animals at risk from telemetry data
{
  v   <- telemetryCounts[,c("CalvesDied", "CowsDied", "BullsDied")]
  n_v <- telemetryCounts[,c("CalvesAtRisk", "CowsAtRisk", "BullsAtRisk")]
  
  rownames(v)   <- c(yearRange)
  rownames(n_v) <- c(yearRange)
}

## Define number of years and age-and-sex classes
{
  Y <- nrow(h1)
  A <- ncol(h1)
  yearRange <- as.vector(stateHarvest[,"Year"])
}

# Set reasonable limits for survival and vulnerability
{
  limitForC <- c(0.001, 0.200)
  limitForS <- c(0.200, 0.999)
}

###############################################################################
########################## Define Objective Function ##########################
###############################################################################

# Define objective function using a multinomial likelihood formulation
objectiveFunction <- function(par) {
  
  ## Import initial (diagonal) cohort values
  {
    N <- matrix(NA, nrow = Y, ncol = A)
    N[1, 1:A] <- par[1:A]
    N[2:Y, 1] <- par[(A + 1):(A + Y - 1)]
  }
  
  ## Import vulnerability and survival estimates
  {
    C1 <- matrix(par[(Y + 2*A - 1)], nrow = Y, ncol = A, byrow = T)
    C2 <- matrix(par[(2*Y + 3*A - 2)], nrow = Y, ncol = A, byrow = T)
    S  <- matrix(par[((3*Y + 3*A - 2) + A - 1)], nrow = Y, ncol = A, byrow = T)
    
    if (separateVulnerabilityByClass == 1) {
      C1[] <- rep(par[(Y + A) : (Y + 2*A - 1)], each = Y)
      C2[] <- rep(par[(2*Y + 2*A - 1):(2*Y + 3*A - 2)], each = Y)
    }
    if (separateVulnerabilityByYears == 1) {
      C1[2:Y,3] <- par[(Y + 2*A) : (2*Y + 2*A - 2)]
      vulnerabilityRatio1 <- C1[,3] / C1[1,3]
      for (i in 1:(A-1)){
        C1[,i] = C1[1,i] * vulnerabilityRatio1
      }
      
      C2[2:Y,3] <- par[(2*Y + 3*A - 1) : ((2*Y + 3*A - 1) + Y - 2)] 
      vulnerabilityRatio2 <- C2[,3] / C2[1,3]
      for (i in 1:(A-1)){
        C2[,i] <- C2[1,i] * vulnerabilityRatio2
      }
    }
    if (separateSurvivabilityByClass == 1) {
      S[,] <- rep(par[(3*Y + 3*A - 2):((3*Y + 3*A - 2) + A - 1)], each = Y)
    }
    if (separateSurvivabilityByYears == 1) {
      S[2:Y,3] <- par[(3*Y + 4*A - 2):((3*Y + 4*A - 2) + Y-2)]
      survivabilityRatio <- S[,3] / S[1,3]
      for (i in 1:(A-1)){
        S[,i] <- S[1,i] * survivabilityRatio
      }
    }
    
    C1 <- C1 / scaleFactor
    C2 <- C2 / scaleFactor
    S  <- S  / scaleFactor
  }
  
  ## Define separate harvest probabilities for state and tribal seasons
  {
    P1 <- (1 - exp(-C1 * f1)) # State
    P2 <- (1 - exp(-C2 * f2)) # Tribal
    
    if (combineStateAndTribalHarvest) {
      P1[] <- (1 - exp(-C1 * f3))
      P2[] <- 0 
      h1[] <- h3 ######
      h2[] <- 0 ############
    }
  }
  
  ## Exclude calves from both state and tribal harvest seasons
  {
    P1[,1] <- 0
    P2[,1] <- 0
  }
  
  # Exclude years of no harvest from parameter optimization
  {
    P1[which(f1 == 0),] <- 0
    P2[which(f2 == 0),] <- 0
  }
  
  ## Define expected population sizes (assuming 50:50 sex ratio among calves)
  {
    for (i in 2:Y) {
      for (j in 2:A) {
        N[i,j] <- N[i-1,j] * (1 - P1[i-1,j]) * (1 - P2[i-1,j]) * S[i-1,j] +
                  N[i-1,1] * (1 - P1[i-1,1]) * (1 - P2[i-1,1]) * S[i-1,1] * 0.5
      }
    }
  }
  
  # Return population size estimate (if requested)
  {
    if (returnPopulationAbundance) return (N)
  }
  
  # Perform initial test of model boundaries
  {
    #if (any(N < h1, 
    #        N < h2))               return (9000009)
    #if (any(C1 < limitForC[1],
    #        C2 < limitForC[1],
    #        C1 > limitForC[2],
    #        C2 > limitForC[2]))    return (9000008)
    #if (any(S  < limitForS[1], 
    #        S  > limitForS[2]))    return (9000007)
    #if (any(P1 + (1 - S) > 0.500,
    #        P2 + (1 - S) > 0.500)) return (9000006)
  }
  
  # Define expected harvest values
  {
    E1 <- N * P1
    E2 <- N * (1 - P1) * P2
  }
  
  ## Return chi-square estimate for inflation factor if requested
  {
    if (returnChiSquareCorrection) {
      return (sum(((h1 - E1) ^ 2 / E1)[-which((h1 - E1) ^ 2 / E1 == Inf)], 
                  na.rm = T) + 
              sum(((h2 - E2) ^ 2 / E2)[-which((h2 - E2) ^ 2 / E2 == Inf)], 
                  na.rm = T))
    }
  }
  
  ## Calculate contributions of each likelihood component
  {
    ### Initialize the age-at-harvest and catch-effort components
    {
      logL_AAH <- matrix(0, nrow = Y, ncol = A)
      logL_ABN <- matrix(0, nrow = Y, ncol = A)
      logL_TEL <- matrix(0, nrow = Y, ncol = A)
    }
    
    ### Contribution of the age-at-harvest component
    {
      for (i in 1:Y) {
        for (j in 1:A) {
          if (P1[i,j] > 0) {
            logL_AAH[i,j] <- logL_AAH[i,j] + 
                             binomial_coeff_log(N[i,j]          , h1[i,j]) +
                             (h1[i,j]) * log(P1[i,j])
          }
          if (P2[i,j] > 0) {
            logL_AAH[i,j] <- logL_AAH[i,j] + 
                             binomial_coeff_log(N[i,j] - h1[i,j], h2[i,j]) +
                             (h2[i,j]) * log((1 - P1[i,j]) * P2[i,j])
          }
          logL_AAH[i,j] <- logL_AAH[i,j] + 
                           (N[i,j] - h1[i,j] - h2[i,j]) * log(1 - (P1[i,j] + (1 - P1[i,j]) * P2[i,j]))
        }
      }
    }
    
    ### Contribution of the annual-abundance-estimate component
    {
      for (i in 1:Y) {
        for (j in 1:A) {
          logL_ABN[i,j] <- -0.5 * (log(2*pi)) + log(a_se[i,j]) +
                           -0.5 * (N[i,j] * (1 - P1[i,j]) * 
                                     (1 - P2[i,j]) - a[i,j]) ^ 2 / 
                                  (a_se[i,j]) ^ 2
        }
      }
      
      #### Account for year when abundance estimates were unavailable
      logL_ABN[17,] <- 0
    }
    
    ### Contribution of the radio-telemetry component
    {
      for (i in 1:Y) {
        for (j in 1:A) {
          logL_TEL[i,j] <- binomial_coeff_log(n_v[i,j], v[i,j]) + 
                           (v[i,j])            * log(1 - S[i,j]) + 
                           (n_v[i,j] - v[i,j]) * log(    S[i,j])
        }
      }
    }
  }
  
  ## Return value of objective function if valid
  {
    logLikelihood <- -c(logL_AAH, logL_ABN, logL_TEL)
    
    if (any(is.na(logLikelihood)))     return (9000003) else
      if (any(logLikelihood == -Inf))  return (9000002) else
        if (any(logLikelihood == Inf)) return (9000001) else
          return(sum(logLikelihood))
  }
}

###############################################################################
############### Perform NUMERICAL OPTIMIZATION of Observed Data ###############
###############################################################################

# Initialize data frames to store results and starting values for optimization
{
  pointEstimates <- data.frame(matrix(NA, nrow = 2 ^ 5, ncol = (A + Y - 1) * 4 + 5 + 1))
  
  pointEstimates <- setNames(pointEstimates, c("N_1_1", "N_1_2", "N_1_3", 
                                               paste(paste("N_", c(2:Y), sep = ""), "_1", sep = ""),
                                               "C_S_1_1", "C_S_1_2", "C_S_1_3", 
                                               paste(paste("C_S_", c(2:Y), sep = ""), "_3", sep = ""),
                                               "C_T_1_1", "C_T_1_2", "C_T_1_3", 
                                               paste(paste("C_T_", c(2:Y), sep = ""), "_3", sep = ""),
                                               "S_1_1", "S_1_2", "S_1_3", 
                                               paste(paste("S_", c(2:Y), sep = ""), "_3", sep = ""),
                                               "M_C","M_C_Class", "M_C_Years", "M_S_Class", "M_S_Years", 
                                               "VAL"))
  
  populationSize <- data.frame(matrix(NA, nrow = 2 ^ 5, ncol = Y))
}

# Initialize starting values for parameterization
{
  scaleFactor <- 1e5 # Helps with optimization to standardize parameter scales
  initialValues <- c(as.numeric(a[1, 1:A]), as.numeric(a[c(2:16,18,18:Y), 1]),
                     rep(mean(h1 / a, na.rm = T), Y + A - 1) * scaleFactor,
                     rep(mean(h2 / a, na.rm = T), Y + A - 1) * scaleFactor,
                     rep(1 - sum(v)/sum(n_v),     Y + A - 1) * scaleFactor)
}

# Optimize parameter space using numerical optimization
{
  ## Initialize model counter
  {
    modelCount <- 1
  }
  
  ## Cycle through possible vulnerability and survival delineations
  for (combineStateAndTribalHarvest in 0:1) {
    for (separateVulnerabilityByClass in c(0:1)) {
      for (separateVulnerabilityByYears in c(0:1)) {
        for (separateSurvivabilityByClass in c(0:1)) {
          for (separateSurvivabilityByYears in c(0:1)) {
            
            ### Optimize objective function using particle swarm optimization
            {
              returnPopulationAbundance <- FALSE
              returnChiSquareCorrection <- FALSE
              optimized <- psoptim(par = initialValues, fn = objectiveFunction,
                                   lower = c(rep(c(00000), A + Y - 1), 
                                             rep(limitForC[1] * scaleFactor, Y + A - 1),
                                             rep(limitForC[1] * scaleFactor, Y + A - 1),
                                             rep(limitForS[1] * scaleFactor, Y + A - 1)),
                                   upper = c(rep(c(10000), A + Y - 1), 
                                             rep(limitForC[2] * scaleFactor, Y + A - 1),
                                             rep(limitForC[2] * scaleFactor, Y + A - 1),
                                             rep(limitForS[2] * scaleFactor, Y + A - 1)),
                                   control = list(maxit = 1000, trace = T, s = 250,
                                                  hybrid = "improved",
                                                  vectorize = TRUE,
                                                  type = "SPSO2011"))
              
              returnPopulationAbundance <- TRUE
              estimatedN <- objectiveFunction(optimized$par)
              pointEstimates[modelCount, ] <- c(optimized$par,
                                                combineStateAndTribalHarvest,
                                                separateVulnerabilityByClass,
                                                separateVulnerabilityByYears,
                                                separateSurvivabilityByClass,
                                                separateSurvivabilityByYears,
                                                optimized$value)
              populationSize[modelCount, ] <- rowSums(estimatedN)
              returnPopulationAbundance <- FALSE
              
              #### Plot intermediate results
              {
                plot(x = yearRange, y = rowSums(estimatedN), type = "b", ylim = c(0,10000),
                     main = paste(combineStateAndTribalHarvest,
                                  separateVulnerabilityByClass,
                                  separateVulnerabilityByYears,
                                  separateSurvivabilityByClass,
                                  separateSurvivabilityByYears))
                points(x = yearRange, y = rowSums(a), type = "l")
              }
            }
            
            ### Increment model counter
            {
              modelCount <- modelCount + 1
            }
          }
        }
      }
    }
  }
  
  ## Compute k for each reconstruction estimate
  {
    pointEstimates$k[pointEstimates$M_C == 0] <- (A + Y - 1) + 
      (pointEstimates$M_C_Class[pointEstimates$M_C == 0] * (A - 1 - 1)) + 1 + 
      (pointEstimates$M_C_Class[pointEstimates$M_C == 0] * (A - 1 - 1)) + 1 +
      (pointEstimates$M_C_Years[pointEstimates$M_C == 0] * (Y - 1 - sum(f1 == 0))) +
      (pointEstimates$M_C_Years[pointEstimates$M_C == 0] * (Y - 1 - sum(f2 == 0))) +
      (pointEstimates$M_S_Class[pointEstimates$M_C == 0] * (A - 1)) + 1 +
      (pointEstimates$M_S_Years[pointEstimates$M_C == 0] * (Y - 1))
    
    pointEstimates$k[pointEstimates$M_C == 1] <- (A + Y - 1) + 
      (pointEstimates$M_C_Class[pointEstimates$M_C == 1] * (A - 1 - 1)) + 1 + 
      (pointEstimates$M_C_Years[pointEstimates$M_C == 1] * (Y - 1 - sum(f3 == 0))) +
      (pointEstimates$M_S_Class[pointEstimates$M_C == 1] * (A - 1)) + 1 +
      (pointEstimates$M_S_Years[pointEstimates$M_C == 1] * (Y - 1))
  }
  
  ## Compute AIC for each reconstruction estimate
  {
     pointEstimates$AIC <- 2 * pointEstimates$VAL + 
                           2 * pointEstimates$k
  }
}

# Extract uncertainty estimates for best-fit model
{
  ## Construct data frame to compile subsequent results
  {
    bestFitModel <- data.frame(Year = yearRange,
                                 Estimate = NA,
                                 Lo95CI = NA, Hi95CI = NA,
                                 Measure = rep(c("Abundance", "Recruitment"),
                                               each = length(yearRange)))
  }

  ## Calculate standard errors for best-fit model parameters
  {
    bestModel <- which.min(pointEstimates$AIC)
    combineStateAndTribalHarvest <- pointEstimates$M_C[bestModel]
    separateVulnerabilityByClass <- pointEstimates$M_C_Class[bestModel]
    separateVulnerabilityByYears <- pointEstimates$M_C_Years[bestModel]
    separateSurvivabilityByClass <- pointEstimates$M_S_Class[bestModel]
    separateSurvivabilityByYears <- pointEstimates$M_S_Years[bestModel]
    bestParameters <- as.numeric(pointEstimates[bestModel,
                                                1:(ncol(pointEstimates) - 5 - 3)])

    returnPopulationAbundance <- FALSE
    hessianMatrix <- abs(hessian(x = bestParameters,
                                 func = objectiveFunction))

    removedIndices <- which(rowSums(hessianMatrix) < 1e-10)
    if (length(removedIndices) > 0) {
      hessianMatrix <- hessianMatrix[-removedIndices, -removedIndices]
    }
    
    standardErrors <- sqrt(abs(diag(solve(hessianMatrix))))
    
    inflatedErrors <- rep(0, length(bestParameters))
    inflatedErrors[removedIndices] <- NA
    
    returnChiSquareCorrection <- TRUE
    inflatedFactor <- sqrt(objectiveFunction(unlist(bestParameters)) / 
                             pointEstimates$k[bestModel])
    returnChiSquareCorrection <- FALSE
    
    inflatedErrors[!is.na(inflatedErrors)] <- standardErrors * inflatedFactor
    inflatedErrors[is.na(inflatedErrors)] <- 0
    
    bestParameters[removedIndices] <- NA
  }

  ## Calculate stochastic abundance estimates using standard errors
  {
    abundanceProjection <- matrix(NA, nrow = 10000, ncol = Y)
    recruiterProjection <- matrix(NA, nrow = 10000, ncol = Y)
    tempParameters <- bestParameters

    returnPopulationAbundance <- TRUE

    for (i in 1:nrow(abundanceProjection)) {
      for (j in 1:(A + Y - 1)) {
        tempParameters[j] <- rtruncnorm(1, a = 0,
                                        mean = bestParameters[j],
                                        sd = inflatedErrors[j])
      }
      for (j in (A + Y - 1 + 1):(4 * (A + Y - 1))) {
        if (!is.na(tempParameters[j])) {
          tempParameters[j] <- rtruncnorm(1, a = 0, b = 1.0 * scaleFactor,
                                          mean = bestParameters[j],
                                          sd = inflatedErrors[j])
        }
      }

      returnPopulationAbundance <- TRUE
      abundanceProjection[i,] <- rowSums(objectiveFunction(tempParameters))
      recruiterProjection[i,] <- (objectiveFunction(tempParameters))[,1]
    }
  }

  ## Compute confidence intervals
  {
    bestFitModel$Estimate[bestFitModel$Measure == "Abundance"] <-
      as.numeric(populationSize[bestModel,])

    bestFitModel$Lo95CI[bestFitModel$Measure == "Abundance"] <-
      apply(abundanceProjection, 2, quantile, prob = 0.025)

    bestFitModel$Hi95CI[bestFitModel$Measure == "Abundance"] <-
      apply(abundanceProjection, 2, quantile, prob = 0.975)

    bestFitModel$Estimate[bestFitModel$Measure == "Recruitment"] <-
      as.numeric(pointEstimates[bestModel,c(1, (A - 1) + c(2:Y))])

    bestFitModel$Lo95CI[bestFitModel$Measure == "Recruitment"] <-
      apply(recruiterProjection, 2, quantile, prob = 0.025)

    bestFitModel$Hi95CI[bestFitModel$Measure == "Recruitment"] <-
      apply(recruiterProjection, 2, quantile, prob = 0.975)
  }
}

###############################################################################
############## Generate VISUALIZATION of Best-Fit Reconstruction ##############
###############################################################################

# Generate figure of pre-harvest abundance and recruitment estimates
{
  tiff("Figure 1.tiff", units = "in", width = 8, height = 5, res = 400)

  ggplot(aes(y = Estimate, x = Year, fill = Measure),
         data = bestFitModel) +
    geom_ribbon(aes(ymin = Lo95CI, ymax = Hi95CI,
                    fill = Measure), alpha = 0.25) +
    geom_point(aes(colour = Measure)) +
    geom_path(aes(colour = Measure, linetype = Measure), size = 1) +
    scale_x_continuous(breaks = seq(min(yearRange), max(yearRange), 2)) +
    scale_y_continuous(breaks = seq(round(min(bestFitModel$Lo95CI), -3),
                                    round(max(bestFitModel$Hi95CI), -3),
                                    2000)) +
    scale_fill_manual(values = c("#F8766D", "#00BA38")) +
    scale_color_manual(values = c("#F8766D", "#00BA38"))

  dev.off()
}

###############################################################################
############### Synthesize ANALYSIS of Reconstruction Estimates ###############
###############################################################################

# Compile summary statistics under the likelihood objective function
{
  # lambda <- ((bestFitModel$Estimate[bestFitModel$Measure == "Abundance"])[Y] /
  #            (bestFitModel$Estimate[bestFitModel$Measure == "Abundance"])[1]) ^
  #           (1 / (Y - 1))
  # 
  # print(paste("      BEST-FIT MODEL USING PARTICLE SWARM OPTIMIZATION      "))
  # print(paste("------------------------------------------------------------"))
  # 
  # print(paste("The model detected", vulnerability, 
  #             "separate vulnerability coefficient(s):", sep = " "))
  # print(paste(round(bestParameters[(A + Y - 1 + 1:A)] / scaleFactor, 3),
  #             " (SD = ",
  #             round(standardErrors[(A + Y - 1 + 1:A)] / scaleFactor, 3),
  #             ")",
  #             sep = ""))
  # 
  # print(paste("These corresponded to harvest rates between ",
  #             round((1 - exp(-mean(bestParameters[(A + Y - 1 + 1:A)] / 
  #                                    scaleFactor, na.rm = TRUE) * 
  #                            min(f))) * 100, 1),
  #             "% and ",
  #             round((1 - exp(-mean(bestParameters[(A + Y - 1 + 1:A)] / 
  #                                    scaleFactor, na.rm = TRUE) * 
  #                            max(f))) * 100, 1),
  #             "%", sep = ""))
  # 
  # print(paste("The model detected", survivability,
  #             "separate survival coefficient(s):", sep = " "))
  # print(paste(round(bestParameters[(A + Y - 1 + A + 1:A)] / scaleFactor, 3),
  #             " (SD = ",
  #             round(standardErrors[(A + Y - 1 + A + 1:A)] / scaleFactor, 3),
  #             ")",
  #             sep = ""))
  # 
  # if (lambda > 1) {
  #   print(paste("The model detected a positive annual growth rate of",
  #               round(lambda, 3), sep = " "))
  # } else {
  #   print(paste("The model detected a negative annual growth rate of",
  #               round(lambda, 3), sep = " "))
  # }
  # 
  # print(paste("------------------------------------------------------------"))
}


###############################################################################
###############################################################################
####################################################################  __     ##
#################################################################### /  \__  ##
#################################################################### \__/  \ ##
#################################################################### /  \__/ ##
#################################################################### \__/    ##
####################################################################         ##
###############################################################################