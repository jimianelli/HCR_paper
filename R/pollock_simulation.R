#' Bering Sea Pollock-like Population Simulation
#' Age-structured population model with Beverton-Holt stock-recruitment
#'
#' @description
#' This script simulates a pollock-like population over 50 years with:
#' - Ages 1-15 (plus group)
#' - Beverton-Holt recruitment (steepness h = 0.6)
#' - Logistic maturity (a50 = 3.6 years)
#' - Logistic fishery selectivity (a50 = 4 years)
#' - Log-normal recruitment variability (CV = 0.7)
#' - Female-only SSB (50:50 sex ratio assumed)

# =============================================================================
# Model Parameters
# =============================================================================

#' Set up simulation parameters
#' @return List of model parameters
set_parameters <- function() {

  params <- list(
    # Simulation dimensions
    n_years = 50,
    n_ages = 15,
    ages = 1:15,

    # Biological parameters
    M = 0.3,                    # Natural mortality

    # Empirical weight-at-age (kg)
    weight_at_age_vec = c(
      0.028664353,   # Age 1
      0.179556496,   # Age 2
      0.371928063,   # Age 3
      0.491302938,   # Age 4
      0.613339000,   # Age 5
      0.749083313,   # Age 6
      0.883110688,   # Age 7
      1.013851875,   # Age 8
      1.132731875,   # Age 9
      1.243061438,   # Age 10
      1.347911625,   # Age 11
      1.388225625,   # Age 12
      1.468646938,   # Age 13
      1.552769000,   # Age 14
      1.790591875    # Age 15+
    ),

    # Sex ratio for SSB calculation (proportion female)
    prop_female = 0.5,

    # Maturity at age (females)
    maturity_at_age_vec = c(
      0.000,  # Age 1
      0.008,  # Age 2
      0.289,  # Age 3
      0.641,  # Age 4
      0.842,  # Age 5
      0.901,  # Age 6
      0.947,  # Age 7
      0.963,  # Age 8
      0.970,  # Age 9
      1.000,  # Age 10
      1.000,  # Age 11
      1.000,  # Age 12
      1.000,  # Age 13
      1.000,  # Age 14
      1.000   # Age 15+
    ),

    # Selectivity ogive (logistic)
    sel_a50 = 4.0,              # Age at 50% selectivity
    sel_slope = 2.0,            # Slope parameter

    # Stock-recruitment (Beverton-Holt)
    h = 0.6,                    # Steepness
    R0 = 13.5e9,                # Virgin recruitment - calibrated for Bmsy ~ 2.3 million t

    # Recruitment variability
    rec_cv = 0.7,               # CV of recruitment

    # Initial fishing mortality
    F_init = 0.3                # Initial F value for simulation
  )

  # Calculate derived parameters
  params$rec_sigma <- sqrt(log(1 + params$rec_cv^2))  # SD on log scale

  return(params)
}

# =============================================================================
# Biological Functions
# =============================================================================

#' Get weight at age from empirical vector
#' @param age Vector of ages
#' @param params Parameter list
#' @return Vector of weights (kg)
weight_at_age <- function(age, params) {
  params$weight_at_age_vec[age]
}

#' Calculate maturity at age (logistic ogive)
#' @param age Vector of ages
#' @param params Parameter list
#' @return Vector of maturity proportions
maturity_at_age <- function(age, params) {
  if (!is.null(params$maturity_at_age_vec)) {
    return(params$maturity_at_age_vec[age])
  }

  1 / (1 + exp(-params$mat_slope * (age - params$mat_a50)))
}

#' Calculate fishery selectivity at age (logistic ogive)
#' @param age Vector of ages
#' @param params Parameter list
#' @return Vector of selectivity proportions
selectivity_at_age <- function(age, params) {
  1 / (1 + exp(-params$sel_slope * (age - params$sel_a50)))
}

#' Calculate spawning stock biomass (females only)
#' @param N Vector of numbers at age (total population)
#' @param params Parameter list
#' @return SSB value (female biomass only)
calc_ssb <- function(N, params) {
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  # SSB is female spawning biomass only (assume 50:50 sex ratio)
  sum(N * wt * mat * params$prop_female)
}

# =============================================================================
# Stock-Recruitment Function
# =============================================================================

#' Beverton-Holt stock-recruitment relationship
#' @param SSB Spawning stock biomass (females only)
#' @param params Parameter list
#' @param SSB0 Virgin spawning stock biomass (females only)
#' @return Expected recruitment (total, both sexes)
beverton_holt <- function(SSB, params, SSB0) {
  h <- params$h
  R0 <- params$R0

  # Standard Beverton-Holt parameterization
  # R = (4*h*R0*SSB) / (SSB0*(1-h) + SSB*(5*h-1))
  alpha <- (4 * h * R0) / (SSB0 * (1 - h))
  beta <- (5 * h - 1) / (SSB0 * (1 - h))

  R <- (alpha * SSB) / (1 + beta * SSB)
  return(R)
}

#' Add log-normal process error to recruitment
#' @param R_det Deterministic recruitment
#' @param params Parameter list
#' @param bias_correct Apply bias correction (default TRUE)
#' @return Recruitment with process error
add_recruitment_error <- function(R_det, params, bias_correct = TRUE) {
  sigma <- params$rec_sigma
  epsilon <- rnorm(1, 0, sigma)

  if (bias_correct) {
    # Bias correction: exp(epsilon - sigma^2/2)
    R <- R_det * exp(epsilon - sigma^2 / 2)
  } else {
    R <- R_det * exp(epsilon)
  }

  return(R)
}

# =============================================================================
# Equilibrium Calculations
# =============================================================================

#' Calculate numbers at age at equilibrium for given F
#' @param F Fishing mortality
#' @param params Parameter list
#' @return Vector of numbers at age per recruit
calc_N_equilibrium <- function(F, params) {
  n_ages <- params$n_ages
  sel <- selectivity_at_age(params$ages, params)
  M <- params$M

  # Start with recruitment = 1 (will scale later)
  N <- numeric(n_ages)
  N[1] <- 1

  for (a in 2:(n_ages - 1)) {
    Z <- M + F * sel[a - 1]
    N[a] <- N[a - 1] * exp(-Z)
  }

  # Plus group: geometric series
  Z_last <- M + F * sel[n_ages - 1]
  Z_plus <- M + F * sel[n_ages]
  N[n_ages] <- N[n_ages - 1] * exp(-Z_last) / (1 - exp(-Z_plus))

  return(N)
}

#' Calculate SSB per recruit at equilibrium (females only)
#' @param F Fishing mortality
#' @param params Parameter list
#' @return SSB per recruit (female biomass)
calc_spr <- function(F, params) {
  N_per_R <- calc_N_equilibrium(F, params)
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  # Female SSB per recruit
  spr <- sum(N_per_R * wt * mat * params$prop_female)
  return(spr)
}

#' Calculate virgin SSB (SSB0) - females only
#' @param params Parameter list
#' @return SSB0
calc_SSB0 <- function(params) {
  spr0 <- calc_spr(0, params)
  SSB0 <- params$R0 * spr0
  return(SSB0)
}

#' Calculate equilibrium recruitment for given F
#' @param F Fishing mortality
#' @param params Parameter list
#' @return Equilibrium recruitment
calc_R_eq <- function(F, params) {
  h <- params$h
  R0 <- params$R0
  spr <- calc_spr(F, params)
  spr0 <- calc_spr(0, params)
  SSB0 <- R0 * spr0

  alpha <- (4 * h * R0) / (SSB0 * (1 - h))
  beta <- (5 * h - 1) / (SSB0 * (1 - h))

  if (alpha * spr > 1) {
    R_eq <- (alpha * spr - 1) / (beta * spr)
  } else {
    R_eq <- 0
  }

  return(R_eq)
}

#' Calculate equilibrium SSB for given F (females only)
#' @param F Fishing mortality
#' @param params Parameter list
#' @return Equilibrium SSB
calc_SSB_eq <- function(F, params) {
  spr <- calc_spr(F, params)
  R_eq <- calc_R_eq(F, params)
  SSB_eq <- R_eq * spr
  return(SSB_eq)
}

#' Calculate equilibrium yield for given F
#' @param F Fishing mortality
#' @param params Parameter list
#' @return Equilibrium yield
calc_yield_eq <- function(F, params) {
  R_eq <- calc_R_eq(F, params)

  if (R_eq == 0) {
    return(0)
  }

  # Calculate yield per recruit
  N_per_R <- calc_N_equilibrium(F, params)
  sel <- selectivity_at_age(params$ages, params)
  wt <- weight_at_age(params$ages, params)
  M <- params$M

  # Baranov catch equation
  ypr <- 0
  for (a in 1:params$n_ages) {
    Z <- M + F * sel[a]
    ypr <- ypr + N_per_R[a] * wt[a] * (F * sel[a] / Z) * (1 - exp(-Z))
  }

  yield <- R_eq * ypr
  return(yield)
}

# =============================================================================
# Reference Points (Fmsy, Bmsy)
# =============================================================================

#' Find Fmsy and Bmsy using numerical optimization
#' @param params Parameter list
#' @return List with Fmsy, Bmsy, MSY
calc_reference_points <- function(params) {

  # Optimize to find Fmsy
  opt <- optimize(
    f = function(F) -calc_yield_eq(F, params),
    interval = c(0.001, 2),
    tol = 1e-6
  )

  Fmsy <- opt$minimum
  MSY <- -opt$objective
  Bmsy <- calc_SSB_eq(Fmsy, params)
  SSB0 <- calc_SSB0(params)
  R_msy <- calc_R_eq(Fmsy, params)

  return(list(
    Fmsy = Fmsy,
    Bmsy = Bmsy,
    MSY = MSY,
    SSB0 = SSB0,
    Bmsy_B0_ratio = Bmsy / SSB0,
    R_msy = R_msy
  ))
}

# =============================================================================
# SPR Reference Points and Control Rules
# =============================================================================

#' Calculate F that achieves a target SPR fraction (e.g., F40%)
#' @param params Parameter list
#' @param spr_target Target SPR fraction (e.g., 0.4 for F40%)
#' @param F_upper_max Maximum F to search
#' @return F value that matches the target SPR, or NA if not found
calc_F_spr <- function(params, spr_target = 0.4, F_upper_max = 10) {
  spr0 <- calc_spr(0, params)
  target <- spr_target * spr0

  f_obj <- function(F) calc_spr(F, params) - target

  lower <- 0
  upper <- 2
  f_upper <- f_obj(upper)

  while (f_upper > 0 && upper < F_upper_max) {
    upper <- upper * 2
    f_upper <- f_obj(upper)
  }

  if (f_upper > 0) {
    warning("Target SPR not reached within search bounds.")
    return(NA_real_)
  }

  uniroot(f_obj, lower = lower, upper = upper)$root
}

#' Sloping harvest control rule based on SSB and B0
#' @param SSB Vector of spawning biomass values
#' @param B0 Unfished spawning biomass
#' @param F_target Target fishing mortality (e.g., F40%)
#' @param b_trigger_frac Biomass fraction where F is F_target
#' @param b_min_frac Biomass fraction where F is zero
#' @return Vector of F values under the control rule
calc_sloping_hcr <- function(SSB, B0, F_target, b_trigger_frac = 0.4, b_min_frac = 0.05) {
  B_trigger <- b_trigger_frac * B0
  B_min <- b_min_frac * B0

  F_vals <- ifelse(
    SSB >= B_trigger,
    F_target,
    ifelse(
      SSB <= B_min,
      0,
      F_target * (SSB - B_min) / (B_trigger - B_min)
    )
  )

  return(F_vals)
}

#' Calculate catch for a given F and numbers-at-age
#' @param N Vector of numbers at age
#' @param F Fishing mortality
#' @param params Parameter list
#' @return Catch (biomass)
calc_catch_at_F <- function(N, F, params) {
  sel <- selectivity_at_age(params$ages, params)
  wt <- weight_at_age(params$ages, params)
  M <- params$M

  catch_total <- 0
  for (a in 1:params$n_ages) {
    Z <- M + F * sel[a]
    catch_total <- catch_total + N[a] * wt[a] * (F * sel[a] / Z) * (1 - exp(-Z))
  }

  return(catch_total)
}

#' Calculate F that satisfies a catch cap for a given N
#' @param N Vector of numbers at age
#' @param params Parameter list
#' @param catch_cap Catch cap (biomass)
#' @param F_max Maximum F to search
#' @param n_grid Number of grid points for bracketing
#' @return F value that matches the catch cap, or F_max if cap is not binding
calc_F_for_catch_cap <- function(N, params, catch_cap, F_max = 5, n_grid = 200) {
  if (catch_cap <= 0) {
    return(0)
  }

  catch_0 <- calc_catch_at_F(N, 0, params)
  if (catch_0 >= catch_cap) {
    return(0)
  }

  F_grid <- seq(0, F_max, length.out = n_grid)
  catch_grid <- sapply(F_grid, function(F) calc_catch_at_F(N, F, params))
  idx <- which(catch_grid >= catch_cap)

  if (length(idx) == 0) {
    return(F_max)
  }

  j <- idx[1]
  if (j == 1) {
    return(0)
  }

  f_root <- function(F) calc_catch_at_F(N, F, params) - catch_cap
  uniroot(f_root, lower = F_grid[j - 1], upper = F_grid[j])$root
}

# =============================================================================
# Population Simulation
# =============================================================================

#' Initialize population at equilibrium for a given F
#' @param params Parameter list
#' @param F_init Fishing mortality used to set the initial equilibrium
#' @return Vector of numbers at age (n_ages x 1)
init_population <- function(params, F_init) {
  N_per_R <- calc_N_equilibrium(F_init, params)
  R_eq <- calc_R_eq(F_init, params)
  N0 <- N_per_R * R_eq
  return(N0)
}

#' Run population simulation
#' @param params Parameter list
#' @param F_series Vector of F values for each year (or single value)
#' @param seed Random seed for reproducibility
#' @param init_F Fishing mortality used to set the initial equilibrium
#' @return List with simulation results
run_simulation <- function(params, F_series = NULL, seed = NULL, init_F = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_years <- params$n_years
  n_ages <- params$n_ages

  # Default F series if not provided
  if (is.null(F_series)) {
    F_series <- rep(params$F_init, n_years)
  } else if (length(F_series) == 1) {
    F_series <- rep(F_series, n_years)
  }

  # Pre-calculate biological parameters
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  sel <- selectivity_at_age(params$ages, params)
  SSB0 <- calc_SSB0(params)

  # Storage matrices
  N <- matrix(0, nrow = n_ages, ncol = n_years)
  SSB <- numeric(n_years)
  Recruitment <- numeric(n_years)
  Catch <- numeric(n_years)
  F_year <- F_series

  if (is.null(init_F)) {
    if (!is.null(params$Fmsy)) {
      init_F <- params$Fmsy
    } else {
      init_F <- calc_reference_points(params)$Fmsy
    }
  }

  # Initialize population at fished equilibrium (Fmsy by default)
  N[, 1] <- init_population(params, init_F)
  SSB[1] <- calc_ssb(N[, 1], params)
  Recruitment[1] <- N[1, 1]

  # Calculate catch in year 1
  F1 <- F_series[1]
  catch_year1 <- 0
  for (a in 1:n_ages) {
    Z <- params$M + F1 * sel[a]
    catch_year1 <- catch_year1 + N[a, 1] * wt[a] * (F1 * sel[a] / Z) * (1 - exp(-Z))
  }
  Catch[1] <- catch_year1

  # Run simulation
  for (t in 2:n_years) {
    F_t <- F_series[t]

    # Calculate recruitment from SSB in previous year
    R_det <- beverton_holt(SSB[t - 1], params, SSB0)
    R_t <- add_recruitment_error(R_det, params)
    Recruitment[t] <- R_t
    N[1, t] <- R_t

    # Survival for ages 2 to n_ages-1
    for (a in 2:(n_ages - 1)) {
      Z <- params$M + F_series[t - 1] * sel[a - 1]
      N[a, t] <- N[a - 1, t - 1] * exp(-Z)
    }

    # Plus group
    Z_last <- params$M + F_series[t - 1] * sel[n_ages - 1]
    Z_plus <- params$M + F_series[t - 1] * sel[n_ages]
    N[n_ages, t] <- N[n_ages - 1, t - 1] * exp(-Z_last) +
                    N[n_ages, t - 1] * exp(-Z_plus)

    # Calculate SSB (females only)
    SSB[t] <- calc_ssb(N[, t], params)

    # Calculate catch (Baranov equation)
    catch_t <- 0
    for (a in 1:n_ages) {
      Z <- params$M + F_t * sel[a]
      catch_t <- catch_t + N[a, t] * wt[a] * (F_t * sel[a] / Z) * (1 - exp(-Z))
    }
    Catch[t] <- catch_t
  }

  # Return results
  results <- list(
    N = N,
    SSB = SSB,
    Recruitment = Recruitment,
    Catch = Catch,
    F = F_series,
    years = 1:n_years,
    ages = params$ages,
    params = params,
    SSB0 = SSB0,
    weight = wt,
    maturity = mat,
    selectivity = sel
  )

  return(results)
}

#' Run population simulation with time-varying steepness and specified F series
#' @param params Parameter list
#' @param h_series Vector of steepness values by year (length >= n_years)
#' @param F_series Vector of F values for each year (or single value)
#' @param seed Random seed for reproducibility
#' @param init_F Fishing mortality used to set the initial equilibrium
#' @return List with simulation results
run_simulation_variable_h <- function(params,
                                      h_series,
                                      F_series = NULL,
                                      seed = NULL,
                                      init_F = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_years <- params$n_years
  n_ages <- params$n_ages

  if (length(h_series) < n_years) {
    stop("h_series must have length >= n_years")
  }

  # Default F series if not provided
  if (is.null(F_series)) {
    F_series <- rep(params$F_init, n_years)
  } else if (length(F_series) == 1) {
    F_series <- rep(F_series, n_years)
  } else if (length(F_series) < n_years) {
    stop("F_series must have length >= n_years")
  }

  beverton_holt_h <- function(SSB, params, SSB0, h) {
    R0 <- params$R0
    alpha <- (4 * h * R0) / (SSB0 * (1 - h))
    beta <- (5 * h - 1) / (SSB0 * (1 - h))
    (alpha * SSB) / (1 + beta * SSB)
  }

  # Pre-calculate biological parameters
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  sel <- selectivity_at_age(params$ages, params)
  SSB0 <- calc_SSB0(params)

  # Storage matrices
  N <- matrix(0, nrow = n_ages, ncol = n_years)
  SSB <- numeric(n_years)
  Recruitment <- numeric(n_years)
  Catch <- numeric(n_years)

  if (is.null(init_F)) {
    if (!is.null(params$Fmsy)) {
      init_F <- params$Fmsy
    } else {
      init_F <- calc_reference_points(params)$Fmsy
    }
  }

  # Initialize population at fished equilibrium (Fmsy by default)
  N[, 1] <- init_population(params, init_F)
  SSB[1] <- calc_ssb(N[, 1], params)
  Recruitment[1] <- N[1, 1]

  # Calculate catch in year 1
  F1 <- F_series[1]
  catch_year1 <- 0
  for (a in 1:n_ages) {
    Z <- params$M + F1 * sel[a]
    catch_year1 <- catch_year1 + N[a, 1] * wt[a] * (F1 * sel[a] / Z) * (1 - exp(-Z))
  }
  Catch[1] <- catch_year1

  # Run simulation
  for (t in 2:n_years) {
    F_t <- F_series[t]

    # Calculate recruitment from SSB in previous year with time-varying h
    h_t <- h_series[t]
    R_det <- beverton_holt_h(SSB[t - 1], params, SSB0, h_t)
    R_t <- add_recruitment_error(R_det, params)
    Recruitment[t] <- R_t
    N[1, t] <- R_t

    # Survival for ages 2 to n_ages-1
    for (a in 2:(n_ages - 1)) {
      Z <- params$M + F_series[t - 1] * sel[a - 1]
      N[a, t] <- N[a - 1, t - 1] * exp(-Z)
    }

    # Plus group
    Z_last <- params$M + F_series[t - 1] * sel[n_ages - 1]
    Z_plus <- params$M + F_series[t - 1] * sel[n_ages]
    N[n_ages, t] <- N[n_ages - 1, t - 1] * exp(-Z_last) +
                    N[n_ages, t - 1] * exp(-Z_plus)

    # Calculate SSB (females only)
    SSB[t] <- calc_ssb(N[, t], params)

    # Calculate catch (Baranov equation)
    catch_t <- 0
    for (a in 1:n_ages) {
      Z <- params$M + F_t * sel[a]
      catch_t <- catch_t + N[a, t] * wt[a] * (F_t * sel[a] / Z) * (1 - exp(-Z))
    }
    Catch[t] <- catch_t
  }

  results <- list(
    N = N,
    SSB = SSB,
    Recruitment = Recruitment,
    Catch = Catch,
    F = F_series,
    years = 1:n_years,
    ages = params$ages,
    params = params,
    SSB0 = SSB0,
    weight = wt,
    maturity = mat,
    selectivity = sel
  )

  return(results)
}

#' Run population simulation with a sloping HCR based on SSB
#' @param params Parameter list
#' @param F_target Target fishing mortality (e.g., F40%)
#' @param B0 Unfished spawning biomass
#' @param catch_cap Catch cap (biomass) applied to annual catch
#' @param b_trigger_frac Biomass fraction where F is F_target
#' @param b_min_frac Biomass fraction where F is zero
#' @param seed Random seed for reproducibility
#' @param init_F Fishing mortality used to set the initial equilibrium
#' @return List with simulation results
run_simulation_hcr <- function(params,
                               F_target,
                               B0,
                               catch_cap = NULL,
                               b_trigger_frac = 0.4,
                               b_min_frac = 0.05,
                               seed = NULL,
                               init_F = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_years <- params$n_years
  n_ages <- params$n_ages

  # Pre-calculate biological parameters
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  sel <- selectivity_at_age(params$ages, params)
  SSB0 <- calc_SSB0(params)

  if (is.null(init_F)) {
    if (!is.null(params$Fmsy)) {
      init_F <- params$Fmsy
    } else {
      init_F <- calc_reference_points(params)$Fmsy
    }
  }

  # Storage matrices
  N <- matrix(0, nrow = n_ages, ncol = n_years)
  SSB <- numeric(n_years)
  Recruitment <- numeric(n_years)
  Catch <- numeric(n_years)
  F_year <- numeric(n_years)

  # Initialize population at fished equilibrium (Fmsy by default)
  N[, 1] <- init_population(params, init_F)
  SSB[1] <- calc_ssb(N[, 1], params)
  Recruitment[1] <- N[1, 1]

  # Fishing mortality in year 1 from HCR (apply catch cap if provided)
  F_base <- calc_sloping_hcr(SSB[1], B0, F_target, b_trigger_frac, b_min_frac)
  if (!is.null(catch_cap)) {
    F_cap <- calc_F_for_catch_cap(N[, 1], params, catch_cap)
    F_year[1] <- min(F_base, F_cap)
  } else {
    F_year[1] <- F_base
  }

  # Calculate catch in year 1
  F1 <- F_year[1]
  catch_year1 <- 0
  for (a in 1:n_ages) {
    Z <- params$M + F1 * sel[a]
    catch_year1 <- catch_year1 + N[a, 1] * wt[a] * (F1 * sel[a] / Z) * (1 - exp(-Z))
  }
  Catch[1] <- catch_year1

  # Run simulation
  for (t in 2:n_years) {
    # Base F is set from previous-year SSB
    F_base <- calc_sloping_hcr(SSB[t - 1], B0, F_target, b_trigger_frac, b_min_frac)

    # Calculate recruitment from SSB in previous year
    R_det <- beverton_holt(SSB[t - 1], params, SSB0)
    R_t <- add_recruitment_error(R_det, params)
    Recruitment[t] <- R_t
    N[1, t] <- R_t

    # Survival for ages 2 to n_ages-1
    for (a in 2:(n_ages - 1)) {
      Z <- params$M + F_year[t - 1] * sel[a - 1]
      N[a, t] <- N[a - 1, t - 1] * exp(-Z)
    }

    # Plus group
    Z_last <- params$M + F_year[t - 1] * sel[n_ages - 1]
    Z_plus <- params$M + F_year[t - 1] * sel[n_ages]
    N[n_ages, t] <- N[n_ages - 1, t - 1] * exp(-Z_last) +
                    N[n_ages, t - 1] * exp(-Z_plus)

    # Calculate SSB (females only)
    SSB[t] <- calc_ssb(N[, t], params)

    if (!is.null(catch_cap)) {
      F_cap <- calc_F_for_catch_cap(N[, t], params, catch_cap)
      F_year[t] <- min(F_base, F_cap)
    } else {
      F_year[t] <- F_base
    }

    # Calculate catch (Baranov equation)
    catch_t <- 0
    for (a in 1:n_ages) {
      Z <- params$M + F_year[t] * sel[a]
      catch_t <- catch_t + N[a, t] * wt[a] * (F_year[t] * sel[a] / Z) * (1 - exp(-Z))
    }
    Catch[t] <- catch_t
  }

  results <- list(
    N = N,
    SSB = SSB,
    Recruitment = Recruitment,
    Catch = Catch,
    F = F_year,
    years = 1:n_years,
    ages = params$ages,
    params = params,
    SSB0 = SSB0,
    weight = wt,
    maturity = mat,
    selectivity = sel
  )

  return(results)
}

#' Run population simulation with a sloping HCR and time-varying steepness
#' @param params Parameter list
#' @param h_series Vector of steepness values by year (length >= n_years)
#' @param F_target Target fishing mortality (e.g., F40%)
#' @param B0 Unfished spawning biomass
#' @param catch_cap Catch cap (biomass) applied to annual catch
#' @param b_trigger_frac Biomass fraction where F is F_target
#' @param b_min_frac Biomass fraction where F is zero
#' @param seed Random seed for reproducibility
#' @param init_F Fishing mortality used to set the initial equilibrium
#' @return List with simulation results
run_simulation_hcr_variable_h <- function(params,
                                          h_series,
                                          F_target,
                                          B0,
                                          catch_cap = NULL,
                                          b_trigger_frac = 0.4,
                                          b_min_frac = 0.05,
                                          seed = NULL,
                                          init_F = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_years <- params$n_years
  n_ages <- params$n_ages

  if (length(h_series) < n_years) {
    stop("h_series must have length >= n_years")
  }

  beverton_holt_h <- function(SSB, params, SSB0, h) {
    R0 <- params$R0
    alpha <- (4 * h * R0) / (SSB0 * (1 - h))
    beta <- (5 * h - 1) / (SSB0 * (1 - h))
    (alpha * SSB) / (1 + beta * SSB)
  }

  # Pre-calculate biological parameters
  wt <- weight_at_age(params$ages, params)
  mat <- maturity_at_age(params$ages, params)
  sel <- selectivity_at_age(params$ages, params)
  SSB0 <- calc_SSB0(params)

  if (is.null(init_F)) {
    if (!is.null(params$Fmsy)) {
      init_F <- params$Fmsy
    } else {
      init_F <- calc_reference_points(params)$Fmsy
    }
  }

  # Storage matrices
  N <- matrix(0, nrow = n_ages, ncol = n_years)
  SSB <- numeric(n_years)
  Recruitment <- numeric(n_years)
  Catch <- numeric(n_years)
  F_year <- numeric(n_years)

  # Initialize population at fished equilibrium (Fmsy by default)
  N[, 1] <- init_population(params, init_F)
  SSB[1] <- calc_ssb(N[, 1], params)
  Recruitment[1] <- N[1, 1]

  # Fishing mortality in year 1 from HCR (apply catch cap if provided)
  F_base <- calc_sloping_hcr(SSB[1], B0, F_target, b_trigger_frac, b_min_frac)
  if (!is.null(catch_cap)) {
    F_cap <- calc_F_for_catch_cap(N[, 1], params, catch_cap)
    F_year[1] <- min(F_base, F_cap)
  } else {
    F_year[1] <- F_base
  }

  # Calculate catch in year 1
  F1 <- F_year[1]
  catch_year1 <- 0
  for (a in 1:n_ages) {
    Z <- params$M + F1 * sel[a]
    catch_year1 <- catch_year1 + N[a, 1] * wt[a] * (F1 * sel[a] / Z) * (1 - exp(-Z))
  }
  Catch[1] <- catch_year1

  # Run simulation
  for (t in 2:n_years) {
    # Base F is set from previous-year SSB
    F_base <- calc_sloping_hcr(SSB[t - 1], B0, F_target, b_trigger_frac, b_min_frac)

    # Calculate recruitment from SSB in previous year with time-varying h
    h_t <- h_series[t]
    R_det <- beverton_holt_h(SSB[t - 1], params, SSB0, h_t)
    R_t <- add_recruitment_error(R_det, params)
    Recruitment[t] <- R_t
    N[1, t] <- R_t

    # Survival for ages 2 to n_ages-1
    for (a in 2:(n_ages - 1)) {
      Z <- params$M + F_year[t - 1] * sel[a - 1]
      N[a, t] <- N[a - 1, t - 1] * exp(-Z)
    }

    # Plus group
    Z_last <- params$M + F_year[t - 1] * sel[n_ages - 1]
    Z_plus <- params$M + F_year[t - 1] * sel[n_ages]
    N[n_ages, t] <- N[n_ages - 1, t - 1] * exp(-Z_last) +
                    N[n_ages, t - 1] * exp(-Z_plus)

    # Calculate SSB (females only)
    SSB[t] <- calc_ssb(N[, t], params)

    if (!is.null(catch_cap)) {
      F_cap <- calc_F_for_catch_cap(N[, t], params, catch_cap)
      F_year[t] <- min(F_base, F_cap)
    } else {
      F_year[t] <- F_base
    }

    # Calculate catch (Baranov equation)
    catch_t <- 0
    for (a in 1:n_ages) {
      Z <- params$M + F_year[t] * sel[a]
      catch_t <- catch_t + N[a, t] * wt[a] * (F_year[t] * sel[a] / Z) * (1 - exp(-Z))
    }
    Catch[t] <- catch_t
  }

  results <- list(
    N = N,
    SSB = SSB,
    Recruitment = Recruitment,
    Catch = Catch,
    F = F_year,
    years = 1:n_years,
    ages = params$ages,
    params = params,
    SSB0 = SSB0,
    weight = wt,
    maturity = mat,
    selectivity = sel
  )

  return(results)
}

# =============================================================================
# Plotting Functions
# =============================================================================

#' Plot simulation results
#' @param sim Simulation results from run_simulation()
#' @param ref_points Reference points from calc_reference_points()
plot_simulation <- function(sim, ref_points = NULL) {

  df_ts <- data.frame(
    Year = sim$years,
    SSB = sim$SSB / 1e9,
    Recruitment = sim$Recruitment / 1e9,
    Catch = sim$Catch / 1e9
  )

  p_ssb <- ggplot2::ggplot(df_ts, ggplot2::aes(x = Year, y = SSB)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(
      x = "Year",
      y = "Female SSB (million tonnes)",
      title = "Spawning Stock Biomass (Females)"
    ) +
    ggthemes::theme_few()

  if (!is.null(ref_points)) {
    p_ssb <- p_ssb +
      ggplot2::geom_hline(
        yintercept = ref_points$Bmsy / 1e9,
        color = "red",
        linetype = "dashed",
        size = 0.8
      ) +
      ggplot2::geom_hline(
        yintercept = ref_points$SSB0 / 1e9,
        color = "blue",
        linetype = "dashed",
        size = 0.8
      )
  }

  p_rec <- ggplot2::ggplot(df_ts, ggplot2::aes(x = Year, y = Recruitment)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(
      x = "Year",
      y = "Recruitment (billions)",
      title = "Recruitment (Age 1)"
    ) +
    ggthemes::theme_few()

  p_catch <- ggplot2::ggplot(df_ts, ggplot2::aes(x = Year, y = Catch)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(
      x = "Year",
      y = "Catch (million tonnes)",
      title = "Catch"
    ) +
    ggthemes::theme_few()

  if (!is.null(ref_points)) {
    p_catch <- p_catch +
      ggplot2::geom_hline(
        yintercept = ref_points$MSY / 1e9,
        color = "red",
        linetype = "dashed",
        size = 0.8
      )
  }

  sr_df <- data.frame(
    SSB = sim$SSB[-length(sim$SSB)] / 1e9,
    Recruitment = sim$Recruitment[-1] / 1e9
  )
  ssb_seq <- seq(0, max(sim$SSB) * 1.1, length.out = 100)
  r_expected <- sapply(ssb_seq, function(s) beverton_holt(s, sim$params, sim$SSB0))
  sr_line <- data.frame(
    SSB = ssb_seq / 1e9,
    Recruitment = r_expected / 1e9
  )

  p_sr <- ggplot2::ggplot(sr_df, ggplot2::aes(x = SSB, y = Recruitment)) +
    ggplot2::geom_point(color = "gray40", size = 1.5) +
    ggplot2::geom_line(data = sr_line, ggplot2::aes(x = SSB, y = Recruitment),
                       color = "red", size = 1) +
    ggplot2::labs(
      x = "Female SSB (million tonnes)",
      y = "Recruitment (billions)",
      title = "Stock-Recruitment"
    ) +
    ggthemes::theme_few()

  plots <- list(
    ssb = p_ssb,
    recruitment = p_rec,
    catch = p_catch,
    stock_recruit = p_sr
  )

  lapply(plots, print)
  invisible(plots)
}

#' Plot biological parameters
#' @param params Parameter list
plot_biology <- function(params) {

  ages <- params$ages

  # Weight at age
  wt <- weight_at_age(ages, params)
  df_wt <- data.frame(Age = ages, Weight = wt)
  p_wt <- ggplot2::ggplot(df_wt, ggplot2::aes(x = Age, y = Weight)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(yintercept = c(0.3, 0.5), linetype = "dashed", color = "gray") +
    ggplot2::annotate("text", x = 2, y = 0.32, label = "Age 3 target", size = 3, color = "gray40") +
    ggplot2::annotate("text", x = 2, y = 0.52, label = "Age 5 target", size = 3, color = "gray40") +
    ggplot2::labs(x = "Age", y = "Weight (kg)", title = "Weight at Age") +
    ggthemes::theme_few()

  # Maturity at age
  mat <- maturity_at_age(ages, params)
  df_mat <- data.frame(Age = ages, Maturity = mat)
  p_mat <- ggplot2::ggplot(df_mat, ggplot2::aes(x = Age, y = Maturity)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = params$mat_a50, linetype = "dashed", color = "gray") +
    ggplot2::labs(x = "Age", y = "Proportion Mature", title = "Maturity at Age (Females)") +
    ggthemes::theme_few()

  # Selectivity at age
  sel <- selectivity_at_age(ages, params)
  df_sel <- data.frame(Age = ages, Selectivity = sel)
  p_sel <- ggplot2::ggplot(df_sel, ggplot2::aes(x = Age, y = Selectivity)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = params$sel_a50, linetype = "dashed", color = "gray") +
    ggplot2::labs(x = "Age", y = "Selectivity", title = "Fishery Selectivity at Age") +
    ggthemes::theme_few()

  # Yield vs F curve
  F_seq <- seq(0, 1.5, by = 0.01)
  yield_seq <- sapply(F_seq, function(f) calc_yield_eq(f, params))
  ref <- calc_reference_points(params)
  df_yield <- data.frame(F = F_seq, Yield = yield_seq / 1e9)
  p_yield <- ggplot2::ggplot(df_yield, ggplot2::aes(x = F, y = Yield)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_vline(xintercept = ref$Fmsy, color = "red", linetype = "dashed", size = 0.8) +
    ggplot2::geom_hline(yintercept = ref$MSY / 1e9, color = "red", linetype = "dashed", size = 0.8) +
    ggplot2::labs(
      x = "Fishing Mortality (F)",
      y = "Equilibrium Yield (million tonnes)",
      title = "Yield vs F"
    ) +
    ggthemes::theme_few()

  plots <- list(
    weight = p_wt,
    maturity = p_mat,
    selectivity = p_sel,
    yield = p_yield
  )

  lapply(plots, print)
  invisible(plots)
}

# =============================================================================
# Scenario Comparison Functions
# =============================================================================

#' Create a scenario with modified steepness
#' @param h Steepness value
#' @param scenario_name Name for the scenario
#' @return List with params, reference points, and scenario name
create_scenario <- function(h, scenario_name = NULL) {
  params <- set_parameters()
  params$h <- h

  if (is.null(scenario_name)) {
    scenario_name <- paste0("h = ", h)
  }

  ref_points <- calc_reference_points(params)
  params$Fmsy <- ref_points$Fmsy

  list(
    name = scenario_name,
    params = params,
    ref_points = ref_points
  )
}

#' Compare reference points across scenarios
#' @param scenarios List of scenarios from create_scenario()
#' @return Data frame comparing reference points
compare_reference_points <- function(scenarios) {
  df <- data.frame(
    Scenario = sapply(scenarios, function(s) s$name),
    Steepness = sapply(scenarios, function(s) s$params$h),
    Fmsy = sapply(scenarios, function(s) round(s$ref_points$Fmsy, 4)),
    Bmsy_mt = sapply(scenarios, function(s) round(s$ref_points$Bmsy / 1e9, 2)),
    MSY_mt = sapply(scenarios, function(s) round(s$ref_points$MSY / 1e9, 2)),
    SSB0_mt = sapply(scenarios, function(s) round(s$ref_points$SSB0 / 1e9, 2)),
    Bmsy_B0 = sapply(scenarios, function(s) round(s$ref_points$Bmsy_B0_ratio, 3)),
    R_msy_B = sapply(scenarios, function(s) round(s$ref_points$R_msy / 1e9, 2))
  )
  colnames(df) <- c("Scenario", "Steepness", "Fmsy", "Bmsy (Mt)",
                    "MSY (Mt)", "SSB0 (Mt)", "Bmsy/SSB0", "R_msy (B)")
  return(df)
}

#' Plot yield curves for multiple scenarios
#' @param scenarios List of scenarios from create_scenario()
#' @param colors Vector of colors for each scenario
plot_yield_comparison <- function(scenarios, colors = NULL) {
  if (is.null(colors)) {
    colors <- c("blue", "red", "darkgreen", "purple")[1:length(scenarios)]
  }

  F_seq <- seq(0, 1.5, by = 0.01)

  # Calculate yield curves
  scenario_names <- sapply(scenarios, function(s) s$name)
  yield_df <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
    data.frame(
      F = F_seq,
      Yield = sapply(F_seq, function(f) calc_yield_eq(f, scenarios[[i]]$params)) / 1e9,
      Scenario = scenario_names[i]
    )
  }))

  fmsy_df <- data.frame(
    Scenario = scenario_names,
    Fmsy = sapply(scenarios, function(s) s$ref_points$Fmsy)
  )

  p <- ggplot2::ggplot(yield_df, ggplot2::aes(x = F, y = Yield, color = Scenario)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_vline(
      data = fmsy_df,
      ggplot2::aes(xintercept = Fmsy, color = Scenario),
      linetype = "dashed",
      size = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = stats::setNames(colors, scenario_names)) +
    ggplot2::labs(
      x = "Fishing Mortality (F)",
      y = "Equilibrium Yield (million tonnes)",
      title = "Yield vs F: Scenario Comparison"
    ) +
    ggthemes::theme_few()

  return(p)
}

#' Plot stock-recruitment curves for multiple scenarios
#' @param scenarios List of scenarios from create_scenario()
#' @param colors Vector of colors for each scenario
plot_sr_comparison <- function(scenarios, colors = NULL) {
  if (is.null(colors)) {
    colors <- c("blue", "red", "darkgreen", "purple")[1:length(scenarios)]
  }

  # Find SSB range (use max SSB0 across scenarios)
  ssb0_max <- max(sapply(scenarios, function(s) s$ref_points$SSB0))
  ssb_seq <- seq(0, ssb0_max * 1.1, length.out = 200)

  # Calculate S-R curves
  scenario_names <- sapply(scenarios, function(s) s$name)
  sr_df <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
    data.frame(
      SSB = ssb_seq / 1e9,
      Recruitment = sapply(ssb_seq, function(ssb) {
        beverton_holt(ssb, scenarios[[i]]$params, scenarios[[i]]$ref_points$SSB0)
      }) / 1e9,
      Scenario = scenario_names[i]
    )
  }))

  bmsy_df <- data.frame(
    Scenario = scenario_names,
    Bmsy = sapply(scenarios, function(s) s$ref_points$Bmsy) / 1e9
  )

  repl_ssb <- c(0, ssb0_max * 1.1) / 1e9
  repl_slope <- scenarios[[1]]$params$R0 / scenarios[[1]]$ref_points$SSB0
  repl_df <- data.frame(
    SSB = repl_ssb,
    Recruitment = repl_slope * repl_ssb
  )

  p <- ggplot2::ggplot(sr_df, ggplot2::aes(x = SSB, y = Recruitment, color = Scenario)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_vline(
      data = bmsy_df,
      ggplot2::aes(xintercept = Bmsy, color = Scenario),
      linetype = "dashed",
      size = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = repl_df,
      ggplot2::aes(x = SSB, y = Recruitment, linetype = "Replacement (R0/SSB0)"),
      color = "gray50",
      size = 0.8
    ) +
    ggplot2::scale_color_manual(values = stats::setNames(colors, scenario_names)) +
    ggplot2::scale_linetype_manual(values = c("Replacement (R0/SSB0)" = "dashed")) +
    ggplot2::labs(
      x = "Female SSB (million tonnes)",
      y = "Recruitment (billions)",
      title = "Stock-Recruitment: Scenario Comparison",
      linetype = ""
    ) +
    ggthemes::theme_few() +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2)
    )

  return(p)
}

#' Run and compare simulations across scenarios
#' @param scenarios List of scenarios
#' @param n_sims Number of Monte Carlo simulations
#' @param seed Base random seed
#' @return List with simulation results for each scenario
run_scenario_simulations <- function(scenarios, n_sims = 100, seed = 42) {
  results <- lapply(scenarios, function(s) {
    params <- s$params
    ref_points <- s$ref_points

    ssb_matrix <- matrix(NA, nrow = params$n_years, ncol = n_sims)
    catch_matrix <- matrix(NA, nrow = params$n_years, ncol = n_sims)
    rec_matrix <- matrix(NA, nrow = params$n_years, ncol = n_sims)

    for (i in 1:n_sims) {
      sim <- run_simulation(params, F_series = ref_points$Fmsy, seed = seed + i)
      ssb_matrix[, i] <- sim$SSB
      catch_matrix[, i] <- sim$Catch
      rec_matrix[, i] <- sim$Recruitment
    }

    list(
      name = s$name,
      params = params,
      ref_points = ref_points,
      ssb = ssb_matrix,
      catch = catch_matrix,
      recruitment = rec_matrix
    )
  })

  names(results) <- sapply(scenarios, function(s) s$name)
  return(results)
}

# =============================================================================
# Intermediate F Analysis Functions
# =============================================================================

#' Calculate intermediate F between two Fmsy values
#' @param scenarios List of scenarios
#' @param method Method for calculating intermediate F: "mean", "geometric", or "median"
#' @return Intermediate F value
calc_intermediate_F <- function(scenarios, method = "mean") {
  fmsy_values <- sapply(scenarios, function(s) s$ref_points$Fmsy)

  if (method == "mean") {
    return(mean(fmsy_values))
  } else if (method == "geometric") {
    return(exp(mean(log(fmsy_values))))
  } else if (method == "median") {
    return(median(fmsy_values))
  } else {
    stop("Unknown method. Use 'mean', 'geometric', or 'median'")
  }
}

#' Run simulations at intermediate F for multiple scenarios
#' @param scenarios List of scenarios
#' @param F_intermediate Intermediate F value to apply
#' @param n_sims Number of Monte Carlo simulations
#' @param seed Base random seed
#' @return List with results for each scenario at intermediate F
run_intermediate_F_simulations <- function(scenarios, F_intermediate, n_sims = 100, seed = 42) {
  results <- lapply(scenarios, function(s) {
    params <- s$params
    ref_points <- s$ref_points
    n_years <- params$n_years

    ssb_matrix <- matrix(NA, nrow = n_years, ncol = n_sims)
    catch_matrix <- matrix(NA, nrow = n_years, ncol = n_sims)
    rec_matrix <- matrix(NA, nrow = n_years, ncol = n_sims)

    for (i in 1:n_sims) {
      sim <- run_simulation(params, F_series = F_intermediate, seed = seed + i)
      ssb_matrix[, i] <- sim$SSB
      catch_matrix[, i] <- sim$Catch
      rec_matrix[, i] <- sim$Recruitment
    }

    list(
      name = s$name,
      params = params,
      ref_points = ref_points,
      F_applied = F_intermediate,
      ssb = ssb_matrix,
      catch = catch_matrix,
      recruitment = rec_matrix
    )
  })

  names(results) <- sapply(scenarios, function(s) s$name)
  return(results)
}

#' Compare performance metrics at intermediate F vs Fmsy
#' @param scenarios List of scenarios
#' @param F_intermediate Intermediate F value
#' @param n_sims Number of simulations
#' @param seed Random seed
#' @return Data frame with comparison metrics
compare_F_strategies <- function(scenarios, F_intermediate, n_sims = 100, seed = 42) {

  results <- lapply(scenarios, function(s) {
    params <- s$params
    ref_points <- s$ref_points
    n_years <- params$n_years

    # Run at Fmsy
    ssb_fmsy <- catch_fmsy <- numeric(n_sims)
    for (i in 1:n_sims) {
      sim <- run_simulation(params, F_series = ref_points$Fmsy, seed = seed + i)
      ssb_fmsy[i] <- mean(sim$SSB)
      catch_fmsy[i] <- mean(sim$Catch)
    }

    # Run at intermediate F
    ssb_int <- catch_int <- numeric(n_sims)
    for (i in 1:n_sims) {
      sim <- run_simulation(params, F_series = F_intermediate, seed = seed + i)
      ssb_int[i] <- mean(sim$SSB)
      catch_int[i] <- mean(sim$Catch)
    }

    data.frame(
      Scenario = s$name,
      Steepness = params$h,
      Fmsy = round(ref_points$Fmsy, 3),
      F_intermediate = round(F_intermediate, 3),
      F_ratio = round(F_intermediate / ref_points$Fmsy, 2),
      Mean_SSB_Fmsy = round(mean(ssb_fmsy) / 1e9, 2),
      Mean_SSB_Fint = round(mean(ssb_int) / 1e9, 2),
      SSB_change_pct = round((mean(ssb_int) / mean(ssb_fmsy) - 1) * 100, 1),
      Mean_Catch_Fmsy = round(mean(catch_fmsy) / 1e9, 2),
      Mean_Catch_Fint = round(mean(catch_int) / 1e9, 2),
      Catch_change_pct = round((mean(catch_int) / mean(catch_fmsy) - 1) * 100, 1)
    )
  })

  do.call(rbind, results)
}

#' Plot simulation trajectories comparing Fmsy vs intermediate F
#' @param scenario Single scenario
#' @param F_intermediate Intermediate F value
#' @param seed Random seed
plot_F_comparison <- function(scenario, F_intermediate, seed = 42) {
  params <- scenario$params
  ref_points <- scenario$ref_points

  sim_fmsy <- run_simulation(params, F_series = ref_points$Fmsy, seed = seed)
  sim_int <- run_simulation(params, F_series = F_intermediate, seed = seed)

  years <- 1:params$n_years

  label_fmsy <- paste0("F=Fmsy (", round(ref_points$Fmsy, 2), ")")
  label_int <- paste0("F=", round(F_intermediate, 2))

  make_df <- function(sim, label) {
    data.frame(
      Year = years,
      Strategy = label,
      SSB = sim$SSB / 1e9,
      Catch = sim$Catch / 1e9,
      Recruitment = sim$Recruitment / 1e9,
      Depletion = sim$SSB / ref_points$Bmsy
    )
  }

  df_fmsy <- make_df(sim_fmsy, label_fmsy)
  df_int <- make_df(sim_int, label_int)
  df_all <- rbind(df_fmsy, df_int)

  p_ssb <- ggplot2::ggplot(df_all, ggplot2::aes(x = Year, y = SSB, color = Strategy)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = ref_points$Bmsy / 1e9, color = "red",
                        linetype = "dashed", size = 0.8) +
    ggplot2::labs(
      x = "Year",
      y = "Female SSB (Mt)",
      title = paste(scenario$name, "- SSB")
    ) +
    ggthemes::theme_few()

  p_catch <- ggplot2::ggplot(df_all, ggplot2::aes(x = Year, y = Catch, color = Strategy)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = ref_points$MSY / 1e9, color = "red",
                        linetype = "dashed", size = 0.8) +
    ggplot2::labs(
      x = "Year",
      y = "Catch (Mt)",
      title = paste(scenario$name, "- Catch")
    ) +
    ggthemes::theme_few()

  p_rec <- ggplot2::ggplot(df_all, ggplot2::aes(x = Year, y = Recruitment, color = Strategy)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(
      x = "Year",
      y = "Recruitment (B)",
      title = paste(scenario$name, "- Recruitment")
    ) +
    ggthemes::theme_few()

  p_dep <- ggplot2::ggplot(df_all, ggplot2::aes(x = Year, y = Depletion, color = Strategy)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 0.8) +
    ggplot2::labs(
      x = "Year",
      y = "SSB / Bmsy",
      title = paste(scenario$name, "- Depletion")
    ) +
    ggthemes::theme_few()

  plots <- list(
    ssb = p_ssb,
    catch = p_catch,
    recruitment = p_rec,
    depletion = p_dep
  )

  lapply(plots, print)
  invisible(plots)
}

# =============================================================================
# Fishing Mortality Sensitivity
# =============================================================================

#' Run simulations across a set of F multipliers
#' @param scenario Single scenario
#' @param F_multipliers Vector of multipliers applied to Fmsy
#' @param n_sims Number of simulations per F multiplier
#' @param seed Base random seed
#' @return Data frame of summary metrics
run_F_sensitivity <- function(scenario, F_multipliers, n_sims = 100, seed = 100) {
  params <- scenario$params
  ref_points <- scenario$ref_points
  n_years <- params$n_years

  results <- lapply(seq_along(F_multipliers), function(i) {
    F_mult <- F_multipliers[i]
    F_val <- F_mult * ref_points$Fmsy

    ssb_matrix <- matrix(NA, nrow = n_years, ncol = n_sims)
    catch_matrix <- matrix(NA, nrow = n_years, ncol = n_sims)

    for (j in 1:n_sims) {
      sim <- run_simulation(params, F_series = F_val, seed = seed + i * 1000 + j)
      ssb_matrix[, j] <- sim$SSB
      catch_matrix[, j] <- sim$Catch
    }

    data.frame(
      Scenario = scenario$name,
      Steepness = params$h,
      F_multiplier = F_mult,
      F_applied = F_val,
      Mean_SSB = mean(ssb_matrix),
      Terminal_SSB = mean(ssb_matrix[n_years, ]),
      P_SSB_lt_Bmsy = mean(ssb_matrix[n_years, ] < ref_points$Bmsy),
      Mean_Catch = mean(catch_matrix)
    )
  })

  do.call(rbind, results)
}

# =============================================================================
# Main Execution
# =============================================================================

if (interactive() || !exists("SOURCED_ONLY")) {

  # Create two scenarios with different steepness
  cat("=== Creating Scenarios ===\n")
  scenario_low <- create_scenario(h = 0.6, "Low steepness (h=0.6)")
  scenario_high <- create_scenario(h = 0.9, "High steepness (h=0.9)")
  scenarios <- list(scenario_low, scenario_high)

  # Compare reference points
  cat("\n=== Reference Point Comparison ===\n")
  comparison <- compare_reference_points(scenarios)
  print(comparison)

  # Calculate intermediate F
  F_int <- calc_intermediate_F(scenarios, method = "mean")
  cat(sprintf("\n=== Intermediate F (mean of Fmsy values) ===\n"))
  cat(sprintf("Fmsy (h=0.6): %.3f\n", scenario_low$ref_points$Fmsy))
  cat(sprintf("Fmsy (h=0.9): %.3f\n", scenario_high$ref_points$Fmsy))
  cat(sprintf("Intermediate F: %.3f\n", F_int))

  # Run simulations at respective Fmsy values
  cat("\n=== Simulations at Fmsy ===\n")
  for (s in scenarios) {
    sim <- run_simulation(s$params, F_series = s$ref_points$Fmsy, seed = 42)
    cat(sprintf("\n%s at F=Fmsy (%.3f):\n", s$name, s$ref_points$Fmsy))
    cat(sprintf("  Mean SSB:   %.2f million tonnes\n", mean(sim$SSB) / 1e9))
    cat(sprintf("  Mean Catch: %.2f million tonnes\n", mean(sim$Catch) / 1e9))
  }

  # Run simulations at intermediate F
  cat("\n=== Simulations at Intermediate F ===\n")
  for (s in scenarios) {
    sim <- run_simulation(s$params, F_series = F_int, seed = 42)
    cat(sprintf("\n%s at F=%.3f (intermediate):\n", s$name, F_int))
    cat(sprintf("  Mean SSB:   %.2f million tonnes\n", mean(sim$SSB) / 1e9))
    cat(sprintf("  Mean Catch: %.2f million tonnes\n", mean(sim$Catch) / 1e9))
    cat(sprintf("  F/Fmsy:     %.2f\n", F_int / s$ref_points$Fmsy))
  }

  # Weight at age check
  cat("\n=== Weight at Age Check ===\n")
  wt <- weight_at_age(scenario_low$params$ages, scenario_low$params)
  cat(sprintf("Age 3 weight: %.3f kg (target: ~0.3 kg)\n", wt[3]))
  cat(sprintf("Age 5 weight: %.3f kg (target: ~0.5 kg)\n", wt[5]))
}
