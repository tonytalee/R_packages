#' @title Count control limits for Xbar-R
#' @param value numeric vector, data to count control limits
#' @param sample_id character vector with same length of value, id of sample
#' @return data.frame contains control limits
#' @export
cl_Xbar_R <- function(value, sample_id) {

    # d2 & d3
    d2 <- c(1.128,1.693,2.059,2.326,2.534,2.704,2.847,2.97,3.078,3.173,3.258,3.336,3.407,
            3.472,3.532,3.588,3.64,3.689,3.735,3.778,3.819,3.858,3.895,3.931)
    d3 <- c(0.8525,0.8884,0.8798,0.8641,0.848,0.8332,0.8198,0.8078,0.7971,0.7873,0.7785,
            0.7704,0.763,0.7562,0.7499,0.7441,0.7386,0.7335,0.7287,0.7242,0.7199,0.7159,
            0.7121,0.7084)

    # Form data frame
    df <- data.frame(value= value, sample_id= sample_id)
    df$sample_id <- as.character(df$sample_id)

    # Reform data frame
    df <- df %>%
        group_by(sample_id) %>%
        summarise(Xbar = mean(value, na.rm = TRUE),
                  R = max(value, na.rm = TRUE) - min(value, na.rm = TRUE), # range
                  size = n()) %>%
        as.data.frame(.)
    # Sample size of subgroup
    if (length(unique(df$size)) == 1) {
        size <-  unique(df$size)
    } else {
        size <- df$size
    }

    # SPC constant
    d2 <- d2[size -1]
    d3 <- d3[size -1]
    sd <- sum(df$R / d2) / length(df$R)
    D3 <- 1 - 3 * d3 / d2
    D4 <- 1 + 3 * d3 / d2
    D3 <- ifelse(D3 < 0, 0, D3)

    # Average of X and R
    Xdbar <- mean(df$Xbar, na.rm = TRUE)
    Rbar <- mean(df$R, na.rm = TRUE)


    # Set control limits
    UCLx <- Xdbar + 3 * sd / sqrt(size) # control limit for Xbar or X chart
    LCLx <- Xdbar - 3 * sd / sqrt(size)
    UCLr <- sum(D4 * df$R) / length(df$R) # control limit for R or s chart
    LCLr <- sum(D3 * df$R) / length(df$R)

    #
    data.frame(UCLx = UCLx, LCLx = LCLx, UCLr = UCLr, LCLr = LCLr, CLx = Xdbar,
               CLr = Rbar, sd= sd)
}


#' @title Count control limits of X-mR chart
#' @param value numeric vector, data to count control limits
#' @return data.frame contains control limits
#' @export
cl_X_mR <- function(value) {

    # Reform data frame
    df <- data.frame(X= value) %>%
        mutate(mR = abs(value - lag(value)))

    # SPC constant
    A2 <- 2.66
    D3 <- 0
    D4 <- 3.27

    # Average of X and R
    Xbar <- mean(df$X, na.rm = TRUE)
    Rbar <- mean(df$mR, na.rm = TRUE)
    sd <- Rbar / 1.128

    # Set control limits
    UCLx <- Xbar + A2 * Rbar # control limit for Xbar or X chart
    LCLx <- Xbar - A2 * Rbar
    UCLmr <- D4 * Rbar # control limit for R or s chart
    LCLmr <- D3 * Rbar

    #
    data.frame(UCLx = UCLx, LCLx = LCLx, UCLmr = UCLmr, LCLmr = LCLmr,
         CLx = Xbar, CLmr = Rbar, sd = sd)
}

#' @title Count control limits for Xbar-mR-R chart
#' @param value numeric vector, data to count control limits
#' @param sample_id character vector with same length of value, id of sample
#' @return data.frame contains control limits
#' @export
cl_Xbar_s <- function(value, sample_id) {
    # Form data frame
    df <- data.frame(value= value, sample_id= sample_id)
    df$sample_id <- as.character(df$sample_id)

    # Reform data frame
    df <- df %>%
        group_by(sample_id) %>%
        summarise(Xbar = mean(value, na.rm = TRUE),
                  s = sd(value, na.rm = TRUE), # standard deviation
                  size = n()) %>%
        as.data.frame(.)

    # Sample size of subgroup
    if (length(unique(df$size)) == 1) {
        size <-  unique(df$size)
    } else {
        size <- df$size
    }

    # SPC constant
    c4<- sqrt(2 / (size-1)) * gamma(size/2) / gamma((size-1) / 2)
    sd <- sum(df$s / c4, na.rm = TRUE) / length(df$s)
    A3 <- round(3 / (c4 * sqrt(size)), 3)
    B3 <- round(1 - 3 * sqrt(1 - c4 **2) / c4, 3)
    B4 <- round(1 + 3 * sqrt(1 - c4 **2) / c4, 3)
    B3 <- ifelse(B3 < 0, 0, B3)

    # Average of X and R
    Xdbar <- mean(df$Xbar, na.rm = TRUE)
    Sbar <- mean(df$s, na.rm = TRUE)

    # Set control limits
    UCLx <- Xdbar + A3* Sbar # control limit for Xbar or X chart
    LCLx <- Xdbar - A3 * Sbar
    UCLs <- B4 * Sbar # control limit for R or s chart
    LCLs <- B3 * Sbar

    #
    data.frame(UCLx = UCLx, LCLx = LCLx, UCLs = UCLs, LCLs = LCLs, CLx = Xdbar,
               CLs = Sbar, sd)
}

#' @title Count control limits for p chart
#' @param value numeric vector, data to count control limits
#' @param size integer vector with same length of value, sample size
#' @param pbar if assigned, pbar will not be counted from raw data
#' @param ruleVarySize rule of varying sample size, by percentage
#' @return data.frame contains control limits of p charts
#' @export
cl_p <- function(value, size, pbar= NULL, ruleVarySize= 0.1) {
    # Form dataframe
    df <- data.frame(value= value, size= size) %>%
        mutate(frac = value / size)

    # pbar
    if (is.null(pbar)) pbar <- sum(df$value) / sum(df$size)

    # large sample size
    nbar <- ceiling(mean(df$size, na.rm = TRUE))
    large_sample <- ifelse((nbar * pbar >= 5 & nbar * (1 - pbar)), TRUE, FALSE)

    # varying sample size
    varying_size <- any(abs(df$size  - mean(df$size, na.rm = TRUE)) /
                            mean(df$size, na.rm = TRUE) > ruleVarySize)

    # Count control limits by the conditions large sample and varing size
    if (varying_size) {
        size <- df$size
    } else {
        size <- nbar
    }

    if (large_sample) {
        UCL <- round(pbar + 3 * sqrt(pbar * (1 - pbar) / size), 3)
        LCL <- round(pbar - 3 * sqrt(pbar * (1 - pbar) / size), 3)
        LCL <- ifelse(LCL <0, 0, LCL)
    } else {
        UCL <- qbinom(0.999, size, pbar) / size
        LCL <- qbinom(0.001, size, pbar) / size
    }

    #
    data.frame(UCL = UCL, LCL = LCL, CL = pbar)
}

#' @title Count control limits for np chart
#' @param value numeric vector, data to count control limits
#' @param size integer vector with same length of value, sample size
#' @param npbar if assigned, npbar will not be counted from raw data
#' @param ruleVarySize rule of varying sample size, by percentage
#' @return data.frame contains control limits of p charts
#' @export
cl_np <- function(value, size, npbar= NULL, ruleVarySize= 0.1) {
    # varying sample size
    varying_size <- any(abs(size  - mean(size, na.rm = TRUE)) /
                            mean(size, na.rm = TRUE) > ruleVarySize)
    if (varying_size) {
        return(paste0("Sample size variation is more than ", ruleVarySize,
                      ", it cannot be regarded as fixed sample size"))
    }

    size= ceiling(mean(size, na.rm = TRUE))
    if (is.null(npbar)) {
        pbar <- NULL
    } else {
        pbar <- npbar / size
    }
    dfa <- cl_p(value, size, pbar, ruleVarySize)

    data.frame(UCL= dfa$UCL * size, LCL= dfa$LCL * size, CL= dfa$CL * size)
}

#' @title Count control limits for u chart
#' @param value numeric vector, data to count control limits
#' @param size integer vector with same length of value, sample size
#' @param iu unit of per inspection
#' @param ubar if assigned, ubar will not be counted from raw data
#' @param ruleVarySize rule of varying sample size, by percentage
#' @return data.frame contains control limits
#' @export
cl_u <- function(value, size, iu, ubar= NULL, ruleVarySize= 0.1) {
    # Form dataframe
    df <- data.frame(value= value, size= size) %>%
        mutate(no_iu = size / iu, u = value / no_iu)

    # ubar
    ubar <- sum(df$value) / sum(df$no_iu)

    # large sample size & varying sample size
    large_sample <- ifelse(ubar >= 2, TRUE, FALSE)
    varying_size <- any(abs(df$size  - mean(df$size, na.rm = TRUE)) /
                            mean(df$size, na.rm = TRUE) > ruleVarySize)

    # Count control limits by the conditions large sample and varing size
    if (varying_size) {
        no_iu <- df$size
    } else {
        no_iu <- round(mean(df$no_iu, na.rm = TRUE))
    }

    if (large_sample) {
        UCL <- round(ubar + 3 * sqrt(ubar / no_iu), 3)
        LCL <- round(ubar - 3 * sqrt(ubar / no_iu), 3)
        if (LCL[1] < 0) { LCL <- 0}
    } else { # small lambda
        UCL <- qpois(0.997, ubar) / no_iu
        LCL <- 0
    }

    #
    data.frame(UCL = UCL, LCL = LCL, CL = ubar)
}

#' @title Count control limits for c chart
#' @param value numeric vector, data to count control limits
#' @param size integer vector with same length of value, sample size
#' @param cbar if assigned, cbar will not be counted from raw data
#' @param ruleVarySize rule of varying sample size, by percentage
#' @return data.frame contains control limits
#' @export
cl_c <- function(value, size, cbar= NULL, ruleVarySize= 0.1) {
    varying_size <- any(abs(size  - mean(size, na.rm = TRUE)) /
                            mean(size, na.rm = TRUE) > ruleVarySize)

    if (varying_size) {
        return(paste0("Sample size variation is more than ", ruleVarySize,
                      ", it cannot be regarded as fixed sample size"))
    }

    size= ceiling(mean(size, na.rm = TRUE))

    # Return
    cl_u(value, size, iu= size, ubar= cbar, ruleVarySize)
}
