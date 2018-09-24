# v0.1.1

# * Function to calculate process capability index -----
#' @title Calculate index for process capability, i.e. Cp, Cpk, Pp, Ppk...
#' @param x numeric vector for calculating process capability index
#' @param subgroup vector with the same length of x to label the subgroup
#' @param USL upper specification limit
#' @param targt spec target
#' @param LSL lower specification limit
#' @export
fun_Cp <- function(x,  subgroup= 1, USL = NULL, target = NULL, LSL = NULL) {

    df <- data.frame(subgroup = subgroup, value = x)
    df <- df[! is.na(df$value), ]

    # d2
    d2 <- c(1.128,1.693,2.059,2.326,2.534,2.704,2.847,2.97,3.078)

    #=== count mean and standard variation of overall data
    sample_size <- nrow(df)
    sample_mean <- round(mean(df$value), 3)
    StDev_overall <- round(sd(df$value), 3)

    # within subgroup analysis
    no_subgroup <- unique(subgroup)
    if (length(no_subgroup) >= 2) {
        df1 <- df %>%
            group_by(subgroup) %>%
            summarise(count = n(),
                      range = max(value) - min(value),
                      sd = sd(value))
        sd_esti <- c()
        for (i in 1:nrow(df1)) {
            n = df1$count[i]
            if (n  >= 2) {
                if (n <= 9) {
                    sd_p <- df1$range[i] / d2[n-1]
                } else {
                    c4<- sqrt(2 / (n-1)) * gamma(n/2) / gamma((n-1) / 2)
                    sd_p <- df1$sd[i] / c4
                }
                sd_esti <- c(sd_esti, sd_p)
            }
        }
        StDev_within <- round(mean(sd_esti), 3)
    } else {
        StDev_within <- NA
    }


    #=== Count Cp, Cpk, Pp, Ppk
    if (is.null(USL) & is.null(LSL)) {
        # Both USL and LSL are null
        Ca <- NA
        Cp <- NA
        Cpk <- NA
        Pp <- NA
        Ppk <- NA
        ppm_within <- NA
        ppm_overall <- NA
        ppm_obs <- NA
    } else {
        if (! is.null(USL) & ! is.null(LSL)) {
            # Both USL and LSL are not null
            if (is.null(target)) target <- (USL + LSL) / 2
            Ca <- round(2 * (sample_mean - target) / (USL - LSL) , 2)
            Cp <- round((USL - LSL) / (6 * StDev_within), 2)
            Pp <- round((USL - LSL) / (6 * StDev_overall), 2)
            CpkU <- round((USL - sample_mean) / (3 * StDev_within), 2)
            PpkU <- round((USL - sample_mean) / (3 * StDev_overall), 2)
            CpkL <- round((sample_mean - LSL) / (3 * StDev_within), 2)
            PpkL <- round((sample_mean - LSL) / (3 * StDev_overall), 2)
            Cpk <- round(min(CpkU, CpkL), 2)
            Ppk <- round(min(PpkU, PpkL), 2)
            ppm_within <- round((pnorm(CpkL * -3) + (1 - pnorm(CpkU * 3))) * 10^6, 2)
            ppm_overall <- round((pnorm(PpkL * -3) + (1 - pnorm(PpkU * 3))) * 10^6, 2)
            ppm_obs <- ((sum(df$value > USL) + sum(df$value < LSL)) / length(df$value)) *
                10^6
        } else {
            if (is.null((LSL))) {
                # LSL is null
                Ca <- NA
                Cp <- NA
                Pp <- NA
                CpkU <- round((USL - sample_mean) / (3 * StDev_within), 2)
                PpkU <- round((USL - sample_mean) / (3 * StDev_overall), 2)
                Cpk <- CpkU
                Ppk <- PpkU
                ppm_within <- round((1 - pnorm(CpkU * 3)) * 10^6, 2)
                ppm_overall <- round((1 - pnorm(PpkU * 3)) * 10^6, 2)
                ppm_obs <- (sum(df$value > USL) / length(df$value)) * 10^6
            } else {
                # USL is null
                Ca <- NA
                Cp <- NA
                Pp <- NA
                CpkL <- round((sample_mean - LSL) / (3 * StDev_within), 2)
                PpkL <- round((sample_mean - LSL) / (3 * StDev_overall), 2)
                Cpk <- CpkL
                Ppk <- PpkL
                ppm_within <- round(pnorm(CpkL * -3) * 10^6, 2)
                ppm_overall <- round(pnorm(PpkL * -3) * 10^6, 2)
                ppm_obs <- (sum(df$value < LSL) / length(df$value)) * 10^6
            }
        }
    }
    data.frame(sample_size = sample_size, sample_mean = sample_mean,
      StDev_overall = StDev_overall,
      StDev_within = StDev_within, Ca = Ca,
      Cp = Cp, Cpk = Cpk, CpkU = CpkU, CpkL = CpkL,
      Pp = Pp, Ppk = Ppk, PpkU = PpkU, PpkL = PpkL,
      ppm_within = ppm_within, ppm_overall = ppm_overall, ppm_obs = ppm_obs)
}


# * Functions for plot -----
#' @title  return a dataframe with x, y coordinates for plotting QQ line
QQ_points <- function(x, distribution = "norm", dparams = list(), qprobs = c(.25, .75),
                      detrend, identity, qtype) {
    # distributional function
    qFunc <- eval(parse(text = paste0("q", distribution)))

    oidx <- order(x)
    smp <- x[oidx]
    n <- length(smp)
    quantiles <- ppoints(n)

    # automatically estimate parameters with MLE, only if no parameters are
    # provided with dparams and there are at least one distributional parameter
    # without a default value
    if(length(dparams) == 0) {
        # equivalence between base R and MASS::fitdistr distribution names
        corresp <- function(distName) {
            switch(
                distName,
                beta = "beta",
                cauchy = "cauchy",
                chisq = "chi-squared",
                exp = "exponential",
                f = "f",
                gamma = "gamma",
                geom = "geometric",
                lnorm = "log-normal",
                logis = "logistic",
                norm = "normal",
                nbinom = "negative binomial",
                pois = "poisson",
                t = dt,
                weibull = "weibull",
                NULL
            )
        }

        # initial value for some distributions
        initVal <- function(distName) {
            switch(
                distName,
                beta = list(shape1 = 1, shape2 = 1),
                chisq = list(df = 1),
                f = list(df1 = 1, df2 = 2),
                t = list(df = 1),
                NULL
            )
        }

        suppressWarnings({
            if(!is.null(corresp(distribution))) {
                if(is.null(initVal(distribution))) {
                    dparams <- MASS::fitdistr(x = smp,
                                              densfun = corresp(distribution))$estimate
                } else {
                    dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution),
                                              start = initVal(distribution))$estimate
                }
            }
        })
    }

    theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))

    if (detrend) {
        if (identity) {
            slope <- 1
            intercept <- 0
        } else {
            xCoords <- do.call(qFunc, c(list(p = qprobs), dparams))
            yCoords <- do.call(quantile, list(x = smp, probs = qprobs, type = qtype))

            slope <- diff(yCoords) / diff(xCoords)
            intercept <- yCoords[1] - slope * xCoords[1]
        }

        # calculate new ys for the detrended sample
        dSmp <- NULL
        for (i in 1:n) {
            lSmp <- slope * theoretical[i] + intercept
            dSmp[i] <- smp[i] - lSmp
        }

        out <- data.frame(sample = dSmp, theoretical = theoretical)
    } else {
        out <- data.frame(sample = smp, theoretical = theoretical)
    }

}

#' @title  return a dataframe with x, y coordinates for plotting QQ line
QQ_line <- function(x, distribution = "norm", dparams = list(), qprobs = c(.25, .75),
                    detrend, identity, qtype) {
    # distributional function
    qFunc <- eval(parse(text = paste0("q", distribution)))

    smp <- sort(x)
    n <- length(smp)
    quantiles <- ppoints(n)

    # automatically estimate parameters with MLE, only if no parameters are
    # provided with dparams and there are at least one distributional parameter
    # without a default value
    if(length(dparams) == 0) {
        # equivalence between base R and MASS::fitdistr distribution names
        corresp <- function(distName) {
            switch(
                distName,
                beta = "beta",
                cauchy = "cauchy",
                chisq = "chi-squared",
                exp = "exponential",
                f = "f",
                gamma = "gamma",
                geom = "geometric",
                lnorm = "log-normal",
                logis = "logistic",
                norm = "normal",
                nbinom = "negative binomial",
                pois = "poisson",
                t = dt,
                weibull = "weibull",
                NULL
            )
        }

        # initial value for some distributions
        initVal <- function(distName) {
            switch(
                distName,
                beta = list(shape1 = 1, shape2 = 1),
                chisq = list(df = 1),
                f = list(df1 = 1, df2 = 2),
                t = list(df = 1),
                NULL
            )
        }

        suppressWarnings({
            if(!is.null(corresp(distribution))) {
                if(is.null(initVal(distribution))) {
                    dparams <- MASS::fitdistr(x = smp,
                                              densfun = corresp(distribution))$estimate
                } else {
                    dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution),
                                              start = initVal(distribution))$estimate
                }
            }
        })
    }

    theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))

    if (detrend) {
        out <- data.frame(xx = c(min(theoretical), max(theoretical)))
        out$yy <- 0
    } else {

        if (identity) {
            slope <- 1
            intercept <- 0
        } else {
            xCoords <- do.call(qFunc, c(list(p = qprobs), dparams))
            yCoords <- do.call(quantile, list(x = smp, probs = qprobs, type = qtype))
            slope <- diff(yCoords) / diff(xCoords)
            intercept <- yCoords[1] - slope * xCoords[1]
        }

        out <- data.frame(xx = c(min(theoretical), max(theoretical)))
        out$yy <- slope * out$xx + intercept
    }

    out
}

#' @title Plot Q-Q plot by plotly
#' @param x vector of data to plot QQ
#' @param df_info dataframe of information shown on hover
#' @param info_names character vector to replace names(info) to show on hover
#' @param title plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param qprobs numeric vector with length 2. Rrepresents the quantitles used to
#' \construct the Q-Q line
#' @export
QQ_plot <- function(x, dist= "norm", dparams= list(), df_info= NULL, info_names= NULL,
                    title= "", xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                    detrend = FALSE, identity= FALSE, qtype= 7) {
    # error handling
    if (!(dist %in% c(
        "beta", "cauchy", "chisq", "exp", "f", "gamma", "geom",  "lnorm", "logis",
        "norm", "nbinom", "pois", "t", "weibull")) &
        length(dparams) == 0 &
        table(sapply(formals(eval(parse(text = paste0("q", dist)))),
                     typeof))["symbol"] > 1) {
        stop(
            "MLE is currently not supported for custom distributions.\n",
            "Please provide all the custom distribution parameters to 'dparams'.",
            call. = FALSE
        )
    }

    # error handling
    if (qtype < 1 | qtype > 9) {
        stop("Please provide a valid quantile type: ",
             "'qtype' must be between 1 and 9.",
             call. = FALSE)
    }
    if (length(qprobs) != 2) {
        stop("'qprobs' must have length two.",
             call = FALSE)
    }
    if (sum(qprobs > 1) + sum(qprobs < 0)) {
        stop("'qprobs' cannot have any elements outside the probability domain [0,1].",
             call = FALSE)
    }

    #... QQ points ...
    # Calculate distribution data
    qq_point <- QQ_points(x, distribution = dist, dparams = dparams, qprobs = qprobs,
                          detrend = detrend, identity = identity,  qtype = qtype)

    df <- data.frame(x= qq_point$theoretical, y= qq_point$sample)

    #... QQ lint ...
    qq_line <- QQ_line(x, distribution = dist, dparams = dparams, qprobs = qprobs,
                       detrend = detrend, identity = identity,  qtype = qtype)

    #--- coordinates limits
    xmi <- min(df$x, na.rm = T)
    xma <-  max(df$x, na.rm = T)
    rg <- (xma - xmi) * 0.1
    xmin= xmi - rg
    xmax = xma + rg

    ymi <- min(df$y, na.rm = T)
    yma <-  max(df$y, na.rm = T)
    rg <- (yma - ymi) * 0.1
    ymin= ymi - rg
    ymax = yma + rg

    #--- Hover text information
    if (! is.null(df_info)) {
        df <- cbind(df, df_info)
        if (is.null(info_names)) info_names <- names(df_info)
        info_names <- c("theoretical", "sample", info_names)
    } else {
        info_names <- c("theoretical", "sample")
    }
    df <- df %>% filter(!is.na(y))
    df$x <- round(df$x, 3)
    df$y <- round(df$y, 3)
    hText <- apply(df, 1, function(x) paste(info_names, x, sep= ": "))
    hText <- apply(hText, 2, function(x) paste(x, collapse = " <br> "))


    #--- Plot
    plot_ly(x= qq_line$xx, y = qq_line$yy, type = 'scatter', mode = 'lines',
            color = I("#708090"), hoverinfo= "skip") %>%
        add_markers(data= df, x= ~x, y= ~y, color = I("steelblue"), hoverinfo= "text",
                    text= hText) %>%
        plotly::layout(showlegend= FALSE, title= title,
                       xaxis= list(title= xlab, zeroline= FALSE, range= c(xmin, xmax)),
                       yaxis= list(title= ylab, zeroline= FALSE, range= c(ymin, ymax)))
}

# * Plot histogram for process capability -----
#' @title Plot histogram for process capability
#' @param x vector of data to plot histogram
#' @param mean numeric, if NUll will be count by x
#' @param StDev_overall numeric, if NUll will be count by x
#' @param StDev_within numeric, if NUll will be count by x
#' @param plot_spec c(USL= 1, target= 0, LSL= -1)
#' @param title plot title
#' @export
Plot_hist_norm <- function(x, mean= NULL, StDev_overall= NULL,
                           StDev_within= NULL, plot_spec= NULL,
                           title= "Sample Distribution", xlab= "", ylab= "") {
    # plot_spec: named vector contains USL, target and LSL

    if (is.null(mean)) {mean <- mean(x, na.rm = T)}
    if (is.null(StDev_overall)) {StDev_overall <- sd(x, na.rm = T)}
    if (! is.na(StDev_overall)) {
        p <- ggplot(data.frame(x=c(mean-4*StDev_overall, mean+4*StDev_overall)), aes(x=x))
    } else {
        if (! is.na(StDev_within)) {
            p <- ggplot(data.frame(x=c(mean-4*StDev_within, mean+4*StDev_within)), aes(x=x))
        } else {
            p <- ggplot()
        }
    }

    p <- p + geom_histogram(data=data.frame(x= x),
                            aes(x=x, y = ..density..), bins =30,
                            color="white", fill="#9FB6CD")
    if (! is.na(StDev_overall)) {
        p <- p + stat_function(fun = dnorm, color="#4682B4",
                               args = list(mean =  mean, sd = StDev_overall),
                               show.legend = T)
    } else {
        if (! is.na(StDev_within)) {
            p <- p + stat_function(fun = dnorm, color="#DAA520", linetype= "dashed",
                                   args = list(mean =  mean, sd = StDev_within),
                                   show.legend = T)
        }
    }

    p <- p + geom_hline(yintercept = 0, color= "slategrey") +
        labs(title= title, x= xlab, y= ylab) +
        coord_cartesian(expand = FALSE) +
        theme(plot.title = element_text(size= 18, hjust = 0.5, vjust = 1),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    if (!is.null(plot_spec)) {
        if (! is.na( plot_spec["target"])) {
            p <- p + geom_vline(xintercept =  plot_spec["target"],
                                color = "#CD2626", linetype = "dashed")
        }

        if (! is.na(plot_spec["USL"])) {
            p <- p + geom_vline(xintercept = plot_spec["USL"],
                                color = "#CD2626", linetype = "dashed")
        }

        if (! is.na(plot_spec["LSL"])) {
            p <- p + geom_vline(xintercept =plot_spec["LSL"],
                                color = "#CD2626", linetype = "dashed")
        }
    }
    p
}

# Plot attribute histogram
#' @title Plot histogram for process capability
#' @param plot_spec c(USL= 1, target= 0, LSL= -1)
#' @param title plot title
#' @export
Plot_hist_attr <- function(x, plot_spec= NULL, title= "Sample Distribution",
                           xlab= "", ylab= "") {
    bins <- round(max(x) - min(x)) + 1
    bins <- ifelse(bins >= 30, 30, bins)
    p <- ggplot(data.frame(x= x), aes(x= x)) +
        geom_histogram(bins = bins, color= "white", fill= "#9FB6CD") +
        geom_hline(yintercept = 0, color= "slategrey") +
        labs(title= title, x= xlab, y= ylab) +
        coord_cartesian(expand = FALSE) +
        theme_min(base_size = 16, xGrid_major = FALSE, border_color= NA)

    if (!is.null(plot_spec)) {
        if (! is.na( plot_spec["target"])) {
            p <- p + geom_vline(xintercept =  plot_spec["target"],
                                color = "#CD2626", linetype = "dashed")
        }

        if (! is.na(plot_spec["USL"])) {
            p <- p + geom_vline(xintercept = plot_spec["USL"],
                                color = "#CD2626", linetype = "dashed")
        }

        if (! is.na(plot_spec["LSL"])) {
            p <- p + geom_vline(xintercept =plot_spec["LSL"],
                                color = "#CD2626", linetype = "dashed")
        }
    }
    p
}


# * Functionn to plot process capability analysis -----
#' @title Normal process capability analysis
#' @param x numeric vector
#' @param subgroup vector with the same length of x to label the subgroup
#' @param USL
#' @param target
#' @param LSL
#' @export
fun_proCap_norm_plot <- function(x, subgroup= 1, USL = NA, target = NA, LSL = NA,
                                 df_info= NULL, info_names= NULL) {
    # Form data frame
    df <- data.frame(subgroup = subgroup, value = x)

    # Get sample mean, standard deviation, Cp, Cpk, Pp, Ppk
    sample_statistic <- fun_Cp(x, subgroup, USL, target, LSL)
    sample_size <- sample_statistic$sample_size
    StDev_overall <- sample_statistic$StDev_overall
    StDev_within <- sample_statistic$StDev_within
    mean <- sample_statistic$sample_mean
    Cp <- sample_statistic$Cp
    Cpk <- sample_statistic$Cpk
    Pp <- sample_statistic$Pp
    Ppk <- sample_statistic$Ppk
    ppm_within <- sample_statistic$ppm_within
    ppm_overall <- sample_statistic$ppm_overall
    ppm_obs <- sample_statistic$ppm_obs

    #--- Plot
    # Q-Q plot
    qq_plot <- QQ_plot(x, dist= "norm", dparams= list(),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_norm(x, mean= mean, StDev_overall= StDev_overall,
                                StDev_within= StDev_within, plot_spec= plot_spec,
                                title= "Sample Distribution", xlab= "", ylab= "")

    # Return
    list(sample_statistic= sample_statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}


# * Function to plot binomial process capability analysis -----
#' @title Binomial process capability analysis
#' @param size vector of (inspection) size
#' @param x vector of data
#' @param USL
#' @param target
#' @param LSL
#' @export
fun_proCap_binom_plot <- function(size, x, USL = NA, target = NA, LSL = NA,
                                  df_info= NULL, info_names= NULL) {
    # Form data frame
    df <- data.frame(Size= size, Value= x)

    #
    base_size <- round(mean(df$Size, na.rm = TRUE))
    df$Size <- round(runif(nrow(df), round(0.85 * base_size), round(1.05 * base_size)))
    base_size <- round(mean(df$Size, na.rm = TRUE))

    #--- Statistics
    df <- df %>%
        mutate(perDefect = Value / Size,
               Value_norm = round(perDefect * base_size))
    # Probability
    prob <- sum(df$Value) / sum(df$Size)
    fit <- fitdist(df$Value_norm, "binom", start = list(prob= prob),
                   fix.arg = list(size= base_size))
    prob <- fit$estimate["prob"]

    # Sample size
    sample_size <- length(x)

    # Defective fraction
    b <- binom.test(round(prob * base_size), base_size)
    # mean, lowr CI, upper CI
    fracDefect <- round(as.vector(c(b$estimate, b$conf.int)), 6)

    # Form dataframe of statistic
    statistic <- data.frame(sample_size= sample_size, USL= USL, target= target,
                            LSL= LSL, fracDefect_mean= fracDefect[1],
                            lower_CI= fracDefect[2], upper_CI= fracDefect[3])

    #--- Plot
    # Q-Q plot
    qq_plot <- QQ_plot(x, dist= "binom", dparams= list(size= base_size, prob= prob),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_attr(x, plot_spec, title = "Sample distribution",
                                xlab= "Defective counts", ylab= "Frequency")

    # Return
    list(sample_statistic= statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}

# * Function to plot poisson process capability analysis -----
#' @title Poisson process capability analysis
#' @param size vector of (inspection) size
#' @param x vector of data
#' @param iu inspection unit, how many samples inspected per unit
#' @param USL
#' @param target
#' @param LSL
#' @export
fun_proCap_pois_plot <- function(size, x, iu, USL = NA, target = NA, LSL = NA,
                                 df_info= NULL, info_names= NULL) {
    # Form data frame
    df <- data.frame(Size= size, Value= x)

    #
    base_size <- round(mean(df$Size, na.rm = TRUE))
    df$Size <- round(runif(nrow(df), round(0.85 * base_size), round(1.05 * base_size)))
    base_size <- round(mean(df$Size, na.rm = TRUE))

    #--- Statistics
    sample_size <- length(x)

    df <- df %>%
        mutate(No_iu = Size / iu, u = round(Value / No_iu))
    ubar <- round(fitdist(df$u, "pois")$estimate)
    #ubar <- round(sum(df$Value) / sum(df$No_iu))

    head <- c("Mean    ", "Lower CI", "Upper CI")
    pois <- poisson.test(ubar)

    dpu <- round(as.vector(c(pois$estimate, pois$conf.int)), 2)

    statistic <- data.frame(sample_size= sample_size, USL= USL, target= target,
                            LSL= LSL, inspect_unit= iu,
                            defect_mean= dpu[1], lower_CI= dpu[2],
                            upper_CI= dpu[3])

    #--- Plot
    # Q-Q plot
    qq_plot <- QQ_plot(x, dist= "pois", dparams= list(lambda= ubar),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_attr(x, plot_spec, title = "Sample distribution",
                                xlab= "Defect per Unit", ylab= "Frequency")

    # Return
    list(sample_statistic= statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}
