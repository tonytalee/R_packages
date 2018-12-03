# v0.1.4
# QQ plot add option of ggplot and plotly

# * Function to calculate process capability index -----
#' @title Calculate index for process capability, i.e. Cp, Cpk, Pp, Ppk...
#' @param data numeric vector for calculating process capability index
#' @param group vector with the same length of data to form group for within and
#'        between group
#' @param USL upper specification limit
#' @param targt spec target
#' @param LSL lower specification limit
#' @return list contains: sample_size, sample_mean, StDev_overall,
#'         StDev_within (within-group), Cp, Cpk, CpkU, CpkL, Pp, Ppk, PpkU, PpkL,
#'         ppm_within, ppm_overall, ppm_obs
#' @export
fun_Cp <- function(data,  group= 1, USL = NULL, target = NULL, LSL = NULL) {

    df <- data.frame(group = group, value = data)
    df <- df[! is.na(df$value), ]

    # d2
    d2 <- c(1.128,1.693,2.059,2.326,2.534,2.704,2.847,2.97,3.078)

    #=== count mean and standard variation of overall data
    sample_size <- nrow(df)
    sample_mean <- round(mean(df$value), 3)
    StDev_overall <- round(sd(df$value), 3)

    # within group analysis
    no_groups <- unique(group)
    if (length(no_groups) >= 2) {
        df1 <- df %>%
            group_by(group) %>%
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
QQ_points <- function(data, distribution = "norm", dparams = list(), qprobs = c(.25, .75),
                      detrend, identity, qtype) {
    # distributional function
    qFunc <- eval(parse(text = paste0("q", distribution)))

    oidx <- order(data)
    smp <- data[oidx]
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

        out <- data.frame(sample = dSmp, theoretical = theoretical, order= oidx)
    } else {
        out <- data.frame(sample = smp, theoretical = theoretical, order= oidx)
    }
}

#' @title  return a dataframe with x, y coordinates for plotting QQ line
QQ_line <- function(data, distribution = "norm", dparams = list(), qprobs = c(.25, .75),
                    detrend, identity, qtype) {
    # distributional function
    qFunc <- eval(parse(text = paste0("q", distribution)))

    smp <- sort(data)
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
#' @param data vector of data to plot QQ
#' @param group vector with the same length as data, to group data
#' @param dist character, theorectical distribution function to use,
#'        e.g. "norm", "binom", "pois"... (refer more to package 'qqplotr')
#' @param dparams list of additional parameters passed on to the previous chosen
#'        distribution funcion.
#' @param df_info dataframe of information shown on hover
#' @param info_names character vector with the same length as data,
#'        to replace names(info) to show on hover. It doesn't include x-variable and
#'        y-varialbe. xlab and ylab will be used for x-variable and y-variable.
#' @param title plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param qprobs numeric vector with length 2. Rrepresents the quantitles used to
#'        construct the Q-Q line
#' @param detrend logic, should the plot objexts be detrended?
#' @param identity logic, should an identity line be used as the reference line used to
#'        construct the confidence bands?
#' @param qtype interger between 1 and 9, Type of quantile algorithm to be used by
#'        the quantile function to construct the Q-Q line.
#' @param colors colors set
#' @param plotly logic, TRUE for plotly, FALSE for ggplot
#' @export
QQ_plot <- function(data, group= 1,  dist= "norm", dparams= list(), df_info= NULL,
                    info_names= NULL, title= "Q-Q Plot", xlab= "Theoretical",
                    ylab= "Sample", qprobs = c(.25, .75), detrend = FALSE,
                    identity= FALSE, qtype= 7, colors= color_set4, plotly= TRUE) {
    # if df_infor is vector
    if (is.vector(df_info)) df_info <- data_frame(info= df_info)

    # remove NA
    data <- data[! is.na(data)]
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

    # match levels of colors to groups
    groups <- unique(group)

    if (length(groups) > length(colors)) {
        n <- ceiling(length(groups) / length(colors))
        colors <- rep(colors, n)
    }

    # QQ points and QQ lines
    if (! is.null(df_info)) df_info_c <- data.frame()
    df <- data.frame()
    qqlines <- list()
    for (gp in groups) {
        x1 <- data[group == gp]
        df_info1 <- df_info[group == gp, ]

        qq_point <- QQ_points(x1, distribution = dist, dparams = dparams, qprobs = qprobs,
                          detrend = detrend, identity = identity,  qtype = qtype)
        # ordering if_info
        if (! is.null(df_info)) {
            df_info1 <- df_info1[qq_point$order, ]
            df_info_c <- rbind(df_info_c, df_info1)
        }

        # form dataframe for qq points
        df1 <- data.frame(x= qq_point$theoretical, y= qq_point$sample, group= gp)
        df <- rbind(df, df1)

        # form list for qq lines
        qqline <- QQ_line(x1, distribution = dist, dparams = dparams, qprobs = qprobs,
                   detrend = detrend, identity = identity,  qtype = qtype)
        qqlines[[gp]] <- qqline
    }
    group <- df$group
    df$group <- NULL

    #--- Plot
    if (plotly) {
        #--- coordinates limits
        xmi <- min(df$x, na.rm = T)
        xma <-  max(df$x, na.rm = T)
        rg <- (xma - xmi) * 0.2
        xmin= xmi - rg
        xmax = xma + rg

        ymi <- min(df$y, na.rm = T)
        yma <-  max(df$y, na.rm = T)
        rg <- (yma - ymi) * 0.2
        ymin= ymi - rg
        ymax = yma + rg

        #--- Hover text information
        if (! is.null(df_info)) {
            df <- cbind(df, df_info_c)
            if (is.null(info_names)) info_names <- names(df_info)
            info_names <- c(xlab, ylab, info_names)
        } else {
            info_names <- c(xlab, ylab)
        }
        df <- df %>% filter(!is.na(y))
        df$x <- round(df$x, 3)
        df$y <- round(df$y, 3)
        hText <- apply(df, 1, function(x) paste(info_names, x, sep= ": "))
        hText <- apply(hText, 2, function(x) paste(x, collapse = " <br> "))

        #--- Plotly
        if (length(groups) == 1) {
            p <- plot_ly(data= df, x= ~x, y= ~y,  type = 'scatter', mode = 'markers',
                     color= I(colors[1]), hoverinfo= "text", text= hText)
            legendx= FALSE
        } else {
            p <- plot_ly(data= df, x= ~x, y= ~y,  type = 'scatter', mode = 'markers',
                     color = ~group, colors = colors, hoverinfo= "text", text= hText)
            legendx= TRUE
        }
        for (i in 1:length(groups)) {
            qq_line <- qqlines[[i]]
            p <- p %>%
                add_segments(x= qq_line$xx[1], xend= qq_line$xx[2],
                             y= qq_line$yy[1], yend= qq_line$yy[2],
                             color= I(colors[i]), showlegend= FALSE)
        }
        p %>%
            plotly::layout(showlegend= legendx, title= title,
                       xaxis= list(title= xlab, zeroline= FALSE, range= c(xmin, xmax)),
                       yaxis= list(title= ylab, zeroline= FALSE, range= c(ymin, ymax)))
    } else {
        df$group <- group

        if (length(groups) == 1) {
            slope <- diff(qq_line$yy)/diff(qq_line$xx)
            int <- qq_line$yy[1L] - slope * qq_line$xx[1L]

            p <- ggplot(df, aes(x= x, y= y)) +
                geom_point(color= "steelblue4", alpha= 0.7, size= 2) +
                geom_abline(slope = slope, intercept = int, color= "#708090")
        } else {
            p <- ggplot(df, aes(x= x, y= y, color= group, group= group)) +
                geom_point(alpha = 0.7, size= 2)

            for (i in 1:length(groups)) {
                qq_line <- qqlines[[i]]
                slope <- diff(qq_line$yy)/diff(qq_line$xx)
                int <- qq_line$yy[1L] - slope * qq_line$xx[1L]
                qq_line <- qqlines[[i]]
                p <- p + geom_abline(slope = slope, intercept = int, color= colors[i])
            }
        }
        p + labs(title= title, x= xlab, y= ylab) +
            scale_color_manual(values = colors) +
            theme(plot.title = element_text(size= 14, hjust = 0.5, vjust = 1),
                  plot.background = element_rect(colour = NA, fill = NA),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.border = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  axis.ticks = element_blank())
    }
}

# * Plot histogram for process capability -----
#' @title Plot histogram for process capability
#' @param data vector of data to plot histogram
#' @param mean numeric, if NUll will be count by data
#' @param StDev_overall numeric, if NUll will be count by data
#' @param StDev_within numeric, if NUll will be count by data
#' @param plot_spec vector with three elements, c(USL= 1, target= 0, LSL= -1)
#' @param title plot title
#' @param xlab x label
#' @param ylab y label
#' @param bar_fill bar filling color
#' @return histogram
#' @export
Plot_hist_norm <- function(data, mean= NULL, StDev_overall= NULL,
                           StDev_within= NULL, plot_spec= NULL,
                           title= "Sample Distribution", xlab= "", ylab= "",
                           bar_fill= "#36648B") {
    # plot_spec: named vector contains USL, target and LSL

    if (is.null(mean)) {mean <- mean(data, na.rm = T)}
    if (is.null(StDev_overall)) {StDev_overall <- sd(data, na.rm = T)}
    if (! is.na(StDev_overall)) {
        p <- ggplot(data.frame(x=c(mean-4*StDev_overall, mean+4*StDev_overall)), aes(x=x))
    } else {
        if (! is.na(StDev_within)) {
            p <- ggplot(data.frame(x=c(mean-4*StDev_within, mean+4*StDev_within)), aes(x=x))
        } else {
            p <- ggplot()
        }
    }

    p <- p + geom_histogram(data=data.frame(x= data),
                            aes(x=x, y = ..density..), bins =30,
                            color="white", fill= bar_fill)
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
        theme(plot.title = element_text(size= 16, hjust = 0.5, vjust = 1),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(size = 12),
              axis.text.x = element_text(size = 12),
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
#' @param data vector of data to plot histogram
#' @param plot_spec vector with three elements, c(USL= 1, target= 0, LSL= -1)
#' @param title plot title
#' @param xlab x label
#' @param ylab y label
#' @param bar_fill bar filling color
#' @return histogram
#' @export
Plot_hist_attr <- function(data, plot_spec= NULL, title= "Sample Distribution",
                           xlab= "", ylab= "", bar_fill= "#36648B") {
    bins <- round(max(data) - min(data)) + 1
    bins <- ifelse(bins >= 30, 30, bins)
    p <- ggplot(data.frame(x= data), aes(x= x)) +
        geom_histogram(bins = bins, color= "white", fill= bar_fill) +
        geom_hline(yintercept = 0, color= "slategrey") +
        labs(title= title, x= xlab, y= ylab) +
        coord_cartesian(expand = FALSE) +
        theme(plot.title = element_text(size= 16, hjust = 0.5, vjust = 1),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(size = 12),
              axis.text.x = element_text(size = 12),
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


# * Functionn to plot process capability analysis -----
#' @title Normal process capability analysis
#' @param data numeric vector
#' @param group vector with the same length of data to form group for within and
#'        between group
#' @param USL upper specification limit of process
#' @param target  target (center) of process
#' @param LSL lower specification limit of process
#' @param df_info dataframe of information shown on hover
#' @param info_names character vector with the same length as data,
#'        to replace names(info) to show on hover. It doesn't include x-variable and
#'        y-varialbe. xlab and ylab will be used for x-variable and y-variable.
#' @param bar_fill bar filling color
#' @return list contains sample statistics, QQ plot and histogram
#' @export
fun_proCap_norm_plot <- function(data, group= 1, USL = NA, target = NA, LSL = NA,
                                 df_info= NULL, info_names= NULL, bar_fill= "#36648B") {
    # Form data frame
    df <- data.frame(group = group, value = data)

    # Get sample mean, standard deviation, Cp, Cpk, Pp, Ppk
    sample_statistic <- fun_Cp(data, group, USL, target, LSL)
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
    qq_plot <- QQ_plot(data, dist= "norm", dparams= list(),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_norm(data, mean= mean, StDev_overall= StDev_overall,
                                StDev_within= StDev_within, plot_spec= plot_spec,
                                title= "Sample Distribution", xlab= "", ylab= "",
                                bar_fill= bar_fill)

    # Return
    list(sample_statistic= sample_statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}


# * Function to plot binomial process capability analysis -----
#' @title Binomial process capability analysis
#' @param size numeric vector of inspection size
#' @param data numeric vector with same length as size, events or defects counts
#' @param USL upper specification limit of process
#' @param target  target (center) of process
#' @param LSL lower specification limit of process
#' @param df_info dataframe of information shown on hover
#' @param info_names character vector with the same length as data,
#'        to replace names(info) to show on hover. It doesn't include x-variable and
#'        y-varialbe. xlab and ylab will be used for x-variable and y-variable.
#' @param bar_fill bar filling color
#' @return list contains sample statistics, QQ plot and histogram
#' @export
fun_proCap_binom_plot <- function(size, data, USL = NA, target = NA, LSL = NA,
                                  df_info= NULL, info_names= NULL, bar_fill= "#36648B") {
    # Form data frame
    df <- data.frame(Size= size, Value= data)

    #
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
    sample_size <- sum(data, na.rm = TRUE)

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
    qq_plot <- QQ_plot(data, dist= "binom", dparams= list(size= base_size, prob= prob),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_attr(data, plot_spec, title = "Sample distribution",
                                xlab= "Defective counts", ylab= "Frequency",
                                bar_fill = bar_fill)

    # Return
    list(sample_statistic= statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}

# * Function to plot poisson process capability analysis -----
#' @title Poisson process capability analysis
#' @param size numeric vector of inspection size
#' @param data numeric vector with same length as size, events or defects counts
#' @param iu, ispection unit
#' @param USL upper specification limit of process
#' @param target  target (center) of process
#' @param LSL lower specification limit of process
#' @param df_info dataframe of information shown on hover
#' @param info_names character vector with the same length as data,
#'        to replace names(info) to show on hover. It doesn't include x-variable and
#'        y-varialbe. xlab and ylab will be used for x-variable and y-variable.
#' @param bar_fill bar filling color
#' @return list contains sample statistics, QQ plot and histogram
#' @export
fun_proCap_pois_plot <- function(size, data, iu, USL = NA, target = NA, LSL = NA,
    df_info= NULL, info_names= NULL, bar_fill= "#36648B") {

    # Form data frame
    df <- data.frame(Size= size, Value= data)

    #--- Statistics
    sample_size <- sum(data, na.rm = TRUE)

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
    qq_plot <- QQ_plot(data, dist= "pois", dparams= list(lambda= ubar),
                       df_info= df_info, info_names= info_names, title= "",
                       xlab= "Theoretical", ylab= "Sample", qprobs = c(.25, .75),
                       detrend = FALSE, identity= FALSE, qtype= 7)

    # Histogram
    plot_spec <- c("USL" = USL, "target" = target, "LSL" = LSL)
    hist_plot <- Plot_hist_attr(data, plot_spec, title = "Sample distribution",
                                xlab= "Defect per Unit", ylab= "Frequency",
                                bar_fill = bar_fill)

    # Return
    list(sample_statistic= statistic, qq_plot= qq_plot, hist_plot= hist_plot)
}
