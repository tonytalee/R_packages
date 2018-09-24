#=== * Sub functions =====
#' @title  Count plot width by numbers of x factors
#' @param x vecter to count the unique value and than calculate plot width
#' according the unique value
#' @export
P_width <- function(x) {
    # Unique x variable
    uni_x <- length(unique(as.character(x)))
    if (uni_x == 1) {
        p_width <- 300
    } else {
        if (uni_x >= 17) {
            p_width <- 1000
        } else {
            p_width <- (uni_x - 2) * 40 + 400
        }
    }
    p_width
}

#' @title Preprocess the input of Plotly_XXX functions and turn them to requirement
#' of parameters for plotting.
#' @param df dataframe contain x and y
#' @param info vector with the same rows as df to show information for hover
#' @param xcolor character vector as the same rows of df for color variable
#' @param xtext vector for showing text on plot
##' @param name_vector character vector to names x-label, y-label and info
#' @param xaxis_stype list to layout for setting x-axis style
#' @param yaxis_stype list to layout for setting y-axis style
#' @param xGrid logical, show x grid
#' @param yGrid logical, show y grid
#' @param xFactor logical, turn x-variable to factor
#' @param x_extend(y_extend) c(lower_extend_ratio, upper_extend_ratio), coordinate extension
#' @param mode character, mode of plotly trace, default "X"
#' @export
P_preProcess <- function(df, info= NULL, xcolor= NULL, xtext= NULL, name_vector= NULL,
                         xaxis_style= NULL, yaxis_style= NULL, xGrid= TRUE, yGrid= TRUE,
                         xFactor= FALSE, x_extend = NULL, y_extend = NULL, mode= "X") {
    # This function is to pre-process the input of Plotly_xxx functions and turn them to
    # requirement of parameters for plotting.

    #--- Get labels
    if (is.null(name_vector)) {
        xlab <- names(df)[1]
        ylab <- names(df)[2]
        if (! is.null(info)) {
            infolabe <- names(info)
        } else {
            infolab <- ""
        }
    } else {
        xlab <- name_vector[1]
        ylab <- name_vector[2]
        infolab <- name_vector[3]
    }

    #--- Rename df and set x variable as factor if required
    names(df) <- c("x", "y")
    if (xFactor) df$x <- factor(df$x, levels = unique(df$x))

    #--- Make sure color variable is character
    if (! is.null(xcolor)) {xcolor <- as.character(xcolor)}

    #--- Count the width of plot
    p_width = P_width(df$x)

    #--- Marker style
    if (mode == "markers") {
        marker_lineSize <- 1
        alpha <- 0.7
    } else {
        marker_lineSize <- 2
        alpha <- 1
    }

    marker_style <- list( size= 10, line= list(color= "white", width= marker_lineSize))

    #--- Hover text
    if (is.numeric(df$x)) {
        xtext <- round(df$x, 4)
        if (! is.null(xaxis_style$tickformat)) {
            if (xaxis_style$tickformat == "%") xtext <- scales::percent(xtext)
        }
    } else {
        xtext <- df$x
    }

    if (is.numeric(df$y)) {
        ytext <- round(df$y, 4)
        if (! is.null(yaxis_style$tickformat)) {
            if (yaxis_style$tickformat == "%") ytext <- scales::percent(ytext)
        }
    } else {
        ytext <- df$y
    }

    hText <- paste(xlab, ": ", xtext, "<br>", ylab, ": ", ytext)
    if (! is.null(info)) hText <- paste(hText, "<br>", infolab, ": ", info)

    #--- Axis range
    if (! is.null(x_extend)) {
        mi <- min(df$x, na.rm = T)
        ma <-  max(df$x, na.rm = T)
        rg <- ma - mi

        if (! is.na(x_extend[1])) {
            Xmin <- mi - rg * x_extend[1]
        } else {
            Xmin <- mi - rg * 0.05
        }

        if (! is.na(x_extend[2])) {
            Xmax <- ma + rg * x_extend[2]
        } else {
            Xmax <- ma + rg * 0.05
        }

        xrange <- c(Xmin, Xmax)
    } else {
        xrange <- NULL
    }

    if (! is.null(y_extend)) {
        mi <- min(df$y, na.rm = T)
        ma <-  max(df$y, na.rm = T)
        rg <- ma - mi
        if (! is.na(y_extend[1])) {
            Ymin <- mi - rg * y_extend[1]
        } else {
            Ymin <- mi - rg * 0.05
        }

        if (! is.na(y_extend[2])) {
            Ymax <- ma + rg * y_extend[2]
        } else {
            Ymax <- ma + rg * 0.05
        }

        yrange <- c(Ymin, Ymax)
    } else {
        yrange <- NULL
    }

    #--- Axis style
    xsty <- list(title= xlab, zeroline= F, showgrid= xGrid)
    if (! is.null(xrange)) {xsty <- c(xsty, list(range= xrange))}
    if (! is.null(xaxis_style)) {xsty <- c(xsty, xaxis_style)}

    ysty <- list(title= ylab, zeroline= F, showgrid= yGrid)
    if (! is.null(yrange)) {ysty <- c(ysty, list(range= yrange))}
    if (! is.null(yaxis_style)) {ysty <- c(ysty, yaxis_style)}


    #--- Return
    list(df= df, xcolor= xcolor, xlab= xlab, ylab= ylab, infolab= infolab,
         p_width= p_width, marker_style= marker_style, hText= hText, xsty= xsty,
         ysty= ysty)
}

#=== * Major functions =====
#' @title Function to plotly scatter-line chart
#' @param df dataframe contain x and y
#' @param info vector with the same rows as df to show information for hover
#' @param xcolor character vector as the same rows of df for color variable
#' @param xtext vector for showing text on plot
#' @param name_vector character vector to names x-label, y-label and info
#' @param margin list(l= 1, r= 1, t= 1, b= 1), plot margin
#' @param xaxis_stype list to layout for setting x-axis style
#' @param yaxis_stype list to layout for setting y-axis style
#' @param showLegend logical, show legend on plot
#' @param xGrid logical, show x grid
#' @param yGrid logical, show y grid
#' @param xFactor logical, turn x-variable to factor
#' @param x_extend(y_extend) c(lower_extend_ratio, upper_extend_ratio), coordinate extension
#' @param title title of plot
#' @param title_font list of plotly titlefont
#' @param height height of plot (width is auto calculated by unique values of x)
#' @param textposition posistion of text to show on plot
#' @param mode plot_ly mode, "lines+markers", "markers" or "lollipop"
#' @export
Plotly_scatter <- function(df, info= NULL, xcolor= NULL, xtext= NULL, name_vector= NULL,
                           margin= NULL, xaxis_style= NULL, yaxis_style= NULL,
                           showLegend= FALSE, xGrid= TRUE, yGrid= TRUE,
                           xFactor= FALSE, x_extend = NULL, y_extend = NULL,
                           title= "", title_font= NULL, height= NULL,
                           textposition= "top", mode= "lines+markers") {
    # df: dataframe
    # info: vector with the same rows as df to show information for hover
    # xcolor: character vector as color variable
    # xtext: vector to show text on plot
    # name_vector: character vector to names x-label, y-label and info
    # margin: list to layout for setting margin
    # xaxis_style & yaxis_style: list to layout for setting axis style
    # showLegend: logical, show legend on plot
    # xGrid & yGrid: logical, show grid
    # xFactor: logical, turn x-variable to factor
    # x_extend, y_extend: c(lower_extend_ratio, upper_extend_ratio), coordinate extension
    # title: title of plot
    # mode: plot_ly mode, "lines+markers", "markers" or "lollipop"

    #--- Pre-process
    # names of df will be changed to ("x", "y") after pre-processing
    li_prePorcess <- P_preProcess(df, info, xcolor, xtext, name_vector, xaxis_style,
                                  yaxis_style, xGrid, yGrid, xFactor, x_extend,
                                  y_extend, mode)
    df = li_prePorcess$df
    xcolor = li_prePorcess$xcolor
    xlab = li_prePorcess$xlab
    ylab = li_prePorcess$ylab
    infolab = li_prePorcess$infolab
    p_width = li_prePorcess$p_width
    marker_style = li_prePorcess$marker_style
    hText = li_prePorcess$hText
    xsty = li_prePorcess$xsty
    ysty = li_prePorcess$ysty
    #---

    #--- plotly
    p <- plot_ly(df, text= ~hText, width= p_width, height = height,
                 showlegend= showLegend)

    if (mode == "lollipop") {
        # Lollipop plot
        ymin <- min(min(df$y, na.rm = T) * 0.95, ysty$range[1], na.rm = T)

        if (is.null(xcolor)) {
            p <- p %>%
                add_markers(x= ~x, y= ~y, color= I("steelblue"), hoverinfo= "text",
                            marker = marker_style)
        } else {
            p <- p %>%
                add_markers(x= ~x, y= ~y, color= xcolor, colors= color_set,
                            hoverinfo= "text", marker = marker_style)
        }

        p <- p %>%
            add_segments(x= ~x, xend= ~x, y= ymin, yend= ~y, color= I("steelblue"),
                         size= I(1), hoverinfo= "skip")
    } else {
        # Scatter plot
        if (is.null(xcolor)) {
            p <- p %>%
                add_trace(x= ~x, y= ~y, type= "scatter", mode= mode, hoverinfo= "text",
                          color= I("steelblue"), marker = marker_style)
        } else {
            p <- p %>%
                add_trace(x= ~x, y= ~y, type= "scatter", mode= mode, hoverinfo= "text",
                          color= xcolor, colors= color_set, marker = marker_style)
        }
    }
    if (! is.null(xtext)) {
        p <-p %>% add_text(x= ~x, y= ~y, text= xtext, textposition= textposition,
                           color= I("black"), hoverinfo= "skip")
    }
    p %>% plotly::layout(title= title, xaxis = xsty, yaxis = ysty, margin= margin,
                         autosize= F, titlefont= title_font)
}

#' @title Function to plotly box or violin chart
#' @param df dataframe contain x and y
#' @param info vector with the same rows as df to show information for hover
#' @param xcolor character vector as the same rows of df for color variable
#' @param xtext vector for showing text on plot
#' @param name_vector character vector to names x-label, y-label and info
#' @param margin list(l= 1, r= 1, t= 1, b= 1), plot margin
#' @param xaxis_stype list to layout for setting x-axis style
#' @param yaxis_stype list to layout for setting y-axis style
#' @param showLegend logical, show legend on plot
#' @param xGrid logical, show x grid
#' @param yGrid logical, show y grid
#' @param xFactor logical, turn x-variable to factor
#' @param x_extend(y_extend) c(lower_extend_ratio, upper_extend_ratio), coordinate extension
#' @param title title of plot
#' @param title_font list of plotly titlefont
#' @param height height of plot (width is auto calculated by unique values of x)
#' @param textposition posistion of text to show on plot
#' @param groupmode logical, boxes/violins are grouped (or overlapped)
#' @param type c("violin", "box")
#' @export
Plotly_box <- function(df, info= NULL, xcolor= NULL, xtext= NULL, name_vector= NULL,
                       margin= NULL, xaxis_style= NULL, yaxis_style= NULL,
                       showLegend= FALSE, xGrid= TRUE, yGrid= TRUE,
                       xFactor= FALSE, x_extend = NULL, y_extend = NULL,
                       title= "", title_font= NULL, height= NULL, textposition= "top",
                       groupmode= FALSE, type= "box") {
    # type: c("violin", "box")
    if (! type %in% c("violin", "box")) stop('type must one of c("violin", "box")')


    #--- Pre-process
    # names of df will be changed to ("x", "y") after pre-processing
    li_prePorcess <- P_preProcess(df, info, xcolor, xtext, name_vector, xaxis_style,
                                  yaxis_style, xGrid, yGrid, xFactor, x_extend,
                                  y_extend)
    df = li_prePorcess$df
    xcolor = li_prePorcess$xcolor
    xlab = li_prePorcess$xlab
    ylab = li_prePorcess$ylab
    infolab = li_prePorcess$infolab
    p_width = li_prePorcess$p_width
    marker_style = li_prePorcess$marker_style
    hText = li_prePorcess$hText
    xsty = li_prePorcess$xsty
    ysty = li_prePorcess$ysty
    #---

    #--- plotly
    p <- plot_ly(df, width= p_width, height = height, showlegend= showLegend)

    if (is.null(xcolor)) {
        p <- p %>%
            add_trace(x= ~x, y= ~y, type= type, color= I("steelblue"))
    } else {
        p <- p %>%
            add_trace(x= ~x, y= ~y, type= type, color= xcolor, colors= color_set)
    }

    if (type == "violin") {
        p <- p  %>% plotly::style(box= list(visible= T), meanline= list(visible= T),
                            marker= list(symbol= 4))
        if (groupmode) {p <- p %>% plotly::layout( violinmode = 'group')}
    } else {
        p <- p %>% plotly::style(marker= list(symbol= 4))
        if (groupmode) {p <- p %>% plotly::layout( boxmode = 'group')}
    }


    p  %>%
        plotly::layout(title= title, xaxis = xsty, yaxis = ysty, margin= margin,
                       autosize= F, titlefont= title_font)
}

#' @title Function to plotly bar chart
#' @param df dataframe contain x and y
#' @param info vector with the same rows as df to show information for hover
#' @param xcolor character vector as the same rows of df for color variable
#' @param xtext vector for showing text on plot
#' @param name_vector character vector to names x-label, y-label and info
#' @param margin list(l= 1, r= 1, t= 1, b= 1), plot margin
#' @param xaxis_stype list to layout for setting x-axis style
#' @param yaxis_stype list to layout for setting y-axis style
#' @param showLegend logical, show legend on plot
#' @param xGrid logical, show x grid
#' @param yGrid logical, show y grid
#' @param xFactor logical, turn x-variable to factor
#' @param x_extend(y_extend) c(lower_extend_ratio, upper_extend_ratio), coordinate extension
#' @param title title of plot
#' @param title_font list of plotly titlefont
#' @param height height of plot (width is auto calculated by unique values of x)
#' @param textposition posistion of text to show on plot
#' @param barstack barmode c("group", "stack"), TRUE for "stack"
#' @export
Plotly_bar <- function(df, info= NULL, xcolor= NULL, xtext= NULL, name_vector= NULL,
                       margin= NULL, xaxis_style= NULL, yaxis_style= NULL,
                       showLegend= FALSE, xGrid= TRUE, yGrid= TRUE,
                       xFactor= FALSE, x_extend = NULL, y_extend = NULL, title= "",
                       title_font= NULL, height= NULL, textposition= "top",
                       barstack= TRUE) {
    # barstack: barmode c("group", "stack"), TRUE for "stack"
    if (barstack) {
        barmode <- "stack"
    } else {
        barmode <- "group"
    }

    #--- Pre-process
    # names of df will be changed to ("x", "y") after pre-processing
    li_prePorcess <- P_preProcess(df, info, xcolor, xtext, name_vector, xaxis_style,
                                  yaxis_style, xGrid, yGrid, xFactor, x_extend,
                                  y_extend)
    df = li_prePorcess$df
    xcolor = li_prePorcess$xcolor
    xlab = li_prePorcess$xlab
    ylab = li_prePorcess$ylab
    infolab = li_prePorcess$infolab
    p_width = li_prePorcess$p_width
    marker_style = li_prePorcess$marker_style
    hText = li_prePorcess$hText
    xsty = li_prePorcess$xsty
    ysty = li_prePorcess$ysty
    #---

    # plotly
    p <- plot_ly(df, width= p_width, height = height, text= ~hText,
                 showlegend= showLegend)

    if (is.null(xcolor)) {
        p <- p %>%
            add_bars(x= ~x, y= ~y, color= I("steelblue"), hoverinfo= "text")
    } else {
        p <- p %>%
            add_bars(x= ~x, y= ~y, color= xcolor, colors= color_set, hoverinfo= "text")
    }
    if (! is.null(xtext)) {
        p <-p %>% add_text(x= ~x, y= ~y, text= xtext, textposition= textposition,
                           color= I("black"), hoverinfo= "skip")
    }

    p %>% plotly::layout(title= title, xaxis = xsty, yaxis = ysty, margin= margin,
                         autosize= F, barmode = barmode, titlefont= title_font)
}

#' @title Function to plot Pareto
#' @param category vector of data for category
#' @param frequency vector of data as the same length of category
#' @param title title of plot
#' @param xlab x label
#' @param ylab y label
#' @param topCap integer, numbers of top categoried to show on plot, 0 for all
#' @export
Plotly_pareto <- function(category, frequency, title= "Pareto Plot",
                          xlab= "", ylab= "Frequency", topCat= 0) {

    df <- data.frame(category= category, frequency= frequency) %>%
        arrange(desc(frequency)) %>%
        mutate(cumfreq = cumsum(frequency),
               cumpercent = percent(round(cumfreq / sum(frequency), 2)))

    df$category <- factor(df$category, levels = df$category)

    if (topCat != 0) df <- df[1:topCat, ]

    # Plot width
    p_width <- P_width(df$category )

    ymax <- max(df$cumfreq) * 1.2
    plot_ly(df, x= ~category, y= ~frequency, type= "bar", name= "Frequency",
            color= I("steelblue"), width= p_width, hoverinfo= "text",
            text= ~paste("Category： ", category,
                         "<br>Frequency: ", frequency)) %>%
        add_trace(y= ~cumfreq, type= "scatter", mode= "lines+markers",
                  hoverinfo= "skip", name= "Ratio", color= I("#CD5C5C"),
                  marker = list( size= 10, line= list(color= "white", width= 2))) %>%
        add_text(y= ~cumfreq, text= ~cumpercent, textposition= "top",
                 color= I("black")) %>%
        plotly::layout(xaxis = list(title= xlab, tickangle=45, zeroline= F),
                       yaxis = list(title= ylab, zeroline= F, range= c(0, ymax)),
                       margin = list(b= 200), showlegend= FALSE)
}
