
#=== Colors

#' @title color_set
#' @export
color_set <- c("#4682B4", "#CD3700", "#458B00", "#DAA520", "#A020F0",
               "#8B4513", "#E9967A", "#708090", "#A2CD5A", "#A2B5CD",
               "#FF1493", "#8B2323", "#00BFFF", "#B4EEB4","#FFB90F")

#' @title color_set4
#' @export
color_set4 <- c("#4477AA", "#CC6677", "#117733", "#DAA520", "#88CCEE", "#AA4499",
                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#AA4466")

#' @title color_pair2
#' @export
color_pair2 <- c("#00008B", "#008B8B", "#EE7600", "#CDAD00", "#006400", "#6E8B3D",
                 "#8B0000", "#E9967A", "#483D8B", "#B23AEE")

#' @title color_light
#' @export
color_light <- c("#FCFCFC", "#D1D1D1", "#B0E2FF", "#A2CD5A", "#EEDC82", "#FFC1C1")

#=== Shapes
#' @title Multiple point shape
#' @export
shapes <- c(16,17,18,15,1,2,5,0,10,9,7)

#=== theme_min() =====
#' @title Create mini theme for ggplot2 object.
#' @param base_size A number for base font size.
#' @param face A string for fontface.
#' @param plotTitle_just Justment for plot title, left: 0, mid: 0.5, right: 1.
#' @param xAngle A number for angle of roation of X-axis text
#' @param yAngle A number for angle of rotation of Y-axis text
#' @param legend A string, oosition of lengend, one of top, bottom, left, right, none.
#' @param border_color A string for color of panel border.
#' @param strip_fill A string for color of filling strip.
#' @param grid_color A string for color of grid color.
#' @param xGrid_major Boolean, with major X-grid line or not.
#' @param yGrid_major Boolean, with major Y-grid line or not.
#' @param xText Boolean, showing X-axis text or not.
#' @param yText Boolean, showing Y-axis text or not.
#' @param tick Boolean, showing X and Y axis ticks or not.
#' @param plotBorder String, color of plot border.
#' @param legend_title Boolean, showing title of lengend or not.
#' @export
theme_min <- function (base_size=14, face='plain', plotTitle_just = 0.5,
                       xAngle=0, yAngle=0, legend = "right", border_color= "grey",
                       strip_fill= "gray90", grid_color= "grey80",
                       xGrid_major=TRUE, yGrid_major=TRUE, xText=TRUE, yText=TRUE,
                       tick=FALSE, plotBorder=NA, legend_title= TRUE) {
    if (xAngle >= 90) {
        ax.title.hjust <- 1
        ax.title.vjust <- 0.5
    } else {
        if (xAngle >= 10) {
            ax.title.hjust <- 1
            ax.title.vjust <- 1
        } else {
            ax.title.hjust <- 0.5
            ax.title.vjust <- 0.5
        }
    }
    theme_i <- theme(
        plot.title =  element_text(size = base_size, hjust = plotTitle_just, vjust = 1,
                                   colour = "black", lineheight = .9, face= 'bold'),
        plot.background = element_rect(colour = plotBorder, fill = NA, size = 1),
        axis.title.x = element_text(size = base_size*.9, angle = 0, hjust = 0.5,
                                    vjust = 1, colour = "black", face = "bold"),
        axis.title.y = element_text(size = base_size*.9, angle = 90, hjust = 0.5,
                                    vjust= 1, colour = "black", face = "bold"),
        legend.title = element_text(size = base_size*.9, angle = 0, colour = "black",
                                    face = "bold"),
        legend.text = element_text(size = base_size*.9, angle = 0, colour = "black"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),

        panel.background =  element_rect(fill = "white", colour = border_color),
        panel.border = element_rect(fill = NA, colour= border_color),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid_color, size = 0.2),

        #panel.margin =      unit(0.25, "lines"),

        strip.text = element_text(size = base_size*.8),
        strip.background = element_rect(fill = strip_fill, colour = NA,
                                        size = 0.3),

        # Legend
        legend.position = legend
    )
    if (! xGrid_major) theme_i <- theme_i + theme(panel.grid.major.x = element_blank())
    if (! yGrid_major) theme_i <- theme_i + theme(panel.grid.major.y = element_blank())
    if (! tick) theme_i <- theme_i + theme(axis.ticks = element_blank())
    if (! legend_title) theme_i <- theme_i + theme(legend.title=element_blank())

    if (xText) {
        theme_i <- theme_i + theme(
            axis.text.x = element_text(size = base_size*.8, angle = xAngle,
                                       hjust = ax.title.hjust,
                                       vjust = ax.title.vjust, colour = "black"))
    } else {
        theme_i <- theme_i + theme(axis.text.x = element_blank())
    }

    if (yText) {
        theme_i <- theme_i + theme(
            axis.text.y = element_text(size = base_size*.8, angle = yAngle,
                                       colour = "black"))
    } else {
        theme_i <- theme_i + theme(axis.text.y = element_blank())
    }

    return(theme_i)
}

#=== Customizing box plot =====
#' @title Customizing box plot
#' @param fill_var String, Variable name for aes coloring.
#' @param fill_col String, color of filling box.
#' @param color_col String, color of box outline and outliers
#' @param alpha Number, alpha of box.
#' @param size Number, size of box outline.
#' @param outlier_size Number, size of outliers.
#' @param outlier_alpha Number, alpha of outlier.
#' @param outlier_stroke Number, stroke of outlier.
#' @param width Number between 0 and 1, relative width of box
#' @param rm_na Logic, remove NA with warning or not
#' @export
boxplot <- function(fill_var=NA, fill_col= "steelblue", color_col= "steelblue4",
                    alpha= 0.5, size= 1, outlier_size= 1.5, outlier_alpha= 1,
                    outlier_stroke= 1, width= 0.6, rm_na = TRUE) {
    if (is.na(fill_var)) {
        geom_boxplot(outlier.size=outlier_size, outlier.shape=4, size=size,
                     fill=fill_col, color=color_col, alpha= alpha,
                     outlier.alpha= outlier_alpha, outlier.stroke= outlier_stroke,
                     width= width, na.rm= rm_na)
    } else {
        geom_boxplot(outlier.size=outlier_size, outlier.shape=4, size=size,
                     alpha= alpha, outlier.alpha= outlier_alpha,
                     outlier.stroke= outlier_stroke, width= width, na.rm= rm_na,
                     aes_string(fill= fill_var, color= fill_var))
    }
}

