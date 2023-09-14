library(tidyverse)

names <- c("IVW", "Penalized IVW", "Robust IVW",
           "Simple median", "Weighted median", "Penalized weighted median", 
            "Penalized robust IVW", "MR-Egger", "Penalized MR-Egger",
            "Robust MR-Egger", "Penalized robust MR-Egger")

# we can update colors here...    
#colors <- c("#F8766D", "#F8766D", "#7CAE00",
#            "#7CAE00", "#C77CFF", "#C77CFF",
#            "#00BFC4", "#00BFC4", "#00BFC4",
#            "#F8766D", "#7CAE00")

# hc_theme
colors <- c("#7CB5EC", "#434348", "#90ED7D",
            "#F7A35C", "#8085E9", "#F15C80",
            "#E4D354", "#8085E8", "#8D4653",
            "#91E8E1", "#A18376")

lines <- c("solid", "dashed", "solid",
           "dashed", "solid",  "dashed",
           "dotted", "solid", "dashed",
           "dotted", "dotted")

df_attr <- tibble(Method = names, colors = colors, lines = lines)

scatter_plot_mr <- function(df_g, df_c, BETA.GWAS, BETA.pQTL, SE.GWAS, SE.pQTL, NAME.pQTL, is.cond=FALSE, plot.error=FALSE) {

    df_c <- df_c %>% select(Method:Estimate) %>% pivot_wider(names_from=Type, values_from=Estimate, values_fill=0.)
    if (!"intercept" %in% colnames(df_c)) {
        df_c <- df_c %>% mutate(intercept = 0.0)
    }

    # subset to same Methods    
    df_attr_f <- df_attr %>% filter(Method %in% df_c$Method)

    # just 0-out the errorbars if we don't want to plot    
    scale <- ifelse(plot.error, 1, 0)

    # if its cond or not
    gwas_str <- ifelse(is.cond,
                       "HG cond GWAS",
                       "HG GWAS")
    pqtl_str <- ifelse(is.cond,
                       paste0(str_to_title(NAME.pQTL), " GDF15 cond pQTL"),
                       paste0(str_to_title(NAME.pQTL), " GDF15 pQTL"))

    ggplot(data = df_g, aes(y = {{BETA.GWAS}}, x = {{BETA.pQTL}})) +
        geom_point() +
        geom_errorbar(aes(ymin = scale * ({{BETA.GWAS}} - qnorm(0.975)*{{SE.GWAS}}), ymax = scale * ({{BETA.GWAS}} + qnorm(0.975)*{{SE.GWAS}})), alpha=0.7) +
        geom_errorbar(aes(xmin = scale * ({{BETA.pQTL}} - qnorm(0.975)*{{SE.pQTL}}), xmax = scale * ({{BETA.pQTL}} + qnorm(0.975)*{{SE.pQTL}})), alpha=0.7) +
        geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
        geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
        geom_abline(data = df_c,
                    #aes(intercept = intercept, slope = slope, color = Method, linetype = Method),
                    aes(intercept = intercept, slope = slope, color = Method),
                    show.legend = TRUE, size = 1) +
		scale_colour_manual(name = "Method", labels = df_attr_f$Method, breaks = df_attr_f$Method, values = df_attr_f$colors) +
        #scale_linetype_manual(name="Method", labels = df_attr_f$Method, breaks = df_attr_f$Method, values = df_attr_f$lines) +
        ylab(bquote(.(gwas_str)~~hat(beta))) +
        xlab(bquote(.(pqtl_str)~~hat(beta))) +
        xlim(c(-0.6, 0.6)) +
        ylim(c(-1.5, 1.5)) +
        theme(
              plot.title = element_text(size = rel(1.5), face = "bold"),
              # Background
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
              panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
              legend.key = element_rect(fill = "white")
        )

}
