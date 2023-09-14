library(tidyverse)
source("MR.plot.R")

rdir <- "/home1/nmancuso/projects/projects/hg_project/GDF15_project/results/mr/"

# file format are
# [PROT-STUDY].GDF15.23andMe.[REF-PANEL]LD.tsv.gz
# [PROT-STUDY].GDF15.23andMe.[REF-PANEL]LD.cond.tsv.gz
# for PROT-STUDY: OLINK, ROCHE
#     REF-PANEL: UKBB, 1000G


df_main <- read_tsv(paste0(rdir, "ROCHE.GDF15.23andMe.UKBBLD.tsv.gz")) %>%
            filter(Type == "slope")

df_attr_f <- df_attr %>% filter(Method %in% df_main$Method)
ivw_color <- df_attr %>% filter(Method == "IVW") %>% pull(colors)

# comparison of mr methods on main data using UKBB LD
forest_methods <- ggplot(data = df_main, aes(y = Method, x = Estimate, xmin = L95, xmax = U95, color = Method)) +
                    geom_point() +
                    geom_errorbar(width = 0.25) +
                    geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                    scale_color_manual(name = "Method", labels = df_attr_f$Method, breaks = df_attr_f$Method, values = df_attr_f$colors) +
                    scale_y_discrete(limits=rev) +
                    xlab("Effect of GDF15 levels on HG risk") +
                    xlim(c(-1.15, 0.1)) +
                    theme(
                          panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                          panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                          legend.key = element_rect(fill = "white"),
                          legend.position="none"
                    )

# comparison ld ref panel
df_okg <- read_tsv(paste0(rdir, "ROCHE.GDF15.23andMe.1000GLD.tsv.gz")) %>%
            filter(Type == "slope")


df_plot <- bind_rows(df_main %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (UKBB LD)")),
                     df_okg  %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (1000G LD)")))

forest_ld <- ggplot(data = df_plot, aes(y = Method, x = Estimate, xmin = L95, xmax = U95)) +
                    geom_point(color = ivw_color) +
                    geom_errorbar(width = 0.25, color = ivw_color) +
                    geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                    xlab("Effect of GDF15 levels on HG risk") +
                    xlim(c(-1.15, 0.1)) +
                    theme(
                          panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                          panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                          legend.key = element_rect(fill = "white"),
                          legend.position="none"
                    )

# comparison ld ref panel
df_olink <- read_tsv(paste0(rdir, "OLINK.GDF15.23andMe.UKBBLD.tsv.gz")) %>%
                filter(Type == "slope")


df_plot <- bind_rows(df_main   %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (Roche)")),
                     df_olink  %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (Olink)")))

forest_assay <- ggplot(data = df_plot, aes(y = Method, x = Estimate, xmin = L95, xmax = U95)) +
                    geom_point(color = ivw_color) +
                    geom_errorbar(width = 0.25, color = ivw_color) +
                    geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                    xlab("Effect of GDF15 levels on HG risk") +
                    xlim(c(-1.15, 0.1)) +
                    theme(
                          panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                          panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                          legend.key = element_rect(fill = "white"),
                          legend.position="none"
                    )

# comparison conditional
df_cond <- read_tsv(paste0(rdir, "ROCHE.GDF15.23andMe.UKBBLD.cond.tsv.gz")) %>%
                filter(Type == "slope")


df_plot <- bind_rows(df_main   %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (marginal)")),
                     df_cond  %>% filter(Method == "IVW") %>% mutate(Method = paste0(Method, " (conditional)")))

forest_cond <- ggplot(data = df_plot, aes(y = Method, x = Estimate, xmin = L95, xmax = U95)) +
                    geom_point(color = ivw_color) +
                    geom_errorbar(width = 0.25, color = ivw_color) +
                    geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                    xlab("Effect of GDF15 levels on HG risk") +
                    xlim(c(-1.15, 0.1)) +
                    theme(
                          panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                          panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                          legend.key = element_rect(fill = "white"),
                          legend.position="none"
                    )

ggsave(paste0(rdir, "GDF15.methods.pdf"), forest_methods, height=3.5, width=5.0, units="in")
ggsave(paste0(rdir, "GDF15.ld.pdf"), forest_ld, height=3.5, width=5.0, units="in")
ggsave(paste0(rdir, "GDF15.assay.pdf"), forest_assay, height=3.5, width=5.0, units="in")
ggsave(paste0(rdir, "GDF15.cond.pdf"), forest_cond, height=3.5, width=5.0, units="in")
