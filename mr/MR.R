source("prep.data.R")
source("MR.plot.R")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    print("Usage: Rscript MR.R SUMSTATS (1000G | UKBB) NAME OUT")
    q()
}

sumstat_path <- args[1]
ld_name <- args[2]
name <- args[3]
out_prefix <- args[4]

df_gwas <- read_tsv(sumstat_path) %>%
            filter(P.pQTL < 5e-8 | SNP.GWAS == "rs1058587") %>%
            filter(BETA.GWAS > -2)
if (toupper(ld_name) == "1000G") {
    data <- get_okg_data(df_gwas, offset=TRUE)
} else if (toupper(ld_name) == "UKBB") {
    data <- get_ukbb_data(df_gwas, offset=TRUE)
} else {
    print("UNKNOWN LD PANEL!")
    q()
}

df_j <- data$df_j
R <- data$R

df_f <- krig_data(df_j, R, out_prefix)
write_tsv(df_f, paste0(out_prefix, ".table.tsv.gz"))

library(MendelianRandomization)

MRInputObject <- mr_input(bx = df_f$BETA.pQTL,
                          bxse = df_f$SE.pQTL,
                          by = df_f$BETA.GWAS,
                          byse = df_f$SE.GWAS,
                          corr = R)

MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")

#cur.plot <- mr_plot(MRAllObject_all) + ylab(bquote("HG GWAS"~~hat(beta))) + xlab(bquote(.(name)~~"GDF15 pQTL"~~hat(beta)))

mr_out <- MRAllObject_all@Values %>%
            mutate(Method = ifelse(grepl("(intercept)", Method), paste0(lag(Method), " (intercept)"), Method)) %>%
            rename(L95 = `95% CI `,
                   U95 = ` `,
                   SE = `Std Error`,
                   P = `P-value`) %>%
            mutate(Type = ifelse(grepl("intercept", Method),
                                 "intercept",
                                 "slope")) %>%
            mutate(Method = gsub(" \\(intercept\\)", "", Method)) %>%
            relocate(Type, .after=Method)

cur.plot <- scatter_plot_mr(df_f, mr_out, BETA.GWAS, BETA.pQTL, SE.GWAS, SE.pQTL, name, is.cond=FALSE, plot.error=FALSE)
cur.plot.ivw <- scatter_plot_mr(df_f, mr_out %>% filter(Method == "IVW"), BETA.GWAS, BETA.pQTL, SE.GWAS, SE.pQTL, name, is.cond=FALSE, plot.error=TRUE)

write_tsv(mr_out, paste0(out_prefix, ".tsv.gz"))
ggsave(paste0(out_prefix,".pdf"), cur.plot, height=3.5, width=5.0, units="in")
ggsave(paste0(out_prefix,".IVW.pdf"), cur.plot.ivw, height=3.5, width=5.0, units="in")
