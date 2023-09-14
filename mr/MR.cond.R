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
            filter(abs(BETA.GWAS) < 2)
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

library(MendelianRandomization)

# conditional
get_resid <- function(bhat, se, N, R, snp = "rs1058587") {

    zscore <- bhat / se
    beta_j <- solve(R, zscore) / sqrt(N)

    oidx <- rownames(R) != snp
    sidx <- rownames(R) == snp
    oname <- rownames(R)[oidx]

    Ro <- R[oname, oname]
    Rc <- R[snp, oname]
    Rresid <- Ro - Rc %*% t(Rc)

    se_o <- se[oidx]
    se_c <- se[sidx]

    vc <- se_c^2
    V = diag(se_o) %*% Ro %*% diag(se_o)
    C = diag(se_o) %*% Rc * se_c

    bhat_cond = (bhat[oidx] - diag(se_o) %*% Rc * zscore[sidx])
    se_cond = sqrt(diag(diag(se_o) %*% Rresid %*% diag(se_o)))
    print(head(bhat_cond))
    print(head(se_cond))

    list(snp = rownames(Ro), bhat.cond = c(bhat_cond), se.cond = c(se_cond), Ro = Rresid)
}

gwas.cond <- get_resid(df_f %>% pull(BETA.GWAS),
                       df_f %>% pull(SE.GWAS),
                       median(df_f$N.GWAS, na.rm=TRUE),
                       R)

pQTL.cond <- get_resid(df_f %>% pull(BETA.pQTL),
                       df_f %>% pull(SE.pQTL),
                       median(df_f$N.pQTL, na.rm=TRUE),
                       R)
df_c <- left_join(tibble(SNP.GWAS = gwas.cond$snp), df_f, by="SNP.GWAS") %>%
            arrange(BP.GRCh38)

df_c <- df_c %>% mutate(BETA.cond.pQTL = pQTL.cond$bhat.cond,
                        BETA.cond.GWAS = gwas.cond$bhat.cond,
                        SE.cond.pQTL = pQTL.cond$se.cond,
                        SE.cond.GWAS = gwas.cond$se.cond)

write_tsv(df_c, paste0(out_prefix, ".table.tsv.gz"))

MRInputObject <- mr_input(bx = pQTL.cond$bhat.cond,
                          bxse = pQTL.cond$se.cond,
                          by = gwas.cond$bhat.cond,
                          byse = gwas.cond$se.cond,
                          corr = gwas.cond$Ro)

MRAllObject_cond <- mr_allmethods(MRInputObject, method = "all")

#cond.plot <- mr_plot(MRAllObject_cond) + ylab(bquote("HG cond GWAS"~~hat(beta))) + xlab(bquote(.(name)~~"GDF15 cond pQTL"~~hat(beta)))

mr_out <- MRAllObject_cond@Values %>%
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

cond.plot <- scatter_plot_mr(df_c, mr_out, BETA.cond.GWAS, BETA.cond.pQTL, SE.cond.GWAS, SE.cond.pQTL, name, is.cond=TRUE, plot.error=FALSE)
cond.plot.ivw <- scatter_plot_mr(df_c, mr_out %>% filter(Method == "IVW"), BETA.cond.GWAS, BETA.cond.pQTL, SE.cond.GWAS, SE.cond.pQTL, name, is.cond=TRUE, plot.error=TRUE)

write_tsv(mr_out, paste0(out_prefix, ".tsv.gz"))
ggsave(paste0(out_prefix,".pdf"), cond.plot, height=3.5, width=5.0, units="in")
ggsave(paste0(out_prefix,".IVW.pdf"), cond.plot.ivw, height=3.5, width=5.0, units="in")
