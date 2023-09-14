source("prep.data.R")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    print("Usage: Rscript coloc.R SUMSTATS (1000G | UKBB) OUT")
    q()
}

sumstat_path <- args[1]
ld_name <- args[2]
out_prefix <- args[3]

df_gwas <- read_tsv(sumstat_path)

if (toupper(ld_name) == "1000G") {
    data <- get_okg_data(df_gwas, offset=FALSE)
} else if (toupper(ld_name) == "UKBB") {
    data <- get_ukbb_data(df_gwas, offset=FALSE)
} else {
    print("UNKNOWN LD PANEL!")
    q()
}

df_j <- data$df_j
R <- data$R

df_j <- krig_data(df_j, R, out_prefix)

pqtl_d <- list(snp = df_j$SNP.GWAS,
               position=df_j$BP.GRCh38, beta=df_j$BETA.pQTL, varbeta=df_j$SE.pQTL^2, type="quant", N=median(df_j$N.pQTL, na.rm=TRUE), sdY=1, LD = R)
gwas_d <- list(snp = df_j$SNP.GWAS,
               position=df_j$BP.GRCh38, beta=df_j$BETA.GWAS, varbeta=df_j$SE.GWAS^2, type="cc", N=median(df_j$N.GWAS, na.rm=TRUE), LD = R)

pqtl_fm <- runsusie(pqtl_d, check_prior=FALSE)
gwas_fm <- runsusie(gwas_d, check_prior=FALSE)

coloc_res <- coloc.susie(pqtl_fm, gwas_fm)
write_tsv(coloc_res$summary %>% as_tibble %>% rename(hit.pQTL = hit1, hit.GWAS = hit2, idx.pQTL = idx1, idx.GWAS = idx2),
          paste0(out_prefix, ".tsv.gz"))
