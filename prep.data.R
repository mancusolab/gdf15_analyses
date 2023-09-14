library(tidyverse)
library(BEDMatrix)
library(susieR)
library(coloc)

allele.qc = function(a1,a2,ref1,ref2) {
    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)

    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip

    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;

    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
    snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

    return(snp)
}

get_ukbb_data <- function(df_gwas, offset=FALSE) {

    ld_snp_path <- "../data/ukbb_ld/gdf15_annotated_221121.snp.dbsnp.txt.gz"
    ld_mat_path <- "../data/ukbb_ld/gdf15_annotated_221121.ld.RDat"

    df_snp <- read_tsv(ld_snp_path) %>%
        mutate(idx = 1:nrow(.))

    load(ld_mat_path)

    # Match summary data to input, record NA where summary data is missing
    df_j <- inner_join(df_gwas, df_snp, by=c("SNP.GWAS"="SNP")) %>%
        filter(!duplicated(SNP.GWAS))

    # QC / allele-flip the input and output
    qc = allele.qc(df_j$ALLELE1 , df_j$ALLELE0 , df_j$A1, df_j$A0)
    df_j <- df_j %>% mutate(flip = qc$flip, keep = qc$keep | SNP.GWAS == "rs1058587") %>%
        filter(keep)

    R <- R[df_j$idx, df_j$idx]
    if ( offset ) {
        res <- eigen(R)
        lambda <- min(res$values) + 1e-3
        R <- R + lambda * diag(ncol(R))
    }
    attr(R, "eigen") = eigen(R, symmetric = TRUE)
    colnames(R) <- df_j$SNP.GWAS
    rownames(R) <- df_j$SNP.GWAS

    # flip effects to match LD
    df_j <- df_j %>% mutate(BETA.GWAS = ifelse(flip, -BETA.GWAS, BETA.GWAS),
                            BETA.pQTL = ifelse(flip, -BETA.pQTL, BETA.pQTL))

    return (list(df_j=df_j, R=R))
}


read_plink_custom <- function(root, impute = c('none', 'avg', 'random')) {
    if(impute == 'random') {
        stop("The 'impute' random option has not been implemented.", call. = FALSE)
    }
    
    ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
    proot <- path.expand(root)
    
    bedfile <- paste(proot, ".bed", sep="")
    famfile <- paste(proot, ".fam", sep="")
    bimfile <- paste(proot, ".bim", sep="")
    
    ## Could change this code to use data.table
    bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    ## Set the dimensions
    geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))
    
    ## Convert to a matrix
    geno <- as.matrix(geno)
    if(impute == 'avg') {
        ## Check if any are missing
        geno_na <- is.na(geno)
        if(any(geno_na)) {
            means <- colMeans(geno, na.rm = TRUE)
            geno[geno_na] <- rep(means, colSums(geno_na))
        }
    }
    colnames(geno) <- bim[,2]
    rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

    list(bed=geno, fam=fam, bim=bim)    
}


get_okg_data <- function(df_gwas, offset=FALSE) {
    geno_path <- "../data/1000G_ld/1000G.gdf15.hg38.dbsnp"

    # Load in reference data
    genos = read_plink_custom(geno_path, impute="avg")

    # Match summary data to input, record NA where summary data is missing
    cur.join <- dplyr::inner_join(genos$bim, df_gwas, by=c("V2"="SNP.GWAS")) %>%
                    dplyr::rename(SNP.GWAS = V2) %>%
                    dplyr::filter(!duplicated(SNP.GWAS))
    cur.genos = scale(genos$bed[,cur.join$SNP.GWAS], center=TRUE, scale=TRUE)

    # QC / allele-flip the input and output
    qc = allele.qc( cur.join$ALLELE1 , cur.join$ALLELE0 , cur.join[,5] , cur.join[,6] )
    cur.join$flip = qc$flip

    # add small amount to diag for PSD reasons...
    R = cor(cur.genos)
    if ( offset ) {
        res <- eigen(R)
        lambda <- min(res$values) + 1e-3
        R <- R + lambda * diag(ncol(R))
    }
    attr(R, "eigen") = eigen(R, symmetric = TRUE)

    # flip effects to match LD
    cur.join <- cur.join %>%
                mutate(BETA.GWAS = ifelse(flip, -BETA.GWAS, BETA.GWAS),
                       BETA.pQTL = ifelse(flip, -BETA.pQTL, BETA.pQTL)) %>%
                as_tibble

    return (list(df_j=cur.join, R=R))
}

krig_data <- function(df_j, R, prefix=NULL) {
    # !! update; two rounds seems to work well
    # run one round of kriging to flip ambig alleles
    # 4 ~ qnorm(0.05 / 5000)
    thold <- 2
    zhat.pQTL <- df_j %>% mutate(Z.pQTL = BETA.pQTL / SE.pQTL) %>% pull(Z.pQTL)
    n.pQTL <- median(df_j$N.pQTL, na.rm=TRUE)
    cond.pQTL <- kriging_rss(zhat.pQTL, R, n = n.pQTL)
    p_flip <- cond.pQTL$conditional_dist %>% as_tibble %>% mutate(FLIP = (abs(z) > thold & logLR > thold)) %>% pull(FLIP)
    if (!is.null(prefix)) {
        ggsave(paste0(prefix, ".pQTL.krig.pdf"), cond.pQTL$plot)
    }

    cond.update.pQTL <- kriging_rss(ifelse(p_flip, -1, 1) * zhat.pQTL, R, n = n.pQTL)
    p_update_flip <- cond.update.pQTL$conditional_dist %>% as_tibble %>% mutate(FLIP = (abs(z) > thold & logLR > thold)) %>% pull(FLIP)
    if (!is.null(prefix)) {
        ggsave(paste0(prefix, ".pQTL.krig.update.pdf"), cond.update.pQTL$plot)
    }
    p_flip_final <- ifelse(p_update_flip, -1, 1) * ifelse(p_flip, -1, 1)

    zhat.GWAS <- df_j %>% mutate(Z.GWAS = BETA.GWAS / SE.GWAS) %>% pull(Z.GWAS)
    n.GWAS <- median(df_j$N.GWAS, na.rm=TRUE)
    cond.GWAS <- kriging_rss(zhat.GWAS, R, n = n.GWAS)
    g_flip <- cond.GWAS$conditional_dist %>% as_tibble %>% mutate(FLIP = (abs(z) > thold & logLR > thold)) %>% pull(FLIP)
    if (!is.null(prefix)) {
        ggsave(paste0(prefix, ".GWAS.krig.pdf"), cond.GWAS$plot)
    }

    cond.update.GWAS <- kriging_rss(ifelse(g_flip, -1, 1) * zhat.GWAS, R, n = n.GWAS)
    g_update_flip <- cond.update.GWAS$conditional_dist %>% as_tibble %>% mutate(FLIP = (abs(z) > thold & logLR > thold)) %>% pull(FLIP)
    if (!is.null(prefix)) {
        ggsave(paste0(prefix, ".GWAS.krig.update.pdf"), cond.update.GWAS$plot)
    }
    g_flip_final <- ifelse(g_update_flip, -1, 1) * ifelse(g_flip, -1, 1)

    df_j <- df_j %>% mutate(BETA.pQTL = BETA.pQTL * p_flip_final,
                            BETA.GWAS = BETA.GWAS * g_flip_final)
    return(df_j)
}
