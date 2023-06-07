PlinkLD_transform <- function(plinkLD, all_keep_snps) {
    wkeep <- which(plinkLD$SNP_A %in% all_keep_snps & plinkLD$SNP_B %in% all_keep_snps)
    plinkLD <- plinkLD[wkeep, ]

    if (length(which(is.na(plinkLD))) > 0) {
        stop("")
    }

    ldJ <- plinkLD[, c("SNP_B", "SNP_A", "R")]
    names(ldJ) <- c("SNP_A", "SNP_B", "R")

    ldJ <- rbind(plinkLD[, c("SNP_A", "SNP_B", "R")], ldJ)

    return(ldJ)
}
