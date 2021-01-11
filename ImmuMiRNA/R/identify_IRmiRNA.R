identify_IRmiRNA <- function(mRNA_exp, lncRNA_exp, adjusted, tum_put, pathways = pathways,
                             k = 0.995){
  par.cor <- function (mRNA_exp, miRNA_exp, tum_pur, adjusted)
  {
    fun_mtx_pcr <- function(x, y, z) {
      r12 = cor(t(x), t(y))
      r13 = cor(t(x), z)
      r23 = cor(z, t(y))
      r123 = r13 %*% r23
      rup = r12 - r123
      rd1 = sqrt(1 - r13 * r13)
      rd2 = sqrt(1 - r23 * r23)
      rd = rd1 %*% rd2
      rrr = rup/rd
      return(rrr)
    }
    inter_samples <- intersect(intersect(colnames(mRNA_exp),
                                         colnames(miRNA_exp)), names(tum_pur))
    inter_tumpur <- tum_pur[inter_samples]
    if (length(inter_tumpur) < 1) {
      stop("no same samples")
    }
    if (adjusted) {
      mRNA_exp_inter <- mRNA_exp[, inter_samples]
      miRNA_exp_inter <- miRNA_exp[, inter_samples]
    }
    else {
      mRNA_exp_inter_ori <- mRNA_exp[, inter_samples]
      miRNA_exp_inter_ori <- miRNA_exp[, inter_samples]
      mRow30 <- which(apply(mRNA_exp_inter_ori, 1, function(v) {
        return((sum(v == 0)/length(v)) >= 0.3)
      }))
      mRemv <- mRow30
      if (length(mRemv) == 0) {
        mRNA_out0 <- mRNA_exp_inter_ori
      }
      else {
        mRNA_out0 <- mRNA_exp_inter_ori[-(mRemv), ]
      }
      mRNA_exp_inter <- log2(mRNA_out0 + 0.001)
      lncRow50 <- which(apply(miRNA_exp_inter_ori, 1, quantile,
                              probs = 0.5) == 0)
      lncRow90 <- which(apply(miRNA_exp_inter_ori, 1, quantile,
                              probs = 0.9) <= 0.1)
      lncRemv <- union(lncRow50, lncRow90)
      if (length(lncRemv) == 0) {
        miRNA_out0 <- miRNA_exp_inter_ori
      }
      else {
        miRNA_out0 <- miRNA_exp_inter_ori[-(lncRemv),
                                            ]
      }
      miRNA_exp_inter <- log2(miRNA_out0 + 0.001)
    }
    n = length(inter_samples)
    gn = 1
    pcor <- fun_mtx_pcr(miRNA_exp_inter, mRNA_exp_inter, inter_tumpur)
    statistic <- pcor * sqrt((n - 2 - gn)/(1 - pcor^2))
    p.value <- 2 * pnorm(-abs(statistic))
    rownames(pcor) <- rownames(miRNA_exp_inter)
    rownames(p.value) <- rownames(miRNA_exp_inter)
    colnames(pcor) <- rownames(mRNA_exp_inter)
    colnames(p.value) <- rownames(mRNA_exp_inter)
    par.cor.res <- list(pcor, p.value)
    names(par.cor.res) <- c("pcor.value", "p.value")
    return(par.cor.res)
  }

  ii <- function (mRNA_exp, miRNA_exp, adjusted, tum_put, pathways = pathways,
            k = 0.995)
  {
    for (pkg in c("fgsea")) {
      if (!requireNamespace(pkg, quietly = T)) {
        stop(paste("The ", pkg, " package needed for this function to work. Please install it.",
                   sep = ""), call. = FALSE)
      }
    }
    library(fgsea, warn.conflicts = F)
    lnc_ori <- as.matrix(miRNA_exp)
    m_ori <- as.matrix(mRNA_exp)
    turpur_ori <- tum_put
    test_res <- par.cor(m_ori, lnc_ori, turpur_ori, adjusted)
    pValue_ori <- test_res$p.value
    pcorValue_ori <- test_res$pcor.value
    RS <- -log10(pValue_ori) * sign(pcorValue_ori)
    fgseaRes_all <- c()
    for (i in 1:nrow(RS)) {
      if (sum(is.infinite(RS[i, ])) != 0) {
        (next)()
      }
      ranks <- RS[i, ]
      fgseaRes <- fgsea(pathways, ranks, minSize = 1, maxSize = 5000,
                        nperm = 1000)
      sigValue <- c()
      for (j in 1:nrow(fgseaRes)) {
        if (fgseaRes$ES[j] > 0) {
          sig_ij <- 1 - 2 * fgseaRes$pval[j]
        }
        else {
          sig_ij <- 2 * fgseaRes$pval[j] - 1
        }
        sigValue <- c(sigValue, sig_ij)
      }
      miRNA <- rownames(RS)[i]
      fgseaRes_i <- cbind(miRNA, fgseaRes, sigValue)
      fgseaRes_all <- rbind(fgseaRes_all, fgseaRes_i)
    }
    sig_ind <- which(abs(fgseaRes_all$sigValue) >= k)
    sig_pairs <- fgseaRes_all[sig_ind, 1:2]
    sig_pairs <- as.matrix(sig_pairs)
    gsea.Res <- list(sig_pairs, fgseaRes_all)
    names(gsea.Res) <- c("sig_pairs", "fgseaRes_all")
    return(gsea.Res)
  }
}
