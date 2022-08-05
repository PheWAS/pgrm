meta.F <- function(b.est, se){
  #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
  b.F = sum(b.est / se^2) / sum(1 / se^2)
  se.F = 1 / sqrt(sum(1 / se^2))
  p.F = pchisq( (b.F / se.F)^2, df = 1, lower = F)
  return(list(b.F = b.F, se.F = se.F, p.F = p.F))
}

forest.plot <- function(x, intervals, labels = NULL, main = NULL, xlab = "Effect size",
                        pchs = rep(19,length(x)), cols = rep("black", length(x)),
                        cexs = rep(1,length(x))){
  K = length(x)
  stopifnot(nrow(intervals) == K)
  plot(0, col="white", xlim = c( min(c(intervals[,1],0) - 0.05), max(c(intervals[,2],0) + 0.05)),
       ylim = c(0, K+1), xlab = xlab, ylab = "", yaxt = "n",main = main)
  axis(2, at = K:1, labels = labels, cex.axis = 0.8)
  arrows(intervals[,1], K:1, intervals[,2], K:1,
         code = 3, angle = 90, length = 0.02, col = cols)
  points(x, K:1, pch = pchs, cex = cexs, col = cols)
  abline(v = 0,lty = 2)
}

b.est = log(c(r1$rOR),log(r2$rOR))
ci= log(matrix(c(r1$L95, r1$rU95,
                 r2$L95, r2$rU95), byrow = T, ncol = 2))
se = (ci[,2] - ci[,1])/(2*1.96)

meta.res = meta.F(b.est, se)
meta.res

cbind(OR = exp(meta.res$b.F), low = exp(meta.res$b.F - 1.96*meta.res$se.F),
      up = exp(meta.res$b.F + 1.96*meta.res$se.F), pval = meta.res$p.F)



meta.F <- function(b.est, se){
  #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
  b.F = sum(b.est / se^2) / sum(1 / se^2)
  se.F = 1 / sqrt(sum(1 / se^2))
  p.F = pchisq( (b.F / se.F)^2, df = 1, lower = F)
  return(list(b.F = b.F, se.F = se.F, p.F = p.F))
}


b.est = log(c(1.50, 1.38, 1.39)) #beta is logOR for case-control data
ci = log(matrix(c(1.25, 1.79,
                  1.17, 1.63,
                  1.15, 1.68), byrow = T, ncol = 2))
se = (ci[,2] - ci[,1])/(2*1.96) #length of 95%CI is 2*1.96*SE
meta.res = meta.F(b.est, se)
meta.res


cbind(OR = exp(meta.res$b.F), low = exp(meta.res$b.F - 1.96*meta.res$se.F),
      up = exp(meta.res$b.F + 1.96*meta.res$se.F), pval = meta.res$p.F)

b.est = c(b.est, meta.res$b.F)
ci = rbind(ci, c(meta.res$b.F + c(-1,1)*1.96*meta.res$se.F))
labs = c("Discv", "Rep1", "Rep2", "Meta")
main.txt = "rs11984041 Stroke/LVD"
forest.plot(b.est, ci, labels = labs, main = main.txt, xlab = "logOR",
            pchs = c(19, , 18), cexs = c(.8, , 1.3), cols = c(1, , 4))



forest.plot <- function(x, intervals, labels = NULL, main = NULL, xlab = "Effect size",
                        pchs = rep(19,length(x)), cols = rep("black", length(x)),
                        cexs = rep(1,length(x))){
  K = length(x)
  stopifnot(nrow(intervals) == K)
  plot(0, col="white", xlim = c( min(c(intervals[,1],0) - 0.05), max(c(intervals[,2],0) + 0.05)),
       ylim = c(0, K+1), xlab = xlab, ylab = "", yaxt = "n",main = main)
  axis(2, at = K:1, labels = labels, cex.axis = 0.8)
  arrows(intervals[,1], K:1, intervals[,2], K:1,
         code = 3, angle = 90, length = 0.02, col = cols)
  points(x, K:1, pch = pchs, cex = cexs, col = cols)
  abline(v = 0,lty = 2)
}
b.est = c(b.est, meta.res$b.F)
ci = rbind(ci, c(meta.res$b.F + c(-1,1)*1.96*meta.res$se.F))
labs = c("Discv", "Rep1", "Rep2", "Meta")
main.txt = "rs11984041 Stroke/LVD"
forest.plot(b.est, ci, labels = labs, main = main.txt, xlab = "logOR",
            pchs = c(19, 19, 19, 18), cexs = c(.8, .8, .8, 1.3), cols = c(1, 1, 1, 4))
