library(mice)
library(dplyr)

vars_y <- paste0("y", 0:12)
vars_R <- paste0("R", 1:11)

vars_imp <- c(
  vars_y,
  vars_R,
  "age",
  "trt"
)
meth <- make.method(dat[vars_imp])

# HbA1c normal imputieren
meth[vars_y] <- "pmm"

# Rescue: passiv, monoton
meth["R1"] <- "logreg"
for (j in 2:11) {
  meth[paste0("R", j)] <- paste0("~ I(pmax(R", j-1, ", R", j, "))")
}

# age, trt nicht imputieren
meth[c("age", "trt")] <- ""

pred <- make.predictorMatrix(dat[vars_imp])

# nichts sagt sich selbst vorher
diag(pred) <- 0

# trt nicht imputieren, aber trt wird NICHT als Prädiktor genutzt
pred[, "trt"] <- 0

# baseline immer Prädiktor
pred[, "y0"] <- 1

imp_list <- vector("list", 2)

for (g in 0:1) {
  dat_g <- dat %>%
    filter(trt == g) %>%
    select(all_of(vars_imp))

  imp_list[[g + 1]] <- mice(
    dat_g,
    m = 10,
    method = meth,
    predictorMatrix = pred,
    maxit = 20,
    seed = 123 + g,
    printFlag = FALSE
  )
}
imp_full <- lapply(1:10, function(m) {
  bind_rows(
    complete(imp_list[[1]], m),
    complete(imp_list[[2]], m)
  )
})
fit_models <- lapply(imp_full, function(d) {
  d$chg12 <- d$y12 - d$y0
  lm(chg12 ~ trt + age + y0, data = d)
})

library(mitools)

pooled <- MIcombine(fit_models)
summary(pooled)
est <- pooled$coefficients
se  <- sqrt(diag(pooled$variance))
df <- pooled$df

tval <- est / se
pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)

cbind(
  estimate = est,
  se = se,
  df = df,
  t = tval,
  p = pval
)


