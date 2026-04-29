
library(microbenchmark)
library(miceFast)

pred_list <- apply(pred, 1, \(x){colnames(pred)[x==1]}) |>
  purrr::imap(\(x, name){setdiff(x, name)}) |>
  purrr::keep(\(x){length(x) > 0})


mb <- microbenchmark(
  old = {
    dat_g <- dat_g_old
    xx <- mice::mice(
      dat_g,
      m = m,
      method = meth_g,
      predictorMatrix = pred,
      maxit = maxit,
      ridge = 1e-5,
      visitSequence = "monotone",
      printFlag = FALSE
    )
    res <- lapply(1:m, \(i){
      mice::complete(xx, i)
    })
  },
  new = {
    dat_g <- dat_g_old
    res <- lapply(1:m, \(i){
      purrr::iwalk(pred_list, \(x, y){
        dat_g[, y] <<- miceFast::fill_NA(
          dat_g,
          model = "lm_bayes",
          posit_y = y,
          posit_x = x,
          ridge = 1e-5
        )
        NULL
      })
    })
  }, times=10)
