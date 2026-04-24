library(purrr)
library(stringr)

results_filenames <- list.files() |>
  str_subset("results_onco.*\\.Rdata")

results_list <- map(results_filenames, \(f){
  tmp_env <- new.env()
  load(f, envir=tmp_env)
  tmp_env$results$file <- f
  message(str_c("read ", nrow(tmp_env$results), " rows from ", f))
  tmp_env$results
})

results <- list_rbind(results_list)

if(!all(sapply(2:length(results_list), \(i){
  identical(attr(results_list[[1]], "design_names"), attr(results_list[[i]], "design_names"))
}))){
  stop("Inconsistent design_names")
}

attr(results, "design_names") <- attr(results_list[[1]], "design_names")
