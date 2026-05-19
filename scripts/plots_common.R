library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(CI.RCT.Sim)
library(purrr)
library(fs)

# global color scale
options(
  ggplot2.discrete.colour = function(...) scale_color_brewer(type="Qualitative", palette = "Set1", ...),
  ggplot2.discrete.fill   = function(...) scale_fill_brewer(type="Qualitative",palette = "Set1", ...)
)

# set global theme
theme_set(
  theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
)

# wrapper to set common graphics device settings
save_plot <- function(gg, filename){
  ggsave(
    filename=filename,
    plot = gg,
    scale=1,
    width=7,
    height=6,
    units="in",
    dpi=600
  )
}

# function to read multiple datasets
read_multiple <- function(results_filenames, check_design_names=TRUE){
  results_list <- map(results_filenames, \(f){
    tmp_env <- new.env()
    load(f, envir=tmp_env)
    tmp_env$results$file <- f
    message(str_c("read ", nrow(tmp_env$results), " rows from ", f))
    tmp_env$results
  })

  results <- list_rbind(results_list)

  if(check_design_names){
    if(!all(sapply(2:length(results_list), \(i){
      identical(attr(results_list[[1]], "design_names"), attr(results_list[[i]], "design_names"))
    }))){
      stop("Inconsistent design_names")
    }
  }

  attr(results, "design_names") <- attr(results_list[[1]], "design_names")
  results
}
