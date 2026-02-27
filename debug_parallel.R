library(crayon)


get_startup <- \(){
  get_rprofile <- \(){
    candidates <- c(
      Sys.getenv("R_PROFILE"),
      file.path(getwd(), ".Rprofile"),
      Sys.getenv("R_PROFILE_USER"),
      file.path(Sys.getenv("HOME"), ".Rprofile"),
      file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site")
    )

    Filter(file.exists, candidates) |>
      head(1)
  }

  list(
    R.home = R.home(),
    Env = Sys.getenv(),
    Libs = .libPaths(),
    Rprofile = get_rprofile(),
    CommandArgs = commandArgs()
  )
}




compare_startup <- \(x, y, show_identical=FALSE){

  diff_dlist <- \(a,b){
    a <- a |>
      as.data.frame() |>
      setNames("x") |>
      tibble::rownames_to_column("names")
    b <- b |>
      as.data.frame() |>
      setNames("y") |>
      tibble::rownames_to_column("names")
    both <- merge(a,b, by="names", all=TRUE) |>
      subset(!((x==y) %in% TRUE)) |>
      by(~names, \(row){
        cat(red(row$name, "\n"))
        cat(row$x, "\n\n")
        cat(row$y, "\n")
      })
  }

  compare_var <- \(var){
    if(identical(x[[var]], y[[var]])){
      cat(green(var, "is identical   ==============================\n"))
      if(show_identical){
        cat(x[[var]], "\n")
      }
    } else {
      cat(red(var, "is different   ================================\n"))
      if(is(x[[var]], "Dlist")){
        diff_dlist(x[[var]], y[[var]])
      } else {
        cat(x[[var]], "\n\n")
        cat(y[[var]], "\n")
      }
    }
  }

  compare_var("R.home")
  compare_var("Env")
  compare_var("Libs")
  compare_var("Rprofile")
  compare_var("CommandArgs")
}

# library(parallel)

# cl <- makeCluster(2)
# clusterExport(cl, "get_startup")

res_main    <- get_startup()
res_cluster <- clusterCall(cl, get_startup)

compare_startup(res_main, res_cluster[[1]])
