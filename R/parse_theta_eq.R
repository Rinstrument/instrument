

parse_regression_eq = function(model, data) {

  if(!str_detect(model, "=")) {
    stop("latent factor equations require a = between factor name and items. E.g. t1 = c(1:25).")
  }

  
  return(list())

}