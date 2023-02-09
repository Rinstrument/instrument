
collapse_categories = function(x) {
  return(apply(x, 2, \(y){ match(y, sort(unique(y))) }))
}
