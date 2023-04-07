
fill_match = function(x, d) {
  y = match(1:d, x)
  y[is.na(y)] = 0
  y[y > 0] = x
  return(y)
}
