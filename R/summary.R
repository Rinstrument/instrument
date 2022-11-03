
# #' Summary of inirt
# #' 
# #' @importFrom rstan summary
# #' @export
# #' 
# summary.inirt = function(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), use_cache = TRUE, ...) {
#     class(object) = "stanfit"
#     attr(object, "package") = "rstan"
#     summary(object = object, pars = pars, probs = probs, use_cache = use_cache, ...)
# }