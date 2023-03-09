#' Traceplot of theta2 output
#' 
#' Produce a traceplot of a theta2 parameter. See below for a full list of 
#' parameters. 
#' 
#' This is currently implemented as a convenient wrapper around 
#' rstan::traceplot. In the future, we may implement our own method.
#' 
#' @param object a theta2Obj model object
#' @param pars parameter names as given by the theta2Obj output
#' @param ... arguments passed into rstan::traceplot. See ?rstan::traceplot for 
#' more information
#' 
#' @return a ggplot2 plot is prited upon return.
#' 
#' @importFrom rstan traceplot
#' @export 
traceplot.theta2Obj = function(object, pars, ...) {

	stanfit = object$stanfit

	return(rstan::traceplot(stanfit, pars = pars, ...))

}