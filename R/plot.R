#' 
#' @import ggplot2
#' @exportS3Method plot hltObj
plot.hltObj = function(x, ...) {
  
  args = list(...)
  
  if("type" %in% names(args)) {
    type = args$type
  }
  
  if("param" %in% names(args)) {
    param = args$param
  }
  
  if("item" %in% names(args)) {
    item = args$item
  }
  
  if("min" %in% names(args)) {
    min = args$min
  } else {
    min = -4
  }
  
  if("max" %in% names(args)) {
    max = args$max
  } else {
    max = 4
  }
  
  post = x$post
  nr = nrow(post)
  
  if(type == "trace") {
    ggplot(data.frame(param = 1:nr, y = post[, param]), aes(param, y)) + geom_line() + 
      xlab("iteration") + ylab("value") + get_theme()
  } else if(type == "icc") {
    plot.hltObj.icc(x, item = item, type = type, min = min, max = max)
  } else if(type == "iic") {
    iic_curve(mod = x, x = item, min = min, max = max)
  } else if(type == "tic") {
    tic_curve(x, min = min, max = max)
  }
  
}

#' 
#' @import ggplot2
#' @exportS3Method plot hltObjList
plot.hltObjList = function(x, ...) {
  args = list(...)
  
  if("type" %in% names(args)) {
    type = args$type
  }
  
  if("param" %in% names(args)) {
    param = args$param
  }
  
  if("item" %in% names(args)) {
    item = args$item
  }
  
  if("min" %in% names(args)) {
    min = args$min
  } else {
    min = -4
  }
  
  if("max" %in% names(args)) {
    max = args$max
  } else {
    max = 4
  }
  
  if(type == "trace") {
    m = merge_chains(x)
    post = m$post
    nr = nrow(post)
    post = cbind.data.frame(post, chain = rep(1:m$nchains, each = nrow(post) / m$nchains))
    ggplot(data.frame(param = rep(1:(nr/m$nchains), m$nchains), y = post[, param], 
                      each = post[, "chain"]), aes(param, y, group = each, color = each)) + 
      geom_line() + xlab("iteration") + ylab("value") + get_theme()
  }
  
}

get_theme = function() {
  theme_bw() + 
  theme(panel.grid.major = element_line(colour="#DDDDDD", size = (.5)),
        panel.grid.minor = element_line(size = (0.2), colour="#DDDDDD"),
        panel.border = element_rect(colour = "black", size = 1.5),
        axis.ticks = element_line(size = 1), axis.ticks.length = unit(.2, "cm"),
        plot.title = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        #legend.key.size = unit(1, 'lines'),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))
}

plot.hltObj.icc = function(mod, item, type, min, max) {
  if(type == "icc") {
    plt = icc_curve(mod, item, min = min, max = max)
  }
  
  return(plt)
}

#' @importFrom tidyr pivot_longer
icc_curve = function(mod, x, min = -4, max = 4) {
  
  if(any(grepl("^[a]", colnames(mod$post)))) {
    alpha = as.vector(summary(mod, param = "alpha", digits = 10, transpose = FALSE)["mean",])
  }
  
  kappa = summary(mod, param = "delta", digits = 10, transpose = FALSE)["mean",]
  kappa_names = names(kappa)
  kappa_id = as.numeric(gsub(".*[d]([^.]+)[_].*", "\\1", kappa_names))
  kappa_list = vector(length = length(unique(kappa_id)), mode = "list")
  
  for(i in 1:length(kappa_list)) {
    kappa_list[[i]] = kappa[kappa_id == i]
  }
  
  pxk = function(theta) {
    #stat = "mean"
    return(as.numeric(sapply(theta, function(theta) {
      cum_sum = exp(cumsum((alpha[x] * theta) - kappa_list[[x]]))
      par_cum_sum = cumsum(cum_sum)
      cum_sum / (1 + par_cum_sum)
    })))
  }
  
  if(length(kappa_list[[x]]) == 1) {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = sapply(seq(min, max, by = 0.1), pxk))
  } else {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = t(sapply(seq(min, max, by = 0.1), pxk)))
  }
  
  kappa_nms = names(kappa_list[[x]])
  names(pdata) = c("x", kappa_nms)
  pdata = tidyr::pivot_longer(pdata, cols = kappa_nms)
  plt = ggplot(pdata, aes(x, value, color = name)) + get_theme() + 
    xlab("theta") + ylab("p") + geom_line()
  
  return(plt)
}

iic_curve = function(mod, x, min, max) {
  
  if(any(grepl("^[a]", colnames(mod$post)))) {
    alpha = as.vector(summary(mod, param = "alpha", digits = 10, transpose = FALSE)["mean",])
  }
  
  kappa = summary(mod, param = "delta", digits = 10, transpose = FALSE)["mean",]
  kappa_names = names(kappa)
  kappa_id = as.numeric(gsub(".*[d]([^.]+)[_].*", "\\1", kappa_names))
  kappa_list = vector(length = length(unique(kappa_id)), mode = "list")
  
  for(i in 1:length(kappa_list)) {
    kappa_list[[i]] = kappa[kappa_id == i]
  }
  
  pxk = function(theta) {
    #stat = "mean"
    return(as.numeric(sapply(theta, function(theta) {
      cum_sum = exp(cumsum((alpha[x] * theta) - kappa_list[[x]]))
      par_cum_sum = cumsum(cum_sum)
      cum_sum / (1 + par_cum_sum)
    })))
  }
  
  if(length(kappa_list[[x]]) == 1) {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = sapply(seq(min, max, by = 0.1), pxk))
  } else {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = t(sapply(seq(min, max, by = 0.1), pxk)))
  }
  
  kappa_nms = names(kappa_list[[x]])
  names(pdata) = c("x", kappa_nms)
  pdata = tidyr::pivot_longer(pdata, cols = kappa_nms)
  
  pdata$inv_value = 1 - pdata$value
  pdata$information = (alpha[x] ^ 2) * pdata$value * pdata$inv_value
  
  plt = ggplot(pdata, aes(x, information, color = name)) + get_theme() + 
    xlab("theta") + ylab("p") + geom_line()
  
  return(plt)
}

tic_curve = function(mod, min, max) {
  
  if(any(grepl("^[a]", colnames(mod$post)))) {
    alpha = as.vector(summary(mod, param = "alpha", digits = 10, transpose = FALSE)["mean",])
  }
  
  nitems = length(alpha)
  
  kappa = summary(mod, param = "delta", digits = 10, transpose = FALSE)["mean",]
  kappa_names = names(kappa)
  kappa_id = as.numeric(gsub(".*[d]([^.]+)[_].*", "\\1", kappa_names))
  kappa_list = vector(length = length(unique(kappa_id)), mode = "list")
  
  for(i in 1:length(kappa_list)) {
    kappa_list[[i]] = kappa[kappa_id == i]
  }
  
  pxk = function(theta) {
    #stat = "mean"
    return(as.numeric(sapply(theta, function(theta) {
      cum_sum = exp(cumsum((alpha[x] * theta) - kappa_list[[x]]))
      par_cum_sum = cumsum(cum_sum)
      cum_sum / (1 + par_cum_sum)
    })))
  }
  
  total_info = vector(mode = "numeric", length = length(seq(min, max, by = 0.1)))
  
  for(i in 1:nitems) {
    
    x = i
    
    if(length(kappa_list[[i]]) == 1) {
      pdata = data.frame(x = seq(min, max, by = 0.1), y = sapply(seq(min, max, by = 0.1), pxk))
    } else {
      pdata = data.frame(x = seq(min, max, by = 0.1), y = t(sapply(seq(min, max, by = 0.1), pxk)))
    }
    
    kappa_nms = names(kappa_list[[i]])
    names(pdata) = c("x", kappa_nms)
    
    inv_value = 1 - pdata[,-1]
    
    information = (alpha[i] ^ 2) * pdata[,-1] * inv_value
    information = rowSums(information)
    
    total_info = total_info + information
  }
  
  dat = data.frame(x = seq(min, max, by = 0.1),
                   information = total_info)
  
  plt = ggplot(dat, aes(x, information)) + get_theme() + 
    xlab("theta") + ylab("Total information") + geom_line()
  
  return(plt)
}
