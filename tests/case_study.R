# ------------------------------------------------------------------------------
# Example data analysis for paper ___
# Michael Kleinsasser

# ------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(data.table)
library(haven)

# ------------------------------------------------------------------------------
# process data

# read sav
dt = read_sav('/Users/mk/Desktop/inirt_project/data/NIAAA_w1to4_and_census_data.sav')

# use data.table
dt = as.data.table(dt)

# remove anything other than w1, w2, w3, w4
all_w14 = setdiff(grep('w1|w2|w3|w4', colnames(dt)), grep('w10|w11|w12', colnames(dt)))

# subset cols
dt4 = dt[, ..all_w14]

# drop missing items
drop_f5 = setdiff(1:ncol(dt4), grep('f5|f4c', colnames(dt4)))

# subset cols
dt4 = dt[, ..drop_f5]

# complete cases analysis
dt4 = na.omit(dt4)

# final dimension
dim(dt4)

# distribution of items
apply(dt4[, -c(1:2)], 2, table)

dt4[, 
  w3g2g := ifelse(w3g2g > 3, 3, w3g2g)][, 
  w4g2g := ifelse(w4g2g > 3, 3, w4g2g)][, 
  w3g2c := ifelse(w3g2c > 3, 3, w3g2c)][, 
  w4g2c := ifelse(w4g2c > 3, 3, w4g2c)][, 
  w1f4n := ifelse(w1f4n > 3, 3, w1f4n)][, 
  w2f4n := ifelse(w2f4n > 3, 3, w2f4n)][, 
  w3f4n := ifelse(w3f4n > 3, 3, w3f4n)][, 
  w4f4n := ifelse(w4f4n > 3, 3, w4f4n)][, 
  w2f4q := ifelse(w2f4q > 3, 3, w2f4q)][, 
  w3f4q := ifelse(w3f4q > 3, 3, w3f4q)][, 
  w4f4q := ifelse(w4f4q > 3, 3, w4f4q)][, 
  w2f4s := ifelse(w2f4s > 3, 3, w2f4s)][, 
  w3f4s := ifelse(w3f4s > 3, 3, w3f4s)][, 
  w4f4s := ifelse(w4f4s > 3, 3, w4f4s)][, 
  w1f4e := ifelse(w1f4e > 4, 4, w1f4e)][, 
  w1f4f := ifelse(w1f4f > 2, 2, w1f4f)][, 
  w2f4f := ifelse(w2f4f > 4, 4, w2f4f)][, 
  w3f4f := ifelse(w3f4f > 4, 4, w3f4f)][, 
  w4f4f := ifelse(w4f4f > 4, 4, w4f4f)][, 
  w1f4j := ifelse(w1f4j > 2, 2, w1f4j)][, 
  w2f4j := ifelse(w2f4j > 2, 2, w2f4j)][, 
  w3f4j := ifelse(w3f4j > 2, 2, w3f4j)][, 
  w4f4j := ifelse(w4f4j > 2, 2, w4f4j)][, 
  w1f4b := ifelse(w1f4b > 3, 3, w1f4b)][, 
  w2f4b := ifelse(w2f4b > 3, 3, w2f4b)][, 
  w3f4b := ifelse(w3f4b > 4, 4, w3f4b)]

# how is the distribution?
apply(dt4[, -1], 2, table)

# how about row sums (naiive scores)
# hist(apply(dt4[, -1], 1, sum))

# transform data to long format

colnames(dt4)

# long_g2a = melt(dt4, id.vars = 'id', measure.vars = grep('g2a', colnames(dt4)))
# long_g2a = long_g2a[, wave := substr(variable, 2, 2)]

items_to_melt = c('g2a', 'g2g', 'g2c', 'g2d', 'g2e', 'f4n', 'f4q', 'f4s', 'f4e',
  'f4f', 'f4j', 'f41a', 'f42a', 'f4b')

setnames(dt4, 
  old = c("id", "w1g2a", "w2g2a", "w3g2a", "w4g2a", "w1g2g", "w2g2g", "w3g2g", 
  "w4g2g", "w1g2c", "w2g2c", "w3g2c", "w4g2c", "w1g2d", "w2g2d", "w3g2d", "w4g2d", 
  "w1g2e", "w2g2e", "w3g2e", "w4g2e", "w1f4n", "w2f4n", "w3f4n", "w4f4n", "w2f4q", 
  "w3f4q", "w4f4q", "w2f4s", "w3f4s", "w4f4s", "w1f4e", "w2f4e", "w3f4e", "w4f4e", 
  "w1f4f", "w2f4f", "w3f4f", "w4f4f", "w1f4j", "w2f4j", "w3f4j", "w4f4j", "w2f4aa", 
  "w3f4aa", "w4f4aa", "w1f4a", "w2f4a", "w3f4a", "w4f4a", "w1f4b", "w2f4b", "w3f4b"),
  new = c("id", "w1g2a", "w2g2a", "w3g2a", "w4g2a", "w1g2g", "w2g2g", "w3g2g", 
  "w4g2g", "w1g2c", "w2g2c", "w3g2c", "w4g2c", "w1g2d", "w2g2d", "w3g2d", "w4g2d", 
  "w1g2e", "w2g2e", "w3g2e", "w4g2e", "w1f4n", "w2f4n", "w3f4n", "w4f4n", "w2f4q", 
  "w3f4q", "w4f4q", "w2f4s", "w3f4s", "w4f4s", "w1f4e", "w2f4e", "w3f4e", "w4f4e", 
  "w1f4f", "w2f4f", "w3f4f", "w4f4f", "w1f4j", "w2f4j", "w3f4j", "w4f4j", "w2f41a", 
  "w3f41a", "w4f41a", "w1f42a", "w2f42a", "w3f42a", "w4f42a", "w1f4b", "w2f4b", "w3f4b"))

dList = 
  lapply(items_to_melt, 
    \(x, d, c_names) {
      y = melt(d, id.vars = 'id', measure.vars = grep(x, c_names))
      y[, wave := substr(variable, 2, 2)]
    }, d = dt4, c_names = colnames(dt4))

# d = merge(dList[[1]], dList[[2]], by = c('id', 'wave'))

# any(duplicated(dList[[2]][, c('id', 'wave')]))

# d = dList[[1]][, 
#   items_to_melt[1] := value][,
#   value := NULL][, 
#   variable := NULL]

dList[[1]] = dList[[1]][, 
    items_to_melt[1] := value][,
    value := NULL][, 
    variable := NULL]

dL = dList[[1]]

distinct(dList[[1]][, 1:2])

for(i in 2:length(dList)) {
  dList[[i]] = dList[[i]][, 
    items_to_melt[i] := value][,
    value := NULL][, 
    variable := NULL]

  dL = merge(dL, dList[[i]], by = c('id', 'wave'), all.x = TRUE)
  cat(dim(dL)[1])
  cat('\n')
}

dL
colnames(dL)

dL = as.data.frame(dL)

# ------------------------------------------------------------------------------
# fit model
model = 'theta1 = g2a + g2g + g2c + g2d + g2e
         theta2 = f4n + f4q + f4s + f4e + f4f
         theta2 = f4j + f41a + f42a + f4b
         theta1 ~ (1 + wave | id)
         theta2 ~ (1 + wave | id)
         theta3 ~ (1 + wave | id)'

fit = theta2::theta2(data = dL, model = model, itype = "2pl", 
  exploratory = TRUE, method = "hmc", iter = 500, chains = 1)

data = dL; model = model; itype = "2pl"; 
exploratory = TRUE; method = "hmc"; iter = 500; chains = 1
library(devtools)
load_all()
# ------------------------------------------------------------------------------
# summarize results


# ------------------------------------------------------------------------------
