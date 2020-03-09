library(tidyverse)

# provide initial code on how the example dataset was generated
tmp <- tempfile()
download.file("https://github.com/kstreet13/slingshot/raw/d6df7c6f4232c2ee8819d93e644754794b738c81/data/slingshotExample.RData", tmp)
load(tmp)

# save data object with usethis::use_data
slingshotExample <- list(rd = rd, cl = cl)
usethis::use_data(slingshotExample, overwrite = TRUE)
