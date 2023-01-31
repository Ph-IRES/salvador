#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(magrittr)
# library(janitor)
# install.packages("fitdistrplus")

# install.packages("rlang")
# install.packages("emmeans")
# library(devtools)
# devtools::install_github("singmann/afex@master")
# devtools::install_github("eclarke/ggbeeswarm")
# install.packages("multcomp")
# install.packages("multcompView")
# install.packages("performance")
# install.packages("optimx")
# install.packages("effects")

library(fitdistrplus)
# library(emmeans)
# library(multcomp)
# library(multcompView)
# library(ggeffects)


# library(rlang)
# library(afex)
# library(ggbeeswarm)
# library(performance)
# library(optimx)
# library(effects)
# library(prediction)

# NOTE: after loading these packages, you may find that tidyverse commands are affected 
#       the solution is to add the appropriate package name before commands that break code
#       such as `dplyr::select` if `select` doesn't work correctly anymore
#       this happens when multiple packages have the same command names. 
#       The last package loaded takes precidence, and your tidyverse commands break.
#       you could load tidyverse last, but invariably, you will load a package after tidyverse
#       so it's impossible to avoid this
#### FUNCTIONS ####
# Functions are commands that you define.  Here we define a function that makes several plots.
# running the code below will only save the function into memory.  Later we will call the function to make the plots
# visualize statistical distributions (see fitdistrplus: An R Package for Fitting Distributions, 2020)
vis_dists <- function(data,
                      response_var){
  
  data_plot <-
    data %>%
    drop_na(!!response_var) %>%
    # log() %>%
    pull(!!response_var)
  
  plotdist(data_plot)
  descdist(data_plot,
           boot=1000)
  
  fw <-
    try(fitdist(data_plot,
                "weibull"))
  fp <-
    try(fitdist(data_plot,
                "pois",
                method="mme"))
  fnb <-
    try(fitdist(data_plot,
                "nbinom",
                method="mme"))
  fg <-
    try(fitdist(data_plot,
                "gamma",
                method="mme"))
  fl <-
    try(fitdist(data_plot,
                "logis"))
  
  fln <-
    try(fitdist(data_plot,
                "lnorm"))
  
  fn <-
    try(fitdist(data_plot,
                "norm"))
  fge <-
    try(fitdist(data_plot,
                "geom",
                method="mme"))
  
  # remove distibutions that choke on neg or zero values if there are neg or zero values
  data_plot_lessorequal_zero = data_plot <= 0
  
  fits_to_plot <- 
    if(TRUE %in% data_plot_lessorequal_zero){
      list(fp, fl, fn, fge)
    } else {
      list(fw, fp, fnb, fg, fl, fln, fn, fge)
    }
  
  plot_legend <- 
    if(TRUE %in% data_plot_lessorequal_zero){
      c("Poisson", "Logis", "Normal", "Geom")
    } else {
      c("Weibull","Poisson","NegBinom","Gamma", "Logis","lognormal", "Normal", "Geom")
    }
  denscomp(fge)
  par(mfrow = c(2, 2))
  try(denscomp(fits_to_plot, legendtext = plot_legend))
  try(qqcomp(fits_to_plot, legendtext = plot_legend))
  try(cdfcomp(fits_to_plot, legendtext = plot_legend))
  try(ppcomp(fits_to_plot, legendtext = plot_legend))
  par(mfrow = c(1, 1))
}


