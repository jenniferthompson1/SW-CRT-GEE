

####
## Code uploaded to Github
####


#Load libraries
library(ggplot2)

# Create global objects

erx <- log(1.3)


#########################
# Generate data:        #
#########################


# Simulate data and analyse
# source("./DataGeneration/DataGeneration_phase1.R")
# source("./DataGeneration/DataGeneration_phase2.R")


#########################
# Phase 1:              #
#########################

source("./Analysis_crstatistics_phase1.R")


#########################
# Phase 2:              #
#########################

source("./Analysis_crstatistics_phase2.R")


#########################
# Post hoc/ KC approx:  #
#########################

# source("./DataGeneration/DataGeneration_phase2_KCapprox.R")

source("./Analysis_crstatistics_kcapprox.R")

rmarkdown::render("./Analysis_kcapprox.Rmd")



#########################
# Paper:                #
#########################

rmarkdown::render("./Paper.Rmd")
rmarkdown::render("./Appendix.Rmd")

