

#### Load data

load("./DataGeneration/geeresults_phase2.RData")

load("./DataGeneration/key.scenario_phase2.RData")

all.results <- as.data.frame(all.results)

#### Neaten up scenarios

key.scenario <- as.data.frame(key.scenario)

# Combine some of the scenarios for graphs
key.scenario$clustersize <- interaction(key.scenario$'mean cluster size', key.scenario$'cv cluster size', drop = TRUE)

levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "24.0"] <- "24 cv=0" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "24.0.4"] <- "24 cv=0.4" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "60.0"] <- "60 cv=0" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "60.0.4"] <- "60 cv=0.4" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "300.0"] <- "300 cv=0" 

all.results <- merge(all.results, key.scenario, by = "scenario")

#### Reshape long

# Datasets too large, splitting into two

all.results <- split(all.results, f = all.results$corrscenario)

l.results <- lapply(all.results,
                    FUN = function(y) reshape(direction = "long",
                                              data = y,
                                              idvar = c("scenario","repid","wc","samplesize"),
                                              times = c("lz","fg","kc"),
                                              varying = c("se_lz","se_fg","se_kc"),
                                              timevar = "sa",
                                              v.names = "se"))

rm(all.results)


# Record bad values that will be set to missing

dta.badse <- lapply(l.results, FUN = function(y) {y[y$se > 5 & !is.na(y$se),]})



l.results <- lapply(l.results, FUN = function(y) {
  
  #Add factor labels
  y$sa <- factor(y$sa, levels = c("lz", "kc", "fg"))
  y$wc <- factor(y$wc, levels = 0:1, labels = c("Independent", "Exchangeable"))
  
  # Remove large coefficients and standard errors
  y$se <- ifelse(y$se > 5, NA, y$se)
  
  y$coef <- ifelse(is.na(y$se), NA,y$coef)
  
  # Add variance
  y$se2 <- y$se^2
  
  # Scenarios with a small number of clusters end up with negative degrees of freedom, Reset these to DF=1
  y$dfcmp <- ifelse(y$dfcmp < 1, 1, y$dfcmp)
  
  
  y$dfdiff <- y$dffgd5 - y$dfcmp
  
  y$dfinf <- Inf
  
  # Create p values and confidence intervals with each df and reshape to long
  y <- lapply(list("dffgd5", "dfcmp","dfinf"), 
              function(x) cbind(y[,c("repid",
                                     "scenario",
                                     "wc",
                                     "sa", 
                                     "coef", 
                                     "se",
                                     "se2",
                                     "dfdiff")],
                                dfname = substr(x,3,nchar(x)),
                                df = y[,x],
                                p = 2 * (1-pt(abs((y$coef)/y$se), y[,x])),
                                lci = y$coef - qt(0.975, y[,x]) * y$se,
                                uci = y$coef + qt(0.975, y[,x]) * y$se))
  
  y <- do.call(rbind, y)
  
  y <- y[(y$sa %in% c("kc","fg") & y$dfname %in% c("fgd5","cmp")) | y$sa == "lz" ,]
  
  #Create type 1 errors  
  
  y$trueinci <- ifelse(erx>y$lci & erx<y$uci, 1, 0)
  
  y$sig <- ifelse(y$p>0.05,0,1)
  
  y$sig <- ifelse(is.na(y$p),NA,y$sig)
  
  return(y)
  
})

#### Create summary statistic dataset

# Create summaries of estimates for bias
dta.stat <- lapply(l.results, 
                   FUN = function(y) {
                     aggregate(cbind(coef, se2, trueinci, sig, df, dfdiff) ~ scenario + wc + sa + dfname , 
                               data = y, 
                               na.action = na.pass,
                               FUN = function(x) c(mean(x, na.rm = TRUE), 
                                                   var(x, na.rm = TRUE), 
                                                   sum(is.na(x))))
                   })

dta.stat <- do.call(rbind, dta.stat)

dta.stat <- merge(dta.stat, key.scenario, by = "scenario")


dta.stat$corrscenario <- factor(dta.stat$corrscenario, levels = 0:3, labels = c("Exchangeable", "AR r=0.6", "AR r=0.8","Reduced intervention"))


dta.stat$sequences <- factor(dta.stat$sequences)


dta.stat$`cv cluster size` <- factor(dta.stat$`cv cluster size`)

dta.stat$`mean cluster size` <- factor(dta.stat$`mean cluster size`)

dta.stat$mdl <- interaction(dta.stat$sa, dta.stat$dfname,  drop = TRUE, lex.order = TRUE)

dta.stat$icc <- factor(dta.stat$icc)

# Create variables that are neater for output

dta.stat$SA <- dta.stat$sa
levels(dta.stat$SA)[levels(dta.stat$SA) == "lz"] <- "Uncorr" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "kc"] <- "KC" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "fg"] <- "FG" 

dta.stat$DF <- dta.stat$dfname
levels(dta.stat$DF)[levels(dta.stat$DF) == "fgd5"] <- "FG"
levels(dta.stat$DF)[levels(dta.stat$DF) == "cmp"] <- "C-P"
levels(dta.stat$DF)[levels(dta.stat$DF) == "inf"] <- "Inf"

dta.stat$Varying_cluster_size <- factor(ifelse(dta.stat$`cv cluster size` == "0.4",1,0),
                                        levels = c(0,1),
                                        labels = c("No", "Yes"))

dta.stat$Sequences <- dta.stat$sequences
dta.stat$Mean_cluster_size <- dta.stat$`mean cluster size`
dta.stat$Correlation_structure <- dta.stat$corrscenario
dta.stat$Working_correlation <- dta.stat$wc
dta.stat$ICC <- dta.stat$icc

dta.stat$Clusters <- factor(dta.stat$clusters, levels = c(54, 48,42, 24, 18, 12, 6))

dta.stat$Model <- interaction(dta.stat$SA, dta.stat$DF, sep = "_")

### 

dta.stat$mean <- dta.stat$coef[,1]
dta.stat$se <- sqrt(dta.stat$coef[,2])

dta.stat$n <- dta.stat$sig[,3]
dta.stat$convergep <- 100 * dta.stat$n/1000

dta.stat$coverage <- 100 * dta.stat$trueinci[,1]
dta.stat$Coverage <- dta.stat$coverage - 95
dta.stat$Coverage_category <- cut(dta.stat$coverage, quantile(dta.stat$coverage, probs = seq(0,1,0.2)))

dta.stat$power <- 100 * dta.stat$sig[,1]

dta.stat$mbvar <- dta.stat$se2[,1]
dta.stat$mbse <- sqrt(dta.stat$se2[,1])
dta.stat$mbsevar <- dta.stat$se2[,2]

dta.stat$lb <- dta.stat$mean - dta.stat$se / 2
dta.stat$ub <- dta.stat$mean + dta.stat$se /2

dta.stat$stdbias <- 100 * (dta.stat$mean - erx) / dta.stat$se

dta.stat$pbias <- 100 * (dta.stat$mean - erx) / erx


dta.stat$relse <- 100 * (dta.stat$mbse / dta.stat$se - 1)

dta.stat$meandf <- dta.stat$df[,1]
dta.stat$vardf <- dta.stat$df[,2]

dta.stat$lbdf <- dta.stat$meandf - sqrt(dta.stat$vardf)/2

dta.stat$ubdf <- dta.stat$meandf + sqrt(dta.stat$vardf)/2

dta.stat <- dta.stat[order(dta.stat$clusters, dta.stat$sequences, dta.stat$corrscenario, dta.stat$clustersize, dta.stat$wc, dta.stat$icc, dta.stat$sa, dta.stat$dfname),]

###  Save datasets

save(dta.stat, file = "./crstatistics_phase_2.RData")

save(dta.badse, file = "./crbadse_phase_2.RData")

save(l.results, file = "./crgeeresults_phase2_long.RData")

### Clear workspace

rm(dta.stat, dta.badse, l.results, key.scenario)
