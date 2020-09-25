

#### Load data

load("./DataGeneration/geeresults_phase2_kcapprox.RData")

load("./DataGeneration/key.scenario_phase2_kcapprox.RData")

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


l.results <- all.results

rm(all.results)


# Record bad values that will be set to missing

l.results$se <- l.results$se_kc

l.results$se_kc <- NULL

dta.badse <- l.results[l.results$se > 5 & !is.na(l.results$se),]

l.results$wc <- factor(l.results$wc, levels = 0:1, labels = c("Independent", "Exchangeable"))


# Remove large coefficients and standard errors
l.results$se <- ifelse(l.results$se > 5, NA, l.results$se)

l.results$coef <- ifelse(is.na(l.results$se), NA,l.results$coef)

# Add variance
l.results$se2 <- l.results$se^2



# Add df and calculate p values and confidence intervals

l.results$dfcmp <- ifelse(l.results$clusters > l.results$sequences,
                          l.results$clusters - l.results$sequences - 1,
                          1)
names(l.results)[names(l.results) == "dffdg5"] <- "dffgd5"

y <- lapply(list("dffgd5", "dfcmp"), 
            function(x) cbind(l.results[,c("repid",
                                   "scenario",
                                   "wc",
                                   "coef", 
                                   "se",
                                   "se2")],
                              dfname = substr(x,3,nchar(x)),
                              df = l.results[,x],
                              p = 2 * (1-pt(abs((l.results$coef)/l.results$se), l.results[,x])),
                              lci = l.results$coef - qt(0.975, l.results[,x]) * l.results$se,
                              uci = l.results$coef + qt(0.975, l.results[,x]) * l.results$se))

l.results <- do.call(rbind, y)
rm(y)


l.results$trueinci <- ifelse(erx > l.results$lci & erx < l.results$uci, 1, 0)

l.results$sig <- ifelse(l.results$p > 0.05, 0, 1)

l.results$sig <- ifelse(is.na(l.results$p), NA, l.results$sig)



######################################################################################################



####Create summary statistic dataset

# Create summaries of estimates for bias
dta.stat <- aggregate(cbind(coef, se2, trueinci, sig) ~ scenario + wc + dfname, 
                      data = l.results, 
                      na.action = na.pass,
                      FUN = function(x) c(mean(x, na.rm = TRUE), 
                                          var(x, na.rm = TRUE), 
                                          sum(is.na(x))))



dta.stat <- merge(dta.stat, key.scenario, by = "scenario")


dta.stat$corrscenario <- factor(dta.stat$corrscenario, levels = 0:3, labels = c("Exchangeable", "AR r=0.6", "AR r=0.8","Reduced intervention"))


dta.stat$sequences <- factor(dta.stat$sequences)


dta.stat$`cv cluster size` <- factor(dta.stat$`cv cluster size`)

dta.stat$`mean cluster size` <- factor(dta.stat$`mean cluster size`)

dta.stat$icc <- factor(dta.stat$icc)

# Create variables that are neater for output

dta.stat$Varying_cluster_size <- factor(ifelse(dta.stat$`cv cluster size` == "0.4",1,0),
                                        levels = c(0,1),
                                        labels = c("No", "Yes"))

dta.stat$Sequences <- dta.stat$sequences
dta.stat$Mean_cluster_size <- dta.stat$`mean cluster size`
dta.stat$Correlation_structure <- dta.stat$corrscenario
dta.stat$ICC <- dta.stat$icc
dta.stat$Working_correlation <- dta.stat$wc

dta.stat$Clusters <- factor(dta.stat$clusters, levels = c(54, 48,42, 24, 18, 12, 6))


### 

dta.stat$mean <- dta.stat$coef[,1]
dta.stat$se <- sqrt(dta.stat$coef[,2])

dta.stat$n <- dta.stat$se2[,3]
dta.stat$convergep <- 100 * dta.stat$n/1000


dta.stat$mbvar <- dta.stat$se2[,1]
dta.stat$mbse <- sqrt(dta.stat$se2[,1])
dta.stat$mbsevar <- dta.stat$se2[,2]

dta.stat$stdbias <- 100 * (dta.stat$mean - erx) / dta.stat$se

dta.stat$pbias <- 100 * (dta.stat$mean - erx) / erx


dta.stat$relse <- 100 * (dta.stat$mbse / dta.stat$se - 1)


dta.stat$coverage <- 100 * dta.stat$trueinci[,1]

dta.stat$power <- 100 * dta.stat$sig[,1]

# Remove variables with multiple columns

dta.stat$coef <- NULL
dta.stat$se2 <-  NULL




###

dta.stat <- dta.stat[order(dta.stat$clusters, dta.stat$sequences, dta.stat$corrscenario, dta.stat$clustersize, dta.stat$icc, dta.stat$wc),]

###  Save datasets

save(dta.stat, file = "./crstatistics_KC_correction.RData")

save(dta.badse, file = "./crbadse_KC_correction.RData")

save(l.results, file = "./crgeeresults_KC_correction_long.RData")

### Clear workspace

rm(dta.stat, dta.badse, l.results, key.scenario)
