

#### Load data

load("./DataGeneration/geeresults_phase1.RData")

load("./DataGeneration/key.scenario_phase1.RData")

all.results <- as.data.frame(all.results)


#### Neaten up scenarios

key.scenario <- as.data.frame(key.scenario)

# Combine some of the scenarios for graphs
key.scenario$clustersize <- interaction(key.scenario$'mean cluster size', key.scenario$'cv cluster size')

levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "24.0"] <- "24 cv=0" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "24.0.4"] <- "24 cv=0.4" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "60.0"] <- "60 cv=0" 
levels(key.scenario$clustersize)[levels(key.scenario$clustersize) == "60.0.4"] <- "60 cv=0.4" 

all.results <- merge(all.results, key.scenario, by = "scenario")

# Rename mk mw
all.results$se_mw <- all.results$se_mk
all.results$dfw_mw <- all.results$dfw_mk
all.results$se_mk <- NULL
all.results$dfw_mk <- NULL

# Add Normal distribution
all.results$dfinf <- Inf

#### Reshape long

l.results <- reshape(direction = "long",
                     data = all.results,
                     idvar = c("scenario","repid","wc","samplesize"),
                     times = c("lz","fg","mw","kc","md","mbn"),
                     varying = list(c("se_lz","se_fg","se_mw","se_kc","se_md","se_mbn"),
                                    c("dfw_lz","dfw_fg","dfw_mw","dfw_kc","dfw_md","dfw_mbn")),
                     timevar = "sa",
                     v.names = c("se","dfw"))



#Add factor labels

l.results$sa <- factor(l.results$sa, levels = c("lz", "kc", "fg", "md", "mbn", "mw"))

l.results$wc <- factor(l.results$wc, levels = 0:1, labels = c("Independent", "Exchangeable"))

# Record bad values that will be set to missing

dta.badse <- l.results[l.results$se > 5 & !is.na(l.results$se),]

# remove large coefficients and standard errors

l.results$se <- ifelse(l.results$se > 5, NA, l.results$se)

l.results$coef <- ifelse(is.na(l.results$se), NA,l.results$coef)

# Wang DF gives crazy high values and some negative values
l.results[,"dfw"] <- ifelse(l.results[,"dfw"] < 1 | l.results[,"dfw"] > 1000, NA, l.results[,"dfw"])

# for (i in c("dfw", "dffgd5", "dfcmp", "dfcpmp", "dfcpmcp")) {
#   l.results[,i] <- ifelse(l.results[,i] < 1, NA, l.results[,i])
#   l.results[,i] <- ifelse(l.results[,i] > 10000, NA, l.results[,i])
# }

### Add vars

l.results$se2 <- l.results$se^2

l.results <- lapply(list("dfw", "dffgd5", "dfcmp", "dfcpmp", "dfcpmcp", "dfinf"), 
                    function(x) cbind(l.results[,c("repid",
                                                   "scenario",
                                                   "wc",
                                                   "sa", 
                                                   "coef", 
                                                   "se",
                                                   "se2")],
                                      dfname = substr(x,3,nchar(x)),
                                      df = l.results[,x],
                                      p = 2 * (1-pt(abs((l.results$coef)/l.results$se), l.results[,x])),
                                      lci = l.results$coef - qt(0.975, l.results[,x]) * l.results$se,
                                      uci = l.results$coef + qt(0.975, l.results[,x]) * l.results$se))
l.results <- do.call(rbind, l.results)

#Sort out dfname as a factor
l.results$dfname <- factor(l.results$dfname, levels = c("inf", "fgd5", "w", "cmp", "cpmp", "cpmcp"))

#Create type 1 errors  

l.results$trueinci <- ifelse(erx>l.results$lci & erx<l.results$uci, 1, 0)

l.results$sig <- ifelse(l.results$p>0.05,0,1)

l.results$sig <- ifelse(is.na(l.results$p),NA,l.results$sig)

#### Create summary statistic dataset

# Create summaries of estimates for bias
dta.stat <- aggregate(cbind(coef, se2, trueinci, sig, df) ~ scenario + wc + sa + dfname , 
                      data = l.results, 
                      na.action = na.pass,
                      FUN = function(x) c(mean(x, na.rm = TRUE), 
                                          var(x, na.rm = TRUE), 
                                          sum(is.na(x))))

dta.stat <- merge(dta.stat, key.scenario, by = "scenario")

dta.stat$icc <- factor(dta.stat$icc)

dta.stat$corrscenario <- factor(dta.stat$corrscenario, levels = 0:3, labels = c("Exchangeable", "AR r=0.6", "AR r=0.8","Reduced intervention"))

dta.stat$sequences <- factor(dta.stat$sequences)

dta.stat$`cv cluster size` <- factor(dta.stat$`cv cluster size`)

dta.stat$`mean cluster size` <- factor(dta.stat$`mean cluster size`)

dta.stat$mdl <- interaction(dta.stat$sa, dta.stat$dfname)

# Create variables that are neater for output

dta.stat$SA <- dta.stat$sa
levels(dta.stat$SA)[levels(dta.stat$SA) == "lz"] <- "Uncorr" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "kc"] <- "KC" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "fg"] <- "FG" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "md"] <- "MD" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "mbn"] <- "MBN" 
levels(dta.stat$SA)[levels(dta.stat$SA) == "mw"] <- "MW" 

dta.stat$DF <- dta.stat$dfname
levels(dta.stat$DF)[levels(dta.stat$DF) == "inf"] <- "Inf"
levels(dta.stat$DF)[levels(dta.stat$DF) == "w"] <- "PW"
levels(dta.stat$DF)[levels(dta.stat$DF) == "fgd5"] <- "FG"
levels(dta.stat$DF)[levels(dta.stat$DF) == "cmp"] <- "C-P"
levels(dta.stat$DF)[levels(dta.stat$DF) == "cpmp"] <- "CP-P"
levels(dta.stat$DF)[levels(dta.stat$DF) == "cpmcp"] <- "CP-C-P"


dta.stat$Varying_cluster_size <- factor(ifelse(dta.stat$`cv cluster size` == "0.4",1,0),
                                        levels = c(0,1),
                                        labels = c("No", "Yes"))

dta.stat$Sequences <- dta.stat$sequences
dta.stat$Mean_cluster_size <- dta.stat$`mean cluster size`
dta.stat$Correlation_structure <- dta.stat$corrscenario
dta.stat$Working_correlation <- dta.stat$wc
dta.stat$ICC <- dta.stat$icc

dta.stat$Model <- interaction(dta.stat$Working_correlation, dta.stat$SA, sep = "_")

### create statistics

dta.stat$mean <- dta.stat$coef[,1]
dta.stat$se <- sqrt(dta.stat$coef[,2])

dta.stat$n <- dta.stat$sig[,3]
dta.stat$convergep <- 100 * dta.stat$n/1000

dta.stat$coverage <- 100 * dta.stat$trueinci[,1]
dta.stat$Coverage <- dta.stat$coverage - 95

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


# Save datasets

save(dta.stat, file = "./crstatistics_phase_1.RData")

save(dta.badse, file = "./crbadse_phase_1.RData")

save(l.results, file = "./crgeeresults_phase1_long.RData")


# Clear up workspace
rm(all.results, dta.badse, dta.stat, key.scenario, l.results)
