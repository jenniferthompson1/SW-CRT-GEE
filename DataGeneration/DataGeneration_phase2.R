
library("Matrix")
library("MASS")
library("matrixcalc")
library("parallel")
library("bindata")
library("geepack")
library("saws")


# Parallel computing set up ---------------------------------

# Set up parallel info
cl <- makeCluster(32)

clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(matrixcalc))
clusterEvalQ(cl, library(MASS))
clusterEvalQ(cl, library(geepack))
clusterEvalQ(cl, library(saws))

clusterSetRNGStream(cl=cl, 285403)

# Parameter set up ---------------------------------

rhoar06 <- 0.6
rhoar08 <- 0.8
rhorx <- 0.5

etime <- log(1.05)
erx <- log(1.3)

iscenario <- 0

reps <- 1000

v.basep <-  0.3
v.sequences <- c(3,6)
v.corrscenario <-  c(0, 1, 2, 3)
v.rhoc <-  c(0.01, 0.05, 0.1)
v.clusters <-  c(6, 12, 18, 24, 42, 48, 54)
v.clustersize <-  c(24, 60, 300)
v.cv <- c(0, 0.4)


# Functions ---------------------------------
## Correction functions =================================

source("./functions.R")

## Data Generation =================================

### Sampling #############################

createdata <- function(x,cps) {
  
  probit <- qnorm(rep(x$mean,x$n))
  
  # Check eigen values are positive for correlation matrix. Random variation and rounding mean some of the big matrices have negative eigen values
  cholStatus <- try(u <- chol(x$corr), silent = FALSE)
  cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  
  newMat <- x$corr
  
  iter <- 0
  
  while (cholError) {
    
    iter <- iter + 1
    cat("iteration ", iter, "\n")
    
    # replace -ve eigen values with small +ve number
    newEig <- eigen(newMat)
    newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
    
    # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
    # eig vectors
    newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
    
    # normalize modified matrix eqn 6 from Brissette et al 2007
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    
    # try chol again
    cholStatus <- try(u <- chol(newMat), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  }
  
  if (cps == 1) {
    norm <- t(matrix(mvrnorm(cps, mu = rep(0, nrow(x$corr)), Sigma = newMat)))
    
  } else {
    norm <- mvrnorm(cps, mu = rep(0, nrow(x$corr)), Sigma = newMat)
  }
  
  y <- matrix(0, cps, nrow(x$corr))
  
  for (i in 1:nrow(x$corr)) y[, i] <- cut(norm[, i] , cbind(-Inf,probit[i], Inf), labels = FALSE)
  y <- 2 - y
  
  colnames(y) <- paste0(paste("t",rep(1:6, x$n),sep="."),rep(1:x$n,each = 6))
  
  y <- cbind(id = seq_len(nrow(y)),
             sequence = x$sequence, 
             y)
  
  y <- cbind(cluster = (y[,"sequence"] - 1) * cps + y[,"id"], 
             y)
  
  y <- reshape(direction = "long",
               data = as.data.frame(y),
               varying = paste0(paste("t",rep(1:6, x$n),sep="."),rep(1:x$n,each = 6)),
               times = paste0(rep(1:6, x$n),rep(1:x$n,each = 6)),
               timevar = "t",
               v.names = "y",
               idvar = "cluster")
  
  y <- cbind(y, 
             period = as.numeric(substr(y[,"t"],1,1)) - 1, 
             subid  = as.numeric(as.numeric(substr(y[,"t"],2,length(y[,"t"])))))
  
  
  
  return(as.matrix(y[,c("y","sequence","cluster","period","subid")]))
  
  
}

####################################################################

### Setting up correlations and analyse #############################

####################################################################

replicatedata <- function(A) {
  
  
  if (cv == 0) {
    n <- rep(clustersize / 6, 6)
  } else {
    s2 <- (cv*clustersize)^2 / 6
    
    cellsize <-  clustersize / 6
    
    actualclustersize <- as.data.frame(cbind(rep(1:6, each = cps * 6),
                                             rep(1:clusters, each = 6),
                                             rep(0:5, clusters),
                                             2 + rnbinom(clusters*6, 
                                                         (cellsize-2)^2/(s2 - (cellsize-2)), 
                                                         (cellsize-2)/s2)))
    names(actualclustersize) <- c("sequence","cluster","period","n")
    
    #set number of observations in each cluster-period
    n <- aggregate(n ~ sequence, data = actualclustersize, FUN = max)[,2]
    
  }
  
  
  # create list of means and correlations for each sequence
  m.list <- list()
  for (i in 1:6) {
    m.list[[i]] <- list(sequence = i,
                        n = n[i],
                        mean = mean[i,],
                        corr = (as.matrix(matrix(rep(1,n[i]*n[i]),ncol=n[i]) %x% moffdiag[[i]]  - 
                                            bdiag(lapply(seq_len(n[i]), function(X) moffdiag[[i]])) + 
                                            bdiag(lapply(seq_len(n[i]), function(X) mdiag[[i]])))))
    
  }
  
  #create dataset
  data <-tryCatch(do.call(rbind,lapply(m.list, 
                                       createdata, cps = cps)),
                  error=function(e) e,
                  warning=function(w) w)
  
  
  if (is(data,"warning") | is(data,"error")) return(list(list(iscenario, .Random.seed), c(iscenario, rep(NA, 8))))
  
  # change cluster size if cv=0.4
  if (cv ==0.4){  
    
    data <- merge(data, actualclustersize, by = c("sequence","cluster","period"))
    
    data <- data[data$subid<=data$n,]
  }
  
  # merge in rx and scenario
  data <- as.matrix(merge(data, design.long, by=c("sequence", "period")))
  
  #Sort data
  data <- data[order(data[,"cluster"], data[,"period"]),]
  
  data1 <- as.data.frame(data)
  
  ##################################
  # Run binary analysis models
  ##################################
  
  results <- c()
  
  formula1 <- formula(y ~  rx + as.factor(period2))
  
  sampsi <- nrow(data1)
  
  for (wc in c("independence","exchangeable")) {
    
    #Run gee model
    geemdl <-  tryCatch(geeglm(formula1,
                               data = data1,
                               id = cluster,
                               family = binomial,
                               corstr = wc),
                        error=function(e) e,
                        warning=function(w) w)
    
    #create code of working correlation matrix
    if (wc == "independence") iwc <- 0 else iwc <- 1
    
    #If gee doesn't work just append a row of missing to results
    if(is(geemdl,"warning") | is(geemdl,"error") )  {
      
      results <- rbind(results,
                       c(iscenario, iwc, sampsi, rep(NA, 7)))
      
      #otherwise, calculate each standard error
    } else {
      
      #Create a row with the scenario etc and intervention effect
      row <- c(iscenario, iwc, sampsi, geemdl$coefficients[2])
      
      
      #Run each correction and append standard error and Wang DF to row
      args <- list( formula = formula1,
                    data = quote(data1),
                    id = "cluster",
                    family = binomial,
                    corstr = wc)
      
      for (sa in c("lz","fg","kc","mk")) {
        
          var <- tryCatch(do.call(paste0("GEEGLM.var.",sa), args),
                          error=function(e) e,
                          warning=function(w) w)
          
          # calculation thows error or doesn't converge append a missing value
          if (is(var,"warning") | is(var,"error")) {
            row <- c(row, NA)
            #Otherwise append the standard error
          } else {
            row <- c(row, sqrt(var$cov.beta[2]))
          }
      }
      
      #Run each DF and append to each row
      
      #fg
      geemdl$scale <- geemdl$geese$gamma
      
      if (wc == "independence") {
        geemdl$working.correlation <- cormax.ind(max(geemdl$geese$clusz))
      } else  geemdl$working.correlation <- cormax.exch(max(geemdl$geese$clusz), geemdl$geese$alpha)
      
      geemdlu <- tryCatch(geeUOmega(geemdl),
                          error=function(e) e,
                          warning=function(w) w)
      
      saws <- tryCatch(saws(geemdlu, method = "d5"),
                       error=function(e) e,
                       warning=function(w) w)
      
      if(is(saws,"warning") | is(saws,"error") )  sdf <- NA else sdf <- saws$df[2]
      
      #Simple:
      ## Clusters minus parameters
      
      row <- c(row,
               sdf,
               clusters - (sequences + 1))
      
      #Append the row with all the standard errors to the results
      results <- rbind(results, row)
      
    }
  }
  
  colnames(results) <- c("scenario", "wc", "samplesize", "coef",
                         "se_lz", "se_fg", "se_kc", "se_mk", 
                         "dffgd5","dfcmp")
  
  return(list(list(iscenario, .Random.seed), as.matrix(results)))
  
}



# Main loop ---------------------------------

ptm <- proc.time()

for (sequences in v.sequences) {
  
  # Create design pattern matrix for 6 sequence design
  if (sequences == 6) {
    
    design <- matrix(1,6,6)
    design[lower.tri(design)]<-0
    
  } else {
    
    design <- matrix(c(rep(1,12), 
                       rep(c(0,0,rep(1,4)),2), 
                       rep(c(rep(0,4),rep(1,2)),2)),
                     6, byrow = TRUE)
    
  }
  
  design.long <- as.vector(t(design))
  design.long <- as.data.frame(cbind(rep(1:6,each = 6),
                                     rep(0:5,6),
                                     design.long))
  
  names(design.long) <- c("sequence","period","rx")
  
  # Create periods for analysis
  if (sequences == 6){
    design.long$period2 <- design.long$period
  } else {
    design.long[design.long$period==0,"period2"] <- 0
    design.long[design.long$period==1,"period2"] <- 0
    design.long[design.long$period==2,"period2"] <- 1
    design.long[design.long$period==3,"period2"] <- 1
    design.long[design.long$period==4,"period2"] <- 2
    design.long[design.long$period==5,"period2"] <- 2
    
  }
  
  for (basep in v.basep) {
    
    # Create matrix of mean log odds in each cell
    meanlogodds <-  log(basep/(1-basep)) + matrix(0:5,6,6, byrow=TRUE)  * etime + erx * design
    
    mean <- exp(meanlogodds) / (1 + exp(meanlogodds))
    
    # create indicator of scenario: 
    # 0 exchangeable
    # 1 autocorrelated 0.6, 
    # 2 autocorrelated 0.8, 
    # 3 intervention drop
    
    for (corrscenario in v.corrscenario) {
      
      if (corrscenario == 0) rhoar <- 1 else if (corrscenario == 1) rhoar <- rhoar06 else rhoar <- rhoar08
      
      #add in each within-period icc value
      for (rhoc in v.rhoc) {
        
        mdiag <- moffdiag <- lat.mdiag <- lat.moffdiag <- list()
        
        #create factors that change with the intervention effect
        for (i in 1:6) {
          
          if (corrscenario == 3) mrx1 <- rhorx ^ abs(cbind(rep(1,6),design[i,]) %*% rbind(design[i,], rep(-1,6))) else mrx1 <- matrix(1,6,6)
          
          lat.mdiag <-  ((matrix(rhoc,6,6) + diag( 1 - rhoc,6)) * 
                           rhoar^abs(row(diag(6))-col(diag(6))) * 
                           mrx1 )
          
          lat.moffdiag <- (rhoc * 
                             rhoar^abs(row(diag(6))-col(diag(6))) * 
                             mrx1) 
          
          latmat <- (as.matrix(matrix(rep(1,2*2),ncol=2) %x% lat.moffdiag  - 
                                 bdiag(lapply(seq_len(2), function(X) lat.moffdiag)) + 
                                 bdiag(lapply(seq_len(2), function(X) lat.mdiag))))
          
          commonprob <- bincorr2commonprob(rep(mean[i,],2), latmat)
          sigma <- commonprob2sigma(commonprob)
          
          mdiag[[i]] <- sigma[1:6,1:6]
          
          moffdiag[[i]] <- sigma[1:6, 7:12]
          
        }
        
        #set number of cluster
        for (clusters in v.clusters) {
          
          cps <- clusters/6
          
          for (clustersize in v.clustersize) {
            
            for (cv in v.cv) {
              
              if (!(clustersize == 300 & (cv == 0.4 | corrscenario %in% c(1,2,3) | 
                    clusters %in% c(6, 18, 48, 54, 60) ))) {
                
                #Export objects to all cores
                clusterExport(cl,c("cv","clustersize","cps","clusters","sequences", "corrscenario", "rhoc",
                                   "iscenario", "design.long", "mdiag", "moffdiag", "mean",
                                   "createdata", "replicatedata",
                                   "GEEGLM.var.lz", "GEEGLM.var.fg", "GEEGLM.var.kc","GEEGLM.var.mk",
                                   "cluster.size", "cormax.ar1", "cormax.exch", "cormax.ind", "mat.sqrt", "mat.prod", "mat.sqrt.inv"))
                #create data
                results <-parSapply(cl, 1:reps, replicatedata, simplify = FALSE)

                randomseeds <- lapply(results, function(l) l[[1]])
                results <- lapply(results, function(l) l[[2]])
                
                randomseeds <-lapply(1:length(randomseeds), function(id) c("repid" = id, randomseeds[[id]]))
                
                results <-lapply(1:length(results), function(id) cbind("repid" = id, results[[id]]))
                
                if (iscenario == 0) {
                  
                  key.scenario <- c(iscenario, sequences, clusters, clustersize, cv, rhoc, corrscenario)
                  all.results <- results
                  all.randomseeds <- randomseeds
                  
                  
                } else {
                  
                  key.scenario <- rbind(key.scenario, c(iscenario, sequences, clusters, clustersize, cv, rhoc, corrscenario))
                  all.results <- c(all.results, results)
                  all.randomseeds <- c(all.randomseeds, randomseeds)
                  
                }
                
                print(c(as.character(Sys.time()) ,iscenario, sequences, clusters, clustersize, cv, rhoc, corrscenario))
                iscenario <- iscenario + 1
              }
            }
          }
        }
      }
    }
  }
}

stopCluster(cl)

all.results <- do.call(rbind, all.results)

# Send to file ---------------------------------

save(all.results, file = "./DataGeneration/geeresults_phase2.RData")

colnames(key.scenario) <- c("scenario", "sequences", "clusters", "mean cluster size", "cv cluster size", "icc", "corrscenario")

save(key.scenario, file = "./DataGeneration/key.scenario_phase2.RData")


save(all.randomseeds, file = "./DataGeneration/allrandomseeds_phase2.RData")
proc.time() - ptm


