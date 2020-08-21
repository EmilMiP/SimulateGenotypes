library(MASS)
library(stringr)
library(doSNOW)
library(foreach)
library(parallel)
library(progress)
library(data.table)
library(tidyverse)
library(dplyr)

#normalising genotypes:
normalize_geno = function(geno_dat, MAF) {
  geno_dat_scaled = sweep(geno_dat, 1, 2*MAF, FUN = "-")
  return(sweep(geno_dat_scaled, 1, sqrt(2*MAF*(1 - MAF)), FUN = "/"))
}

#function to generate binary phenotype for parents and offspring
#and unscaled genotype for the offspring.
generateOffspring = function(parents.halved, nsib, M, MAF, environ, lia.T, lia.beta) {
  offspring.geno.unscaled = replicate(1 + nsib, rowSums(round(parents.halved + rnorm(2*M, sd = 0.02))))
  offspring.geno.scaled = normalize_geno(offspring.geno.unscaled, MAF = MAF)
  offspring.geno.lia.first = lia.beta %*% offspring.geno.scaled
  offspring.lia.first = offspring.geno.lia.first  + environ[-(2:3)]
  return(list(
    "offspring.cc" = as.numeric( offspring.lia.first > lia.T ) ,
    "offspring.geno.unscaled" = offspring.geno.unscaled[,1], # I do not need the scaled offsprings genotype after this point.
    "offspring.geno.lia" = offspring.geno.lia.first,
    "offspring.lia" = offspring.lia.first
  ))
}

generateOffspring_geno = function(parents.halved, nsib, M, MAF) {
  return("offspring.geno.unscaled" = replicate(1 + nsib, rowSums(round(parents.halved + rnorm(2*M, sd = 0.02)))) # I do not need the scaled offsprings genotype after this point.
  )
}

#function to generate binary phenotype for parents and offspring
#and unscaled genotype for the offspring.
generateOffspring_pheno = function(offspring.geno.unscaled, MAF, environ, lia.T, lia.beta) {
  offspring.geno.scaled = normalize_geno(offspring.geno.unscaled, MAF = MAF)
  offspring.geno.lia.first = lia.beta %*% offspring.geno.scaled
  offspring.lia.first = offspring.geno.lia.first  + environ[-(2:3)]
  return(list(
    "offspring.cc" = as.numeric( offspring.lia.first > lia.T ) ,
    "offspring.geno.lia" = offspring.geno.lia.first,
    "offspring.lia" = offspring.lia.first
  ))
}


generateChildren = function(
  NoChildren = 10000,
  M = 10000,
  MAF,
  C = M*0.1,
  h2_1 = 0.5,
  h2_2 = 0.5,
  gen_cor = 0.5,
  K = 0.05,
  nsib = 2,
  nthreads = 6,
  allele.mat = matrix(c("A","A","A","G", "G","G"), nrow = 3, ncol = 2, byrow = T),
  out = "D:/Work/Project1/simulatedData/ph"
) {
  #liability threshold:
  lia.T = qnorm(1 - K) 
  #liability scale effect sizes:
  #construct cov_mat:
  cov_mat = diag(c(h2_1, h2_2))
  cov_mat[1,2] <- cov_mat[2,1] <- gen_cor * sqrt(h2_1 * h2_2)
  #get betas
  nonzero.betas = mvrnorm(n = C, mu = rep(0,2), Sigma = cov_mat) * 1/sqrt(C)

  #vector of liability scale betas:
  lia.beta = matrix(0, nrow = 2, ncol = M)
  #positions of nonzero betas:
  beta.pos = sample(1:M, size = C, replace = FALSE)
  
  #assigning nonzero betas to their positions:
  lia.beta[,beta.pos] = t(nonzero.betas)
  
  #covariance matrix for environmen:  
  covarMat = diag(1 - lia.h2, 3 + nsib) #this can be extended by eq (4) from the article to include environmental correlation.
  
  #mean value vector for the environment variable:
  muVec = rep(0, 3 + nsib)
  
  #environment for parents and offspring:
  environ = mvrnorm(n = 2, mu = muVec, Sigma = covarMat) 
  
  ##### doing first indiv outside parallel loop to start the output file:  
  #generates a matrix containing genotype [0/1/2] for each potential parent:
  parents = replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))

  parents.scaled = normalize_geno(parents, MAF = MAF)
  
  #assign parental liabilities:
  parents.geno.lia.first = lia.beta %*% parents.scaled #each row is a phenotype
  
  parents.lia = parents.geno.lia.first + environ[,2:3]
  #assign parental phenotypes:
  parents.cc = parents.lia > lia.T
  ##generating child geno and phenotype:
  #Other Method:
  parents.halved = parents/2 #same parents for all siblings, no need to recalculate this

  offspring_genotype = generateOffspring_geno(parents.halved = parents.halved,
                                              nsib = nsib,
                                              M = M,
                                              MAF = MAF)
  
  offspring = lapply(1:2, FUN = function(n) {
    generateOffspring_pheno(offspring.geno.unscaled = offspring_genotype, MAF = MAF, environ[n,], lia.T, lia.beta[n,])
  }) 
  

  lock = tempfile()
  
  fwrite(as.data.table(matrix(c(1,1,0,0,0,-9, #FID, IID, Father, Mother, Sex, dummy phenotype
                                "offspring_geno" = t(allele.mat[offspring_genotype[,1] + 1,]) ), nrow = 1)),
         file = paste(out, ".ped", sep = ""), sep = " ",
         row.names = F, col.names = F, quote = F)
  
  
  initial_offspring_matrix = lapply(1:2, FUN = function(n) {
    res = matrix(c(unlist(offspring[[n]]), parents.cc[n,], parents.geno.lia.first[n,], parents.lia[n,]), ncol = 9 + 3*nsib, byrow = T)
    colnames(res) =   c("offspring_cc",
                        if (nsib  > 0) paste("siblings_cc", 1:nsib, sep = "_"),
                        "offspring_geno_lia",
                        if (nsib  > 0) paste("siblings_geno_lia", 1:nsib, sep = "_"),
                        "offspring_lia",
                        if (nsib  > 0) paste("siblings_lia", 1:nsib, sep = "_"),
                        paste("parent_cc", 1:2, sep = "_"),
                        paste("parent_geno_lia", 1:2, sep = "_"),
                        paste("parent_lia", 1:2, sep = "_"))
    res
  })
  #here "offspring.geno" is technically a 2xM matrix, however when it is being written to the file it is collapsed
  #into a vetor where each column comes one after the other.
  
  
  cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
  cl = makeCluster(nthreads, type = "SOCK")
  registerDoSNOW(cl)
  iterations = NoChildren
  
  pb = progress_bar$new(
    format = "[:bar] :percent",
    total = iterations,
    width = 100)
  
  progress_num = 1:iterations
  progress = function(n){
    pb$tick(tokens = list(letter = progress_num[n]))
  }
  
  opts = list(progress = progress)
  
  ph = foreach(i = 2:iterations, .packages = c("MASS", "stringr", "flock", "data.table"), .options.snow = opts, .export = c("generateOffspring_geno", "generateOffspring_pheno", "normalize_geno"), .inorder = TRUE) %dopar% {
        
    #environment for parents and offspring:
    environ = mvrnorm(n = 2, mu = muVec, Sigma = covarMat) 
    
    
    #generates a matrix containing genotype [0/1/2] for each potential parent:
    parents = replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))
    
    parents.scaled = normalize_geno(parents, MAF = MAF)
    
    #assign parental liabilities:
    parents.geno.lia = lia.beta %*% parents.scaled #each row is a phenotype
    
    parents.lia = parents.geno.lia + environ[,2:3]
    #assign parental phenotypes:
    parents.cc = parents.lia > lia.T
    ##generating child geno and phenotype:
    #Other Method:
    parents.halved = parents/2 #same parents for all siblings, no need to recalculate this
    
    #generate genotypes for offspring
    offspring_genotype = generateOffspring_geno(parents.halved = parents.halved,
                                                nsib = nsib,
                                                M = M,
                                                MAF = MAF)
    #looping over different sets of liability betas:
    offspring = lapply(1:2, FUN = function(n) {
      generateOffspring_pheno(offspring.geno.unscaled = offspring_genotype, MAF = MAF, environ[n,], lia.T, lia.beta[n,])
    }) 
    
    mat.ph = as.data.table(matrix(c(i,i,0,0,0,0, #FID, IID, Father, Mother, Sex, dummy phenotype
                                    "offspring_geno" = t(allele.mat[offspring_genotype[,1] + 1,]) ), nrow = 1))

    out_dist = paste(out,".ped", sep = "")
    locked = flock::lock(lock)
    fwrite(mat.ph,
           file = out_dist, sep = " ",
           append = T)
    flock::unlock(locked)
    
    offspring_out = lapply(1:2, FUN = function(n){
      res = offspring[[n]]
      res$parents.cc = parents.cc[n,]
      res$parents.geno.lia = parents.geno.lia[n,]
      res$parents.lia = parents.lia[n,]
      res
    })
    offspring_out

  }
  # close(pb)
  stopCluster(cl)
  
 ## pheno_distinct_list = list()
#
 # for (i in 1:2) {
 #   pheno_distinct_list[[i]] = lapply(1:length(ph), FUN = function(n) {
 #     ph[[n]][[i]]
 #   })
 # }
 # ph2 = matrix(unlist(pheno_distinct_list[[1]]), ncol = 9 + 3*nsib, byrow = T)
 # pheno_distinct_matrix = lapply(1:2, FUN = function(n) {
 #   res = matrix(unlist(pheno_distinct_list[[n]]), ncol = 9 + 3*nsib, byrow = T)
 #   res = rbind(initial_offspring_matrix[[n]], res)
 #   res
 # })
 # pheno_distinct_matrix = lapply(1:2, FUN = function(n) as_tibble(pheno_distinct_matrix[[n]]))
  phenotypes = list()
  for (i in 1:2) {
    phenotypes[[i]] = tibble(
      "offspring_cc"   = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.cc[1]),
      "offspring_geno_lia"      = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.geno.lia[1]),
      "offspring_lia"      = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.lia[1]),
      "parent_cc_2" = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.cc[1]),
      "parent_geno_lia_2"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.geno.lia[1]),
      "parent_lia_2"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.lia[1]),
      "parent_cc_1" = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.cc[2]),
      "parent_geno_lia_1"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.geno.lia[2]),
      "parent_lia_1"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.lia[2]),
    )
    if (nsib > 0) {
      for (ii in 1:nsib) {
        cur_sib_stat = paste("siblings_cc", ii, sep = "_")
        cur_sib_gen = paste("siblings_geno_lia", ii, sep = "_")
        cur_sib_lia = paste("siblings_lia", ii, sep = "_")
        phenotypes[[i]][[cur_sib_stat]] = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.cc[-1][[ii]])
        phenotypes[[i]][[cur_sib_gen]]  = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.geno.lia[-1][[ii]])
        phenotypes[[i]][[cur_sib_lia]]  = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.lia[-1][[ii]])
      }
    }
  }
  phenotypes = lapply(1:2, FUN = function(n) {
    initial_offspring_matrix[[n]] = as_tibble(initial_offspring_matrix[[n]])[colnames(phenotypes[[n]])]
    rbind(initial_offspring_matrix[[n]],phenotypes[[n]])
  })

  
  return(list("phenotypes" = phenotypes,
              "lia_betas" = t(lia.beta),
              "lia_threshold" = lia.T))
}


#No. snps:
M = 10000
#No. Children:
NoChildren = 100000
#simulated maf values:
MAF = runif(n = M, min = 0.01, max = 0.5)
#No. causal snps:
C = M*0.1
#liability scale heritability:
lia.h2 = 0.5
#prevalence in parents:
K = .05
#out = "C:\\Users\\FIUN7293\\CommandStationBeta\\EphemeralFolder\\Results\\sim100kx100k_v10"
nsib = 2



for (out in paste("D:\\Work\\Project1\\simulatedData\\sibs2_100kx10k_2traits_C",C,"_v", 1:10,"_maf001", sep = "")) {
  cat("\n:-================================================================================-:\n")
  cat("\nworking on:\n", out, "\n")
  children = generateChildren(MAF = MAF,
                              h2_1 = .5,
                              h2_2 = .5,
                              gen_cor = .5,
                              M = M, 
                              C = C, 
                              nthreads = 10, 
                              NoChildren = NoChildren,
                              K = K,
                              out = out,
                              nsib = nsib)
  
  #save df.map:
  out_file_map = paste(out, ".map", sep = "")
  fwrite(list("CHR" = rep(1,M),
              "SNP" = paste("rs", 1:M, sep = ""),
              "cM" = rep(0, M),
              "BP" = 1:M), 
         file = out_file_map, sep = " ", col.names = F)
  
  # sNP info file
  out_file_snpinfo = paste(out, ".snpinfo", sep = "")
  ph_snpinfo = tibble("lia_betas_pheno1" = children$lia_betas[,1],
                      "lia_betas_pheno2" = children$lia_betas[,2],
                      "maf" = MAF)
  fwrite(ph_snpinfo, file = out_file_snpinfo, sep = " ")
  
  
  for (n in 1:2) {
    ph = tibble("FID" = 1:NoChildren,
                "IID" = 1:NoChildren,
                "CHILD_STATUS" = children[[1]][[n]]$offspring_cc,
                "P1_STATUS" = children[[1]][[n]]$parent_cc_1,
                "P2_STATUS" = children[[1]][[n]]$parent_cc_2,
                "NUM_SIBS" = rep(nsib, NoChildren),
                "SIB_STATUS" = ifelse(rowSums(children[[1]][[n]][,c("siblings_cc_1", "siblings_cc_2")]) > 0, 1, 0))
    #save input file for LTFH, i.e., case/ctrl status for offspring and parents. sibs not included here:
    out_file_phen = paste(out,"_pheno", n, ".phen", sep = "")
    fwrite(ph, file = out_file_phen, sep = " ")
    
    
    
    
    
    
    
    ph = data.frame(
                "FID" = 1:NoChildren,
                "IID" = 1:NoChildren,
                "offspring_lia" = children$phenotypes[[n]]$offspring_lia,
                "offspringgeno_lia" = children$phenotypes[[n]]$offspring_geno_lia,
                "parents_lia_1" = children$phenotypes[[n]]$parent_lia_1,
                "parents_lia_2" = children$phenotypes[[n]]$parent_lia_2,
                "parents_geno_lia_1" = children$phenotypes[[n]]$parent_geno_lia_1,
                "parents_geno_lia_2" = children$phenotypes[[n]]$parent_geno_lia_2,
                "lia_threshold" = children$lia_threshold,
                "prevalens" = K,
                "No_Causal_SNPs" = C,
                "sex" = as.numeric(runif(NoChildren) >= .5))
    for (ii in 1:nsib) {
      sib_lia_name = paste("siblings_lia", ii, sep = "_")
      ph[[sib_lia_name]] = children$phenotypes[[n]][[sib_lia_name]]
      sib_gen_lia_name = paste("siblings_geno_lia", ii, sep = "_")
      ph[[sib_gen_lia_name]] = children$phenotypes[[n]][[sib_gen_lia_name]]
      ph[[paste("sibling_sex", ii, sep = "_")]] = as.numeric(runif(NoChildren) >= .5)
    }
    
    
    out_file_true = paste(out,"_pheno", n, ".true", sep = "")
    fwrite(ph, file = out_file_true,  sep = " ", col.names = T)
    
  }
  
}


##### check if simulated behaves as expected ####
#betas = children$lia_betas[which(children$lia_betas[,1] != 0),]
#cov(betas)
#h2_1 <- h2_2 <- gen_cor <- .5
#cov_mat = diag(c(h2_1, h2_2))
#cov_mat[1,2] <- cov_mat[2,1] <- gen_cor * sqrt(h2_1 * h2_2)
#1/sqrt(C) * diag(2) %*% cov_mat %*% diag(2) * 1/sqrt(C)
#
#true = as_tibble(fread("D:/Work/Project1/simulatedData/sibs2_10kx10k_2traits_C1000_v1_maf002_pheno1.true"))
#
#
#qqnorm(true$offspringgeno_lia)
#qqline(true$offspringgeno_lia)
#qqnorm(true$offspring_lia)
#qqline(true$offspring_lia)
#
#par(mfrow = c(2,5))
#for(n in 1:10) {
#  qqnorm(true[,str_detect(colnames(true), "lia")][,-7][[n]])
#  qqline(true[,str_detect(colnames(true), "lia")][,-7][[n]])
#  #print(p)
#}
#par(mfrow = c(1,1))
#
#
#plot(children$phenotypes[[1]]$offspring_geno_lia, children$phenotypes[[2]]$offspring_geno_lia)
#cov(children$phenotypes[[1]]$offspring_geno_lia, children$phenotypes[[2]]$offspring_geno_lia)
#cor(children$phenotypes[[1]]$offspring_geno_lia, children$phenotypes[[2]]$offspring_geno_lia)
#table(children$phenotypes[[1]]$offspring_cc, children$phenotypes[[2]]$offspring_cc)
#cor(children$phenotypes[[1]]$offspring_cc, children$phenotypes[[2]]$offspring_cc)
#cov(children$phenotypes[[1]]$offspring_cc, children$phenotypes[[2]]$offspring_cc)
#
#
#