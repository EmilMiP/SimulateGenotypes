library(MASS)
library(stringr)
library(future)
library(future.apply)
library(progressr)
library(data.table)
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
generateOffspring_pheno = function(offspring.geno.unscaled, MAF, environ, lia.beta) {
  offspring.geno.scaled = normalize_geno(offspring.geno.unscaled, MAF = MAF)
  offspring.geno.lia.first = lia.beta %*% offspring.geno.scaled
  offspring.lia.first = offspring.geno.lia.first  + environ[-(2:3)]
  return(list(
    "offspring.geno.lia" = offspring.geno.lia.first,
    "offspring.lia" = offspring.lia.first
  ))
}


generateChildren = function(
  NoChildren = 10000,
  M = 10000,
  MAF,
  C = M*0.1,
  gen_cov_mat = matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2, nrow = 2),
  nsib = 2,
  allele.mat = matrix(c("A","A","A","G", "G","G"), nrow = 3, ncol = 2, byrow = T),
  out = "D:/Work/Project1/simulatedData/ph",
  overlap = "none"
) {
  
  #get betas
  nonzero.betas = mvrnorm(n = C, mu = rep(0, nrow(gen_cov_mat)), Sigma = gen_cov_mat) * 1/sqrt(C)
  
  #vector of liability scale betas:
  lia.beta = matrix(0, nrow = 2, ncol = M)
  #positions of nonzero betas:
  if(overlap == "none") { #no overlap between causal alleles
    ctr = 1
    positions = 1:M
    while (ctr <= nrow(gen_cov_mat)) {
      beta.pos = sample(positions, size = C, replace = FALSE)
      lia.beta[ctr, beta.pos] = nonzero.betas[,ctr]
      positions = positions[-beta.pos]
      ctr = ctr + 1
    }
    
  } else if (overlap == "complete") { #complete overlap 
    lia.beta[, sample(1:M, size = C, replace = FALSE)] = t(nonzero.betas)
  } else { # partial overlap
    #assigning first positions, common positions
    overlap_freq = as.numeric(overlap)
    ctr = 1
    all_positions = 1:M
    removed_positions = c()
    beta.pos = sample(all_positions, size = round(overlap_freq * C), replace = FALSE)
    lia.beta[, beta.pos] = t(nonzero.betas[1:length(beta.pos),])
    
    removed_positions = c(removed_positions, beta.pos) # slow, cba to find fast solution for now
    positions = all_positions[-removed_positions]
    
    #assigning the rest individually
    while (ctr <= nrow(gen_cov_mat)) {
      beta.pos = sample(positions, size = round( (1 - overlap_freq) * C), replace = FALSE)
      lia.beta[ctr, beta.pos] = t(nonzero.betas[-(1:length(beta.pos)),ctr])
      removed_positions = c(removed_positions, beta.pos) # slow, cba to find fast solution for now
      positions = all_positions[-removed_positions]
      ctr = ctr + 1
    }
    
  }

  #covariance matrix for environment:  
  covarMat = lapply(diag(gen_cov_mat), function(x) diag(1 - x, 3 + nsib)) 
  #this can be extended by eq (4) from the article to include environmental correlation.
  
  #mean value vector for the environment variable:
  muVec = rep(0, 3 + nsib)
  
  #environment for parents and offspring:
  environ = sapply(1:length(covarMat), function(n) mvrnorm(n = 1, mu = muVec, Sigma = covarMat[[n]])) 
  
  ##### doing first indiv outside parallel loop to start the output file:  
  #generates a matrix containing genotype [0/1/2] for each potential parent:
  parents = replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))

  parents.scaled = normalize_geno(parents, MAF = MAF)
  
  #assign parental liabilities:
  parents.geno.lia.first = lia.beta %*% parents.scaled #each row is a phenotype
  
  parents.lia = parents.geno.lia.first + t(environ[2:3,])

  ##generating child geno and liabilities:
  parents.halved = parents/2 #same parents for all siblings, no need to recalculate this
  #we are considering the "average" passed on. variation for siblings are in the ties, i.e.
  #0.5 and 1.5 values being rounded with equal proba up or down.
  offspring_genotype = generateOffspring_geno(parents.halved = parents.halved,
                                              nsib = nsib,
                                              M = M,
                                              MAF = MAF)
  
  offspring = lapply(1:nrow(gen_cov_mat), FUN = function(n) {
    generateOffspring_pheno(offspring.geno.unscaled = offspring_genotype, MAF = MAF, environ[,n], lia.beta[n,])
  }) 
  

  lock = tempfile()
  
  fwrite(as.data.table(matrix(c(1,1,0,0,0,-9, #FID, IID, Father, Mother, Sex, dummy phenotype
                                "offspring_geno" = t(allele.mat[offspring_genotype[,1] + 1,]) ), nrow = 1)),
         file = paste(out, ".ped", sep = ""), sep = " ",
         row.names = F, col.names = F, quote = F)
  
  
  initial_offspring_matrix = lapply(1:nrow(gen_cov_mat), FUN = function(n) {
    res = matrix(c(unlist(offspring[[n]]), parents.geno.lia.first[n,], parents.lia[n,]), ncol = 6 + 2*nsib, byrow = T)
    colnames(res) =   c("offspring_geno_lia",
                        if (nsib  > 0) paste("siblings_geno_lia", 1:nsib, sep = "_"),
                        "offspring_lia",
                        if (nsib  > 0) paste("siblings_lia", 1:nsib, sep = "_"),
                        paste("parent_geno_lia", 1:2, sep = "_"),
                        paste("parent_lia", 1:2, sep = "_"))
    res
  })
  #here "offspring.geno" is technically a 2xM matrix, however when it is being written to the file it is collapsed
  #into a vetor where each column comes one after the other.
  
  progress_bar_n = 1:NoChildren
  pbn = 1
  p <- progressr::progressor(along = progress_bar_n)
  
  # cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
  # cl = makeCluster(nthreads, type = "SOCK")
  # registerDoSNOW(cl)
  # iterations = NoChildren
  
  # pb = progress_bar$new(
  #   format = "[:bar] :percent",
  #   total = iterations,
  #   width = 100)
  # 
  # progress_num = 1:iterations
  # progress = function(n){
  #   pb$tick(tokens = list(letter = progress_num[n]))
  # }
  # 
  # opts = list(progress = progress)
  
  ph = future.apply::future_lapply(X = 2:NoChildren, FUN = function(i){
        
    #environment for parents and offspring:
    environ = sapply(1:length(covarMat), function(n) mvrnorm(n = 1, mu = muVec, Sigma = covarMat[[n]])) 
    
    
    #generates a matrix containing genotype [0/1/2] for each potential parent:
    parents = replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))
    
    parents.scaled = normalize_geno(parents, MAF = MAF)
    
    #assign parental liabilities:
    parents.geno.lia = lia.beta %*% parents.scaled #each row is a phenotype
    
    parents.lia = parents.geno.lia + t(environ[2:3,])
    ##generating child geno and phenotype:
    parents.halved = parents/2 #same parents for all siblings, no need to recalculate this
    
    #generate genotypes for offspring
    offspring_genotype = generateOffspring_geno(parents.halved = parents.halved,
                                                nsib = nsib,
                                                M = M,
                                                MAF = MAF)
    #looping over different sets of liability betas:
    offspring = lapply(1:nrow(gen_cov_mat), FUN = function(n) {
      generateOffspring_pheno(offspring.geno.unscaled = offspring_genotype, MAF = MAF, environ[, n], lia.beta[n,])
    }) 
    
    mat.ph = as.data.table(matrix(c(i,i,0,0,0,-9, #FID, IID, Father, Mother, Sex, dummy phenotype
                                    "offspring_geno" = t(allele.mat[offspring_genotype[,1] + 1,]) ), nrow = 1))

    out_dist = paste(out,".ped", sep = "")
    locked = flock::lock(lock)
    fwrite(mat.ph,
           file = out_dist, sep = " ",
           append = T)
    flock::unlock(locked)
    
    offspring_out = lapply(1:nrow(gen_cov_mat), FUN = function(n){
      res = offspring[[n]]
      res$parents.geno.lia = parents.geno.lia[n,]
      res$parents.lia = parents.lia[n,]
      res
    })
    offspring_out
    

    p(sprintf("%g", pbn))
    pbn = pbn + 1

  }, future.seed = T)

  phenotypes = list()
  for (i in 1:nrow(gen_cov_mat)) {
    phenotypes[[i]] = tibble(
      "offspring_geno_lia"      = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.geno.lia[1]),
      "offspring_lia"      = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.lia[1]),
      "parent_geno_lia_2"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.geno.lia[1]),
      "parent_lia_2"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.lia[1]),
      "parent_geno_lia_1"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.geno.lia[2]),
      "parent_lia_1"    = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$parents.lia[2]),
    )
    if (nsib > 0) {
      for (ii in 1:nsib) {
        cur_sib_gen = paste("siblings_geno_lia", ii, sep = "_")
        cur_sib_lia = paste("siblings_lia", ii, sep = "_")
        phenotypes[[i]][[cur_sib_gen]]  = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.geno.lia[-1][[ii]])
        phenotypes[[i]][[cur_sib_lia]]  = sapply(1:length(ph), FUN = function(n) ph[[n]][[i]]$offspring.lia[-1][[ii]])
      }
    }
  }
  phenotypes = lapply(1:nrow(gen_cov_mat), FUN = function(n) {
    initial_offspring_matrix[[n]] = as_tibble(initial_offspring_matrix[[n]])[colnames(phenotypes[[n]])]
    rbind(initial_offspring_matrix[[n]],phenotypes[[n]])
  })

  
  return(list("phenotypes" = phenotypes,
              "lia_betas" = as_tibble(t(lia.beta))))
}


#No. snps:
M = 10000
#No. Children:
NoChildren = 10000
#simulated maf values:
MAF = runif(n = M, min = 0.01, max = 0.5)
#No. causal snps:
C = M*0.01
#liability scale heritability:
lia.h2 = 0.5
#prevalence in parents:
K = .05
#out = "C:\\Users\\FIUN7293\\CommandStationBeta\\EphemeralFolder\\Results\\sim100kx100k_v10"
nsib = 2

overlap = "none"
gen_cov_mat = matrix(c(0.5, 0.5 * sqrt(0.5*0.5), 
                       0.5 * sqrt(0.5*0.5), 0.5), ncol = 2, nrow = 2)

nthreads = 3
plan(multisession)

threshold = qnorm( 0.05, lower.tail = F)
handlers("progress")
for (out in paste("C:/Users/emp/simulatedData/ph_C",C,"_v", 1:1,"_maf001", sep = "")) {
  cat("\n:-================================================================================-:\n")
  cat("\nworking on:\n", out, "\n")
  with_progress({children <- generateChildren(MAF = MAF,
                              M = M, 
                              C = C,
                              gen_cov_mat = gen_cov_mat,
                              NoChildren = NoChildren,
                              out = out,
                              nsib = nsib,
                              overlap = overlap)})
  
  #save df.map:
  out_file_map = paste(out, ".map", sep = "")
  fwrite(list("CHR" = rep(1,M),
              "SNP" = paste("rs", 1:M, sep = ""),
              "cM" = rep(0, M),
              "BP" = 1:M), 
         file = out_file_map, sep = " ", col.names = F)
  
  # sNP info file
  out_file_snpinfo = paste(out, ".snpinfo", sep = "")
  ph_snpinfo = tibble(as_tibble(children$lia_betas),
                      "maf" = MAF)
  colnames(ph_snpinfo)[1:nrow(gen_cov_mat)] = paste0("lia_betas_pheno", 1:nrow(gen_cov_mat))
  fwrite(ph_snpinfo, file = out_file_snpinfo, sep = " ")
  
  
  for (n in 1:nrow(gen_cov_mat)) {
    ph = tibble("FID" = 1:NoChildren,
                "IID" = 1:NoChildren,
                "CHILD_STATUS" = (children[[1]][[n]]$offspring_lia > threshold) + 0L,
                "P1_STATUS" = (children[[1]][[n]]$parent_lia_1 > threshold) + 0L,
                "P2_STATUS" = (children[[1]][[n]]$parent_lia_2 > threshold) + 0L,
                "NUM_SIBS" = rep(nsib, NoChildren),
                "SIB_STATUS" = ifelse(rowSums(children[[1]][[n]][, paste0("siblings_lia_", 1:nsib)]) > threshold, 1, 0))
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

