library(doParallel)
library(foreach)

# RunNetMHCpan ----------
# A user-friendly launcher for single-thread/parallel netMHCpan runs.
# __inputs__
# alleles_loc: path to a file containing alleles to predict binding of epitopes to. Allele names in netMHCpan-compatible format.
# pep_loc: path to .pep file, containing one epitope per line.
# results_out_loc: path to an output directory, where the netMHCpan result files should appear
# parallel threads: an integer number, containing the number of parallel threads to perform predictions. If FALSE, single thread will work.
# netmhcpan_loc: path to the software.
# __outputs__
# No in-R output. Result files are being transferred into path defined in results_out_loc.

RunNetMHCpan <- function(alleles_loc, pep_loc, results_out_loc, parallel_threads = FALSE, netmhcipan_loc = "~/Programok/netMHCpan-4.0/netMHCpan") {
  # reading allele file
  alleles = readLines(alleles_loc)
  
  # producing commands to run
  strings = sapply(alleles, function(allele) {paste0(netmhcipan_loc, " -inptype 1 -f ", pep_loc  ," -BA -a ", allele, " > ", results_out_loc, allele, ".txt")}) 
  
  # creating directory for results
  dir.create(results_out_loc, recursive = TRUE)
  
  # running commands
  if (parallel_threads == FALSE) {
    for (i in 1:length(strings)) {system(strings[i])}
  } else {
    c3 <- makeCluster(parallel_threads)
    registerDoParallel(c3)
    
    foreach (i = 1:length(strings)) %dopar% {system(strings[i])}
    stopCluster(c3)
  }
}

# CreateRecognMatrix
# Loads and processes netMHCpan 4.0 result files into recogn matrices
# __inputs__
# results_out_loc: path to directory containing netMHCpan result files.
# __outputs__
# A list, containing two matrices (one for aff and one for rp data).
# In each matrix the rows represent alleles, the cols represent epitopes, and the numerics in the cells represent binding values

CreateRecognMatrix <- function(results_out_loc) {
  # collecting result file paths
  resfiles_paths <- list.files(results_out_loc, full.names = T)
  
  # loading result files into a simplified list
  pblapply(resfiles_paths, function(current_path) {
    system(paste0("tail -n +51 ", current_path,  " | head -n -5 | tr -s ' '"), intern = T) %>% 
      strsplit(" ") %>% sapply(function(x) {x[c(3, 4, 14, 15)]}) %>% t()
  }) -> resdata
  
  # collecting alleles and epitopes
  alleles <- sapply(resdata, function(x) {x[1, 1]})
  epitopes <- resdata[[1]][, 2]
  
  # transforming everything into matrix format
  recogn_matrix_aff <- t(sapply(resdata, function(x) {as.numeric(x[, 3])})) %>% set_rownames(alleles) %>% set_colnames(epitopes)
  recogn_matrix_rp <- t(sapply(resdata, function(x) {as.numeric(x[, 4])})) %>% set_rownames(alleles) %>% set_colnames(epitopes)
  
  # returning matrices
  list(aff = recogn_matrix_aff, rp = recogn_matrix_rp)
}
