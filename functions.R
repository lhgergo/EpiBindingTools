library(doParallel)
library(foreach)

# RunNetMHCpan ----------
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
