#' Extract the haplotypes for one sample.
#'
#' @param sex character string that is either 'feamles' or 'males'.
#' @param index index of the sample to extract.
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#'
#' @return list containing one element per chromosome. Each element is a numeric matrix wtih markers in rows and two haplotypes in columns.
#' @export
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
#' smpl = extract_sample('females', index = 5, pop)
extract_sample = function(sex, index, pop) {
  
  # Verify that the index is withing bounds.
  if(index > dim(pop[[1]][[sex]])[3]) {
    
    stop(paste('Sample index', index, 'is greater than the number of smaples (',
               dim(pop[[1]][[sex]])[3],')'))
    
  } # if(index > dim(pop[[1]][[sex]])[3])
  
  retval = setNames(as.list(names(pop)), names(pop))
  
  for(chr in seq_along(pop)) {
    
    retval[[chr]] = pop[[chr]][[sex]][,,index]
    
  } # for(chr)
  
  return(retval)
  
} # extract_sample()
