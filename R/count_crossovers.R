#' Count the number of crossovers in each sample on each chromosome.
#'
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#'
#' @return numeric matrix containing the number of crossovers for each sample on each chromosome.
#' @export
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
#' xo = count_crossovers(pop)
count_crossovers = function(pop) {
  
  # Get the autosomes and X Chr.
  auto = which(names(pop) != 'X')
  x    = which(names(pop) == 'X')
  
  # Get the sample names.
  f_names = dimnames(pop[[1]][['females']])[[3]]
  m_names = dimnames(pop[[1]][['males']])[[3]]
  
  # Get the number of samples for each sex.
  n_females = length(f_names)
  n_males   = length(m_names)
  
  # Create crossover output object.
  xo = matrix(0, nrow = c(n_females + n_males), ncol = length(pop),
              dimnames = list(c(f_names, m_names), names(pop)))
  
  # Get the crossovers on the autosomes.
  for(chr in auto) {
    
    # Count crossovers by getting a diff of the haplotypes on each chromosome.
    # Females
    h1 = apply(diff(pop[[chr]][['females']][,'hap1',]) != 0, 2, which)
    h2 = apply(diff(pop[[chr]][['females']][,'hap2',]) != 0, 2, which)
    xo[f_names, chr] = sapply(h1, length) + sapply(h2, length)
    
    
    # Males
    h1 = apply(diff(pop[[chr]][['males']][,'hap1',]) != 0, 2, which)
    h2 = apply(diff(pop[[chr]][['males']][,'hap2',]) != 0, 2, which)
    xo[m_names, chr] = sapply(h1, length) + sapply(h2, length)
    
  } # for(chr)
  
  # Females: Chr X
  h1 = apply(diff(pop[[x]][['females']][,'hap1',]) != 0, 2, which)
  h2 = apply(diff(pop[[x]][['females']][,'hap2',]) != 0, 2, which)
  xo[f_names, x] = sapply(h1, length) + sapply(h2, length)
  
  # Males: Chr X
  if(dim(pop[[x]][['males']])[3] > 1) {
    
    # Approach for multiple samples.
    h1 = apply(diff(pop[[x]][['males']][,'hap1',]) != 0, 2, which)
    xo[m_names,x] = sapply(h1, length)
    
  } else {
    
    # Approach for one sample.
    h1 = which(diff(pop[[x]][['males']][,'hap1',]) != 0)
    xo[m_names,x] = length(h1)
    
  } # else
  
  return(xo)
  
} # count_crossovers()
