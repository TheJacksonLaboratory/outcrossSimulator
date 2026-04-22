#' Calculate the recombination block size for each sample on each chromosome.
#'
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#'
#' @return list containing recombination block sizes for each sample on each chromosome.
#' @export
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
#' blocks = calc_block_size(pop, markers)
calc_block_size = function(pop, markers) {
  
  # Get the autosomes and X Chr.
  auto = which(names(pop) != 'X')
  x    = which(names(pop) == 'X')
  
  # Get the sample names.
  f_names = dimnames(pop[[1]][['females']])[[3]]
  m_names = dimnames(pop[[1]][['males']])[[3]]
  
  # Get the number of samples for each sex.
  n_females = length(f_names)
  n_males   = length(m_names)
  
  # Create recombination block output object.
  blocks = setNames(vector('list', length(pop)), names(pop))
  for(chr in seq_along(pop)) {
    
    blocks[[chr]] = setNames(vector('list', length(pop[[chr]])), names(pop[[chr]]))
    
    blocks[[chr]][['females']] = setNames(as.list(1:n_females), f_names)
    blocks[[chr]][['males']]   = setNames(as.list(1:n_males),   m_names)
    
  } # for(chr)
  
  # Get the recombination blocks on the autosomes.
  for(chr in auto) {
    
    # Get recombination blocks by getting a diff of the haplotypes on each
    # chromosome.
    # Females
    fem = pop[[chr]][['females']]
    
    # Place -1 in first and last marker to get start and end of chromosome.
    tmp = array(-1, dim = c(dim(fem) + c(1,0,0)), 
                dimnames = list(c('start', rownames(fem)), colnames(fem), dimnames(fem)[[3]]))
    tmp[2:(nrow(tmp) - 1),,] = fem[-nrow(fem),,]
    
    # Take the diff of each haplotype and get the marker names at the crossover points.
    h1 = apply(diff(tmp[,'hap1',]) != 0, 2, which)
    h2 = apply(diff(tmp[,'hap2',]) != 0, 2, which)
    
    # Map the crossovers to the genetic map.
    h1 = lapply(h1, function(z) { diff(markers[[chr]][names(z)]) })
    h2 = lapply(h2, function(z) { diff(markers[[chr]][names(z)]) })
    
    # Put the results in the output object.
    blocks[[chr]][['females']] = mapply(c, h1, h2)
    
    # Males
    mal = pop[[chr]][['males']]
    
    # Place -1 in first and last marker to get start and end of chromosome.
    tmp = array(-1, dim = c(dim(mal) + c(1,0,0)), 
                dimnames = list(c('start', rownames(mal)), colnames(mal), dimnames(mal)[[3]]))
    tmp[2:(nrow(tmp) - 1),,] = mal[-nrow(mal),,]
    
    # Take the diff of each haplotype and get the marker names at the crossover points.
    h1 = apply(diff(tmp[,'hap1',]) != 0, 2, which)
    h2 = apply(diff(tmp[,'hap2',]) != 0, 2, which)
    
    # Map the crossovers to the genetic map.
    h1 = lapply(h1, function(z) { diff(markers[[chr]][names(z)]) })
    h2 = lapply(h2, function(z) { diff(markers[[chr]][names(z)]) })
    
    # Put the results in the output object.
    blocks[[chr]][['males']] = mapply(c, h1, h2)
    
  } # for(chr)
  
  # Females: Chr X
  fem = pop[[x]][['females']]
  
  # Place -1 in first and last marker to get start and end of chromosome.
  tmp = array(-1, dim = c(dim(fem) + c(1,0,0)), 
              dimnames = list(c('start', rownames(fem)), colnames(fem), dimnames(fem)[[3]]))
  tmp[2:(nrow(tmp) - 1),,] = fem[-nrow(fem),,]
  
  # Take the diff of each haplotype and get the marker names at the crossover points.
  h1 = apply(diff(tmp[,'hap1',]) != 0, 2, which)
  h2 = apply(diff(tmp[,'hap2',]) != 0, 2, which)
  
  # Map the crossovers to the genetic map.
  h1 = lapply(h1, function(z) { diff(markers[[chr]][names(z)]) })
  h2 = lapply(h2, function(z) { diff(markers[[chr]][names(z)]) })
  
  # Put the results in the output object.
  blocks[[x]][['females']] = mapply(c, h1, h2)
  
  # Males: Chr X
  if(dim(pop[[x]][['males']])[3] > 1) {
    
    # Approach for multiple samples.
    mal = pop[[x]][['males']]
    
    # Place -1 in first and last marker to get start and end of chromosome.
    tmp = array(-1, dim = c(dim(mal) + c(1,0,0)), 
                dimnames = list(c('start', rownames(mal)), colnames(mal), dimnames(mal)[[3]]))
    tmp[2:(nrow(tmp) - 1),,] = mal[-nrow(mal),,]
    
    # Take the diff of each haplotype and get the marker names at the crossover points.
    h1 = apply(diff(tmp[,'hap1',]) != 0, 2, which)
    
    # Map the crossovers to the genetic map.
    h1 = lapply(h1, function(z) { diff(markers[[x]][names(z)]) })
    
    # Put the results in the output object.
    blocks[[x]][['males']] = h1
    
  } else {
    
    # Approach for one sample.
    mal = pop[[x]][['males']]
    
    # Place -1 in first and last marker to get start and end of chromosome.
    tmp = array(-1, dim = c(dim(mal) + c(1,0,0)), 
                dimnames = list(c('start', rownames(mal)), colnames(mal), dimnames(mal)[[3]]))
    tmp[2:(nrow(tmp) - 1),,] = mal[-nrow(mal),,]
    
    # Take the diff of each haplotype and get the marker names at the crossover points.
    h1 = which(diff(tmp[,'hap1',]) != 0)
    
    # Map the crossovers to the genetic map.
    h1 = diff(markers[[x]][names(h1)])
    
    # Put the results in the output object.
    blocks[[x]][['males']] = h1
    
  } # else
  
  return(blocks)
  
} # calc_block_size()
