#' Simulate a first generation outcross.
#' 
#' Given a number of breeding pairs, number of founders, and the genetic markers, create a population of breeders derived directly from the breeding funnels.
#'
#' @param n_pairs integer that is the number of breeding pairs. Total animals will be 2 * n_pairs.
#' @param n_founders integer that is the number of founder animals. Must be a power of 2.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#'
#' @return data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#' @export
#' @importFrom stats setNames
#' 
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
simulate_founders = function(n_pairs = 175, n_founders = 8, markers) {
  
  if(n_founders %% 2 != 0) {
    
    stop(paste('n_founders =', n_founders, 'Number of founders must be a power of 2.'))
  
  } # if(n_founders %% 2 != 0)
  
  # Create founder numbers.
  founder_numbers = 1:n_founders
  
  # Get number of F1 pairs.
  num_f1_pairs = max(n_founders / 4, 1)
  
  # Create the data structure to hold the F1 breeders.
  pop = create_pop_list(markers, n_pairs) 
  
  # Get autosomes and Chr X.
  auto = which(names(pop) != 'X')
  x    = which(names(pop) == 'X')
  
  # Create n_pairs of first generation crosses by creating the F1s from the founders
  # and breeding them.
  for(pair in 1:n_pairs) {
    
    funnel = create_breeding_funnel(founder_numbers, markers, num_f1_pairs)
    
    # At this point, we have the F1 breeders for the funnel.
    # Next we need to breed them until we have two first generation mice and then
    # put them in the population.
    first_gen = breed_funnel(funnel, markers)
    
    for(chr in seq_along(first_gen)) {
      
      pop[[chr]][['females']][,,pair] = first_gen[[chr]][['females']][,,1]
      pop[[chr]][['males']][,,pair]   = first_gen[[chr]][['males']][,,1]
      
    } # for(chr)
    
  } # for(pair)
  
  return(pop)
  
} # simulate_founders()
