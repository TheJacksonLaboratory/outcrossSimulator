#' Breed a funnel of founders to a single mating pair.
#' 
#' This function accpets a breeding funnel and continues to call itself until there are only two animals, one female and one male, left.
#'
#' @param funnel list containing the founders set up to breed.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#' @param nu floating point number to pass to the xoi::simStahl function. Interference parameter.
#' @param p floating point number to pass to the xoi::simStahl function. Proportion of crossovers with no crossover interference.
#'
#' @return population list containing the haplotypes for a pair of animals. 
#' @export
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' founders = 1:4
#' funnel = create_breeding_funnel(sample(founders), markers = markers, num_f1_pairs = 2)
#' mating_pair = breed_funnel(funnel, markers)
breed_funnel = function(funnel, markers, nu = 11, p = 0.05) {
  
  n_pairs = dim(funnel[[1]][[1]])[3]
  
  # Get autosomes and Chr X.
  auto = which(names(funnel) != 'X')
  x    = which(names(funnel) == 'X')
  
  # Create a return object of offspring.
  next_gen = create_pop_list(markers, n_pairs)
  
  # For each chromosome, simulate founder haplotypes.
  for(chr in auto) {
    
    next_gen[[chr]] = breed_one_chr(funnel[[chr]], markers[[chr]], nu, p)
    
  } # for(chr)
  
  next_gen[[x]] = breed_one_x_chr(funnel[[x]], markers[[x]], nu, p)
  
  # If there is more than one female, then we need to keep breeding the funnel.
  # Otherwise, we return the next generation that we just bred.
  n_females = dim(next_gen[[1]][[1]])[3]
  
  if(n_females > 1) {
    
    # Take the first half of the offspring. We only need one female and one male
    # from this round of breeding.
    for(chr in seq_along(next_gen)) {
      
      # Subset the females and males.
      first_half  = 1:(n_females / 2)
      second_half = (n_females / 2 + 1):n_females
      
      next_gen[[chr]][['females']] = next_gen[[chr]][['females']][,,first_half, drop = FALSE]
      next_gen[[chr]][['males']]   = next_gen[[chr]][['males']][,,second_half,  drop = FALSE]
      
    } # for(chr)
    
    next_gen = breed_funnel(funnel = next_gen, markers = markers)
    
  } # else
  
  return(next_gen)
  
} # breed_funnel()
