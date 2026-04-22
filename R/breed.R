#' Given a population of anaimals, breed them for one generation.
#' 
#' We keep the order of the females constant and randomizes teh males before mating. We return a population of animals of the same size as the input population.
#'
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#'
#' @return a population of animals which have been bred.
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 4, markers)
breed = function(pop, markers) {
  
  # Randomly reorder the males and keep the females in place.
  new_order = sample(1:dim(pop[[1]][['males']])[3])
  
  for(chr in seq_along(pop)) {
    
    pop[[chr]][['males']] = pop[[chr]][['males']][,,new_order] 
    
  } # for(i)
  
  auto = which(names(pop) != 'X')
  x    = which(names(pop) == 'X')
  
  new_pop = pop
  
  for(chr in auto) {
    
    new_pop[[chr]] = breed_one_chr(pop[[chr]], markers[[chr]])
    
  } # for(chr)
  
  new_pop[[x]] = breed_one_x_chr(pop[[x]], markers[[x]])
  
  return(new_pop)
  
} # breed()           
