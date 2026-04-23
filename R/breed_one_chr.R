#' Breed the females and males on one chromosome.
#'
#' @param breeders list containing the haplotype for a population of animals on one chromosome.
#' @param mkr named numeric vector containing the markers for the current chromosome.
#' @param nu floating point number to pass to the xoi::simStahl function. Interference parameter.
#' @param p floating point number to pass to the xoi::simStahl function. Proportion of crossovers with no crossover interference.
#'
#' @return list containing the haplotype for a population of animals on one chromosome.
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' founders = 1:4
#' breeders = create_breeding_funnel(hap_order = 1:4, markers)
#' one_chr = breed_one_chr(breeders[[1]], markers[[1]])
breed_one_chr = function(breeders, mkr, nu = 11, p = 0.05) {
  
  n_pairs = dim(breeders[[1]])[3]
  
  # Create new breeder data structure.
  new_breeders = create_one_chr(mkr, n_pairs)
  
  for(i in 1:n_pairs) {
    
    # Get the haplotypes for the breeding pair.
    f = breeders[['females']][,,i]
    m = breeders[['males']][,,i]
    
    # Create new female.
    new_breeders[['females']][,1,i] = get_transmitted_chr(haps = f, mkr = mkr, nu = nu, p = p)
    new_breeders[['females']][,2,i] = get_transmitted_chr(haps = m, mkr = mkr, nu = nu, p = p) 
    
    # Create new male.
    new_breeders[['males']][,1,i] = get_transmitted_chr(haps = f, mkr = mkr, nu = nu, p = p)
    new_breeders[['males']][,2,i] = get_transmitted_chr(haps = m, mkr = mkr, nu = nu, p = p)
    
  } # for(i)
  
  return(new_breeders)
  
} # breed_one_chr()
