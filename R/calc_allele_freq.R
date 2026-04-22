#' Tabulate the frequency of each founder allele.
#' 
#' This function will count a founder allele even when the allele count is zero for one allele.
#' This is used to count the frequency of each founder letter in z. 
#'
#' @param z numeric matrix containing haplotypes for one animal.
#' @param founder_numbers integers representing the founder alleles.
#'
#' @return numeric vector containing the frequency of each founder allele.
allele_count_fxn = function(z, founder_numbers) {
  
  colSums(sapply(founder_numbers, '==', z))
  
} # allele_count_fxn()


#' Get the frequency of each founder allele on each chromosome.
#'
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#' @param founder_numbers integers representing the founder alleles.
#'
#' @return list with one element per chromosome. Each element is a numeric matrix containing the proportion of each sample that is derived from each founder.
#' @export
#'
#' @examples
#' n_mkr = 100
#' n_founders = 8
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = n_founders, markers)
#' allele_freq = calc_allele_freq(pop, 1:n_founders)
calc_allele_freq = function(pop, founder_numbers) {
  
  auto = which(names(pop) != 'X')
  x    = which(names(pop) == 'X')
  
  freq = setNames(as.list(names(pop)), names(pop))
  
  for(i in auto) {
    
    # When alleles are not present at a given marker, the apply()
    # returns a list. We'd like a 2D numeric matrix.
    f = apply(pop[[i]][['females']], 1, allele_count_fxn, founder_numbers)
    m = apply(pop[[i]][['males']],   1, allele_count_fxn, founder_numbers)
    
    freq[[i]] = (f + m) / 
      (4 * dim(pop[[i]][['females']])[3])
    
  } # for(i)
  
  f = apply(pop[[x]][['females']],    1, allele_count_fxn, founder_numbers)
  m = apply(pop[[x]][['males']][,1,], 1, allele_count_fxn, founder_numbers)
  
  freq[[x]] = (f + m) / 
    (4 * dim(pop[[i]][['females']])[3])
  
  return(freq)
  
} # calc_allele_freq()
