#' Calculate the proportion of homozygosity for each sample on each chromosome.
#'
#' @param pop data structure containing breeding population haplotypes. List containing one element per chromosome. Each element is a list with two elements named: females and males. Each element is a 3D numeric array with markers in rows, haplotypes in columns, and samples in slices.
#'
#' @return numeric matrix containing the proportion of homozygous markers for each sample.
#' @export
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
#' hom = calc_hom(pop)
calc_hom = function(pop) {
  
  samples   = sapply(pop[[1]], function(z) { dimnames(z)[[3]] })
  n_samples = length(samples)
  n_markers = sum(sapply(pop, function(z) { nrow(z[['females']]) }))
  
  p_hom = matrix(0, nrow = n_samples, ncol = 1, dimnames = list(samples, 'p_hom'))
  
  for(chr in seq_along(pop)) {
    
    # Females.
    pop_chr = pop[[chr]][['females']]
    hom     = colSums(pop_chr[,1,] == pop_chr[,2,])
    p_hom[names(hom),1] = p_hom[names(hom),1] + hom
    
    # Males.
    pop_chr = pop[[chr]][['males']]
    hom     = colSums(pop_chr[,1,] == pop_chr[,2,])
    p_hom[names(hom),1] = p_hom[names(hom),1] + hom
    
  } # for(chr)
  
  
  return(p_hom / n_markers)
  
} # calc_hom()
