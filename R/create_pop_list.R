#' Create data structure for one chromosome. 
#'
#' @param mkr named numeric vector containing marker positions in cM and marker names.
#' @param n_pairs integer that is the number of pairs of breeding pairs to create.
#'
#' @return list contining two elements: females and males. Each list element is a named, three dimensional matrix with markers in rows, haplotypes in columns, and breeders in slices.
#' @export
#' @importFrom stats setNames
#' 
#' @examples
#' markers = lapply(1:2, function(z) { setNames(1:10, paste0('mkr', 1:10)) })
#' names(markers) = c(1, 'X')
#' pop4chr = create_one_chr(mkr = markers, n_pairs = 10)
create_one_chr = function(mkr, n_pairs) {
  
  return(list(females = array(0L, dim = c(length(mkr), 2, n_pairs), 
                              dimnames = list(names(mkr), c('hap1', 'hap2'), 
                                              paste0('F', 1:n_pairs))),
              males   = array(0L, dim = c(length(mkr), 2, n_pairs), 
                              dimnames = list(names(mkr), c('hap1', 'hap2'), 
                                              paste0('M', 1:n_pairs)))))
  
} # create_one_chr()


#' Create a list to hold the population haplotypes.
#'
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#' @param n_pairs integer that is the number of pairs of breeders to create.
#'
#' @return list containing one element per chromosome. Each list element contains a list with two elements: female and male. The female and male elements are 3D integer arrays with markers in rows, haplotypes in columns, and sample in rows.
#' @export
#'
#' @examples
#' markers = lapply(1:2, function(z) { setNames(1:10, paste0('mkr', 1:10)) })
#' names(markers) = c(1, 'X')
#' pop = create_pop_list(markers = markers, n_pairs = 10)
create_pop_list = function(markers, n_pairs) {
  
  # Create the founder population as a list.
  pop = setNames(as.list(1:length(markers)), names(markers))
  
  # For each chromosome, simulate founder haplotypes.
  for(chr in seq_along(markers)) {
    
    # Get the markers for the current chromosome.
    curr_mkr = markers[[chr]]
    
    # Create a list for the population.
    # There will be two list elements; one for females and one for males.
    # Each list element will be a 3D array with markers in rows, hap1 & 2 in columns, and
    # samples in slices.
    pop[[chr]] = create_one_chr(curr_mkr, n_pairs)
    
  } # for(chr)
  
  return(pop)
  
} # create_pop_list()
