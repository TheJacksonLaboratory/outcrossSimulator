#' Create one breeding funnel which combines all founders in a random order.
#' 
#' This function creates pairs of breeders based on the hap_order. The females will have haplotypes from the odd values in hap_order and the males will have haplotypes from the even numbers. There will be n_founders / 2 breeding pairs.
#'
#' @param hap_order character vector with founder haplotypes permuted for this funnel.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#'
#' @return list containing one breeding funnel.
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' markers = lapply(1:2, function(z) { setNames(1:10, paste0('mkr', 1:10)) })
#' names(markers) = c(1, 'X')
#' founders = 1:4
#' funnel = create_breeding_funnel(sample(founders), markers = markers)
create_breeding_funnel = function(hap_order, markers) {
  
  # Get the number of founders.
  n_pairs = length(hap_order) / 2
  
  # Create the breeding funnel at the F1 level.
  funnel = create_pop_list(markers, n_pairs)
  
  # Get autosomes and Chr X.
  auto = which(names(markers) != 'X')
  x    = which(names(markers) == 'X')
  
  # Populate the autosomes with F1 haplotypes.
  for(chr in auto) {
    
    # Get the markers for the current chromosome.
    curr_mkr = markers[[chr]]
    
    # Create inbred founders in the order found in 'hap_order'.
    for(i in 1:n_pairs) {
      
      # Populate female.
      funnel[[chr]][['females']][,,i] = hap_order[2 * i - 1]
      
      # Populate male.
      funnel[[chr]][['males']][,,i]   = hap_order[2 * i]
      
    } # for(i)
    
  } # for(chr)
  
  # Populate Chr X with F1 haplotypes.
  curr_mkr = markers[[x]]
  
  # Populate Chr X with F1 haplotypes.
  for(i in 1:n_pairs) {
    
    # Populate female.
    funnel[[x]][['females']][,,i] = hap_order[2 * i - 1]
    
    # Populate male.
    funnel[[x]][['males']][,'hap1',i]   = hap_order[2 * 1]
    
  } # for(i)
  
  return(funnel)
  
} # create_breeding_funnel()
