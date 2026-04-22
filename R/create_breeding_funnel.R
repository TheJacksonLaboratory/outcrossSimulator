#' Create one breeding funnel which combines all founders in a random order and produces one pair of breeders.
#'
#' @param hap_order character vector with founder haplotypes permuted for this funnel.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#' @param num_f1_pairs integer that is the number of F1 breeding pairs to create.
#'
#' @return list containing one breeding funnel.
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' markers = lapply(1:2, function(z) { setNames(1:10, paste0('mkr', 1:10)) })
#' names(markers) = c(1, 'X')
#' founders = 1:4
#' funnel = create_breeding_funnel(sample(founders), markers = markers, num_f1_pairs = 2)
create_breeding_funnel = function(hap_order, markers, num_f1_pairs) {
  
  # Create the breeding funnel at the F1 level.
  funnel = create_pop_list(markers, num_f1_pairs)
  
  # Get autosomes and Chr X.
  auto = which(names(markers) != 'X')
  x    = which(names(markers) == 'X')
  
  # Populate the autosomes with F1 haplotypes.
  for(chr in auto) {
    
    # Get the markers for the current chromosome.
    curr_mkr = markers[[chr]]
    
    # Create F1 hybrids in the order found in 'hap_order'.
    for(i in 1:num_f1_pairs) {
      
      # Populate female.
      funnel[[chr]][['females']][,'hap1',i] = hap_order[(4 * (i-1)) + 1]
      funnel[[chr]][['females']][,'hap2',i] = hap_order[(4 * (i-1)) + 2]
      
      # Populate male.
      funnel[[chr]][['males']][,'hap1',i]   = hap_order[(4 * (i-1)) + 3]
      funnel[[chr]][['males']][,'hap2',i]   = hap_order[(4 * (i-1)) + 4]
      
    } # for(i)
    
  } # for(chr)
  
  # Populate Chr X with F1 haplotypes.
  curr_mkr = markers[[x]]
  
  # Populate Chr X with F1 haplotypes.
  for(i in 1:num_f1_pairs) {
    
    # Populate female.
    funnel[[x]][['females']][,'hap1',i] = hap_order[(4 * (i-1)) + 1]
    funnel[[x]][['females']][,'hap2',i] = hap_order[(4 * (i-1)) + 2]
    
    # Populate male.
    funnel[[x]][['males']][,'hap1',i]   = hap_order[(4 * (i-1)) + 3]
    
  } # for(i)
  
  return(funnel)
  
} # create_breeding_funnel()
