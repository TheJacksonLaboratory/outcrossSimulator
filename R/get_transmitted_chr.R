#' Get the transmitted chromosome from one meiosis.
#'
#' Given a pair of haplotypes for one animal, simulate crossovers and return a single haplotype, possibly with crossovers.
#'
#' @param haps integer matrix with markers in rows and two haplotypes in columns.
#' @param mkr named numeric vector containing the cM positions for each marker on the current chromosome.
#' @param nu floating point number to pass to the xoi::simStahl function. Interference parameter.
#' @param p floating point number to pass to the xoi::simStahl function. Proportion of crossovers with no crossover interference.
#'
#' @return numeric vector containing the sequence on one haplotype.
#' @export
#' @importFrom stats setNames
#' @importFrom xoi simStahl
#'
#' @examples
#' n_mkr = 51
#' haps = matrix(rep(1:2, each = n_mkr), nrow = n_mkr, ncol = 2, 
#' dimnames = list(paste0('mkr', 1:n_mkr), paste0('hap', 1:2)))
#' mkr = setNames(1:n_mkr, paste0('mkr', 1:n_mkr))
#' trans_hap = get_transmitted_chr(haps, mkr, nu = 11, p = 0.05)
get_transmitted_chr = function(haps, mkr, nu, p) {
  
  # Simulate crossover positions along the transmitted chromatid.
  xo_pos = xoi::simStahl(n.sim            = 1L,
                         nu               = nu,
                         p                = p,
                         L                = max(mkr),
                         obligate_chiasma = TRUE)[[1]]
  
  # Randomly choose which parental haplotype the gamete starts from.
  current_hap = sample(1:2, 1)
  
  # Start with the selected haplotype.
  result = haps[, current_hap]
  
  # Apply crossovers if any are present.
  if(length(xo_pos) > 0) {
    
    # Sort positions.
    xo_pos = sort(xo_pos)
    
    # Retain crossovers within the marker range.
    xo_pos = xo_pos[xo_pos >= mkr[1] & xo_pos <= mkr[length(mkr)]]
    
    if(length(xo_pos) > 0) {
      
      for(j in seq_along(xo_pos)) {
        
        # Get nearest markers around crossover site.
        proximal = max(which(mkr <= xo_pos[j]))
        distal   = min(proximal + 1, length(mkr))
        
        # Switch to the other parental haplotype and copy the distal segment.
        current_hap                   = 3 - current_hap
        result[distal:length(result)] = haps[distal:length(result), current_hap]
        
      } # for(j)
      
    } # if(length(xo_pos) > 0)
    
  } # if(length(xo_pos) > 0)
  
  return(result)
  
} # get_transmitted_chr()
