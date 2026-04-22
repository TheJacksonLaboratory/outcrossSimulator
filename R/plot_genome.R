#' Extract the haplotypes for one sample.
#'
#' @param haps haplotypes for one sample from extract_sample.
#' @param markers list containing one chromosome per element. Each element contains a named numeric vector containing marker positions in cM and marker names.
#' @param colors vector of colors, one per founder.
#' @param ... other argument to pass to plot.
#'
#' @return list containing one element per chromosome. Each element is a numeric matrix wtih markers in rows and two haplotypes in columns.
#' @export
#' @importFrom graphics axis par rect
#'
#' @examples
#' n_mkr = 100
#' markers = lapply(1:2, function(z) { setNames(1:n_mkr, paste0('mkr', 1:n_mkr)) })
#' names(markers) = c(1, 'X')
#' pop = simulate_founders(n_pairs = 10, n_founders = 8, markers)
#' smpl = extract_sample('females', index = 5, pop)
#' plot_genome(smpl, markers, colors = 1:4)
plot_genome = function(haps, markers, colors, ...) {
  
  chr_width = 0.3
  chr_names = names(haps)
  
  par(mgp = c(1.5, 0.5, 0))
  plot(1, 1, col = 'white', xlim = c(0, max(sapply(markers, max))), 
       ylim = c(1, length(markers)), las = 1,
       xlab = 'cM', ylab = 'Chromosome', yaxt = 'n')
  axis(side = 2, at = 1:length(chr_names), labels = chr_names, las = 1)
  
  for(chr in seq_along(haps)) {
    
    # Get the haplotype blocks on each chromosome.
    blocks = apply(haps[[chr]], 2, rle)
    
    # Haplotype 1
    len    = markers[[chr]][cumsum(blocks$hap1$lengths)]
    fndr   = blocks$hap1$values
    rect(c(0, len[-length(len)]), chr - chr_width, len, chr, col = colors[fndr])
    
    # Haplotype 2
    len    = markers[[chr]][cumsum(blocks$hap2$lengths)]
    fndr   = blocks$hap2$values
    rect(c(0, len[-length(len)]), chr, len, chr + chr_width, col = colors[fndr])
    
  } # for(chr)
  
} # plot_genome()