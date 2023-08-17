# Turn the typical mass header in the pmsf (for example "99.08079_99.0806600670743_24.75") into numerical masses.
## function

extract_mean_mass <- function(x,decimals = 6){
    bounds <- strsplit(x, "_")[[1]]
    return(round(mean(as.numeric(bounds[1:2])),decimals))
}
