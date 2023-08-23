# Turn the typical mass header in the pmsf (for example "99.08079_99.0806600670743_24.75") into numerical masses.
## function

extract_mean_mass <- function(x,decimals = 6){
    bounds <- strsplit(sub("^.", "", x), "_")[[1]]
    return(round(mean(as.numeric(bounds[1:2])),decimals))
}


# Extract time from the mass header in the pmsf data (for example "99.08079_99.0806600670743_24.75") into numerical masses.
extract_pmsf_time <- function(x){
    bounds <- strsplit(x, "_")[[1]]
    return(as.numeric(bounds[3])) # working under the assumption that the third digit is the time value
}

# count the number of rows in a data frame like object with no zeros
# The way for this to work is to take in an arbritary number of vector-like objects of the same lenght and make a data table like object out of them and count the number of rows with non zero objects in them

# TODO add a way to make sure that all the vectors are of the same length
count_zeroes <- function(...,x = 0){
    current_table <- as.data.table(cbind(...))
    count_of_zeros <- rowSums(current_table == x)
    unique_table <- table(count_of_zeros)
    return(unique_table)
}