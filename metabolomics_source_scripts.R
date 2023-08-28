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

# the setup to calculate the most likely formulas for a given mass

myElements <- c("C","H","N","O","S","P")

best_fit_mf <- function(mean_mass){
  df <- calcMF(mz = mean_mass, z = 1, ppm = 5, top = NULL,
  elements = Rdisop::initializeElements(myElements), maxCounts = TRUE,
  SeniorRule = TRUE, HCratio = TRUE, moreRatios = TRUE,
  elementHeuristic = TRUE, Filters = list(DBErange = c(-5, 40)
  , minElements = "C0H0N0S0O0", maxElements = "C40H60N20S20O20", parity = "e", maxCounts =
  TRUE, SENIOR3 = 0, HCratio = TRUE, moreRatios = TRUE, elementHeuristic =
  TRUE), summarize = FALSE, BPPARAM = NULL)

  if(is.null(df)){
    return('NA')
  }else{
    return(df[order(abs(df$ppm)),][1,2])
  }

}

# heritability for a given data set, please pay attention to the data that needs to be supplied

anova_vals <- function(daltons = list(),treatment = list(),genotype = list(), drop_zeros = FALSE ){

    data_table <- as.data.table(cbind(daltons,treatment,genotype))

    data_table$daltons <- as.numeric(data_table$daltons)

    if(drop_zeros){
        data_table <- data_table[data_table$daltons != 0,]
    }

    model <- lm(daltons~genotype + treatment + genotype*treatment, data = data_table)

    anov_data <- anova(model)
    
    variance_data <- sum(anov_data[["Sum Sq"]][1:3]) / sum(anov_data[["Sum Sq"]][1:4])
    return(variance_data)
}

# The function to extract the specific date and time for a feature

# this funciton depends on the user supplying very specific params for the 

meta_data <- function(feature = "", data_table = feature_metadata,feature_array = "mass_features",time_array = "feature_time", mass_array = "mean_masses"){
    
    metadata <- list()
    metadata$mean_mass <- data_table[data_table[[feature_array]] == feature,][[mass_array]]
    metadata$time <- data_table[data_table[[feature_array]] == feature,][[time_array]]
    
    return(metadata)
    
}