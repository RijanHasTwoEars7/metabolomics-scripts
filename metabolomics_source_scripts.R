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
tabulate_zeroes <- function(...,x = 0){
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

meta_data <- function(feature = "", data_table = data.table(),feature_array = "mass_features",time_array = "feature_time", mass_array = "mean_masses"){
    
    metadata <- list()
    metadata$mean_mass <- data_table[data_table[[feature_array]] == feature,][[mass_array]]
    metadata$time <- data_table[data_table[[feature_array]] == feature,][[time_array]]
    
    return(metadata)
    
}

# individual heritability of two arrays worth of data against the pmsf genotype and treatment data and their mixed heritability
# please pay attention to the data that need to be supplied

best_heritability_augumentation <- function(feature_one = list(),feature_two = list(),count_zeroes = TRUE,treatment_data = list(),genotype_data = list()){
    
    proto_mixed_set <- as.data.table(cbind(feature_one,feature_two)) 

    mixed_data <- apply(proto_mixed_set,1,function(x) {
        if(length(which(x == 0)) == 1){
            return(sum(x))
        }
        else{
            return(mean(unlist(x)))
        }
    })

    values <- list()

    if(count_zeroes){
        values$feature_one_zeroes <- length(which(feature_one == 0))
        values$feature_two_zeroes <- length(which(feature_two == 0))
        values$combined_zeroes <- length(which(mixed_data == 0))
    }

    values$feature_one_heritability <- anova_vals(daltons = feature_one,treatment = treatment_data,genotype = genotype_data)
    values$feature_two_heritability <- anova_vals(daltons = feature_two,treatment = treatment_data,genotype = genotype_data)
    values$mixed_heritability <- anova_vals(daltons = mixed_data,treatment = treatment_data,genotype = genotype_data)

    return(values)
}

# Is.salt() function implementation
# This funciton will use the formula discussed in the following paper  
# https://www.researchgate.net/publication/307894872_Post-acquisition_filtering_of_salt_cluster_artefacts_for_LC-MS_based_human_metabolomic_studies

is.salt = function(x){
    x = as.numeric(x)

    integer = floor(x)

    diff = abs(x - integer)
    
    litmus = 0.00112*x + 0.01953

    if (diff > litmus){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}

# Function to count zeroes in a array style object

count_items = function(data_list = list(), item = NULL){
    
    counts = length(which(data_list == item))

    return(counts)
}

# QOL funciton for ease of use
range_of_all = function(...){
    return(range(cbind(...)))
}

# The funciton to expand formulas

formula_expander <- function(formula){
  expanded_formula <- list()
  molecules <- list()
  functional_groups <- strsplit(formula,"[()]")
  for(functional_group in functional_groups){
    typeof(unlist(strsplit(functional_group, "(?<=[a-z])(?=[A-Z])", perl=TRUE)))
    molecules <- c(molecules,unlist(strsplit(gsub("(?<=.)(?=[A-Z])", " ", functional_group, perl=TRUE), " ")))
  }
  for(i in molecules){
    if(grepl("[0-9]",i)){
      multiplier <- as.numeric(gsub("[^0-9]", "", i))
      element <- gsub("[^[:alpha:]]", "", i)
      expanded_formula <- c(expanded_formula,rep(element, multiplier))
    } else{
      expanded_formula <- c(expanded_formula, i)
    }
  }
  return(paste(expanded_formula, collapse = ""))
}

# The funciton that takes in the molecular formulas and makes matrices

make_matrix <- function(received_molecule) {

    if(received_molecule == 'NA') {
        return('NA')
    } else{
        current_molecule_list <- list()
        for (i in myElements){
            current_molecule <- table(strsplit(received_molecule,""))
            current_molecule_list[i] <- tryCatch(current_molecule[[i]], error = function(e) 0)

        }
        current_molecules_matrix <- matrix(unlist(current_molecule_list), ncol=1, byrow=TRUE)
        rownames(current_molecules_matrix) <- c(myElements)
    }
    return(as.vector(current_molecules_matrix))
}