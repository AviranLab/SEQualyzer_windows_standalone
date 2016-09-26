#1. Define ranges for low, medium and high reactivities
#2. Calculate requirement factor for reactivities in all three ranges
#3. In every range, find 95 percentile as the predicted index
#4. If a single index is desired, use the index for medium range of reactivities.

CQI <- function(reads, num_replicates, confidence, percent_var) { #, single_end= FALSE) {
  low_range <- c(0.1, 0.3)
  medium_range <- c(0.3, 0.7)
  high_range <- 0.7
  #what percentile to take for giving indices
  index_criteria <- 0.95
  min_reactivity <- 0.10
  z_value <- qnorm(1- ((1-confidence)/2))
  
  n <- nrow(reads)
  indices <- matrix(, num_replicates, 4)
  colnames(indices) <- c("Replicate", "Low", "Medium", "High")
  
  for (j in 1:num_replicates) {
    eval(parse(text= paste("minus.counts <- matrix(reads$minus.counts.", j,", nrow=n, ncol=1)", sep="")))
    eval(parse(text= paste("minus.pass <- matrix(reads$minus.pass.", j,", nrow=n, ncol=1)", sep="")))
    eval(parse(text= paste("plus.counts <- matrix(reads$plus.counts.", j,", nrow=n, ncol=1)", sep="")))
    eval(parse(text= paste("plus.pass <- matrix(reads$plus.pass.", j,", nrow=n, ncol=1)", sep="")))
    
    supply_minus <- sum(minus.counts, na.rm=T)
    supply_plus <- sum(plus.counts, na.rm=T)
    size_factor <- supply_minus/supply_plus
    
    #Minus channel noise
    g <- matrix(minus.counts/(minus.counts + minus.pass), nrow=n, ncol=1)
    
    #Signal
    b <- matrix(local.coverage.react.calc(minus.counts, minus.pass, plus.counts, plus.pass, 1, 0, 0, 1),
                nrow= n, ncol=1)
    
    per_residue_info <- matrix(plus.counts + plus.pass, nrow=n, ncol=1)
    per_residue_req <- matrix(ceiling((((b + g - (b * g)) * (1- (b + g - (b * g)))) + (g*(1-g)/size_factor))/((percent_var*b/z_value)^2 )),
                              nrow= n, ncol= 1)
    
    norm <- get.normalizer(b, 1)
    index_low_reactivity <- which(b/norm >= low_range[1] & b/norm < low_range[2])
    index_medium_reactivity <- which(b/norm >= medium_range[1] & b/norm < medium_range[2])
    index_high_reactivity <- which(b/norm >= high_range)
    
    #req_factor
    rf <- per_residue_req/per_residue_info
    
    rf_low <- rf[index_low_reactivity]
    rf_medium <- rf[index_medium_reactivity]
    rf_high <- rf[index_high_reactivity]
    
    sorted_low <- rf_low[order(rf_low)]
    sorted_medium <- rf_medium[order(rf_medium)]
    sorted_high <- rf_high[order(rf_high)]
    
    req_factor_low <- quantile(sorted_low, index_criteria)
    req_factor_medium <- quantile(sorted_medium, index_criteria)
    req_factor_high <- quantile(sorted_high, index_criteria)
    
    indices[j, ] <- c(j, req_factor_low, req_factor_medium, req_factor_high)
  }
  
  return(indices)
}
