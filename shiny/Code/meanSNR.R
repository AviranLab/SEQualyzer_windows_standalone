#TW mean SNR
meanSNR <- function(dat, session_name, est.flag, nor.flag, ln.flag, zero.flag, reps,
                    loc_cov) {
  if (!.Platform$OS.type == "windows") registerDoParallel(cores=detectCores(all.tests=TRUE))
  
  b.normalize <- function(raw.estimates) {
    outlier.limit <- 0
    non.zero.only <- 0
    
    sorted.estimates <- raw.estimates[order(raw.estimates)]
    box.data <- boxplot(sorted.estimates, plot = FALSE)
    trimmed <- sorted.estimates[! sorted.estimates %in% box.data$out]
    trimmed <- trimmed[!is.na(trimmed) & !is.infinite(trimmed)]
    if (non.zero.only == 1) {
      box.data.pos <- boxplot(sorted.estimates[sorted.estimates > 0], plot = F)
      trimmed <- sorted.estimates[! sorted.estimates %in% box.data.pos$out]
      trimmed <- trimmed[!is.na(trimmed) & !is.infinite(trimmed)] 
    }
    
    #optional 10% max threshold for outliers
    outlier.threshold <- .12
    if (outlier.limit == 1) {
      if (length(box.data$out) > (length(raw.estimates) * outlier.threshold)) {
        trimmed <- sorted.estimates[1:round(length(sorted.estimates) * 
                                              (1 - outlier.threshold))]
      }
    }
    normalize.range <- c(round(length(trimmed) * .9), length(trimmed) * 1)
    normalizer <- mean(trimmed[normalize.range[1]:normalize.range[2]], na.rm = TRUE)
    normalized.reactivity <- raw.estimates / normalizer
    return(normalized.reactivity)
  }
  
  t.e.normalize <- function(raw.estimates) {
    sorted <- raw.estimates[order(raw.estimates)]
    if (any(is.na(sorted))) {
      normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))    
    } else {
      normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
    }
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
    return(raw.estimates/normalizer)
  }
  
  log_m <- function(x) {
    x[x== 0] <- 1
    return(log(x))
  }
  
  if (is.null(dat$dat_tw$mean_SNR)) {
    data <- dat$dat
    registerDoParallel(cores=detectCores(all.tests=TRUE))
    
    # calculate replicate SNR
    mean.SNR <- foreach(i =1:length(data), .combine= c, .multicombine= TRUE) %dopar% {
      r <- matrix(rep(0, nrow(data[[i]])*reps), nrow(data[[i]]), reps)
      if (loc_cov) {
        for (rep in 1:reps) {
          
          minus.counts <- data[[i]][, paste0("minus.counts.", rep)]
          plus.counts <- data[[i]][, paste0("plus.counts.", rep)] 
          minus.pass <- data[[i]][, paste0("minus.pass.", rep)]
          plus.pass <- data[[i]][, paste0("plus.pass.", rep)]
          
          if (ln.flag== 1) {
            minus.counts[minus.counts == 0] <- 1
            plus.counts[plus.counts == 0] <- 1
            minus.pass[minus.pass == 0] <- 1
            plus.pass[plus.pass == 0] <- 1
            
            minus.counts <- log(minus.counts)
            plus.counts <- log(plus.counts)
            minus.pass <- log(minus.pass)
            plus.pass <- log(plus.pass)
          }
          
          minus.rate <- minus.counts/minus.pass
          plus.rate <- plus.counts/plus.pass
          
          if (ln.flag == 2) {
            minus.rate <- log_m(minus.rate)
            plus.rate <- log_(plus.rate)
          }
          
          if (est.flag == 1) {react <- plus.rate - minus.rate}
          if (est.flag == 2) {react <- plus.rate / minus.rate - 1}
          
          #weird bug, sometimes output is a 1 column df. This fixes it
          if (length(dim(react)) > 0) {react <- as.vector(react[, 1])}
          
          if(zero.flag == 1) {react[react < 0] <- 0}
          
          if (nor.flag == 1) {react <- t.e.normalize(react)}
          if (nor.flag == 2) {react <- b.normalize(react)}
          r[, rep] <- react
        }
        
      } else {
        for (rep in 1:reps) {
          
          minus.counts <- if (ln.flag == 1) log_m(data[[i]][, paste0("minus.counts.", rep)]) else data[[i]][, paste0("minus.counts.", rep)]
          plus.counts <- if (ln.flag==1) log_m(data[[i]][, paste0("plus.counts.", rep)]) else data[[i]][, paste0("plus.counts.", rep)]
          
          minus.rate <- minus.counts/mean(minus.counts, na.rm=T)
          plus.rate <- plus.counts/mean(plus.counts, na.rm=T)
          
          if (ln.flag == 2) {
            minus.rate <- log_m(minus.rate)
            plus.rate <- log_(plus.rate)
          }
          
          if (est.flag == 1) {react <- plus.rate - minus.rate}
          if (est.flag == 2) {react <- plus.rate / minus.rate - 1}
          
          #weird bug, sometimes output is a 1 column df. This fixes it
          if (length(dim(react)) > 0) {react <- as.vector(react[, 1])}
          
          if(zero.flag == 1) {react[react < 0] <- 0}
          
          if (nor.flag == 1) {react <- t.e.normalize(react)}
          if (nor.flag == 2) {react <- b.normalize(react)}
          
          r[, rep] <- react
          
        }
      }
      
      mean(apply(r, 1, function(x) min(mean(x)/sd(x), 35)), na.rm=T)
    }
    
    dat$dat_tw$mean_SNR <- mean.SNR
    save(dat, file= file.path(getwd(), "Saved_results", paste("data_", session_name, 
                                                              "_with_mean_SNR.RData", sep="") ))
  } else mean.SNR = dat$dat_tw$mean_SNR
  
  return(mean.SNR)
}