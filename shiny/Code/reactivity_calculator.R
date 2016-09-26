local.coverage.react.calc <- function(minus.stops, minus.coverage, plus.stops, plus.coverage,
                                      estimate.flag, norm.flag, log.flag, zero.flag) {
  
  # estimate flags: 1 = difference, 2 = ratio
  # norm.flag: 1 = 2%-8% normalization, 2 = boxplot normalization
  # log.flag: 1 = log counts and pass individually (before summing), 2 = log rates
  # zero.flag: 0 = keep negative reactivities, 1 = set negative reactivities to zero
  
  # log flag = 1 counts become log pseudocounts
  if (log.flag == 1) {
    minus.stops <- log.maker(minus.stops)
    minus.coverage <- log.maker(minus.coverage)
    plus.stops <- log.maker(plus.stops)
    plus.coverage <- log.maker(plus.coverage)
  }
  
  #  for SHAPE-Seq
  minus.stoprate <- minus.stops / minus.coverage
  plus.stoprate <- plus.stops / plus.coverage
  
  
  
  # single-end coverage normalization
  
  if (log.flag == 2) {
    minus.stoprate <- log.maker(minus.stoprate)
    plus.stoprate <- log.maker(plus.stoprate)
  }
  
  if (estimate.flag == 1) {react <- plus.stoprate - minus.stoprate}
  if (estimate.flag == 2) {react <- plus.stoprate / minus.stoprate - 1}
  
  #weird bug, sometimes output is a 1 column df. This fixes it
  if (length(dim(react)) > 0) {react <- as.vector(react[, 1])}
  
  if(zero.flag == 1) {react[react < 0] <- 0}
  
  if (norm.flag == 1) {react <- two.eight.normalize(react)}
  if (norm.flag == 2) {react <- boxplot.normalize(react)}
  
  return(react)
  
}

log.postprocess <- function(react) {
  real.react <- react[!is.infinite(react) & !is.na(react)]
  b <- min(real.react)
  new.react <- react - b
  new.react[is.infinite(react) | is.na(react)] <- 0
  return(new.react)
}


log.coverage.creater <- function(counts.df) {
  # just sum the stops of every stop after each desired residue, but log first
  counts.df$minus.counts <- log.maker(counts.df$minus.counts)
  counts.df$plus.counts <- log.maker(counts.df$plus.counts)
  
  counts.df$minus.pass <- 0
  counts.df$plus.pass <- 0
  r <- length(counts.df$minus.pass)
  for (i in 2:r) {
    counts.df$minus.pass[i] <- sum(counts.df$minus.counts[1:i - 1])
    counts.df$plus.pass[i] <- sum(counts.df$plus.counts[1: i - 1])
  }
  return(counts.df)
}



log.maker <- function(counts) {
  mod.counts <- counts
  mod.counts[mod.counts == 0] <- 1
  counts.output <- log(mod.counts)
  return(counts.output)
}


react.calculator <- function(minus, plus, estimate.flag, norm.flag, log.flag, zero.flag) {
  
  # Flags: 
  # log.flag: 0 = no log, 1 = log 
  # estimate.flag: 1 = Ding Formula, 2 = Talkish Formula
  # norm.flag: 1 = 2% 8% norm, 2 = boxplot norm
  # zero.flag: 1 = Shift to zero
  
  data.set <- data.frame(minus, plus)
  
  if (log.flag == 1) {
    data.set$minus <- log.maker(data.set$minus)
    data.set$plus <- log.maker(data.set$plus)
  }
  
  if (estimate.flag == 1) {
    estimate <- diff.MR(data.set$minus, data.set$plus)
  }
  
  if (estimate.flag == 2) {
    estimate <- ratio.MR(data.set$minus, data.set$plus)
  }
  
  if (zero.flag == 1) {estimate[estimate < 0] <- 0}
  
  if (log.flag == 2) { 
    estimate <- log.postprocess(estimate)
  }
  
  if (norm.flag == 1) {estimate<- two.eight.normalize(estimate)}
  if (norm.flag == 2)  {estimate <- boxplot.normalize(estimate)}  
  
  return(estimate)
}

diff.MR <- function(minus, plus) {
  # input needs to be either raw counts or already logged
  minus.scaled <- minus / (mean(minus, na.rm=T))
  plus.scaled <- plus / (mean(plus, na.rm=T))
  reactivity <- plus.scaled - minus.scaled
  return(reactivity)
}

ratio.MR <- function(minus, plus) {
  
  minus.scaled <- minus / sum(minus, na.rm=T)
  plus.scaled <- plus / sum(plus, na.rm=T)
  reactivity <- plus.scaled /minus.scaled # beware of log(1) == 0
  return(reactivity)
}

two.eight.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  if (any(is.na(sorted))) {
    normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))    
  } else {
    normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
  }
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  return(raw.estimates/normalizer)
}

boxplot.normalize <- function(raw.estimates) {
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



get.normalizer <- function(raw.estimates, flag) {
  if (flag ==0) normalizer <- 1
  
  if (flag == 1) {
    sorted <- raw.estimates[order(raw.estimates)]
    if (any(is.na(sorted))) {
      normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))    
    } else {
      normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
    }
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  if (flag == 2) {
    sorted.estimates <- raw.estimates[order(raw.estimates)]
    box.data <- boxplot(sorted.estimates, plot = FALSE)
    trimmed <- sorted.estimates[! sorted.estimates %in% box.data$out]
    
    #optional 10% max threshold for outliers
    outlier.limit <- 0
    outlier.threshold <- .1
    if (outlier.limit == 1) {
      if (length(box.data$out) > (length(raw.estimates) * outlier.threshold)) {
        trimmed <- sorted.estimates[1:round(length(sorted.estimates) * .9)]
      }
    }
    normalize.range <- c(round(length(trimmed) * .9), length(trimmed) * 1)
    normalizer <- mean(trimmed[normalize.range[1]:normalize.range[2]], na.rm = TRUE)
  }
  
  if (flag == 3) {
    sorted <- raw.estimates[order(raw.estimates)]
    normalize.range <- c(round(length(sorted) * .85), round(length(sorted) * .95))
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 4) {
    sorted <- raw.estimates[order(raw.estimates)]
    normalize.range <- c(round(length(sorted) * .80), round(length(sorted) * .90))
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 5) {
    sorted <- raw.estimates[order(raw.estimates)]
    normalize.range <- c(round(length(sorted) * .90), round(length(sorted) * 1))
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 6) {
    sorted <- raw.estimates[order(raw.estimates)]
    normalize.range <- c(round(length(sorted) * .60), round(length(sorted) * .7))
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  return(normalizer)
}