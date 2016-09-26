#----------------------------------------------------------------------------------------
#A function to generate minus channel and plus channel fragment counts
#Requires frequency of k-fragments in minus and plus channel and target sample sizes for paired end reads.
#Requires stops and passes for minus and plus channel and sample sizes for single end reads.
#----------------------------------------------------------------------------------------
SampleGenerator <- function (freq, sampleSize) {
  
  n <- length(freq$minus.counts)
  freq$minus
  sample_gen <- list(minus.counts="", plus.counts="")
  sample_gen$minus.counts <- tabulate(sample(n, sampleSize$minus.counts, rep= TRUE, 
                                             prob= replace(freq$minus.counts, which(is.na(freq$minus.counts)), 0 )), nbins=n)
  sample_gen$plus.counts <- tabulate(sample(n, sampleSize$plus.counts, rep= TRUE, 
                                            prob= replace(freq$plus.counts, which(is.na(freq$plus.counts)), 0 )), nbins=n)
  return(sample_gen)
}


#----------------------------------------------------------------------------------------
#A function for bootstrapping data.
#----------------------------------------------------------------------------------------
bootstrap <- function(num_bootstrapSamples,num_biological_replicates, reads, single_end= FALSE) {
  
  registerDoParallel(cores=detectCores(all.tests=TRUE))
  #a function to tell foreach how to combine lists of kind data as they are generated in parallel
  combine <- function(...) {
    
    all_lists <- list(...)
    list.combined <- list()
    
    list.combined$minus.counts <- sapply(all_lists, '[[', 'minus.counts')
    list.combined$plus.counts <- sapply(all_lists, '[[', 'plus.counts')
    return(list.combined)
  }
  
  SampleGenerator <- function (freq, sampleSize) {
    
    n <- length(freq$minus.counts)
    freq$minus
    sample_gen <- list(minus.counts="", plus.counts="")
    sample_gen$minus.counts <- tabulate(sample(n, sampleSize$minus.counts, rep= TRUE, 
                                               prob= replace(freq$minus.counts, which(is.na(freq$minus.counts)), 0 )), nbins=n)
    sample_gen$plus.counts <- tabulate(sample(n, sampleSize$plus.counts, rep= TRUE, 
                                              prob= replace(freq$plus.counts, which(is.na(freq$plus.counts)), 0 )), nbins=n)
    return(sample_gen)
  }
  
  boot <- list()
  
  for (j in 1:num_biological_replicates) {
    
    freq <- list(minus.counts="", plus.counts="")
    sampleSize <- list(minus.counts="", plus.counts="")
    
    sampleSize$minus.counts <- round(sum(reads[, paste0("minus.counts.", j)], na.rm=T))
    sampleSize$plus.counts <- round(sum(reads[, paste0("plus.counts.", j)], na.rm=T))
    
    freq$minus.counts <- reads[, paste0("minus.counts.", j)]/sampleSize$minus.counts
    freq$plus.counts <- reads[, paste0("plus.counts.", j)]/sampleSize$plus.counts
    
    boot_current <- foreach(i=1:num_bootstrapSamples , .combine = combine, .multicombine=TRUE) %dopar% SampleGenerator(freq, sampleSize)
    eq <- paste("boot","$","rep", j, " <- ", "boot_current",sep="" )
    eval(parse(text= eq))
    
    if (!single_end) {
      boot[[paste0("rep", j)]]$minus.pass <- reads[[paste0("minus.pass.", j)]]
      boot[[paste0("rep", j)]]$plus.pass <- reads[[paste0("plus.pass.", j)]]
    }
    
    rm(boot_current)
  }
  return(boot)
}


#----------------------------------------------------------------------------------------
#Get reactivities from a bootstrapped sample
#----------------------------------------------------------------------------------------
getReactivities <- function(boot, est.flag, norm.flag, log.flag, zero.flag, single_end= FALSE) {
  
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
  
  num_bio_reps <- length(names(boot))
  reactivities <- list()
  num_resamples <- ncol(boot[[1]]$plus.counts)
  
  for(j in 1:num_bio_reps) {
    eval(parse(text= paste("reactivities$rep", j, " <- cbind()", sep="" )))
  }
  
  for (j in 1:num_bio_reps) {
    if (single_end) {
      eval(parse(text= paste("reactivities$rep", j, " <- foreach(r= 1:num_resamples, .combine=cbind, .multicombine= TRUE) %dopar% react.calculator(boot[[j]]$minus.counts[, r], boot[[j]]$plus.counts[, r],
                           est.flag, norm.flag, log.flag, zero.flag) ", sep="")))
      
    } else {
      eval(parse(text= paste("reactivities$rep", j, " <- foreach(r= 1:num_resamples, .combine=cbind, .multicombine= TRUE) %dopar% local.coverage.react.calc(boot[[j]]$minus.counts[, r], boot[[j]]$minus.pass, boot[[j]]$plus.counts[, r], boot[[j]]$plus.pass,
                           est.flag, norm.flag, log.flag, zero.flag) ", sep="")))
    }
  }
  
  return(reactivities)
}

#----------------------------------------------------------------------------------------
#A function to get nucleotide level SNR (bootstrap).
#----------------------------------------------------------------------------------------
SNRcalc <- function(reactivities, clip_at =35) {
  
  num_bio_reps <- length(reactivities)
  SNR <- matrix(, nrow=nrow(reactivities[[1]]), 1+num_bio_reps)
  colnames(SNR) <- c("Nucleotide", paste("Replicate", c(1: num_bio_reps), sep="_"))
  SNR[,1] <- c(1:nrow(reactivities[[1]]))
  
  for (j in 1:num_bio_reps) {
    SNR[,j+1] <- apply(reactivities[[j]], 1, mean, na.rm=TRUE)/ apply(reactivities[[j]], 1, sd, na.rm=TRUE)
    SNR[which(SNR[, j+1] > clip_at) ,j+1] <- clip_at
  }
  
  return(SNR)
}

#----------------------------------------------------------------------------------------
#A function to get nucleotide level mean reactivity (bootstrap).
#----------------------------------------------------------------------------------------
mean.calc <- function(reactivities) {
  
  num_bio_reps <- length(reactivities)
  mean_rea <- matrix(, nrow=nrow(reactivities[[1]]), 1+num_bio_reps)
  colnames(mean_rea) <- c("Nucleotide", paste("Replicate", c(1: num_bio_reps), sep="_"))
  mean_rea[,1] <- c(1:nrow(reactivities[[1]]))
  
  for (j in 1:num_bio_reps) {
    mean_rea[,j+1] <- apply(reactivities[[j]], 1, mean, na.rm=TRUE)
  }
  
  return(mean_rea)
}

#----------------------------------------------------------------------------------------
#A function to get nucleotide level mean reactivity (bootstrap).
#----------------------------------------------------------------------------------------
sd.calc <- function(reactivities, as_numeric= TRUE) {
  
  num_bio_reps <- length(reactivities)
  if (as_numeric) {
    mean_rea <- c()
    for (j in 1:num_bio_reps) {
      mean_rea<- c(mean_rea, apply(reactivities[[j]], 1, sd, na.rm=TRUE))
    }
    mean_rea <- data.frame(sdev= mean_rea)
  } else {
    mean_rea <- matrix(, nrow=nrow(reactivities[[1]]), 1+num_bio_reps)
    colnames(mean_rea) <- c("Nucleotide", paste("Replicate", c(1: num_bio_reps), sep="_"))
    mean_rea[,1] <- c(1:nrow(reactivities[[1]]))
    
    for (j in 1:num_bio_reps) {
      mean_rea[,j+1] <- apply(reactivities[[j]], 1, sd, na.rm=TRUE)
    }
  }
  
  return(mean_rea)
}

#----------------------------------------------------------------------------------------
#A function to get window level SNR (bootstrap).
#----------------------------------------------------------------------------------------
SNR.scanner <- function(ind.snr, window.size) {
  num_bio_reps <- ncol(ind.snr)-1
  n <- nrow(ind.snr)
  half.window <- round(window.size / 2)
  
  window.snr <- matrix(, nrow= n,ncol= 1+num_bio_reps)
  colnames(window.snr) <- c("Nucleotide", paste("Replicate", 1:num_bio_reps, sep="_"))
  window.snr[,1] <- 1:n
  for (i in 1: num_bio_reps) {
    #change na to clip_at
    ind.snr[is.infinite(ind.snr[,i+1]), i+1] <- NA
    for (j in window.snr[,1]) {
      window.snr[j, i+1] <- if (j < half.window | j > n-half.window) NA else mean(ind.snr[(j - half.window) : (j + half.window), i+1], na.rm= TRUE)
    }
  }
  return(data.frame(window.snr))
}

#----------------------------------------------------------------------------------------
#A function to get SNR (formula based).
#----------------------------------------------------------------------------------------
formula_sd_calc <- function(reactivity, noise, norm, data, est) {
  num_bio_reps <- ncol(noise)-1
  sd <- matrix(, nrow= nrow(data),ncol= 1+num_bio_reps)
  colnames(sd) <- c("Nucleotide", paste("Replicate", 1:num_bio_reps, sep="_"))
  sd[,1] <- 1:nrow(data)
  
  for (j in 1:num_bio_reps) {
    eval(parse(text= paste("plus.pass <- data$plus.pass.",j, sep="" )))
    eval(parse(text= paste("minus.pass <- data$minus.pass.",j, sep="" )))
    b <- reactivity[, j+1]
    g <- noise[, j+1]
    
    if (est == 1) {
      #formula for difference
      sd[, j+1] <- (((b + (g*(1-b)))*(1-(b + (g*(1-b))))/plus.pass) + ((g*(1-g)/minus.pass)))^0.5
    } else {
      #formula for ratio
      sd[,j+1] <- (((b + (g*(1-b)))*(1-b-(g*(1-b)))/(g*g*plus.pass)) + (((b + (g*(1-b)))^2)*(1-g)/(g*g*g*minus.pass)))^0.5
    }
    
    sd[, j+1] <- sd[, j+1]/norm[1, j]
    
  }
  
  return(data.frame(sd))
}

#----------------------------------------------------------------------------------------
#A function to get bootstrap mean SNR.
#----------------------------------------------------------------------------------------
bootSNR <- function(data, loc_cov, est.flag, nor.flag, log.flag, zero.flag, clip_at=35) {
  boot <- bootstrap(100, 1, data, !loc_cov)
  r <- getReactivities(boot, est.flag, nor.flag, log.flag, zero.flag, !loc_cov)$rep1
  snr <- round(mean(apply(r, 1, function(x) min(mean(x, na.rm=T)/sd(x, na.rm=T), clip_at)), na.rm=T), 1)
  return(snr())
}