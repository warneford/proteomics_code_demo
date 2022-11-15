# RBRID Functions 
# Author: Robert Warneford-Thomson
# 
# Update 04Aug2016 < Corrected RBR-ID score to square p-value. Added Flag normalization. Added Peptideplot to compare normalization
# Update 17Aug2016 <- Added code to plot ideal peptide coverage in Residuescoreplot
# Update 18Aug2016 <- Added function pval to calculate p values
# Update 19Aug2016 <- Added IPR domains and custom domain plotting functionality to residuescoreplot
# Update 22Aug2016 <- Added F_statistic function
# Update 30Aug2016 <- Added peptide highlight functionality to Residuescoreplot
# Update 18Sep2016 <- Updated/debugged NormRBR
# Update 03Oct2016 <- Added highlight range functionality to ResiduescorePlot.
# Update 19Sep2018 <- Updated coverage.lapply to sync with Bioconductor 3.7
# Update 15Jan2021 <- reordered factors in OverlayScorePlot to ensure experiment labels aren't switched
# Update 03Feb2022 <- added median normalization to NormRBR
# Update 31Mar2022 <- updated NormRBR to perform fisher test instead of t-test for peptides with no observations in either -/+ 4SU condition

## NormRBR

### temp code from debugging, disregard
# rbrid <- df
# Normalization = "median"
# convert.zero.means = T
# p.adjust = T
# paired = F
# var.equal = T
# control_string = "minus"
# test_string = "plus"
#####

NormRBR <-function(rbrid = Prot_table,
                   Normalization = "mean",
                   paired = F,
                   var.equal = T,
                   control_string = "minus",
                   test_string = "plus",
                   convert.zero.means = F,
                   p.adjust = F)
  # last updated 14 Apr 2022
  ## Function to normalize and analyze RBR peptide intensities. -4su and +4sU replicates are determined by finding columns containing "minus" and "plus" respectively 
  # Normalization method "none" performs no normalization, when you wish to normalize manually, minus and plus 4su columns must contain "minus.*Norm" and "plus.*Norm" strings
  # Normalization method "signal" normalizes all peptide intensities to total signal peptide, then scales by mean
  # Normalization method "mean" normalizes all peptide intensities to mean peptide
  # Normalization method "median" normalizes all peptide intensities to median peptide
  # Normalization method "global_median" normalizes all peptide intensities to median peptide, then multiplies by median intensity observed across all replicates
  # Normalization method "mol.cell.2016" normalizes all peptide intensities by total signal, do not remove contaminants
  # Normalization method "protein" normalizes all peptide intensities by dividing by sum of peptide intensities for each protein
  # Normalization method "silac" normalizes peptide intensities by median, then calculates ratios for each sample light and heavy
  # Normalization method "silac_pseudo" uses same normalization as above, but inserts pseudocounts to mitigate sparse data
  # intensity in that run
  # silac_pseudocount refers to what values to replace NA values with when no other peptide intensities are available
  # "col_minimum" will replace NA with minimum in column, "silac_pseudocount = "protein" will replace with minimum value for same protein
  ## If a paired test is desired all vectors must have same number of test and control values (i.e. no NA's).
  # p.adjust flag specifies whether to correct for FDR via B-H method (only corrects non-missing peptides)
{
  if(Normalization == "silac" | Normalization == "silac_pseudo") {
    silac_flag = T
  } else { silac_flag = F}
  
  pacman::p_load(tidyverse,
                  purrr)
  
  #data columns must contain control or test string and end with "t\\d" eg. 't1'
  if (control_string == "minus" & test_string == "plus") {
  minusCols <- grep(paste0(control_string, ".*[tb]\\d+.*"), names(rbrid), value = T)
  plusCols <- grep(paste0(test_string, ".*[tb]\\d+.*"), names(rbrid), value = T)
  } else { #use arbitrary grep strings to match test and control columns
    minusCols <- grep(control_string, names(rbrid), value = T)
    plusCols <- grep(test_string, names(rbrid), value = T)
}
  
  Samplecols <- c(minusCols, plusCols)
  
  #for silac samples, check that minus and plus cols are paired properly
  if (silac_flag == T) {
    mc <- map_chr(minusCols, ~gsub("(light|heavy)_(minus|plus)", "", .))
    pc <- map_chr(plusCols, ~gsub("(light|heavy)_(minus|plus)", "", .))
    if(any(mc != pc)) {
      stop(paste0("Minus and Plus columns cannot be paired properly for silac ratio calculation\n",
                  "Minus columns: \n",
                  paste0(minusCols, collapse = " "),
                  "\nPlus columns: \n",
                  paste0(plusCols, collapse = " ")))
                  
    }
  }
  
  #ungroup
  rbrid <- rbrid %>%
    ungroup() %>%
    as_tibble()
  
  #convert zero to NA (unless you need the zero's for paired t-test)
  if (Normalization != "mol.cell.2016" | paired != T) {
    rbrid[,Samplecols][rbrid[,Samplecols] == 0] <- NA
  }
  
  # Convert Contaminants to logical vector
  if (any(grepl("contaminant", colnames(rbrid))) ==T) {  
    rbrid$contaminant<-as.logical(rbrid$contaminant)
  } else {
    rbrid$contaminant <- F
    
  }
  
  if (Normalization == "mol.cell.2016") {
    # Normalize peptide intensity by Sum of total intensity (not removing contaminants)
    rbrid <- rbrid %>% 
      mutate(
        across(all_of(Samplecols), ~ .x / sum(.x, na.rm = T),
               .names = "{.col}_Norm")
      )
  }
  
  if (Normalization == "signal") {
    # Normalize peptide intensity by Sum of total intensity (removing contaminants)
    # then multiple all values by mean intensity (to make numbers more human friendly)
    rbrid <- rbrid %>% 
      mutate(
        across(all_of(Samplecols), 
               ~ (.x * mean(.x[!contaminant], na.rm = T) ) / sum(.x[!contaminant], na.rm = T),
               .names = "{.col}_Norm")
      )
  }
  
  if (Normalization == "median" | silac_flag == T) {
    # Normalize peptide intensity by median intensity (removing contaminants) observed in replicate
    rbrid <- rbrid %>% 
      mutate(
        across(all_of(Samplecols), ~ .x / median(.x[!contaminant], na.rm = T),
               .names = "{.col}_Norm")
      )
  }
  
  if (Normalization == "mean") {
    # Normalize peptide intensity by mean intensity (removing contaminants) observed in replicate
    rbrid <- rbrid %>% 
      mutate(
        across(all_of(Samplecols), ~ .x / mean(.x[!contaminant], na.rm = T),
               .names = "{.col}_Norm")
      )
  }
  
  # Update sample names
  if (Normalization == "none") {
    MinusSamples_Norm <- minusCols
    PlusSamples_Norm <- plusCols
  } else {
    MinusSamples_Norm <- paste0(minusCols, "_Norm")
    PlusSamples_Norm <- paste0(plusCols, "_Norm")
  }
  
  if (silac_flag == F) {
  # Compute mean of biological replicates
  rbrid$minus_mean <- rowMeans(rbrid[,MinusSamples_Norm], na.rm = T)
  rbrid$plus_mean <- rowMeans(rbrid[,PlusSamples_Norm], na.rm = T)
  
  rbrid$minus_mean[is.nan(rbrid$minus_mean)] <- NA
  rbrid$plus_mean[is.nan(rbrid$plus_mean)] <- NA
  } else {
    
    #for silac data, split data into minus and plus
    minus_data <- rbrid %>%
      dplyr::select(all_of(MinusSamples_Norm))
    plus_data <- rbrid %>%
      dplyr::select(all_of(PlusSamples_Norm))
    
    #identify how many replicates show peptide only in either heavy or light, or are missing in both
    tic("Count minus and plus 4SU peptide missing values")
    rbrid <- rbrid %>%
      mutate(both_present = rowSums(!is.na(minus_data) & !is.na(plus_data)),
             minusOnly = rowSums(!is.na(minus_data) & is.na(plus_data)),
             plusOnly = rowSums(is.na(minus_data) & !is.na(plus_data)),
             both_missing = rowSums(is.na(minus_data) & is.na(plus_data))
      ) %>%
      mutate(empty = both_missing == length(MinusSamples_Norm))
    #find peptides that are only present in minus or plus across all replicates
    rbrid <- rbrid %>%
      mutate(all_minusOnly = minusOnly > 0 & plusOnly == 0 & both_present == 0,
             all_plusOnly = plusOnly  > 0 & minusOnly == 0 & both_present == 0,
             minusObs = rowSums(!is.na(rbrid[,MinusSamples_Norm]) | rbrid[,MinusSamples_Norm] == 0, na.rm = T),
             plusObs =  rowSums(!is.na(rbrid[,PlusSamples_Norm]) | rbrid[,PlusSamples_Norm] == 0 , na.rm = T)) %>%
      mutate(empty = minusObs + plusObs == 0)
    
    # add pseudocounts
    if (Normalization == "silac_pseudo") {
      
      #for NA values on peptides where there is another measurement for the same 4su condition for the same peptide, 
      # replace NA value with the other measurement
        #replace missing minus 4su values
      #calculate minimum intensity before adding pseudocounts in each column (for later)
      minus_col_mins <- map_dbl(MinusSamples_Norm, ~ min(rbrid[[.x]], na.rm = T))
      plus_col_mins <- map_dbl(PlusSamples_Norm, ~ min(rbrid[[.x]], na.rm = T))
      
      groups <- c("minus", "plus")
      opposite_groups <- c("plus", "minus")
      group_cols <- c("MinusSamples_Norm", "PlusSamples_Norm")
      opposite_cols <- c("PlusSamples_Norm", "MinusSamples_Norm")
      
      for (i in seq_along(groups)) {
        
        group = groups[i]
        opp_group = opposite_groups[i]
        Columns = eval(str2lang(group_cols[i]))
        opposite = eval(str2lang(opposite_cols[i]))
      
        idx <- which(!rbrid[[paste0("all_",opp_group, "Only")]] & rbrid[[paste0(opp_group, "Only")]] > 0)
      # calculate minimum values
      pseudo_min <- map_dbl(idx,
                            ~ min(rbrid[.x , grepl(paste0(group, ".*Norm$"),names(rbrid))], na.rm = T))
      
      #replace minus NA values (that have matching observation in opposite isotope) with row-wise minimums
      rbrid[idx,Columns] <- 
        map_dfr(seq_along(idx), 
                function(index) {
                  row_index <- idx[[index]]
                  unformatted <- rbrid[row_index,Columns] %>% unlist()
                  #find columns where sample was detected in plus 4su isotope
                  col_index <- map_lgl(seq_along(Columns), ~ !is.na(rbrid[row_index,opposite[.x]]))
                  #replace NA values
                  unformatted[is.na(unformatted) & col_index] <- pseudo_min[[index]]
                  
                  #if only one observation is present in minus and plus 4su samples, merge them into single observation instead of duplicating values
                  present <- !is.na(rbrid[row_index, Columns])
                  singleObs <- sum(present) == 1 &  sum(col_index) == 1
                  if (singleObs) { unformatted[present] <- NA_real_ }
                  
                  return(unformatted)
                  })
    
      }
    }
      #update missing values for fisher test determination
    rbrid <- rbrid %>%
      mutate(both_present = rowSums(!is.na(rbrid[,MinusSamples_Norm]) & !is.na(rbrid[,PlusSamples_Norm])))
    
    #for silac data, split data into minus and plus (updated)
    minus_data <- rbrid %>%
      dplyr::select(all_of(MinusSamples_Norm))
    plus_data <- rbrid %>%
      dplyr::select(all_of(PlusSamples_Norm))
    
    ratios <- plus_data / minus_data
    ratios <- ratios %>%
      dplyr::rename_with(. , ~gsub(pattern = "(.*)_(light|heavy)_(minus|plus)_(.*)", replacement = "\\1_\\4_plus_minus_ratio", .))
    
    rbrid <- bind_cols(rbrid, ratios)
    
    #calculate mean ratio of plus / minus intensities for each peptide
      pseudo_minus <- mean(apply(minus_data, 2, min, na.rm=T))
      pseudo_plus <- mean(apply(plus_data, 2, min, na.rm=T))
      
      rbrid <- rbrid %>%
        #for peptides with at least one replicate with both isotopes, use average of ratios
        mutate(mean_plus_minus_ratio = rowMeans(rbrid[,grep("Norm_plus_minus_ratio", names(rbrid))], na.rm = T))
      
      #for completely missing values assign NA
      rbrid$mean_plus_minus_ratio[rbrid$empty] <- NA
      
      # add pseudo count to calculate average log fold changes for silac peptides with missing values
      
      # for peptides only present in minus or plus, divide by mean minimum intensity in the samples with the opposite treatment
      rbrid$mean_plus_minus_ratio[rbrid$all_minusOnly] <- pseudo_plus / rowMeans(rbrid[rbrid$all_minusOnly,MinusSamples_Norm], na.rm = T)
      rbrid$mean_plus_minus_ratio[rbrid$all_plusOnly] <- rowMeans(rbrid[rbrid$all_plusOnly,PlusSamples_Norm], na.rm = T) / pseudo_minus
      
      #if not using pseudocounts
      if (Normalization == "silac") {
      #for peptides where no replicate detects both isotopes, calculate ratio by merging plus and minus replicates first
      idx <- rbrid$both_present == 0 & rbrid$minusOnly > 0 & rbrid$plusOnly > 0
      rbrid$mean_plus_minus_ratio[idx] <- rowMeans(rbrid[idx,PlusSamples_Norm], na.rm = T) / rowMeans(rbrid[idx,MinusSamples_Norm], na.rm = T)
      }
      
      rbrid$log2fold <- log2(rbrid$mean_plus_minus_ratio)
      
  }
  
  if (convert.zero.means == T & silac_flag == F) 
  {
    #replace missing mean values with average of the lowest non-zero intensities across all samples
    mins <- rbrid %>%
      dplyr::select(all_of(c(MinusSamples_Norm, PlusSamples_Norm))) %>%
      summarise(across(everything(), ~ min(.x[.x != 0], na.rm =T))) %>%
      unlist()
    
    mean_min <- mean(mins)
    rbrid$minus_mean[rbrid$minus_mean==0 |is.na(rbrid$minus_mean)]<- mean_min
    rbrid$plus_mean[rbrid$plus_mean==0 | is.na(rbrid$plus_mean)]<- mean_min
  }

  if (silac_flag == F) {

  # Calculate +4su/-4sU log2-fold change
  rbrid$log2fold <- log2(rbrid$plus_mean/rbrid$minus_mean)
  rbrid$log2fold[is.infinite(rbrid$log2fold)] <- max(rbrid$log2fold)
  
  #count number of observations in minus and plus 4SU samples
  rbrid <- rbrid %>%
    mutate(minusObs =  rowSums(!is.na(rbrid[,MinusSamples_Norm]) | rbrid[,MinusSamples_Norm] == 0, na.rm = T),
           plusObs =  rowSums(!is.na(rbrid[,PlusSamples_Norm]) | rbrid[,PlusSamples_Norm] == 0 , na.rm = T)) %>%
    mutate(empty = minusObs + plusObs == 0)
  }
  
  
  #calculate statistical significance
  
  # calculate statistical signficance of depletion in 2 ways:
  # 1: for peptides with observations in both - and + 4SU samples, use a standard
  #    2 sample Student's t-test with equal varianc
  #   (or for silac a one-sample t-test comparing ratio against 1
  
  # 2: For peptides that only have observations in plus or minus 4su samples, instead perform a
  # fisher exact test where the test evaluates the presence/absence of a peptide between conditions
  
  #divide peptides into subsets for t-test or fisher test
  
  
  if (silac_flag == T) {
    #for silac identify peptides that do not have a single replicate with both heavy and light values for fisher test
    fisher_index <- rbrid$both_present == 0 & !rbrid$empty
    
  } else {
  #peptides where all peptides are missing in either condition
  fisher_index <- rbrid$minusObs ==0 | rbrid$plusObs == 0 
  }
  
  #subset peptides for t-test
  student_peps <- rbrid %>%
    dplyr::filter(!fisher_index)
  
  if (silac_flag == T) {
    # Calculate p values from one-sample two-sided t-test
    student_peps <- student_peps %>%
      mutate(test = "t_test")
    student_peps$pval <- apply(student_peps, 1, function(row) {
      pval(test_values =  row[grep("Norm_plus_minus_ratio",names(row))],
           type = "one.sample",
           paired = paired,
           var.equal = var.equal,
           mu = 1,
           alternative = "two.sided")
      })
  } else {
  # Calculate p values from two-sided t-test for non-silac data
  student_peps <- student_peps %>%
    mutate(test = "t_test")
  student_peps$pval = apply(student_peps, 1, FUN = function(row) {
    pval(control_values = row[MinusSamples_Norm],
         test_values =  row[PlusSamples_Norm],
         type = "two.sample",
         paired = paired,
         var.equal = var.equal,
         alternative = "two.sided")})
  
  }
  
  if (sum(fisher_index)> 0) {
    # Calculate p-values using fisher exact test
    fisher_peps <- rbrid %>%
      dplyr::filter(fisher_index)
    
    fisher_peps <- fisher_peps %>%
      mutate(test = "fisher_exact")
    fisher_peps$pval = apply(fisher_peps, 1, function(row) {
      if(row["empty"]) {return(1)}
      fisher_pval(minusVals = row[MinusSamples_Norm],
                  plusVals = row[PlusSamples_Norm])})
                                           
                                           
    #merge t-test and fisher exact test results
    
    rbrid <- bind_rows(student_peps,
                   fisher_peps)
  } else {
    #if no fisher peptides present, default to t-test
    rbrid <- student_peps
  }
  
  #adjust p values, exclude empty peptide rows (if any)
    rbrid$padj <- 1
    rbrid$padj[!rbrid$empty] <- stats::p.adjust(p = rbrid$pval[!rbrid$empty],
                                                             method = "BH")

  # Calculate RBR-ID score
  if (p.adjust) {
    rbrid$score_uses = "padj"
    rbrid$score <- -(log10(rbrid$padj)^2)*rbrid$log2fold
    
  } else {
    rbrid$score_uses = "pval"
    rbrid$score <- -(log10(rbrid$pval)^2)*rbrid$log2fold
  }

  rbrid$score[is.na(rbrid$score)] <- 0
  
  #drop silac columns
  if(silac_flag == T) {
    rbrid <- rbrid %>%
      dplyr::select(-all_minusOnly,
                    -all_plusOnly,
                    -both_missing,
                    -both_present)  %>%
      dplyr::select(1:empty,mean_plus_minus_ratio:score, everything())
  }

#Re-organize columns
leftcols <- names(rbrid)[!names(rbrid) %in% c(Samplecols, MinusSamples_Norm, PlusSamples_Norm)]

rbrid <- rbrid %>%
  dplyr::select(all_of(leftcols), everything()) %>%
  #Sort by descending score
  arrange(-score)

return(rbrid)
  }
  
pval <- function(control_values, test_values, type = "two.sample", paired = F, alternative="two.sided", var.equal = F, mu = 0, minObservations = 2, ...)
    #updated 06 Jun 2018
  {
    if (type == "two.sample") 
      ## If a paired test is desired all vectors must have same number of test and control values (i.e. no NA's).
    {
      x <- as.numeric(control_values)
      y <- as.numeric(test_values)
      
      # Returns p=1 if there are less than 2 non-zero/NA values in  minus 4SU observations
      # or if all values are zero
      if ( sum(!is.na(x))< minObservations | sum(!is.na(y)) < minObservations | 
           (sum(x, na.rm = T)==0 & sum(y, na.rm = T)==0)
      ){ return(pval = 1)}
      
      #remove NA values from comparison
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) != length(y) & paired == T)
        stop("Columns do not have same # of observations for paired t-test")
      
      if (((all(x == 0) & all(y == 0)) | ( all(is.na(x)) | all(is.na(y))))) {return(pval = 1)}
      
      # if all values are 1 return p = 1  (for protein normalization)
      if ( all(x == 1) | all (y == 1) ) {
        return(pval = 1)
      }
      
      else {pval <- t.test(x, y, paired = paired, var.equal=var.equal, alternative = alternative, ...)$p.value}
      
      return(pval) }
    
    if (type == "one.sample")
    {
      x <- as.numeric(test_values)
      
      
      # Returns p=1 if there are less than 2 values in any one distribution or there is no difference b/t replicates
      # Returns p = 1 if there are no observations any one of the replicates
      if (sum(!is.na(x))<2 ) # sd(x, na.rm = T) == 0 )
      {pval = 1}
      
      else { pval <- t.test(x, mu = mu, alternative = alternative, var.equal = var.equal, ...)$p.value}
      
      return(pval) }
  }
  
fisher_pval <- function(minusVals, plusVals)
    #updated 11 Apr 2022
{
    # count observations present or missing in minus / plus 4su samples
    Obs <- c(minus_missing = sum(minusVals %in% c(NA, 0)), 
             minus_present = sum(!minusVals %in% c(NA, 0)), 
             plus_missing = sum(plusVals %in% c(NA, 0)),
             plus_present = sum(!plusVals %in% c(NA, 0))
    )
    mat <- matrix(Obs,nrow = 2)
    p.value <- fisher.test(mat)$p.value
    return(p.value)
}
  
  
  ##Peptideplot
  
  Peptideplot <- function(method, samples, Treatment = "minus", pep_percent = 0.1, experiment = "Your Favorite Experiment", ...) 
    # Function to plot normalized peptide intensities for all signal peptides. "method" delivers normalization method to function RBRnorm, accepted values are "Signal", "Flag", or "RNase". 
    # Treatment = "plus" or "minus" indicates whether minus or plus 4sU samples are to be shown.
    # pep_percent specifies what percent of top peptides after sorting to show on plot
    # Additional arguments can be passed as options to RBRnorm
  {
    library(dplyr)
    Prot_rbrid <- NormRBR(rbrid = samples, Normalization = method, ...)
    
    ordermethod <- paste0(Treatment, "_mean")
    
    signalpeps <- Prot_rbrid[!Prot_rbrid$contaminant,]
    
    
    # reorder rbrid for -/+ 4sU mean intensities
    signalpeps <- signalpeps %>%
      arrange(!is.na(signalpeps[[paste0(ordermethod)]]), signalpeps[[paste0(ordermethod)]])
    
    #Filter top % of peptides after sorting
    signalpeps <- signalpeps[round(nrow(signalpeps)*(1-pep_percent)):nrow(signalpeps),]
    
    # reorder peptide ID factor
    signalpeps$Pept_ID <- factor(signalpeps$Pept_ID, levels = signalpeps$Pept_ID)  
    
    # identify number of normalized replicates
    reps <- names(signalpeps)[grepl(pattern = c(paste0(Treatment, ".*Norm")) , names(signalpeps))]
    
    # Plot peptides
    plot(x = signalpeps$Pept_ID, y = signalpeps[[paste0(reps[1])]], 
         xlab = paste0("Top ", pep_percent*100, "% Peptides ranked by ", ordermethod), 
         ylab = paste0("Normalized intensity (", method, ")"), 
         main = paste0(experiment, " ", Treatment, " 4sU ", method, " Normalization"), 
         xaxt = 'n', pch = 16, ylim=c(0, max(signalpeps[reps], na.rm = T)), cex= 0.5)
    
    for (i in 2:length(reps)) {
      samplerep <- signalpeps[[paste0(reps[i])]]
      points(signalpeps$Pept_ID, samplerep, col=i, pch = 16, cex = 0.5)
    }
    
    legend(x = length(signalpeps$Pept_ID)*0.1, y = max(signalpeps[reps], na.rm = T), legend =  reps, col = c(1:length(reps)), pch = 16)
    
  }
  
  
  ##Coverage.apply
  
  calculateCoverage.apply<-function(range,peps, index = NA)
  {
    if (is.na(index)) {index <- 1:length(range)}
    hits<-subsetByOverlaps(peps,range[index])
    covered.regions<- GenomicRanges::reduce(hits)
    coverage <- sum(width(covered.regions))/end(range[index])
    
    return(coverage)
  }
  
  coverage.lapply <- function(pep.gr.list, annot.gr, threads = 5, annotate.experiment = F, ...) {
    #Last updated 15oct2018
    #calculates coverage for multiple pep.gr files (e.g. different experiments)
    #requires genomic ranges
    #Updated to be compatible with Bioconductor 3.7
    #initialize parallelization
    cl<-makeCluster(threads)
    
    prot.gr.coverage <- lapply(1:length(pep.gr.list), function(i, ...) {
      pep.gr <- pep.gr.list[[i]]
      
      clusterExport(cl, list("subsetByOverlaps","reduce","width","end", "calculateCoverage.apply")) # must pass non-conventional functions to worker nodes  
      coverage <- parSapply(cl=cl, 1:length(annot.gr), FUN = function(index, ...) {
        temp.range = annot.gr[index]
        calculateCoverage.apply(range = temp.range, peps = pep.gr)})
      
      # this takes about 1'15" with 30 nodes
      
      annot.gr$coverage<-round(coverage*100, digits = 2)
      if(annotate.experiment) {annot.gr$Experiment <- unique(pep.gr$Experiment)}
      
      return(annot.gr)})
    
    stopCluster(cl) # close parallel cluster
    return(prot.gr.coverage)}
  
  
  ##residueScorePlot
  # 06Jun2018 - added code to plot lines with gaps for NA, removed glitchy whitespace overlay
  # 11Jun2018 - added code to export missing values as gaps with return.gaps argument
  # 05Dec2018 - changed code to avoid overlaying domain labels for clarity
  
  residueScorePlot<-function(protein, uniprot = "", range, peps, IPR_annotations = NA, Manual_annotations = NA, highlight_peptide = NA, highlight_range = NA, idealpeps = NA, smooth=T, no.negs=F,span=0.0001, return.values=F, title = "", color = "Greys", enzym = "Trypsin", excludeDomain = "", ONLYdomain = "", return.gaps = F, ...)
    ## LAST UPDATED 05Dec018
    ## Draws score plot for indicated protein symbol, defaults to uniprot ID with highest coverage unless uniprot is specified using uniprot argument
    #  smoothing and white space on NA peptides can be controlled by options
    ## Shows ideal peptide coverage 
    ## Domains can now be plotted using the IPR_annotations and Manual_annotations arguments.
    # IPR_annotations takes as input a rbrid.frame of IPR annotations with the following columns: 
    # "uniprotID","iprID","IPR_name","IPR_signature","IPR_start","IPR_end"
    # If domains as shown appear cluttered, specific domains can be excluded via the excludeDomain argument, or only desired domains can be shown via ONLYdomain
    # Both exludeDomain and ONLYdomain recognize regex
    # If instead you want to input the domains manually, feed a rbrid.frame to Manual_annotations argument.
    # Manual annotations - 1st col (Domain name), 2nd (start aa), 3rd (end aa)
  # Can highlight desired peptide if specified using highlight_peptide
  # Can highlight desired amino acid range via highlight_range argument
  # Updated code to correctly show domain annotation (y value plot object) for no.negs = F
  {
    library(ggplot2)
    library(GenomicRanges)
    if (nchar(uniprot) == 0)
      
    {
      # If there are multiple uniprot IDs present for a given protein, will default to one with highest coverage
      
      pos <- which(range$symbol == protein) # ID's matching uniprots
      if (length(pos) == 0) return(paste0("Protein ", protein, " is not identified in range. Check Uniprot IDs"))
      
      pos <- pos[which.max(range[pos]$coverage)] # selects max coverage ID
      uniprot <- as.character(seqnames(range[pos]))} 
    
    # ID's matching uniprots
    range <- range[seqnames(range) == uniprot]
    peps <- peps[seqnames(peps) == uniprot]
    
    # Check if protein is present
    check.uniprot <- which(seqnames(range) == uniprot)
    if (length(check.uniprot) ==0) { return(paste0("Protein ", protein, " / Uniprot ID: ", uniprot, " is not identified in range. Check Uniprot IDs")) }
    
    if (!any(is.na(idealpeps))) {idealpeps <- idealpeps[[uniprot]]}
    
    score<-rep(0,end(range))
    hits<-subsetByOverlaps(peps,range)
    
    #export missing values
    if(return.gaps == T) {
      values <- list(x =1:end(range), missing = rep(0, end(range)))
      gaps <- GenomicRanges::setdiff(range, GenomicRanges::reduce(hits))
      
      for (i in 1:length(gaps)) {     
        gap.pos<-GenomicRanges::start(gaps[i]):GenomicRanges::end(gaps[i])
        values$missing[gap.pos]<- 1}
      return(values) }
    
    for(i in 1:length(hits)){
      pep.pos<-start(hits[i]):end(hits[i])
      pep.pos<-pep.pos[pep.pos > 0]
      score[pep.pos]=score[pep.pos]+hits$score[i]
    }
    
    score.mod<- score
    if(no.negs) score.mod[score.mod<0]<-0
    if(smooth){
      values<-supsmu(1:length(score.mod),score.mod,span=span)
    } else {
      values<-(list(x=1:length(score.mod),y=score.mod))
    }
    
    gaps <- GenomicRanges::setdiff(range, GenomicRanges::reduce(hits))
    for(i in 1:length(gaps)){
      gap.pos<-start(gaps[i]):end(gaps[i])
      values$y[gap.pos]<-NA
      
      # Assign NA values adjacent to non-missing values (at edges of gap) 
      # as 0 to avoid discontinuous peaks during plotting
      values$y[c(start(gaps[i]), end(gaps[i]))] <- 0
    }
    
    rbrid <- rbrid.frame(residue = values$x, score = values$y)
    
    # Output values
    if(return.values) return(values)
    
    # Only show specific amino acid range on plot
    # e.g. highlight_range = c(0,50) will show amino acids 0 to 50
    if (!any(is.na(highlight_range))) {
      rbrid <- rbrid[highlight_range[1]:highlight_range[2],]}
    
    #generate groups for plotting NA values as gaps, have to plot separate lines grouped by factor to avoid interpolating gaps
    #remove missing values
    rbrid <- rbrid[!is.na(rbrid$score),]
    # generate the groups automatically 
    idx <- c(1, diff(rbrid$residue))
    i2 <- c(1,which(idx != 1), nrow(rbrid)+1)
    rbrid$grp <- rep(1:length(diff(i2)), diff(i2))
    
    # generate plot object
    {
      legend.vector <- c(F, T)
      names(legend.vector) <- c("colour", "linetype")
      v <- ggplot(rbrid, mapping = aes(x = residue, y = score)) +
        theme_classic() +
        theme_bw()+
        geom_line(aes(group = grp, col = Experiment, linetype = Experiment), show.legend = legend.vector) +
        xlab("residue #") +
        ylab("RBR-ID score") +
        scale_color_manual(labels = legend.labs, values = colors) +
        if (!any(is.na(idealpeps))) { ggtitle(paste0(range$symbol, " (",seqnames(range),")- Adj. Cov. ", round(range$Adjusted_Cov, 1), "%/ Total Cov. ", round(range$coverage, 1), " % ", title))} else 
        { ggtitle(paste0(range$symbol, " (",seqnames(range),")- Coverage ", round(range$coverage, 1), "% " , title)) }
    }
    # Set lower limit for y axis for annotations
    if(no.negs) 
    {v <- v + ylim(-0.05*max(rbrid$score), layer_scales(v)$y$range$range[[2]])
    ymin <- layer_scales(v)$y$limits[[1]]}
    else if (abs(min(rbrid$score)) < 0.05*max(rbrid$score)) 
    {v <- v + ylim(-0.05*max(rbrid$score), layer_scales(v)$y$range$range[[2]])
    ymin <- layer_scales(v)$y$limits[[1]]}
    else {ymin <- (min(rbrid$score)-0.1*abs(min(rbrid$score)))}
    
    if (!any(is.na(idealpeps)))
      # Show ideal coverage overlaid on plot
    {idpep <- GenomicRanges::reduce(idealpeps)
    vertices <- lapply(1:length(idpep), function(x) {
      pos <- rbrid.frame(xx = c(rep(start(idpep[x]), 2), rep(end(idpep[x]), 2)),
                         yy = c(ymin, rep(0, 2), ymin),
                         group = x)
    })
    
    
    vertices <- do.call(rbind.rbrid.frame,vertices)
    
    v <- v + geom_polygon(rbrid = vertices, mapping = aes(x = xx, y = yy, group=group, fill = "peptide"), alpha = 0.4)
    
    v <- v + scale_fill_manual(
      name   = paste0('Ideal Coverage - ', round(range$Idealcoverage, 1), " %"),
      breaks = c("peptide"), # <<< corresponds to fill aesthetic labels
      values = c("light green"),
      labels = paste0(enzym, " (", idealpeps$lowMW[1], " Da to ", idealpeps$hiMW[1], " Da)")
    )
    }  
    
    # Highlight specific peptide and scores on plot
    if (!is.na(highlight_peptide)) {
      
      extracted_range <-regexpr(highlight_peptide, range$sequence)
      
      if (extracted_range[1] == -1) {
        return(paste0("Highlighted peptide, ", highlight_peptide, " is not identified in given range"))}
      
      
      
      highlight <- rbrid.frame(Sequence = highlight_peptide, 
                               start.position = extracted_range[1], 
                               end.position = extracted_range[1] + attr(extracted_range, which = "match.length") -1)
      
      # Convert peptide query to Granges
      highlight <- GRanges(seqnames=Rle(seqnames(range)), ranges=IRanges(start=highlight$start.position, end=highlight$end.position))
      highlight$sequence <- highlight_peptide
      
      
      highlight_hits <- subsetByOverlaps(query = peps, subject = highlight) # find matching peptides in rbrid
      if (length(highlight_hits)==0) {
        return(paste0("Highlighted peptide,", highlight_peptide, " is not identified in given range"))}
      
      highlight_table <- rbrid.frame(
        Domain = paste0(highlight_hits$sequence, " - Score ", round(highlight_hits$score, 2)),
        start = start(highlight_hits),
        end = end(highlight_hits),
        score = round(highlight_hits$score, 2))
      
      domain_indices <- lapply(1:nrow(highlight_table), function(x) {
        pos <- rbrid.frame(
          xx = c(rep(highlight_table[x,"start"], 2), rep(highlight_table[x,"end"], 2)),
          yy = c(ymin, rep(0, 2), ymin),
          Domain = highlight_table[x,1])
      })
      
      
      domain_indices <- do.call(rbind.rbrid.frame,domain_indices)
      # Reorder domain factors so legend displays domains N-to-C terminal
      order <- domain_indices$Domain[order(domain_indices$xx)]
      order <- order[!duplicated(order)]
      domain_indices$Domain <- factor(domain_indices$Domain, levels = order)
      
      v <- v + geom_polygon(rbrid = domain_indices, mapping = aes(x = xx, y = yy, group=Domain, fill =
                                                                    Domain), alpha = 0.4)
      # Adjust x axes to focus on highlights
      min <- floor((min(domain_indices$xx)/length(values$x))*length(values$x)*0.9)
      max <- ceiling((max(domain_indices$xx)/length(values$x))*length(values$x)*1.1)
      v <- v + xlim(min, max)
      return(v)
    }
    
    if (any(!is.na(Manual_annotations))){
      names(Manual_annotations)[2:3] <- c("start", "end")
      
      # If desired, exclude specified domains 
      if (nchar(excludeDomain) > 0) 
      {Manual_annotations <- Manual_annotations[!grepl(excludeDomain, Manual_annotations[,1]),] }
      
      # If desired, display only specific domains
      if (nchar(ONLYdomain) > 0) 
      {Manual_annotations <- Manual_annotations[grepl(ONLYdomain, Manual_annotations[,1]),]}
      
      domain_indices <- lapply(1:nrow(Manual_annotations), function(x) {
        # Split domain y-position to avoid overlaying domain labels
        y.unit <- (min(rbrid$score) - ymin)/nrow(Manual_annotations)
        y.upper <- (min(rbrid$score)-y.unit*(x-1))
        y.lower <- (min(rbrid$score)-y.unit*(x))
        
        pos <- rbrid.frame(
          xx = c(rep(Manual_annotations[x,"start"], 2), rep(Manual_annotations[x,"end"], 2)),
          yy = c(y.lower, rep(y.upper, 2), y.lower),
          Domain = Manual_annotations[x,1]) })
      
      domain_indices <- do.call(rbind.rbrid.frame,domain_indices)
      
      # Reorder domain factors so legend displays domains N-to-C terminal
      order <- domain_indices$Domain[order(domain_indices$xx)]
      order <- order[!duplicated(order)]
      domain_indices$Domain <- factor(domain_indices$Domain, levels = order)
      
      v <- v + geom_polygon(rbrid = domain_indices, mapping = aes(x = xx, y = yy, group=Domain, fill = Domain), alpha = 0.4)
    }
    
    ## To plot IPR domains onto Score plot
    if (!any(is.na(IPR_annotations))) {
      
      # Find matching domains for protein
      hits<-which(IPR_annotations$uniprotID==uniprot)
      
      if (length(hits) != 0) {
        Prot_domains <- IPR_annotations[hits,]
        
        
        # If desired, exclude specified domains 
        if (nchar(excludeDomain) > 0) 
        {Prot_domains <- Prot_domains[!grepl(excludeDomain, Prot_domains$IPR_name),] }
        
        # If desired, display only specific domains
        if (nchar(ONLYdomain) > 0) 
        {Prot_domains <- Prot_domains[grepl(ONLYdomain, Prot_domains$IPR_name),]}
        
        
        # Filter for domains within range (highlight_range option)
        if (!any(is.na(highlight_range))) {
          Prot_domains <- Prot_domains[Prot_domains$IPR_end >= highlight_range[1] & Prot_domains$IPR_start <= highlight_range[2],]
          
          Prot_domains$IPR_start[Prot_domains$IPR_start <= highlight_range[1]] <- highlight_range[1]
          Prot_domains$IPR_end[Prot_domains$IPR_end >= highlight_range[2]] <- highlight_range[2]
        }
        
        
        # Convert to genomic ranges
        domains.gr<- GRanges(seqnames=Rle(Prot_domains$uniprotID), ranges=IRanges(start=Prot_domains$IPR_start, end=Prot_domains$IPR_end))
        domains.gr$Domain <- Prot_domains$IPR_name
        domains.gr$iprID <- Prot_domains$iprID
        
        # Collapse domains with slightly different names
        Prot_domains$IPR_name <- Collapse_Domains(Domains = Prot_domains$IPR_name, breakwords = 2)
        
        temp_domain <- Prot_domains$IPR_name[!duplicated(Prot_domains$IPR_name)]
        
        # Remove redundant overlapping domains
        domain_list <- sapply(temp_domain, function(x) {
          filtered <- GenomicRanges::reduce(domains.gr[domains.gr$Domain == x])
          filtered$Domain <- x 
          filtered$iprID <- domains.gr[domains.gr$Domain == x]$iprID[1]
          return(filtered)
        })
        
        # For clarity, remove domains that span entire range of protein
        domain_filtered <- lapply(domain_list, function(x) {
          x <- x[width(x) != end(range.list[[1]])] 
        })
        
        empty <- which(lapply(domain_filtered, length) == 0)
        if (length(empty) != 0) {domain_filtered <- domain_filtered[-empty]}
        
        # If no domains are left after filtering, output plot
        if (length(domain_filtered) == 0) {return(v)}
        
        domain_indices <- lapply(1:length(domain_filtered), function(index) {
          domain_type <- domain_filtered[[index]]
          pos <- lapply(1:length(domain_type), function(x) {
            # Split domain y-position to avoid overlaying domain labels
            y.unit <- (min(rbrid$score) - ymin)/length(domain_filtered)
            y.upper <- (min(rbrid$score)-y.unit*(index-1))
            y.lower <- (min(rbrid$score)-y.unit*(index))
            
            Num <- c(" ", 2:length(domain_type)) # Domain number of same type
            rbrid.frame(
              xx = c(rep(start(domain_type)[x], 2), rep(end(domain_type)[x], 2)),
              yy = c(y.lower, rep(y.upper, 2), y.lower),
              Domain = paste0(domain_type$Domain[x], " ", Num[x]))
          })
          comb_index <- do.call(rbind.rbrid.frame, pos)
        })       
        
        if (length(domain_indices) == 1) {domain_indices <- domain_indices[[1]]} else {domain_indices <- do.call(rbind.rbrid.frame,domain_indices)}
        
        
        # Reorder domain factors so legend displays domains N-to-C terminal
        order <- domain_indices$Domain[order(domain_indices$xx)]
        order <- order[!duplicated(order)]
        domain_indices$Domain <- factor(domain_indices$Domain, levels = order)
        domain_indices$xx <- as.numeric(domain_indices$xx)
        domain_indices$yy <- as.numeric(domain_indices$yy)
        domain_indices$Domain <- as.character(domain_indices$Domain)
        # Plot IPR domains
        v <- v + geom_polygon(rbrid = domain_indices, mapping = aes(x = xx, y = yy, group= Domain, fill = Domain), alpha = 0.4)
      }}
    
    # Crop x axis range for highlight_range = T
    if (!any(is.na(highlight_range))) {
      v <- v + xlim(highlight_range[1], highlight_range[2]) }
    
    return(v)
  }
  
  
  ##OverlayScorePlot
  # Updated 30Jan2020
  
  OverlayresidueScorePlot<-function(protein, uniprot = "", range.list, pep.list, IPR_annotations = NA, Manual_annotations = NA, highlight_peptide = NA, highlight_range = NA, smooth=T, no.negs=F,span=0.0001, return.values=F, title = "", palette = "Greys", man.color = NA, enzym = "Trypsin", excludeDomain = "", ONLYdomain = "", ymax = NA, ...)
    ## LAST UPDATED 15Jan2021
    ## Draws score plot for indicated protein symbol across multiple experiments, defaults to uniprot ID with highest coverage in first experiment unless uniprot is specified using uniprot argument
    #  smoothing and white space on NA peptides can be controlled by options
    ## Domains can now be plotted using the IPR_annotations and Manual_annotations arguments.
    # IPR_annotations takes as input a rbrid.frame of IPR annotations with the following columns: 
    # "uniprotID","iprID","IPR_name","IPR_signature","IPR_start","IPR_end"
    # If domains as shown appear cluttered, specific domains can be excluded via the excludeDomain argument, or only desired domains can be shown via ONLYdomain
    # Both exludeDomain and ONLYdomain recognize regex
    # If instead you want to input the domains manually, feed a rbrid.frame to Manual_annotations argument.
    # Manual annotations - 1st col (Domain name), 2nd (start aa), 3rd (end aa)
    # Can highlight desired peptide if specified using highlight_peptide
  # Can highlight desired amino acid range via highlight_range argument
  # Updated code to correctly show domain annotation (y value plot object) for no.negs = F
  # Updated code to correctly show gaps missing in at least one experiment
  # 30Jan2020 Updated code to remove dodge plotting of manual annotations which was buggy
  # 15Jan2021 Reordered factors to remove bug where experiment labels would get switched
  {
    # Create a new ggplot theme
    theme_custom <- function (base_size = 11, base_family = "") {
      theme_bw() %+replace% 
        theme(
          panel.grid.major  = element_line(color = "grey", size= 0.2, linetype = 5),
          panel.grid.minor = element_blank())}
    
    library(ggplot2)
    library(GenomicRanges)
    library(RColorBrewer)
    
    # Extract uniprotID if none specified
    if (nchar(uniprot) == 0)
      
    {
      # If there are multiple uniprot IDs present for a given protein, will default to one with highest coverage
      exp.hits <- sapply(range.list, function(gr) {index <- any(grepl(protein, gr$symbol))})
      temp.range <- range.list[exp.hits]
      protein.hits <- lapply(temp.range, function(gr) {filter <- gr[gr$symbol == protein]
      df <- rbrid.frame(symbol = filter$symbol, uniprot = as.character(seqnames(filter)), coverage = filter$coverage)
      return(df)})
      protein.hits <- do.call("rbind", protein.hits)
      
      if (is.null(protein.hits)) return(paste0("Protein ", protein, " is not identified in selected experiments"))
      
      
      uniprot <- as.character(protein.hits$uniprot[which.max(protein.hits$coverage)]) # selects max coverage ID
      
    } 
    
    # extract range and peps matching uniprotID
    range.list <- lapply(range.list, function(prot.gr) {range.list <- prot.gr[seqnames(prot.gr) == uniprot]})
    pep.list <- lapply(pep.list, function(pep.gr) {peps <- pep.gr[seqnames(pep.gr) == uniprot]})
    Experiments <- names(range.list)
    
    # Check if protein is present in each experiment
    index <- sapply(Experiments, function(Exp) {
      peps <- pep.list[[Exp]]
      range <- range.list[[Exp]]
      if (length(peps) == 0) {print(paste0("Protein ", protein, " / Uniprot ID: ", uniprot, " is not identified in experiment ", Exp, ". Check Uniprot IDs"))  
        return(F)} else {return(T)} })
    
    # remove ranges where uniprot ID is not detected
    range.list <- range.list[index]
    Experiments <- Experiments[index]
    
    ## generate plotting values ###
    rbrid.values <- lapply(Experiments, function(Exp) { 
      # Extract ranges & peps from each experiment
      range <- range.list[[Exp]]
      peps <- pep.list[[Exp]]
      
      # Calculate residue level scores
      score<-rep(0,end(range))
      hits<-suppressWarnings(subsetByOverlaps(peps,range))
      for(i in 1:length(hits)){
        pep.pos<-start(hits[i]):end(hits[i])
        pep.pos<-pep.pos[pep.pos > 0]
        score[pep.pos]=score[pep.pos]+hits$score[i]
      }
      
      score.mod<- score
      if(no.negs) score.mod[score.mod<0]<-0
      if(smooth){
        values<-supsmu(1:length(score.mod),score.mod,span=span)
      } else {
        values<-(list(x=1:length(score.mod),y=score.mod))}
      
      gaps<- suppressWarnings(GenomicRanges::setdiff(range, GenomicRanges::reduce(hits)))
      for(i in 1:length(gaps)){
        gap.pos<-start(gaps[i]):end(gaps[i])
        values$y[gap.pos]<-NA
        
        # Assign NA values adjacent to non-missing values (at edges of gap) 
        # as 0 to avoid discontinuous peaks during plotting
        values$y[c(start(gaps[i]), end(gaps[i]))] <- 0
      }
      
      rbrid <- rbrid.frame(residue = values$x, 
                           score = values$y, 
                           Experiment = as.character(range$Experiment))
      #Anchor end of protein peaks at zero
      if (rbrid[nrow(rbrid),"score"] !=0) {
        rbrid[nrow(rbrid)+1,] <- rbrid[nrow(rbrid),]
        rbrid[nrow(rbrid),c("residue", "score")] <- c(rbrid$residue[nrow(rbrid)-1]+1, 0) 
      }
      return(rbrid) })
    
    #merge rbrid
    rbrid <- do.call("rbind", rbrid.values)
    
    # Output values
    if(return.values) return(rbrid)
    
    #Extract coverage from each experiment
    cov.df <- rbrid.frame(uniprot = rep(uniprot, length(Experiments)),
                          Coverage = sapply(range.list, function(range) {range$coverage}),
                          Experiment = sapply(range.list, function(range) {range$Experiment}))
    # legend label text
    exp.labs <- sapply(range.list, function(list) {list$Experiment[1]})
    legend.labs <- paste(exp.labs,"/", round(cov.df$Coverage, 1), "% cov.")
    # Only show specific amino acid range on plot
    # e.g. highlight_range = c(0,50) will show amino acids 0 to 50
    if (!any(is.na(highlight_range))) {
      start <- highlight_range[1]
      end <- highlight_range[2]
      rbrid <- rbrid[rbrid$residue >= start & rbrid$residue <= end,]}
    
    #generate groups for plotting NA values as gaps, have to plot separate lines grouped by factor to avoid interpolating gaps
    #remove missing values
    rbrid <- rbrid[!is.na(rbrid$score),]
    # generate the groups automatically 
    idx <- c(1, diff(rbrid$residue))
    i2 <- c(1,which(idx != 1), nrow(rbrid)+1)
    rbrid$grp <- rep(1:length(diff(i2)), diff(i2))
    
    # generate color palette
    colors <- brewer.pal(n = ifelse(length(Experiments) < 3, yes = 3, no = length(Experiments)), name = palette) #palette
    if (length(Experiments) == 1) {colors <- "black"}
    if (!is.na(man.color[1])) {colors <- man.color}
    # generate linetype palette (alternate solid and dotted)
    
    
    # reorder factors
    rbrid$Experiment <- factor(rbrid$Experiment, levels = exp.labs)
    
    # generate plot object
    v <- ggplot(rbrid, mapping = aes(x = residue, y = score)) +
      theme_custom() +
      geom_line(aes(group = grp, col = Experiment, linetype = Experiment)) +
      xlab("residue #") +
      ylab("RBR-ID score") +
      ggtitle(paste0(protein, " (uniprotID: ", uniprot,") RBR-ID residue scores "), title) +
      scale_color_manual(name = "Experiment", labels = legend.labs, values = colors) +
      scale_linetype_manual(name = "Experiment", labels = legend.labs, values = rep(c("solid", "dotted"), length.out = length(Experiments))) 
    
    if(!is.na(ymax)) {v <- v + ylim(min(rbrid$score), ymax)}
    # Set lower limit for y axis for annotations
    if(no.negs) 
    { if(!is.na(ymax)) {v <- v + ylim(-0.15*ymax, ymax)
    ymin <- -0.15*ymax} else {
      v <- v + ylim(-0.1*max(rbrid$score), layer_scales(v)$y$range$range[[2]])
      ymin <- layer_scales(v)$y$limits[[1]]} }
    
    # Plot domain annotations onto protein
    if (any(!is.na(Manual_annotations))){
      names(Manual_annotations)[2:3] <- c("start", "end")
      
      # If desired, exclude specified domains 
      if (nchar(excludeDomain) > 0) 
      {Manual_annotations <- Manual_annotations[!grepl(excludeDomain, Manual_annotations[,1]),] }
      
      # If desired, display only specific domains
      if (nchar(ONLYdomain) > 0) 
      {Manual_annotations <- Manual_annotations[grepl(ONLYdomain, Manual_annotations[,1]),]}
    }
    
    ## To plot IPR domains onto Score plot
    if (!any(is.na(IPR_annotations))) {
      
      # Find matching domains for protein
      hits<-which(IPR_annotations$uniprotID==uniprot)
      
      if (length(hits) != 0) {
        Prot_domains <- IPR_annotations[hits,]
        
        # If desired, exclude specified domains 
        if (nchar(excludeDomain) > 0) 
        {Prot_domains <- Prot_domains[!grepl(excludeDomain, Prot_domains$IPR_name),] }
        
        # If desired, display only specific domains
        if (nchar(ONLYdomain) > 0) 
        {Prot_domains <- Prot_domains[grepl(ONLYdomain, Prot_domains$IPR_name),]}
        
        # Filter for domains within range (highlight_range option)
        if (!any(is.na(highlight_range))) {
          Prot_domains <- Prot_domains[Prot_domains$IPR_end >= highlight_range[1] & Prot_domains$IPR_start <= highlight_range[2],]
          
          Prot_domains$IPR_start[Prot_domains$IPR_start <= highlight_range[1]] <- highlight_range[1]
          Prot_domains$IPR_end[Prot_domains$IPR_end >= highlight_range[2]] <- highlight_range[2]
        }
        
        # Convert to genomic ranges
        domains.gr<- GRanges(seqnames=Rle(Prot_domains$uniprotID), ranges=IRanges(start=Prot_domains$IPR_start, end=Prot_domains$IPR_end))
        domains.gr$Domain <- Prot_domains$IPR_name
        domains.gr$iprID <- Prot_domains$iprID
        
        # Collapse domains with slightly different names
        Prot_domains$IPR_name <- Collapse_Domains(Domains = Prot_domains$IPR_name, breakwords = 2)
        temp_domain <- Prot_domains$IPR_name[!duplicated(Prot_domains$IPR_name)]
        
        # Remove redundant overlapping domains
        domain_list <- sapply(temp_domain, function(x) {
          filtered <- GenomicRanges::reduce(domains.gr[domains.gr$Domain == x])
          filtered$Domain <- x 
          filtered$iprID <- domains.gr[domains.gr$Domain == x]$iprID[1]
          return(filtered)
        })
        
        # For clarity, remove domains that span entire range of protein
        domain_filtered <- lapply(domain_list, function(x) {x <- x[width(x) != end(range.list[[1]])]  
        })
        
        empty <- which(lapply(domain_filtered, length) == 0)
        if (length(empty) != 0) {domain_filtered <- domain_filtered[-empty]}
        
        # If no domains are left after filtering, output plot
        if (length(domain_filtered) == 0) {return(v)}
        
        #convert to df
        domain_filtered <- GRangesList(domain_filtered)
        domain_table <- as.rbrid.frame(GRangesList(domain_filtered)) %>% mutate(uniprotID = seqnames) %>% select(Domain, start, end, uniprotID) 
        
        if (!any(is.na(Manual_annotations))) { domain_table <- rbind(domain_table, Manual_annotations) }
      }
    }
    
    if (!any(is.na(Manual_annotations)) | !any(is.na(IPR_annotations))) {
      #prepare annotations for plotting
      dom.names <- unique(domain_table$Domain) %>% as.character()
      domain_indices <- lapply(dom.names, function(dom) {
        #Split domain y-position to avoid overlaying domain labels
        x=which(dom.names == dom)
        y.unit <- (min(rbrid$score) - ymin)/length(dom.names)
        y.upper <- (min(rbrid$score)-y.unit*(x-1))
        y.lower <- (min(rbrid$score)-y.unit*(x))
        
        temp.df <- domain_table %>% filter(Domain == dom)
        dom.merge <- rbrid.frame()
        pos.list <- lapply(1:nrow(temp.df), function(i) {
          pos <- rbrid.frame(
            xx = rep(as.numeric(temp.df[i,2:3]), each=2),
            yy = c(y.lower, y.upper, y.upper, y.lower),
            Domain = dom)
          #for domains spanning single residue increase length +1 (so they're not invisible)
          if(var(pos$xx) ==0) {pos$xx[3:4] <- pos$xx[3:4]+1}
          return(pos)
        })
        dom.merge <- do.call("rbind", pos.list)
        return(dom.merge)
      })
      
      domain_indices <- do.call(rbind.rbrid.frame,domain_indices)
      
      # Reorder domain factors so legend matches position of domains
      order <- domain_indices$Domain[order(domain_indices$yy, decreasing=T)]
      order <- order[!duplicated(order)]
      domain_indices$Domain <- factor(domain_indices$Domain, levels = order)
      
      # Plot domains
      v <- v + geom_polygon(rbrid = domain_indices, mapping = aes(x = xx, y = yy, group= Domain, fill = Domain), alpha = 0.4)
      # Crop x axis range for highlight_range = T
      if (!any(is.na(highlight_range))) {
        v <- v + xlim(highlight_range[1], highlight_range[2]) }
      
      return(v)
    }
  }
  
  
  ##Collapse_domains
  Collapse_Domains <- function(Domains, breakwords = 2) 
    # Function to collapse to merge domains that have similar but non-identical names.
    # breakwords specifies how many of the starting words in each name to look for redundancies. I.E. for breakwords = 2 Collapse_Domains will merge all names that are identical up to the second word
    # e.g. overlapping domains named "Zinc finger" and "Zinc finger C2H2" would be collapsed to union named "Zinc finger"
  {for (i in 1:length(Domains)) {
    foo <- unlist(strsplit(Domains[i], "\\W+"))
    pattern <- paste(foo[1:breakwords], collapse = " ")
    hits <- grepl(pattern, x = Domains)
    Domains[hits] <- Domains[hits][1] # Change name of all matching domains to that of first hit
  }
    return(Domains)}
  
  
  ##Listscoreplot
  
  listScoreplot <- function (Proteins, Ties = "", PDB = NA, range, peps, wd = "~/", Pval = "0.2", Experiment = "", Descriptor = "", Date = "01Jan2000", ...)
    # Last updated 04Oct2016
    # This is a wrapper function for residueScorePlot which will take a list of proteins, plot their RBRID scores, and output the scores
    # as an attribute file for use with the Chimera protein visualization software for a desired .pdb file
    # If there are multiple uniprot IDs present in rbridset for a single protein, 
    # Ties argument can break ties and specify which isoform to use for plot. 
    # Ties are given as uniprot IDs in character vector, matching corresponding index in Proteins argument.
    # If Ties is left blank, plot defaults to isoform with highest coverage
    # If PDB = NA, listscoreplot will output plots using residue score plot with specified parameters. 
    # If PDB = PDB code, function will instead output attribute files for each protein that can be used to render attributes in Chimera
    # Wd specifies working directory to output attribute files for Chimera
  {
    library(GenomicRanges)
    
    if (length(Proteins) != length(Ties)) {Ties <- c(Ties, rep("", length(Proteins)- length(Ties)))} # Adds empty rows to Ties argument to match lengths
    
    lapply(1:length(Proteins), function(pos) 
    {
      candidate <- Proteins[pos]
      tiebreaker <- Ties[pos]
      
      if (is.na(PDB)) {  
        if (!any(grepl(candidate, range$symbol))) {return(paste0("Candidate ", candidate, " is not identified. Check Uniprot IDs"))} 
        else {
          
          
          plot <- residueScorePlot(candidate, uniprot = tiebreaker, range = range, peps = peps, ...)  
          return(plot)
        }}
      
      else if  (!is.na(PDB))  
      {
        # to make score file that gets inserted in PDB for coloring
        t <- residueScorePlot(candidate, range = range, peps = peps, return.values = T, no.negs = T)  
        t$label<-paste0(":",t$x)
        scoretable<- cbind(rep("", length(t$x)), t$label, round(t$y*10)/10)
        
        fileName <- file.path(wd, paste0(Descriptor, "_", Date,"_", candidate, "scoreSmooth_", PDB, ".txt"))
        write.table(scoretable,sep="\t",quote=F, file= fileName, row.names = F, col.names = F)
        
        header <- c(paste0("#  PDB entry ", PDB, " - " ,candidate, " RBR-ID RNA binding prediction score"),
                    "#  (4sU crosslinked peptide mass spectrometry depletion)",
                    "#  Method credit to Chongsheng He & Roberto Bonasio",
                    paste0("#", "  Experiment: ", Experiment),
                    paste0("#", "  Details: ", Descriptor, ". P-value = ", Pval),
                    "#  Use this file to assign the attribute to crystal structures in in Chimera with the",
                    "#  Define Attribute tool or the command defattr.",
                    paste0("#  Make sure to specify attribute only to ", candidate, " chain"),
                    "#",
                    paste0("attribute: score", candidate),
                    "match mode: 1-to-1",
                    "recipient: residues")
        
        fConn <- file(fileName)
        Lines <- readLines(fConn)
        writeLines(c(header, Lines), sep = "\n" , con = fConn)
        close(fConn)
        return(paste0("Attribute Output file: ", fileName))
      }
    })
  }
  
  
  
  ## ListAttrib.plot
  # 07Jun2018 - Added code to allow for assigning scores to specific models in chimera
  # 11Jun2019 - Added code to allow exporting missing residues as logical attribute file with return.gaps argument
  listAttrib.plot <- function ( DBREF.df, PDB = NA, range, peps, Pval = "1", Experiment = "", Descriptor = "", Date = "01Jan2000", out.file = "", return.gaps = F, ...)
    # Last updated 11 Jun 2018
    # This is a wrapper function for residueScorePlot which will take a list of proteins, plot their RBRID scores, and output the scores
    # as an attribute file for use with the Chimera protein visualization software for a desired .pdb file
    # Note, only use alphanumeric characters or underscores in experiment value, special characters interfere with Chimera attribute names
    # If PDB = PDB code, function will instead output attribute files for each protein that can be used to render attributes in Chimera
    # DBREF.df should contain mappings for which uniprotID is associated with which pdb chain, must be a rbridframe with variables "uniprotID" and "chainID", can also specify model number for each chain with "model" if there are conflicting chain identities
    # out.file specifies attribute file name + prefix
  {
    library(GenomicRanges)
    
    if (!any(names(DBREF.df) %in% "model")) {DBREF.df$model <- ""
    range$model <- ""}
    
    # Remove duplicate chains within models in DBREF.df
    DBREF.df$model.chain <- paste0(DBREF.df$model, ".", DBREF.df$chainId)
    dedup.index <- which(!duplicated(DBREF.df$model.chain))
    DBREF.df <- DBREF.df[dedup.index,]
    
    Score.list <- lapply(DBREF.df$model.chain, function(model.chain) {
      model <- str_split(model.chain, "[.]")[[1]][1]
      chain <- str_split(model.chain, "[.]")[[1]][2]
      
      temp.df <- DBREF.df[DBREF.df$model == model ,]
      uniprot <- temp.df$uniprotID[temp.df$chainId == chain]
      symbol <- as.character(temp.df$symbol[temp.df$chainId == chain])
      model.range <- range[range$model == model & range$chain == chain]
      
      # to make score file that gets inserted in PDB for coloring
      t <- residueScorePlot(uniprot = uniprot, range = model.range, peps = peps, return.values = T, no.negs = T, smooth = T, return.gaps = return.gaps)  
      if(return.gaps == T) {
        #remove NA values (residues with coverage)
        t <- as.rbrid.frame(t)
        t$boolean <- as.logical(t$missing)
        t$label<-paste0("#", model, ":",t$x, ".", chain)
        
        #assign blue transparent values for missing residues
        missing.table <- cbind(rep("", length(t$x)), t$label, t$boolean)
        return(missing.table)}
      
      t$label<-paste0("#", model, ":",t$x, ".", chain)
      # convert NA values to zero (Chimera attributes do not accept NA)
      t$y[is.na(t$y)] <- 0
      
      scoretable<- cbind(rep("", length(t$x)), t$label, round(t$y*10)/10) 
      
      return(scoretable)})
    
    scoretable <- do.call("rbind", Score.list)
    
    fileName <- out.file
    
    write.table(scoretable,sep="\t",quote=F, file= fileName, row.names = F, col.names = F)
    
    header <- c(paste0("#  PDB entry ", PDB, " - (all chains) RBR-ID RNA binding prediction score"),
                "#  (4sU crosslinked peptide mass spectrometry depletion)",
                paste0("#", "  Experiment: ", Experiment),
                paste0("#", "  Details: ", Descriptor, ". P-value = ", Pval),
                "#  Use this file to assign the attribute to crystal structures in Chimera with the",
                "#  Define Attribute tool or the command defattr.",
                "# Chain mapping key:")
    
    mapping <- paste("#","PDB", "model", "chain", "symbol", collapse = "\t")
    for (i in 1:nrow(DBREF.df))
    { mapping <- c(mapping, paste("#", DBREF.df[i,"structureId"], DBREF.df[i, "model"], DBREF.df[i,"chainId"], DBREF.df[i,"symbol"], collapse = "\t"))}
    gap.label <- ifelse(return.gaps, "gaps", "")
    header <- c(header, mapping, paste0("attribute: score_", PDB, "_", Experiment, gap.label),
                "match mode: 1-to-1",
                "recipient: residues")
    
    fConn <- file(fileName)
    Lines <- readLines(fConn)
    writeLines(c(header, Lines), sep = "\n" , con = fConn)
    close(fConn)
    return(paste0("Attribute Output file: ", fileName))
  } 
  
  
  ##Annotate peptides
  
  annotatePeptides<-function(index, rbrid, IPR.info)
  {
    # Roberto Bonasio / updated 04 Apr 2018
    ## Annotates overlapping IPRs for protein and peptide
    ## Computationally intensive, use parallelization if possible
    ## Takes 1 min for 5,000 proteins in mouse if no parallelization
    ## Takes 1.5 min for ~80,000 proteins with 30 nodes
    peptide<-rbrid[index,]
    hits<-which(IPR.info$uniprotID==peptide$uniprotID)  # finds hits in IPR table
    prot_IPR_hits<-paste(unique(IPR.info$iprID[hits]),collapse="|") # annotates with all IPR hits for the protein (regardless of position)
    
    # the following decides whether the actual peptide overlaps the start, end, or is contained within the various IPRs
    overlap_start<-peptide$start<=IPR.info$IPR_start[hits] & peptide$end>=IPR.info$IPR_start[hits]
    overlap_end<-peptide$start<=IPR.info$IPR_end[hits] & peptide$end>=IPR.info$IPR_end[hits]
    inside<-peptide$start>=IPR.info$IPR_start[hits] & peptide$end<=IPR.info$IPR_end[hits]
    
    pept_IPR_hits<-paste(unique(IPR.info$iprID[hits][overlap_start | overlap_end | inside]),collapse="|") # annotates IPRs that overlap
    output <- c(interpro_protein = prot_IPR_hits,
                interpro_peptide = pept_IPR_hits)
    return(output)
  }
  
  ## WhatSpecies
  
  Whatspecies <- function(string)
    # Function that takes a string and matches it to a species based on what cell lines 
    # contained within its substrings
  {
    # extract species information from dataset name
    human_cells <- "human|K562|HCT116|HepG2|Hela"
    mouse_cells <- "mouse|mESC|NIH3T3|NPC"
    
    if (grepl(human_cells, string, ignore.case = T)) {
      return("human")
    } else if (grepl(mouse_cells, string, ignore.case = T)) {
      return("mouse")
    } else stop("Cannot identify species")
    
  }
  
  ##F_statistic
  
  F_statistic <- function(rbrid, groups)
    # Function to calculate F-statistic for ANOVA test. between-group variability / within-group variability.
    # rbrid = rbrid.frame, groups = list of columns for each group
  { 
    allSamples <- as.numeric(unlist(groups))
    Inter <- sapply(groups, function(x) {
      num <- length(x)*(rowMeans(rbrid[x], na.rm = T)-rowMeans(rbrid[allSamples], na.rm = T))^2
      denom <- length(groups)-1
      var <- num/denom
      final <- sum(var, na.rm = T)
      return(final)})
    
    inter_sum <- sum(as.numeric(unlist(Inter)), na.rm = T) # sum variance across all groups
    
    
    Within <- sapply(groups, function(group) {
      
      
      foo <- sapply(1:length(group), function(rep) {
        num <- (rbrid[group[rep]] - rowMeans(rbrid[group], na.rm = T))^2
        denom <- length(group) - length(groups) # # of replicates - # of groups
        var <- num/denom
        return(var)})
      final <- sum(foo, na.rm = T)
    })
    Within_sum <- sum(as.numeric(Within), na.rm = T) # Sum variance within groups
    
    Fstat <- inter_sum/Within_sum
    return(Fstat)
  }
  