##' @title Biomass for age groups
##' @param age.grp Age groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t Miligram to tons scalar
##' @param x.cn Carbon to Nitrogen transformation (Value/scalar)
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo
bio.age <- function(age.grp, nc.out, ctg, mg2t, x.cn){
    grp.bio <- NULL
    for(age in 1 : nrow(age.grp)){
        cohort <- NULL
        for(coh in 1 : age.grp[age, 'numcohorts']){
            name.fg <- paste0(age.grp$name[age], coh)
            b.coh   <- (ncdf4::ncvar_get(nc.out, paste0(name.fg, '_ResN'))  +
                        ncdf4::ncvar_get(nc.out, paste0(name.fg, '_StructN')))  *
                ncdf4::ncvar_get(nc.out, paste0(name.fg, '_Nums')) * mg2t * x.cn
            b.coh   <- apply(b.coh, 3, sum, na.rm = TRUE)
            cohort  <- cbind(cohort, b.coh)
        }
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, nrow(cohort)), FG  = as.character(age.grp$code[age]),
                                             Biomass  = rowSums(cohort, na.rm  = TRUE), Simulation = ctg))
    }
    return(grp.bio)
}

##' @title Biomass for age groups
##' @param pol.grp Biomass pool groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t Miligrams to tons scalar
##' @param x.cn Ration from nitrogen to carbon
##' @param box.info Information by box and layer
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo
bio.pool <- function(pol.grp, nc.out, ctg, mg2t, x.cn, box.info){
    grp.bio <- NULL
    for(pool in 1 : nrow(pol.grp)){
        name.fg <- paste0(pol.grp$name[pool], '_N')
        biom    <- ncdf4::ncvar_get(nc.out, name.fg)
        if(length(dim(biom))== 3){
           if(pol.grp$grouptype[pool] == 'LG_INF'){
               biom <- apply(biom, 3, '*', box.info$VolInf)
           }else{
               biom <- apply(biom, 3, '*', box.info$Vol)
           }
            biom <- apply(biom, 2, sum, na.rm = TRUE)
        } else {
            biom <- apply(biom, 2, function(x) x * box.info$info$Area)
            biom <- apply(biom, 2, sum, na.rm = TRUE)
        }
        biom    <- biom * mg2t * x.cn
        biom[1] <- biom[2]
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, length(biom)), FG  = as.character(pol.grp$code[pool]), Biomass  = biom, Simulation = ctg))
    }
    return(grp.bio)
}

##' @title Biomass for age groups
##' @param pwn.grp Biomass pools with age classes
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t mg C converted to wet weight in tonnes == 20 / 1000000000
##' @param x.cn Redfield ratio of C:N 5.7
##' @param box.info Information for each box and layer
##' @return a dataframe with the biomass for all the functional groups with Biomass pool age classes
##' @author Demiurgo
bio.pwn <- function(pwn.grp, nc.out, ctg, mg2t, x.cn, box.info){
    grp.bio <- NULL
    for(pwn in 1 : nrow(pwn.grp)){
        cohort <- NULL
        for(coh in 1 : pwn.grp[pwn, 'numcohorts']){
            name.fg <- paste0(pwn.grp$name[pwn], '_N', coh)
            b.coh   <- ncdf4::ncvar_get(nc.out, name.fg) * mg2t * x.cn
            if(length(dim(b.coh)) > 2){
                ## for Volumen
                b.coh   <- apply(b.coh, 3, '*', box.info$Vol)
            } else {
                ## for area
                b.coh   <- apply(b.coh, 2, '*', box.info$info$Area)
            }
            b.coh   <- apply(b.coh, 2, sum, na.rm = TRUE)
            cohort  <- cbind(cohort, b.coh)
        }
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, nrow(cohort)), FG  = as.character(pwn.grp$code[pwn]), Biomass  = rowSums(cohort, na.rm  = TRUE), Simulation = ctg))
    }
    return(grp.bio)
}

##' @title Calculate the value of a time serie based on the first value
##' @param df Data frame of biomass by FG
##' @param biomass Is biomass included
##' @param Vec True is the information in vector
##' @return A vector of relative values
##' @author Demiurgo
relative <- function(df, biomass = TRUE, Vec = NULL){
    if(biomass)       vector <- df$Biomass
    if(!is.null(Vec)) vector <- df[Vec]
    df$Relative  <- vector / vector[1]
    return(df)
}




.calc.mort <- function(SpeMort, PredMort, FG, Sums = FALSE){
    F.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'F', sep = '-')
    N.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'M1', sep = '-')
    P.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'M2', sep = '-')
    ## Calculating Fishing mortality
    Fmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% F.mort)), with = FALSE]
    colnames(Fmort) <- c('Time', paste0('Age0', 1 : (ncol(Fmort) - 1)))
    ## Calculating Natural mortality
    Nmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% N.mort)), with = FALSE]
    colnames(Nmort) <- c('Time', paste0('Age0', 1 : (ncol(Nmort) - 1)))
    ## Calculating Predation mortality
    Pmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% P.mort)), with = FALSE]
    colnames(Pmort) <- c('Time', paste0('Age0', 1 : (ncol(Pmort) - 1)))
    Fmort$Mtype <- 'Fishing Mortality'
    Nmort$Mtype <- 'Natural Mortality'
    Pmort$Mtype <- 'Predation Mortality'
    out <- rbind(Fmort, Nmort, Pmort)
    if(Sums){
        out <- out[, lapply(.SD, sum, na.rm=TRUE), by= Mtype]
    }
    return(out)
}



read.age <- function(grp, nc.file){
    n.coh <- grp$numcohorts
    init   <- sec <- num  <- end <- NULL
    r.init <- r.end  <- resn <- NULL
    for(coh in 1 : n.coh){
        name.fg <- paste0(grp$name, coh)
        nums    <- ncdf4::ncvar_get(nc.file, paste0(name.fg, '_Nums'))
        resN    <- ncdf4::ncvar_get(nc.file, paste0(name.fg, '_ResN'))
        nums.t  <- apply(nums, 3, sum, na.rm = TRUE)
        init    <- c(init, nums.t[1])
        sec     <- c(sec, nums.t[1 + which(diff(nums.t) > 0)[1]])
        end     <- c(end, tail(nums.t, 1))
        nums.t  <- log(mean(abs(nums.t - nums.t[1]) / nums.t[1], na.rm = TRUE), 10)
        resN.t  <- apply(resN, 3, sum, na.rm = TRUE)
        r.init  <- c(r.init,  resN.t[1])
        r.end   <- c(r.end, tail(resN.t, 1))
        resN.t  <- log(mean(abs(resN.t - resN.t[1]) /  resN.t[1], na.rm = TRUE), 10)
        num     <- c(num, nums.t)
        resn    <- c(resn, resN.t)
    }
    age <- 1 : coh
    out <- cbind(num, resn, age, init, sec, end, r.init, r.end)
    colnames(out) <- c('numbers',  'reserve', 'age', 'ninit', 'nsec', 'nend', 'rinit', 'rend')
    return(out)
}

##' @title Calculation of mortality level
##' @param mort.file Mortality File
##' @param target Target functional group
##' @return A list with two data frames with level of mortality and Predators and their predation level
##' @author Javier Porobic
morality.lvl <- function(mort.file, target){
    lines  <- readLines(mort.file)
    data1  <- c(1 : length(lines))[-c(1, grep('Time', lines)[ -1])] ## lines with wierd end
    header <- c(unlist(strsplit(lines[2], ' ')))
    vals   <- rle(1 : length(lines) %in% data1)
    idx    <- c(1, cumsum(vals$lengths))[which(vals$values)]
    ## indx: start = starting row of sequence, length = length of sequence (compare to s)
    indx   <- data.frame(start=idx, length=vals$length[which(vals$values)])[ - 1, ]
    result <- do.call(rbind, apply(indx, 1, function(x){
        return(fread(mort.file, nrows=x[2], skip = x[1] ))}))
    colnames(result) <- header
    tar.sp  <- result[which(result$Prey %in% target),  - c(1, 2)]
    mort.sp <- colMeans(tail(tar.sp, 10), na.rm = TRUE) / colMeans(tail(tar.sp, 10), na.rm = TRUE)[1]
    Pred.lvl <- ifelse(mort.sp['TotalPredMort'] > 0.8,  'Too high',
                ifelse(mort.sp['TotalPredMort'] > 0.6, 'High',
                ifelse(mort.sp['TotalPredMort'] > 0.3, 'Normal',
                ifelse(mort.sp['TotalPredMort'] > 0.1, 'Low', 'Too Low'))))
    ml.eff <- mort.sp['mL'] / mort.sp['TotalPredMort']
    mq.eff <- mort.sp['mQ'] / mort.sp['TotalPredMort']
    ml.lvl<- ifelse(ml.eff > 0.8,  'Too high',
             ifelse(ml.eff > 0.6, 'High',
             ifelse(ml.eff > 0.3, 'Normal',
             ifelse(ml.eff > 0.1, 'Low', 'Too Low'))))
    mq.lvl<- ifelse(mq.eff > 0.8,  'Too high',
             ifelse(mq.eff > 0.6, 'High',
             ifelse(mq.eff > 0.3, 'Normal',
             ifelse(mq.eff > 0.1, 'Low', 'Too Low'))))
    pred      <- names(which(mort.sp[- c(1 : 4)]>0))
    Predators <- data.frame(Predators = pred, Proportion = mort.sp[which(names(mort.sp) %in% pred)])
    factors   <- 0.1 / c(mort.sp['TotalPredMort'], ml.eff, mq.eff)
    Lvl.Mortality <- data.frame(Mortality = c('Pred', 'mL', 'mQ'), Level = c(Pred.lvl, ml.lvl, mq.lvl),
                                Factor = factors)
    return(list(Predators, Lvl.Mortality))
}
