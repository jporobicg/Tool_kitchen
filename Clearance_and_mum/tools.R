
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Estimation of groth rate 'mum' for the holling type II equation implemented in Atlantis for each functional groups
##' @param length Vecto of length in cm
##' @param weight Vector of weight in gr
##' @param metric The default metric for the wet weight is 'mg' (miligrams), the other options are:  'g' for grams and 'Kg' for Kilograms
##' @param wet Logical vector for the type weight, the default the is wet weight
##' @param spw.rate Spawning rate, vector of values that represent the increase in waeight due to reproduction
##' @param mature Logical, vector of TRUE or FALSE to descrive if the FG at that weight/length its mature
##' @param AgeClass Age classes by cohort
##' @return A vector with the values of mum for each length (mgNd-1)
##' @author Demiurgo
mum.f <- function(length, weight, metric  = 'mg', wet = TRUE, spw.rate, mature, AgeClass){
    we  <- c(0, weight)
    ## transformation to weight in nitrogen
    we  <- weight2N(we, metric = metric)
    len <- c(0, length)
    d.w <- vector('numeric')
    for(i in 2 : length(we)){
        wet        <- ifelse(isTRUE(mature[i - 1]), we[i] * spw.rate, we[i]) ## increase in weight when mature
        d.w[i - 1] <- wet - we[i - 1] # change in weight.
    }
    d.l         <- diff(len)
    growth.rate <- d.w / d.l
    growth.rate <- growth.rate / (AgeClass * 365)
    return(growth.rate)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Estimation of the allometric relationship
##' @param length Vector of values of length
##' @param alfa Parameter
##' @param beta Parameter
##' @return a vector with the weight
##' @author Demiurgo
alometric <- function(length, alfa, beta){
    weight <- alfa  * length ^ beta
    return(weight)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Transform from conventional wet weight to gram of nitrogen
##' @param weight Vector of weight
##' @param metric The default metric for the wet weight is 'mg' (miligrams), the other options are:  'g' for grams and 'Kg' for Kilograms
##' @param wet Logical vector for the type weight, the default the is wet weight
##' @return a vector witht the weight in miligrams of nitrogen (mgN)
##' @author Demiurgo
weight2N <- function(weight, metric = 'mg', wet = TRUE, reserve = TRUE){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~            transform weight to Ng        ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if(metric == 'g')  weight <- weight * 1000
    if(metric == 'Kg') weight <- weight * 1000000
    ## Wet to dry weight
    if(wet) weight <- (weight / 20)
    ## From dry weight to nitrogen
    weight  <- weight / 5.7
    ## Structural weight
    weight  <- weight / 3.65
    ## If is reserve
    if(isTRUE(reserve)){
        weight  <- weight * 2.65
    }
    return(weight)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' \description{This function estimates clearance based on biological parameters.}
##' @title Clearance
##' @param fg Vector with the names of the functional groups.
##' @param speed Vector with the speed in m/h. If the speed is not available (NA), the speed will be calculated using the length-at-age class of the functional group.
##' @param len Vector with the length at age class of the functional groups.
##' @param height Vector of the ratio of height compared to the length at age of the functional group. If not provided, the estimation will use 1/5 of the total length.
##' @param ratio The ratio between the height and the width of the individual. If not provided, the width will be the same as the height.
##' @param time.l Vector proportion of time in a day invested by the species searching for food.
##' @param max.speed Maximum speed reported for the functional group. This sets the upper boundary for the calculation of the speed.
##' @param alfa Alometric coefficient for calculating mass based on length.
##' @param beta Alometric coefficient for calculating mass based on length.
##' @param by.group Arrange the output for easy manipulation.
##' @return A vector with the values of clearances for each functional group in mg/Nm³/day.
##' @author Demiurgo
estimateClearance <- function(fg = NULL, speed, len, height = NA, ratio = NULL, time.l = NULL, max.speed = NULL, alfa = NULL, beta = NULL, by.group = TRUE){
    ## General assumption, If I dont have the speed I will use the mass to calculate the speed
    mass   <- alometric((len * 100), alfa, beta) / 1000  ## in kilograms
    max.age <- max(table(fg))
    ## Assumption based on Sato et al 2007 Swimming speed
    speed <- ifelse(is.na(speed), mass ^ 0.27 * 3600, speed) # speed mh-1
    ## Adjusting the speed based on the functional group size
    name.fg <- unique(fg)
    for(g in 1 : length(unique(fg))){
        tmp.fg      <- which(fg %in% name.fg[g])
        max.spd.tmp <- max.speed[tmp.fg]
        speed.tmp   <- speed[tmp.fg]
        if(!all(is.na(max.speed[tmp.fg])) && any(max.speed > speed)){
            speed[tmp.fg] <- (speed.tmp/max(speed.tmp, na.rm=TRUE)) * max.spd.tmp
        }
    }
    ## Assume that the height is at least 1/5 of the length
    height <- ifelse(is.na(height), len / 5, len * height)
    ## The width can be the same as the height, but it can have a deformation (flat fish)
    if(!is.null(ratio)){
        width  <- height * ratio
    } else {
        width  <- height
    }
    ## The time that the fish is looking for food
    if(is.null(time.l)) time.l <- 0.5
    ## all the time is in hours, the idea is to have it per day and per volumen (d/m3)
    out <- speed * pi() * (height / 2) * (width / 2) * 24 * time.l
    if(!is.null(fg)){
        out <- data.frame(FG = fg, Clearance = out)
        if(isTRUE(by.group)){
            temp.out <- split(out, out$FG)
            m        <- as.vector(table(fg))
            out      <- matrix(NA, nrow = max(m), ncol = length(m))
            for(i in 1 : length(temp.out)){
                out[, i] <- c(temp.out[[i]][, 2], rep(NA, max.age - length(temp.out[[i]][, 1])))
            }
            out           <- t(out)
            rownames(out) <- names(temp.out)
            colnames(out) <- c(1 : (ncol(out)))
        }
    }
    return(out)
}
