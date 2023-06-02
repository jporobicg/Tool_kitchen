library(data.table)
library(ncdf4)
source('estimation/tools.R')
source('/home/por07g/Documents/Code_Tools/Package/ReactiveAtlantis/R/Util-text.R')
#-Run model
biom_indx <- fread('estimation/output01/outputNorthSeaBiomIndx.txt',header=T)

prm        <- readLines('estimation/01NorthSea_biol_fishing.prm', warn = FALSE)
grp <- fread('estimation/functionalGroups.csv')
colnames(grp) <- tolower(colnames(grp))
SpeMort       <- fread('estimation/output01/outputNorthSeaSpecificMort.txt', sep = ' ')
PredMort      <- fread('estimation/output01/outputNorthSeaSpecificPredMort.txt', sep = ' ')
Pred.level.f  <- 'estimation/output01/outputNorthSeaMortPerPred.txt'
nc.file       <- nc_open('estimation/output01/outputNorthSea.nc')


biom_indx <- biom_indx[ , -c(which(colnames(biom_indx) %in% paste0('Rel', c('DL', 'DR', 'DC', 'PB')))), with = FALSE]
## getting the biggest differences
##- Centered in 0
absolute_change <- biom_indx[, which(colnames(biom_indx) %in% paste0('Rel',grp$code)), with = FALSE]
absolute_change[which(absolute_change == 0, arr.ind = TRUE)] <- 1e-16
change_direction = biom_indx[, which(colnames(biom_indx) %in% paste0('Rel',grp$code)), with = FALSE]
change_direction <- colMeans(change_direction, na.rm = TRUE)
big.c   <- colMeans(abs(log(absolute_change, 10)), na.rm = TRUE)
score.c <- sum(big.c, na.rm = TRUE)
Changes <- list()
names.chng <- NULL
changes=1
for( changes in 1 :  5){
    grp.c   <- big.c[which.max(big.c)]
    ## Zoo or phyto
    target <- gsub('Rel', '', names(grp.c))
    cat(paste0(' -- Doint the change for : ', target, '\n'))

    ## Predation mortality
    mort.sp <- .calc.mort(SpeMort, PredMort, FG = target, Sums = TRUE)
#    debug(text2num)

    if(inf.grp$ispredator){
        ## Predator
        if(inf.grp$numcohorts > 1){
            ## AgeClass
            data    <-  read.age(grp=as.data.frame(inf.grp), nc.file)
            chnge   <- as.matrix(abs(data[, 1 : 2]))
            chnge   <- which(matrix(chnge %in% head(sort(chnge, TRUE), 2), nr = nrow(chnge)), arr.ind = TRUE)
            mq      <- text2num(text=prm, pattern = paste0(target, '_mQ'), FG='look', Vector=TRUE)
            ml      <- text2num(text=prm, pattern = paste0(target, '_mL'), FG='look', Vector=TRUE)
            mature  <- text2num(text=prm, pattern = paste0(target, '_age_mat'), FG='look')
            mum     <- text2num(text=prm, pattern = paste0('^mum_', target), FG='look', Vector=TRUE)
            c       <- text2num(text=prm, pattern = paste0('C_', target), FG='look', Vector=TRUE)
            stage   <- ifelse(c(chnge[, 2]) > mature$Value, 'adult', 'juvenile')
            st.pos  <- ifelse(c(chnge[, 2]) > mature$Value, 2, 1)
            new.par <- list(ml = c(NA, NA), mq  = c(NA, NA),
                            mum = rep(NA, inf.grp$numcohorts),
                            c   = rep(NA, inf.grp$numcohorts))
            ## NUMBERS or Reserve?
            for(focus in 1 : 2){
                if(chnge[focus, 2] == 1){ ## Numbers
                    ## Reproduction
                    if(chnge[focus, 1] == 1){
                        cat('It\'s necessary to reduce the reproduction')
                    }
                    ## Change predation
                    cat('It\'s mortality what we want to change first\n')
                    if(sum(mort.sp[3, -c(1, 2)]) == 0){
                        cat('WARNING!: \n The functional group ', inf.grp$name, ' has 0 predation mortality\n')
                    }
                    ## mortality
                    if(inf.grp$grouptype == 'MAMMAL' && stage[focus] == 'adult'){
                        ## increase Cuadratic mortality due to densodependency
                        new.par$mq[st.pos[focus]] <- (data[chnge[focus, 1], 'nsec'] - data[chnge[focus, 1], 'ninit']) / (data[chnge[focus, 1], 'nend']) ^ 2
                    } else {
                        ## in fish you alwyas spect to have more ML that MQ
                        if(mq[, st.pos[focus]] >=  ml[, st.pos[focus]]) {
                            ml.tmp <- (data[chnge[focus, 1], 'nsec'] - data[chnge[focus, 1], 'ninit']) / data[chnge[focus, 1], 'nend']
                            if(ml.tmp < 0) {
                                ## making sure that the mortality value decrease
                                ml.tmp <- ml[st.pos[focus]] + ml.tmp
                            }
                            new.par$ml[st.pos[focus]] <- ml.tmp
                        } else {
                            mq.tmp <- (data[chnge[focus, 1], 'nsec'] - data[chnge[focus, 1], 'ninit']) / (data[chnge[focus, 1], 'nend']) ^ 2
                            if(mq.tmp < 0) {
                                ## making sure that the mortality value decrease
                                mq.tmp <- mq[st.pos[focus]] + mq.tmp
                            }
                            new.par$mq[st.pos[focus]] <- mq.tmp
                        }
                    }
                } else {
                    ## RESERVE
                    cat('changing the reserve')
                }


            }


            ## Increasing mortality



            prey <- PredMort[which(PredMort$Group %in% target), ]
            prey[, which(PredMort$Group %in% target)]
            #which(target == Predators)
            ## Change mortality

        }
        if(inf.grp$numcohorts == 1){
            ## Variables
            mq      <- text2num(text=prm, pattern = paste0(target, '_mQ'), FG='look', Vector = FALSE)
            ml      <- text2num(text=prm, pattern = paste0(target, '_mL'), FG='look', Vector = FALSE)
            ## level of predation pressure
            debug(morality.lvl)
            lvl.mort <- morality.lvl(Pred.level.f, target=inf.grp$target)


        }
    }

    Changes[[changes]] <- new.par
    names.chng <- c(names.chng, target)
    big.c   <- big.c[ - which.max(big.c)]
}

aggregate(data[, 1:2], list(data$age), mean)

by(data, 2, mean)
Type of FG [Age Class] or [Biomass pool]

if age class
 is increase [Num] [Reserve]

        -  - if [Number]
        -  -  - - Mortality
        -  -  - - Predatio
        -  -  -  - Reproduction

-  - if[Reserve]
-  -  - [mum]
-  -  - [Predation]



out <- read.age(inf.grp,  nc.file)


        scalrs <- function(Vector){
            for(i in 1 : length(Vector)){

            if(
grp <- as.data.frame(inf.grp)
nc.file

















par01 <- read /home/por07g/Documents/2019/Oil_spill/Salish_Sea_Atlantis_model/01SS_biology.prm

## Pproducer
grp.phy <- grp$code[c(grep('PHY', grp$grouptype), which(grp$grouptype %in% c('SEAGRASS')))]

b.phy <- colMeans(absolute_change[, which(colnames(absolute_change) %in% paste0('Rel', grp.phy))] - 1)

colMeans(absolute_change[, c(46, 47)])
score.phy <- colSums(abs((absolute_change[-1, which(colnames(absolute_change) %in% paste0('Rel', grp.phy))] - 1)))

t.phy <- b.phy[which.max(abs(b.phy))]


grp.zoo <- grp$code[c(grep('ZOO', grp$grouptype), which(grp$grouptype %in% c('PWN', 'LG_INF', 'SM_INF')))]



b.zoo <- colMeans(absolute_change[-1, which(colnames(absolute_change) %in% paste0('Rel', grp.zoo))] - 1)

score.zoo <- colSums(abs((absolute_change[-1, which(colnames(absolute_change) %in% paste0('Rel', grp.zoo))] - 1)))

t.zoo <- b.zoo[which.max(abs(b.zoo))]


## Zoo or phyto
focus  <- c(t.zoo, t.phy)[which.max(c(t.zoo, t.phy))]
target <- gsub('Rel', '', names(focus))

## Predation mortality
mort.sp <- .calc.mort(SpeMort, PredMort, FG = target)


definition_data_frame = grp
target_group = target
##if(focus > 0){
## Not enough mortality
    info_group <- function(target_group, definition_data_frame){
        inf.grp <- cbind(target_group, definition_data_frame[which(definition_data_frame$code %in% target), c('name','grouptype', 'ispredator', 'numcohorts')])
        ## dis function identified what classification have the target functional group
        ## This classification is just for calibration porpouse
        inf.grp$group <- if( length(grep('ZOO', inf.grp$grouptype)) !=0 || inf.grp$grouptype %in% c('CEP', 'PWN', 'LG_INF', 'SM_INF')){
                       'route')
                   } else{
                       cat('npe')
                   }
    grp.phy <- grp$code[c(grep('PHY', grp$grouptype), which(grp$grouptype %in% c('SEAGRASS')))]
