#This script will contain all the functions to perform multiple alignment and circularization
#using dynamic time warping


#this function makes global and local alignments from two univariate series, decides which type of alignment is best
#and returns a dataframe with both query sequences and the result sequence
#it has 5 arguments, the first 3 are required:
    # - df1, df2: two dataframes with the series data to be aligned
    # - dataColumns: a vector of length 2, indicating which columns from each df to be used for the alignment
    # - stepPattern: a list of length 4 with the different step patterns that can be used for the alignment, for more 
    #   information see the dtw package
    # - showDistPlot: a boolean, if TRUE then it will also plot the normalized distance in the alignment over the 
    #   sliding for the local alignment and the global alignment
#Returns: a dataframe with a time, query sequence, ref sequence and aligned sequence columns
pairwiseUnivariateAlignment <- function(df1, df2, dataColumns, stepPattern = NULL, showDistPlot = FALSE){
    
    #argument check
    if(is.null(stepPattern)){
        stepPattern <- list(asymmetricP05, asymmetricP05, asymmetricP05, symmetricP2)
    }else if(length(stepPattern) != 4){
        stop(print('Please provide 4 step patterns as a vector or leave it empty for default settings.'))
    }
    
    if(missing(df1) |  missing(df2) |  missing(dataColumns)){
        stop(print('You are missing the minimum required arguments: df1, df2 and/or dataColumns'))
    }
    if(!is.data.frame(df1) | !is.data.frame(df2)){
        stop(print('df1 and df2 need to be dataframes objects'))
    }
    if(!is.vector(dataColumns)){
        stop(print('dataColumns needs to be a vector object'))
    }
    if(length(dataColumns) != 2){
        stop(print('dataColumns needs to be a vector of length 2, each number corresponds to a column of df1 and df2'))
    }
    
    
    #list where we store all the alignment objects
    aligList <- list()
    #list counter
    k <- 0
    
    #query enters the reference from the left--------------------------------------------
    i <- nrow(df1)
    while(i >= 1){
        #dataframe with a fraction of the whole df data
        tempdf1 <- df1[i:nrow(df1), ]
        
        #alignment dtw function call
        alignment <- try(dtw(tempdf1[,dataColumns[1]], 
                             df2[,dataColumns[2]], #reference
                             step.pattern = stepPattern[[1]], 
                             keep.internals  = TRUE,  
                             open.end = TRUE, 
                             open.begin = FALSE), silent = TRUE)
        
        #if the alignment fails (sometimes sequences are too short, we create a fake alignment object (list of 12),
        #being the 12th element a very big number that will count as normalized distance)
        if(class(alignment) == 'try-error'){
            errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
            #increase the size of the window
            i <- i - 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- errorList
        }else{
            #increase the size of the window
            i <- i - 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- alignment
        }
        
        #just in case an extra break loop condition
        if(length(aligList) == nrow(df1)){
            break
        }
    }
    
    #query enters the reference from the right--------------------------------------------
    i <- 1
    while(i <= nrow(df1)){
        #dataframe with the fraction of the whole df data
        tempdf1 <- df1[1:i, ]
        
        #alignment dtw call
        alignment <- try(dtw(tempdf1[,dataColumns[1]], 
                             df2[,dataColumns[2]], 
                             step.pattern = stepPattern[[2]], 
                             keep.internals  = TRUE,  
                             open.end = FALSE, 
                             open.begin = TRUE), silent = TRUE)
        
        #if the alignment fails (sometimes sequences are too short, we create a fake alignment object (list of 12),
        #being the 12th element a very big number that will count as normalized distance)
        if(class(alignment) == 'try-error'){
            errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
            #increase the size of the window
            i <- i + 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- errorList
        }else{
            #increase the size of the window
            i <- i + 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- alignment
        }
        
        #just in case an extra break loop condition
        if(length(aligList) == (nrow(df1)*2)){
            break
        }
    }
    
    #open local alignment---------------------------------------------------------------
    alignment <- try(dtw(df1[,dataColumns[1]], 
                         df2[,dataColumns[2]], 
                         step.pattern = asymmetricP05, 
                         keep.internals  = TRUE,  
                         open.end = TRUE, 
                         open.begin = TRUE), silent = TRUE)
    
    #if the alignment fails (sometimes sequences are too short, we create a fake alignment object (list of 12),
    #being the 12th element a very big number that will count as normalized distance)
    if(class(alignment) == 'try-error'){
        errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
        #counter for the list of alignments
        k <- k + 1
        aligList[[k]] <- errorList
    }else{
        #counter for the list of alignments
        k <- k + 1
        aligList[[k]] <- alignment
    }
    
    #global alignment---------------------------------------------------------------------
    alignment <- try(dtw(df1[,dataColumns[1]], 
                         df2[,dataColumns[2]], #reference
                         step.pattern = stepPattern[[4]], 
                         keep.internals  = TRUE,  
                         open.end = FALSE, 
                         open.begin = FALSE), silent = TRUE)
    
    #if the alignment fails (sometimes sequences are too short, we create a fake alignment object (list of 12),
    #being the 12th element a very big number that will count as normalized distance)
    if(class(alignment) == 'try-error'){
        errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
        #counter for the list of alignments
        k <- k + 1
        aligList[[k]] <- errorList
    }else{
        #counter for the list of alignments
        k <- k + 1
        aligList[[k]] <- alignment
    }
    
    #model evaluation as function of sliding---------------------------------------------
    if(showDistPlot){
        require(ggplot2)
        df <- data.frame(mod = 1:length(aligList), dist = sapply(aligList, '[[', 12))
        
        p <- ggplot(data = df, aes(x = mod, y = dist)) + 
            geom_point() + 
            coord_cartesian() + 
            theme(axis.title.x = element_text(size= 14, face = NULL), 
                  axis.title.y = element_text(size= 14, face = NULL), 
                  axis.text.x = element_text(size= 14),
                  axis.text.y = element_text(size= 14),
                  axis.line.x = element_line(color="black", size = 1),
                  axis.line.y = element_line(color="black", size = 1),
                  panel.grid.major = element_blank(),                                                 
                  panel.grid.minor = element_blank(),                                                  
                  panel.background = element_blank()) + 
            coord_cartesian(ylim = c(0,1))
        plot(p)
    }
    
    #check best bitting alignment ----------------------------------------------------
    #defining the indexes of the number of models per type of model
    modelLeft <- c(1:nrow(df1))
    modelRight <- c((nrow(df1) + 1):(nrow(df1)*2))
    modelLocal <- length(aligList) - 1
    modelGlobal <- length(aligList) 
    
    #finding the best alignment
    smallestDistMod <- which(sapply(aligList, '[[', 12) == min(sapply(aligList, '[[', 12)))
    
    if(length(smallestDistMod > 1)){
        smallestDistMod <- smallestDistMod[length(smallestDistMod)]
    }
    
    #creating the profile from the alignment of the two dfs in different ways depending
    #on which type of model is best
    #LOCAL -----------
    if(smallestDistMod == modelLocal){
        print('LOCAL alignment')
        profileDf <- data.frame(time = c(1:max(nrow(df1), nrow(df2))), 
                                dfref = NA, 
                                dfquery = NA)
        #get alignment with smallest normalized distance
        model <- aligList[[smallestDistMod]]
        
        #match the indices from both sequences
        for(ind in unique(model$index2)){
            dfQueryData <- df1[model$index1[which(model$index2 == ind)], dataColumns[1]]
            dfRefData <- df2[ind, dataColumns[2]]
            
            profileDf[ind, 2:3] <- c(mean(dfRefData), mean(dfQueryData))
        }
        
        #add missing points from beginning
        if(!(1 %in% unique(model$index2))){
            begDf <- data.frame(time = NA, 
                                dfref = df2[c(1:c(unique(model$index2)[1]-1)), dataColumns[2]], 
                                dfquery = NA)
            profileDf <- rbind(begDf, profileDf)
        }
        
        #add missing points from the end
        if(!(nrow(df2) %in% unique(model$index2))){
            endDf <- data.frame(time = NA, 
                                dfref = df2[c((model$index2[length(model$index2)]+1):nrow(df2)), dataColumns[2]], 
                                dfquery = NA)
            profileDf <- rbind(profileDf, endDf)
        }
        
        #remove double NAs
        profileDf <- profileDf[rowSums(is.na(profileDf[2:3])) != 2, ]
        profileDf$time <- 1:nrow(profileDf)
    
    # from the LEFT -------------
    }else if(smallestDistMod %in% modelLeft){
        print('LEFT coming alignment')
        #how many timepoints are left out
        hanglength <- nrow(df1) - which(modelLeft == smallestDistMod)
        #create our profile df
        profileDf <- data.frame(time = c(1:(max(nrow(df1), nrow(df2)) + hanglength)), 
                                dfref = NA, 
                                dfquery = NA)
        #extract the best model
        model <- aligList[[smallestDistMod]]
        profileDf[1:(nrow(df1) - length(unique(model$index1))), 'dfquery'] <- df1[1:(nrow(df1) - length(unique(model$index1))), dataColumns[1]]
        
        continuation <- which(is.na(profileDf$dfquery))[1]-1
        #do the index alignment
        for(ind in unique(model$index2)){
            dfQueryData <- df1[model$index1[which(model$index2 == ind)] + continuation, dataColumns[1]]
            dfRefData <- df2[ind, dataColumns[2]]
            
            profileDf[ind + continuation, 2:3] <- c(mean(dfRefData), mean(dfQueryData))
        }
        
        hangLeftDf <- data.frame(time = NA, 
                                 dfref = df2[c((unique(model$index2)[length(unique(model$index2))]+1):nrow(df2)),dataColumns[2]],
                                 dfquery = NA)
        profileDf <- rbind(profileDf, hangLeftDf)
        
        profileDf <- profileDf[rowSums(is.na(profileDf[2:3])) != 2, ]
        
        profileDf$time <- c(1:nrow(profileDf))
        
    # from the RIGHT ----------------------
    }else if(smallestDistMod %in% modelRight){
        print('RIGHT coming alignment')
        #how many timepoints are left out
        hanglength <- nrow(df1) - which(modelRight == smallestDistMod)
        #create our profile df
        profileDf <- data.frame(time = c(1:(max(nrow(df1), nrow(df2)) + hanglength)), 
                                dfref = NA, 
                                dfquery = NA)
        #extract the best model
        model <- aligList[[smallestDistMod]]
        
        profileDf[1:(nrow(df2)-length(unique(model$index2))), 'dfref'] <- df2[1:(nrow(df2)-length(unique(model$index2))),dataColumns[2]]
        
        continuation <- which(is.na(profileDf$dfref))[1]
        
        #do the index alignment
        for(ind in unique(model$index2)){
            dfQueryData <- df1[model$index1[which(model$index2 == ind)], dataColumns[1]]
            dfRefData <- df2[ind, dataColumns[2]]
            
            profileDf[ind , 2:3] <- c(mean(dfRefData), mean(dfQueryData))
        }
        
        hangRightDf <- data.frame(time = NA, 
                                  dfref = NA,
                                  dfquery = df1[c((which(modelRight == smallestDistMod)+1):nrow(df1)),dataColumns[1]])
        profileDf <- rbind(profileDf, hangRightDf)
        profileDf <- profileDf[rowSums(is.na(profileDf[2:3])) != 2, ]
        
        profileDf$time <- c(1:nrow(profileDf))
        
    #GLOBAL --------------------------
    }else if(smallestDistMod == modelGlobal){
        print('GLOBAL alignment')
        profileDf <- data.frame(time = c(1:max(nrow(df1), nrow(df2))), 
                                dfref = NA, 
                                dfquery = NA)
        model <- aligList[[smallestDistMod]]
        
        
        for(ind in unique(model$index2)){
            dfQueryData <- df1[model$index1[which(model$index2 == ind)], dataColumns[1]]
            dfRefData <- df2[ind, dataColumns[2]]
            
            profileDf[ind, 2:3] <- c(mean(dfRefData), mean(dfQueryData))
        }
        profileDf$time <- 1:nrow(profileDf)
    }
    
    profileDf$uni <- apply(profileDf[,2:3], 1, mean, na.rm = T)
    return(profileDf)
}


#This function takes a list of time-series sequences and creates a sequence profile by multi-aligning all of them. This
#profile can be later used to align the same or other sequences to it to create a multialignment. This profile is built
#by calculating the distance matrix of the given sequence list, pair-wise aligning the two most similar sequences 
#calling the previous function and using the mean sequence between both. Then the distance matrix is calculated again
#with the new pair-wise alignment in it. This is done until one sequence (the profile is left). Since this is just
#a fancy while loop with the arguments for the pair-wise alignment function you might consider coding this part yourself
#to fit your multi-alignment needs
#it has 6 arguments, the first 2 are required:
# - dfList: list of dataframes that contain the sequences to be aligned
# - dataColumns: a vector of length 2, indicating which columns from each df to be used for the alignment
# - stepPattern: a list of length 4 with the different step patterns that can be used for the alignment, for more 
#   information see the dtw package
# - showDistPlot: a boolean, if TRUE then it will also plot the normalized distance in the alignment over the 
#   sliding for the local alignment and the global alignment
# - showDendrogram = plots the dendrogram of the distance matrix
# - showAlignmentPlot = plots the alignment between the two sequences and the new aligned-sequence
#Returns: a dataframe with a time, query sequence, ref sequence and aligned sequence columns
multiAlignmentUnivariateProfile <- function(dfList, dataColumns, stepPattern = NULL, 
                                            showDistPlot = FALSE, showDendrogram = FALSE, showAlignmentPlot = FALSE){
    
    
    #argument check
    if(is.null(stepPattern)){
        stepPattern <- list(asymmetricP05, asymmetricP05, asymmetricP05, symmetricP2)
    }else if(length(stepPattern) != 4){
        stop(print('Please provide 4 step patterns as a vector or leave it empty for default settings.'))
    }
    
    if(missing(dfList) |  missing(dataColumns)){
        stop(print('You are missing the minimum required arguments: dfList and/or dataColumns'))
    }else if(length(dfList) < 2){
        stop(print('You need at least 2 dataframes with data in the list'))
    }
    
    if(!is.list(dfList)){
        stop(print('df1 and df2 need to be dataframes objects'))
    }
    if(!is.vector(dataColumns)){
        stop(print('dataColumns needs to be a vector object'))
    }
    if(length(dataColumns) != 2){
        stop(print('dataColumns needs to be a vector of length 2, 
                   each number corresponds to a column of df1 and df2 for pairwise alignment'))
    }
    
    require(dtw)
    
    while(length(dfList) > 1){
        for(i in 1:length(dfList)){
            if(i == 1){
                distMatrix <- matrix(data = NA, nrow = length(dfList), ncol = length(dfList))
                colnames(distMatrix) <- names(dfList)
                rownames(distMatrix) <- names(dfList)
            }
            for(j in 1:length(dfList)){
                if(j > i){
                    next
                }
                if('uni' %in% colnames(dfList[[i]])){
                    distMatrix[i,j] <- dtw(x = na.omit(dfList[[i]][,which(colnames(dfList[[i]]) == 'uni')]), 
                                           y = na.omit(dfList[[j]][,dataColumns[2]]), 
                                           distance.only = T)$distance
                }else if('uni' %in% colnames(dfList[[j]])){
                    distMatrix[i,j] <- dtw(x = na.omit(dfList[[i]][,dataColumns[1]]), 
                                           y = na.omit(dfList[[j]][,which(colnames(dfList[[j]]) == 'uni')]), 
                                           distance.only = T)$distance
                }else if('uni' %in% colnames(dfList[[i]]) & 'uni' %in% colnames(dfList[[j]])){
                    distMatrix[i,j] <- dtw(x = na.omit(dfList[[i]][,which(colnames(dfList[[i]]) == 'uni')]), 
                                           y = na.omit(dfList[[j]][,which(colnames(dfList[[j]]) == 'uni')]), 
                                           distance.only = T)$distance
                }else{
                    distMatrix[i,j] <- dtw(x = na.omit(dfList[[i]][,dataColumns[1]]), 
                                           y = na.omit(dfList[[j]][,dataColumns[2]]), 
                                           distance.only = T)$distance
                }
            }
        }
        diag(distMatrix) <- NA
        
        #find the most similar cells
        indexes <- which(distMatrix == min(distMatrix, na.rm = T), arr.ind = T)
        dfname1 <- rownames(distMatrix)[indexes[,1]]
        dfname2 <- colnames(distMatrix)[indexes[,2]]
        
        if(showDendrogram){
            if(length(dfList) > 2){
                cluster <- hclust(d = as.dist(distMatrix))
                plot(cluster)
            }
        }
        
        
        #create a temporary datalist from which the dfs that we use for alignment are removed
        #and the profiles of those are added
        #dfs that we will use for the next alignment
        df1 <- dfList[[which(names(dfList) == dfname1)]]
        df2 <- dfList[[which(names(dfList) == dfname2)]]
        #remove them from the datalist
        dfList <- dfList[-c(which(names(dfList) == dfname1), which(names(dfList) == dfname2))]
        
        #swap positions if df2 is shorter than df1
        if(nrow(df1) > nrow(df2)){
            tdf1 <- df1
            tdf2 <- df2
            df1 <- tdf2
            df2 <- tdf1
        }
        
        #add alignment to the temporary datalist
        print(paste(dfname1, 'and' ,dfname2))
        
        if('uni' %in% colnames(df1) & 'uni' %in% colnames(df2)){
            profileDf <- pairwiseAlignment(df1 = df1, 
                                           df2 = df2, 
                                           dataColumns = c(which(colnames(df1) == 'uni'), 
                                                           which(colnames(df2) == 'uni')), 
                                           stepPattern = stepPattern,
                                           showDistPlot = showDistPlot)
        }else if('uni' %in% colnames(df1)){
            profileDf <- pairwiseAlignment(df1 = df1, 
                                           df2 = df2, 
                                           dataColumns = c(which(colnames(df1) == 'uni'), 
                                                           dataColumns[2]), 
                                           stepPattern = stepPattern,
                                           showDistPlot = showDistPlot)
        }else if('uni' %in% colnames(df2)){
            profileDf <- pairwiseAlignment(df1 = df1, 
                                           df2 = df2, 
                                           dataColumns = c(dataColumns[1], 
                                                           which(colnames(df2) == 'uni')), 
                                           stepPattern = stepPattern,
                                           showDistPlot = showDistPlot)
        }else{
            profileDf <- pairwiseAlignment(df1 = df1, 
                                           df2 = df2, 
                                           dataColumns = dataColumns, 
                                           stepPattern = stepPattern,
                                           showDistPlot = showDistPlot)
        }
        
        if(showAlignmentPlot){
            pp <- ggplot(data = profileDf, aes_string(x = names(profileDf)[1], y = names(profileDf)[4])) + 
                geom_point() + 
                coord_cartesian(ylim = c(-1,1)) + 
                ggtitle('Unified sequences')
            p1 <- ggplot(data = profileDf, aes_string(x = names(profileDf)[1], y = names(profileDf)[3])) + 
                geom_point() + 
                coord_cartesian(ylim = c(-1,1)) + 
                ggtitle('Sequence1')
            p2 <- ggplot(data = profileDf, aes_string(x = names(profileDf)[1], y = names(profileDf)[2])) + 
                geom_point() + 
                coord_cartesian(ylim = c(-1,1))+ 
                ggtitle('Sequence2')
            
            plotList <- list(pp, p2, p1)
            gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)
            plot(gridplot)
        }
        
        dfList[[length(dfList) + 1]] <- profileDf
        names(dfList)[length(dfList)] <- paste0(dfname1, dfname2)
    }
    return(dfList[[1]])
}



circularizeSequenceUnivariate <- function(sequenceDf, dataColumn, stepPattern = NULL, showDistPlot = FALSE){
    
    
    #argument check
    if(is.null(stepPattern)){
        stepPattern <- list(asymmetricP05, asymmetricP05)
    }else if(length(stepPattern) != 2){
        stop(print('Please provide 2 step patterns as a vector or leave it empty for default settings.'))
    }
    
    if(missing(sequenceDf) | missing(dataColumn)){
        stop(print('You are missing the minimum required arguments: sequenceDf and/or dataColumn'))
    }
    if(!is.data.frame(sequenceDf)){
        stop(print('seqeunceDf needs to be dataframes objects'))
    }
    if(!is.integer(dataColumn)){
        stop(print('dataColumn needs to be an integer object'))
    }
    if(length(dataColumn) > 1){
        stop(print('dataColumns needs to be a integer (length = 1)'))
    }
    
    
    require(dtw)
    aligList <- list()
    k <- 0
    
    #cut the profile in half and try to fit it from the left to the complete profile
    cutDfLeft <- sequenceDf[round(nrow(sequenceDf)/2):nrow(sequenceDf),]
    
    #query enters the reference from the left
    i <- nrow(cutDfLeft)
    while(i >= 1){
        #dataframe with the fraction of the whole profile sequence
        tempcutDfLeft <- cutDfLeft[i:nrow(cutDfLeft), ]
        
        #alignment dtw call
        alignment <- try(dtw(tempcutDfLeft[,dataColumn], 
                             sequenceDf[,dataColumn], #reference
                             step.pattern = stepPattern[[1]], 
                             keep.internals  = TRUE,  
                             open.end = TRUE, 
                             open.begin = FALSE), silent = TRUE)
        
        if(class(alignment) == 'try-error'){
            errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
            #increase the size of the window
            i <- i - 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- errorList
        }else{
            #increase the size of the window
            i <- i - 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- alignment
        }
        
        #just in case an extra break loop condition
        if(length(aligList) == nrow(sequenceDf)){
            break
        }
    }
    
    #cut the profile in half and try to fit it from the right to the complete profile
    cutDfRight <- sequenceDf[1:round(nrow(sequenceDf)/2),]
    
    #query enters the reference from the right
    i <- 1
    while(i <= nrow(cutDfRight)){
        #dataframe with the fraction of the whole cell data
        tempcutDfRight <- cutDfRight[1:i, ]
        
        #alignment dtw call
        alignment <- try(dtw(tempcutDfRight[,dataColumn], 
                             sequenceDf[,dataColumn], #reference
                             step.pattern = stepPattern[[2]], 
                             keep.internals  = TRUE,  
                             open.end = FALSE, 
                             open.begin = TRUE), silent = TRUE)
        
        if(class(alignment) == 'try-error'){
            errorList <- list(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000)
            #increase the size of the window
            i <- i + 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- errorList
        }else{
            #increase the size of the window
            i <- i + 1
            #counter for the list of alignments
            k <- k + 1
            aligList[[k]] <- alignment
        }
        
        #just in case an extra break loop condition
        if(length(aligList) == (nrow(sequenceDf)*2)){
            break
        }
    }
    
    if(showDistPlot){
        require(ggplot2)
        df <- data.frame(mod = 1:length(aligList), dist = sapply(aligList, '[[', 12))
        
        p <- ggplot(data = df, aes(x = mod, y = dist)) + 
            geom_point() + 
            coord_cartesian() + 
            theme(axis.title.x = element_text(size= 14, face = NULL), 
                  axis.title.y = element_text(size= 14, face = NULL), 
                  axis.text.x = element_text(size= 14),
                  axis.text.y = element_text(size= 14),
                  axis.line.x = element_line(color="black", size = 1),
                  axis.line.y = element_line(color="black", size = 1),
                  panel.grid.major = element_blank(),                                                 
                  panel.grid.minor = element_blank(),                                                  
                  panel.background = element_blank()) + 
            coord_cartesian(ylim = c(0,1))
        plot(p)
    }
    
    modelLeft <- c(1:nrow(cutDfLeft))
    modelRight <- c((nrow(cutDfRight) + 1):(nrow(sequenceDf)))
    
    #finding the alignment with smalles normalized distance
    smallestDistMod <- which(sapply(aligList, '[[', 12) == min(sapply(aligList, '[[', 12)))
    
    if(length(smallestDistMod) > 1){
        smallestDistMod <- smallestDistMod[length(smallestDistMod)]
    }
    
    if(smallestDistMod %in% modelLeft){
        print('LEFT coming alignment')
        #how many timepoints are left out
        hanglength <- nrow(cutDfLeft) - which(modelLeft == smallestDistMod)
        #create our profile df
        profileDf <- data.frame(time = c(1:nrow(profileDf)), 
                                cellref = profileDf[,dataColumn], 
                                cellquery = NA)
        model <- aligList[[smallestDistMod]]
        
        #do the index alignment
        for(ind in unique(model$index2)){
            cellQueryData <- cutDfLeft[(model$index1[which(model$index2 == ind)] + 
                                                    nrow(cutDfLeft) - 
                                                    length(unique(model$index1))), dataColumn]
            
            cellRefData <- sequenceDf[ind, dataColumn]
            profileDf[ind, 2:3] <- c(mean(cellRefData), mean(cellQueryData))
        }
        
        #remove the cut timepoints from the end
        profileDf <- profileDf[-rev(1:nrow(profileDf))[1:length(unique(model$index1))],]
        
        profileDf$uni <- apply(profileDf[,2:3], 1, mean, na.rm = T)
        profileDf$time <- c(1:nrow(profileDf))
        
        return(profileDf)
        
    }else if(smallestDistMod %in% modelRight){
        print('RIGHT coming alignment')
        #how many timepoints are left out
        hanglength <- nrow(cutDfRight) - which(modelRight == smallestDistMod)
        #create our profile df
        profileDf <- data.frame(time = c(1:nrow(sequenceDf)), 
                                cellref = sequenceDf[,dataColumn], 
                                cellquery = NA)
        #select the model
        model <- aligList[[smallestDistMod]]
        
        for(ind in unique(model$index2)){
            cellQueryData <- cutDfRight[(model$index1[which(model$index2 == ind)]), dataColumn]
            cellRefData <- sequenceDf[ind, dataColumn]
            
            profileDf[ind, 2:3] <- c(mean(cellRefData), mean(cellQueryData))
        }
        
        #remove the cut timepoints from the begining
        profileDf <- profileDf[-c(1:length(unique(model$index1))),]
        profileDf$uni <- apply(profileDf[,2:3], 1, mean, na.rm = T)
        
        profileDf$time <- c(1:nrow(profileDf))
        
        return(profileDf)
        
    }
}
