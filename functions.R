#This script will contain all the functions to perform multiple alignment and circularization
#using dynamic time warping


#this function takes as arguments 2 dataframes that will be aligned, a vector that contains the column numbers
#for each dataframe that indicate which data should be used for the alignment, and a vector with the step
#patterns that want to be used for each type of alignment (from the left, from the right, local and global)
#if the patterns are NULL (default), asymmetrixP05 will be used for the first 3 alignments and symmetricP2 will
#be used for the global alignment.
#showDistPlot, if TRUE it will plot a representation of distance versus the sliding window
localUnivariateAlignment <- function(df1, df2, dataColumns, stepPattern = NULL, showDistPlot = FALSE){
    
    #argument check
    if(is.null(stepPattern)){
        stepPattern <- c(asymmetricP05, asymmetricP05, asymmetricP05, symmetricP2)
    }
    if(length(stepPattern) != 4){
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
        #dataframe with the fraction of the whole df data
        tempdf1 <- df1[i:nrow(df1), ]
        
        #alignment dtw function call
        alignment <- try(dtw(tempdf1[,dataColumns[1]], 
                             df2[,dataColumns[2]], #reference
                             step.pattern = stepPattern[1], 
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
                             step.pattern = stepPattern[2], 
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
        if(length(aligList) == (nrow(df1)*2)){
            break
        }
    }
    
    #open local alignment---------------------------------------------------------------
    alignment <- try(dtw(df1[,dataColumns[1]], 
                         df2[,dataColumns[2]], 
                         step.pattern = stepPattern[3], 
                         keep.internals  = TRUE,  
                         open.end = TRUE, 
                         open.begin = TRUE), silent = TRUE)
    
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
                         step.pattern = stepPattern[4], 
                         keep.internals  = TRUE,  
                         open.end = FALSE, 
                         open.begin = FALSE), silent = TRUE)
    
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
    
    #model evaluation
    if(showDistPlot){
        df <- data.frame(mod = 1:length(aligList), dist = sapply(aligList, '[[', 12))
        
        p <- ggplot(data = df, aes(x = mod, y = dist)) + 
            geom_point() + 
            loesTheme + 
            coord_cartesian()
        plot(p)
    }
    
    #defining the indexes of the number of models per type of model
    modelLeft <- c(1:nrow(df1))
    modelRight <- c((nrow(df1) + 1):(nrow(df1)*2))
    modelLocal <- length(aligList) - 1
    modelGlobal <- length(aligList) 
    
    #finding the best model
    smallestDistMod <- which(sapply(aligList, '[[', 12) == min(sapply(aligList, '[[', 12)))
    
    #creating the profile from the alignment of the two dfs in different ways depending
    #on which type of model is best
    if(smallestDistMod == modelLocal){
        print('LOCAL alignment')
        profileDf <- data.frame(time = c(1:max(nrow(df1), nrow(df2))), 
                                dfref = NA, 
                                dfquery = NA)
        model <- aligList[[smallestDistMod]]
        
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
        profileDf <- profileDf[-which(apply(profileDf[,2:3], 1, sum, na.rm = T) == 0),]
        profileDf$time <- 1:nrow(profileDf)
        
    }else if(smallestDistMod %in% modelLeft){
        print('LEFT coming alignment')
        #how many timepoints are left out
        hanglength <- nrow(df1) - which(modelLeft == smallestDistMod)
        #create our profile df
        profileDf <- data.frame(time = c(1:(max(nrow(df1), nrow(df2)) + hanglength)), 
                                dfref = NA, 
                                dfquery = NA)
        model <- aligList[[smallestDistMod]]
        profileDf[1:(nrow(df1) - length(unique(model$index1))), 'dfquery'] <- df1[1:(nrow(df1) - length(unique(model$index1))), dataColumns[1]]
        #extract the best model
        
        continuation <- which(is.na(profileDf$dfquery))[1]
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
        
        profileDf <- profileDf[-which(apply(profileDf[,2:3], 1, sum, na.rm = T) == 0),]
        
        profileDf$time <- c(1:nrow(profileDf))
        
        
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
            
            profileDf[ind + continuation, 2:3] <- c(mean(dfRefData), mean(dfQueryData))
        }
        
        hangRightDf <- data.frame(time = NA, 
                                  dfref = NA,
                                  dfquery = df1[c((which(modelRight == smallestDistMod)+1):nrow(df1)),dataColumns[1]])
        profileDf <- rbind(profileDf, hangRightDf)
        profileDf <- profileDf[-which(apply(profileDf[,2:3], 1, sum, na.rm = T) == 0),]
        
        profileDf$time <- c(1:nrow(profileDf))
        
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
    loess <- loess(formula = uni ~ time, data = profileDf, span = 0.3)
    profileDf$uni <- predict(loess)
    
    p1 <- ggplot() + 
        geom_rect(data = profileDf, aes(xmin = time, xmax = time, ymin = -1, ymax = Inf, 
                                        color = dfquery, 
                                        fill = dfquery), size = 1, show.legend = F) + 
        geom_line(data = profileDf, aes(x = time, y = dfquery), color = 'black', size = 2, na.rm = T) + 
        scale_color_gradient2(low = 'green', 
                              mid = 'yellow', 
                              high = 'red', 
                              na.value = 'white', 
                              limits = c(-1, 1)) + 
        coord_cartesian(ylim = c(-1,1))
    
    p2 <- ggplot() + 
        geom_rect(data = profileDf, aes(xmin = time, xmax = time, ymin = -1, ymax = Inf, 
                                        color = dfref, 
                                        fill = dfref), size = 1, show.legend = F) + 
        geom_line(data = profileDf, aes(x = time, y = dfref), color = 'black', size = 2, na.rm = T) + 
        scale_color_gradient2(low = 'green', 
                              mid = 'yellow', 
                              high = 'red', 
                              na.value = 'white', 
                              limits = c(-1, 1)) + 
        coord_cartesian(ylim = c(-1,1))
    
    pp <- ggplot() + 
        geom_rect(data = profileDf, aes(xmin = time, xmax = time, ymin = -1, ymax = Inf, 
                                        color = uni, 
                                        fill = uni), size = 1, show.legend = F) + 
        geom_line(data = profileDf, aes(x = time, y = uni), color = 'black', size = 2, na.rm = T) + 
        scale_color_gradient2(low = 'green', 
                              mid = 'yellow', 
                              high = 'red', 
                              na.value = 'white', 
                              limits = c(-1, 1)) + 
        coord_cartesian(ylim = c(-1,1))
    
    
    plotList <- list(pp, p2, p1)
    gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)
    plot(gridplot)
    
    return(profileDf)
}
