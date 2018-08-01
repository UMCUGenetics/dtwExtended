#In this R script we will generate some demo data using the Sinus function and then we will
#the functions in the function.R script to do multiple dtw alignment, local alignment and 
#circular alignment

#We will generate a sinus sequence of 2 frequencies length and subset portions of the sequence 
#of different length that have some overlap between them

#representation of the data
y <- round(sin(seq(from = 0, to = 2.5*pi, length.out = 420)), digits = 5)
t <- c(1:420)

plot(t, y)

#subset of the data
#we will make 10 subsets of the data that span the whole length and have different lengths in 
#the overlap between them
fractionsList <- list(s1 = c(1:120),
                      s2 = c(80:140),
                      s3 = c(110:210),
                      s4 = c(160:250),
                      s5 = c(230:320),
                      s6 = c(290:420),
                      s7 = c(270:390),
                      s8 = c(310:400),
                      s9 = c(300:420),
                      s10 = c(300:420))


#create dataframes of the portions and put them in a list
dataList <- lapply(X = fractionsList, function(x){
    
    return(data.frame(t = c(1:length(x)), y = y[x]))
    
})

#test some univariate alignments

#LEFT ---------
profile <- pairwiseUnivariateAlignment(df1 = dataList[[1]], 
                                       df2 = dataList[[2]], 
                                       dataColumns = c(2, 2), 
                                       showDistPlot = F)

#visual representation of the alignment
require(ggplot2)
require(gridExtra)

pp <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[4])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Unified sequences')
p1 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[3])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Sequence1')
p2 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[2])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1))+ 
    ggtitle('Sequence2')

plotList <- list(pp, p2, p1)
gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)


#RIGHT ---------
profile <- localUnivariateAlignment(df1 = dataList[[3]], 
                                    df2 = dataList[[2]], 
                                    dataColumns = c(2, 2), 
                                    showDistPlot = T)

#visual representation of the alignment
require(ggplot2)
require(gridExtra)

pp <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[4])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Unified sequences')
p1 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[3])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Sequence1')
p2 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[2])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1))+ 
    ggtitle('Sequence2')

plotList <- list(pp, p2, p1)
gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)


#LOCAL ---------
profile <- localUnivariateAlignment(df1 = dataList[[7]], 
                                    df2 = dataList[[8]], 
                                    dataColumns = c(2, 2), 
                                    showDistPlot = T)

#visual representation of the alignment
require(ggplot2)
require(gridExtra)

pp <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[4])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Unified sequences')
p1 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[3])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Sequence1')
p2 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[2])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1))+ 
    ggtitle('Sequence2')

plotList <- list(pp, p2, p1)
gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)


#GLOBAL ---------
profile <- localUnivariateAlignment(df1 = dataList[[9]], 
                                    df2 = dataList[[10]], 
                                    dataColumns = c(2, 2), 
                                    showDistPlot = T)

#visual representation of the alignment
require(ggplot2)
require(gridExtra)

pp <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[4])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Unified sequences')
p1 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[3])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) + 
    ggtitle('Sequence1')
p2 <- ggplot(data = profile, aes_string(x = names(profile)[1], y = names(profile)[2])) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1))+ 
    ggtitle('Sequence2')

plotList <- list(pp, p2, p1)
gridplot <- grid.arrange(grobs = plotList, ncol = 1, as.table = FALSE)
plot(gridplot)


#next we will make a profile from doing a progressive pair-wise alignment of all the portions

profileDf <- multiAlignmentUnivariateProfile(dfList = dataList, dataColumns = c(2,2), showDendrogram = F, showAlignmentPlot = T)


ggplot(data = profileDf, aes(x = time, y = uni)) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) 


#circularization of the profile ------------------------
#since we know that our data is a loop over time. We will use the next function to reduce it to a single frequency 

circularDf <- circularizeSequenceUnivariate(sequenceDf = profileDf, dataColumn = 4)

ggplot(data = circularDf, aes(x = time, y = uni)) + 
    geom_point() + 
    coord_cartesian(ylim = c(-1,1)) 
