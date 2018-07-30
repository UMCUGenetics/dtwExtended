#In this R script we will generate some demo data using the Sinus function

#We will generate a sinus sequence of 2 frequencies length and subset portions of the sequence 
#of different length that have some overlap between them

#representation of the data
y <- sin(seq(from = 0, to = 4*pi, length.out = 720))
t <- c(1:720)

plot(t, y)

#subset of the data
#we will make 10 subsets of the data that span the whole length and have different lengths in 
#the overlap between them
fractionsList <- list(s1 = c(1:120),
                      s2 = c(80:140),
                      s3 = c(110:210),
                      s4 = c(160:250),
                      s5 = c(230:320),
                      s6 = c(290:400),
                      s7 = c(360:430),
                      s8 = c(400:640),
                      s9 = c(510:710),
                      s10 = c(600:720))


#create dataframes of the portions and put them in a list
dataList <- lapply(X = fractionsList, function(x){
    
    return(data.frame(t = c(1:length(x)), y = y[x]))
    
})

#save the object as an RData file that can be loaded in the test environment
save(dataList, file = 'Demo files/demoData.RData')



