library(splines)
make_data<-function(start = '2015-08-01', end = '2015-08-31', fill = c()){
    water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
    inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
    vars<-read.csv('data/processed/dronninglund/variables.csv')
    ambientTemp<-vars$temp_dry

    inputs <- cbind(inputs,ambientTemp)
    idx <- complete.cases(inputs)
    
    data <- cbind(water_sensors[idx,],inputs[idx,])
    # data <- cbind(water_sensors, inputs)
    data$dFtopIn <- c(0,diff(data$FtopIn))
    data$dTtop <- c(0,diff(data$Ttop))
    data$dTbot <- c(0,diff(data$Tbot))
    data$dFbotIn <- c(0,diff(data$FbotIn))
    data$dFbotOut <- c(0,diff(data$FbotOut))

    data$FbotInk1 <- c(0,data$FbotIn[1:(nrow(data)-1)])
    data$Tbotk1 <- c(0,data$Tbot[1:(nrow(data)-1)])

    data$FbotOutk1 <- c(0,data$FbotOut[1:(nrow(data)-1)])
    data$FtopInk1 <- c(0,data$FtopIn[1:(nrow(data)-1)])
    data$FtopInk2 <- c(0,data$FtopIn[1:(nrow(data)-1)])


    data <- data[, !duplicated(colnames(data))]
    # dropna
    # data <- data[complete.cases(data),]
    data <- tidyr::fill(data, fill, .direction = "downup")

    data <- data[which(data$X >= start & data$X <= end),]
    


    return(data.frame(data))
}

