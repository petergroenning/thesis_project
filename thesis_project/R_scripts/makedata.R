library(splines)
make_data<-function(start = '2015-08-01', end = '2015-08-31', splines = 5){
    water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
    inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
    vars<-read.csv('data/processed/dronninglund/variables.csv')
    ambientTemp<-vars$temp_dry


    data <- cbind(water_sensors, inputs,ambientTemp)
    data <- data[, !duplicated(colnames(data))]
    # dropna
    data <- data[complete.cases(data),]
    data <- data[which(data$X >= start & data$X <= end),]


    if (splines)
        bs<-bs(1:16,df=splines,intercept = TRUE)
        for (i in 1:splines){
            for (j in 0:15){
                data[,paste0('b',i,'x',j)]<-bs[(j+1),i]
            }
        }

    return(data.frame(data))
}

