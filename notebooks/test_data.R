library(dplyr)
library(ggplot2)


# path <- 'data/raw/Damvarmelager temperatur målinger 2023-09.csv'
path <- 'data/processed/data.csv'
data <- read.csv(path)
head(data)

A_COLS <- c('A1TT03', 'A2TT03', 'A1TT04', 'A2TT04', 'A1TT05',
            'A2TT05', 'A1TT06', 'A2TT06', 'A1TT07', 'A2TT07',
            'A1TT08', 'A2TT08', 'A1TT09', 'A2TT09', 'A1TT10',
            'A2TT10', 'A1TT11', 'A2TT11', 'A1TT12', 'A2TT12',
            'A1TT13', 'A2TT13', 'A1TT14', 'A2TT14')


boxplot(data[A_COLS])








