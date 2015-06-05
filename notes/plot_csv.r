data <- read.csv("data.csv")
channel0 = data[data$channel == 0,];

data2 <- read.csv("data2.csv")
channel02 = data2[data2$channel == 0,]

plot(channel0$iteration, channel0$objective)
lines(channel02$iteration, channel02$objective)

min(channel0$objective)
min(channel02$objective)
