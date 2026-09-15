
men <- c(102,87,101,96,107,101,91,85,108,67,85,82) 
women <- c(73,81,111,109,143,95,92,120,93,89,119,79,90,126,62,92,77,106, 105,111)

mean(men)
sd(men)
mean(women)
sd(women)

sd(women)/sd(men)
t.test(x=men, y=women, alternative="two.sided", conf.level=0.95, var.equal=TRUE)