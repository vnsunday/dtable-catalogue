

# 
# Interested:
#       Mean Net-weight
# 
# H0: Mean = 80
# HA: Mean != 80

snacks <- c(87.7,80.01,77.28,78.76,81.52,74.2,80.71,79.5,77.87,81.94,80.7,
82.32,75.78,80.19,83.91,79.4,77.52,77.62,81.4,74.89,82.95,
73.59,77.92,77.18,79.83,81.23,79.28,78.44,79.01,80.47,76.23,
78.89,77.14,69.94,78.54,79.7,82.45,77.29,75.52,77.21,75.99,
81.94,80.41,77.7)

n <- length(snacks)
snack.mean <- mean(snacks)
snack.sd <- sd(snacks)
snack.se <- snack.sd/sqrt(n)
snack.T <- (snack.mean - 80)/snack.se 
pt(snack.T, df=n-1)
snack.mean+c(-1,1)*qt(0.975,n-1)*snack.se


