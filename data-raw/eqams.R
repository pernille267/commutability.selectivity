set.seed(2)
library(data.table)
tau <- runif(5, min = 5, max = 20)
tau1 <- rep(tau*1.05, each = 3)
tau2 <- tau1*1.01 - 0.02
tau3 <- tau2*0.99 + 0.02
tau4 <- tau2*0.97 + 0.01
tau5 <- tau4*1.12 - 0.03
sd1 <- 0.03 * mean(tau)
sd2 <- 0.04 * mean(tau)
sd3 <- 0.03 * mean(tau)
sd4 <- 0.05 * mean(tau)
sd5 <- 0.06 * mean(tau)

taus <- data.table(ID1 = rep(1:5, each = 3), ID2 = rep(1:3, times = 5),
                   tau1, tau2, tau3, tau4, tau5)
taus <- taus[,.(ID2=ID2, MP1 = round(tau1 + rnorm(3, mean = 0, sd = sd1), 2),
                MP2 = round(tau2 + rnorm(3, mean = 0, sd = sd2), 2),
                MP3 = round(tau3 + rnorm(3, mean = 0, sd = sd3), 2),
                MP4 = round(tau4 + rnorm(3, mean = 0, sd = sd4), 2),
                MP5 = round(tau4 + rnorm(3, mean = 0, sd = sd5), 2)),by=ID1]
setnames(x = taus, old = c("ID1", "ID2"), new = c("SampleID", "ReplicateID"), skip_absent = TRUE)
sampled_eqam_measurements <- taus

usethis::use_data(sampled_eqam_measurements, overwrite = TRUE, compress = "xz")

usethis::use_r("documentation_sampled_eqam_measurements")
