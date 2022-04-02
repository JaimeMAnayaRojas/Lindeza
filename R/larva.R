####  
# library(brms)
library(rethinking)
library(data.table)
library(dplyr)
library(readxl)
library(tseries)

rm(list=ls(all=TRUE))
###########################################################################################################
# First get the data
getwd()
setwd("~/Dropbox/Jaime M/Projects_JM/Muenster/Lindeza/")
list.files()
data <- read_excel("data/Larvae density_Larvae and adults.xlsx")
names(data)[c(4,5)] <- c("larvae","adults")

head(data)
data$Replicate <- factor(data$Replicate)
ggplot(data, aes(y=larvae, x=Generation, colour = Replicate)) + 
  geom_point()

data$Regime <- factor(data$Regime)
levels(data$Regime)

data$Treat <- factor(paste(data$Regime, data$Replicate, sep = '-'))

length(levels(data$Treat))


op<-par(mfrow=c(1,1), mar = c(4, 5, 2, 1), oma = c(1, 1, 1, 1), cex = 0.75)

for(i in 1:15){
  
  pacf(data$larvae[which(data$Treat == levels(data$Treat)[i]  )], 9,  ylim=range(-1,1))
  
}



for(i in 1:15){
  
  pacf(data$adults[which(data$Treat == levels(data$Treat)[i]  )], 9,  ylim=range(-1,1))
  
}


sum_data <- data %>% 
  group_by(Regime, Generation) %>%
  summarise(larvae_m = mean(larvae), larvae_sd = sd(larvae), adults_m = mean(adults), adults_sd = sd(adults))

sum_data
graphics.off()

#
names(data)
data = data.frame(data)
data$Match <- ifelse(data$Regime == "Matched", 1,0)
data$Unmatch <- ifelse(data$Regime == "Unmatched", 1,0)
data$Replicate <- as.numeric(as.character(data$Replicate))
data$A2 = 0
data$A1 = 0


i = 31

for(i in 1:dim(data)[1]){
  print(i)
  if(data$Generation[i] > 3 & data$Generation[i] < 9){
    id1 = which(data$Regime == data$Regime[i] & data$Replicate == data$Replicate[i] & data$Generation == data$Generation[i]-1 )
    data$A1[i] = data$larvae[id1]
    id2 = which(data$Regime == data$Regime[i] & data$Replicate == data$Replicate[i] & (data$Generation == data$Generation[i]+1)  )
    data$A2[i] = data$adults[id2]
  }else if(data$Generation[i] == 9){
    
    id1 = which(data$Regime == data$Regime[i] & data$Replicate == data$Replicate[i] & data$Generation == data$Generation[i]-1 )
    data$A1[i] = data$larvae[id1]
    
    data$A2[i] = NA
    
  }
  
}

df = subset(data, A1 > 0)

head(df)

library(brms)


plot(larvae ~ Generation, col = Treat, df)

# Exploratory model, I am using the treatment as random effects

df$G1 = df$Generation - 4
df$G2 = df$Generation^2 - 4^2 

# Exploratory model, I am using the treatment as random effects

levels(df$Treat)
# 
# f1 <- bf(larvae ~  G1*Regime + G2  +  A1 + (1|Treat), family = poisson())
# f2 <- bf()
# mod4 <- brm(f1 + f2,  data = df, 
#             control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
# summary(mod4)
# 
# 
# mod5 <- brm(f2,  data = df, 
#             control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
# summary(mod5)
# 

mod6 <- brm(A2 ~  adults * Regime  + (1|Treat), family = poisson(),  data = df, 
            control = list(max_treedepth = 10, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
summary(mod6)

plot((A2) ~ (larvae), col = Regime, df)

names(df)
cor(df[, c(4,9:12)])
plot(mod6)
p_means = apply(posterior_predict(mod6), 2, mean)
p_ci = t(apply(posterior_predict(mod6), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)



#


data$Gen = factor(paste("G-",data$Generation, sep=""))

levels(data$Gen) = c("G3", "G4", "G5", "G6", "G7", "G8", "G9")

mod0 <- brm(larvae ~  Gen*Regime  + (1|Treat), family = poisson() ,data = data, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)

summary(mod0)

conditional_effects(mod0, effects = "Gen:Regime")



mod1 <- brm(larvae ~  G1*Regime + G2  +  A1 + (1|Treat), family = poisson() ,data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)

plot(mod1)
p_means = apply(posterior_predict(mod1), 2, mean)
p_ci = t(apply(posterior_predict(mod1), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)



mod2 <- brm(larvae ~  G1*Regime + G2  +  (1|Treat), family = poisson() ,data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)

plot(mod2)
p_means = apply(posterior_predict(mod2), 2, mean)
p_ci = t(apply(posterior_predict(mod2), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)


df$G3 <- df$Generation^3 - 4^3
mod3 <- brm(larvae ~  G1*Regime + G2 + G3  +  (1|Treat), family = poisson() ,data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)

summary(mod3)
plot(mod3)
p_means = apply(posterior_predict(mod3), 2, mean)
p_ci = t(apply(posterior_predict(mod3), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)



plot(p_means ~ df$Generation)

conditional_effects(mod3, effects = "G1:Regime")


post = data.frame(posterior_samples(mod3))

sum = as.data.frame(precis(post, prob = .95, digits = 3, depth = 2 ))

write.csv(sum, "Parameter estimations.csv")

##

p_link <- function(post, Generations, Matched, Unmatched){
  
  G1 = Generations - 4
  G2 =  Generations^2 - 4^2
  G3 =  Generations^3 - 4^3
  U = Unmatched
  M = Matched
  out <- with(post, b_Intercept + b_G1*G1 + b_RegimeMatched* M + b_RegimeUnmatched * U + b_G2*G2 + b_G3*G3 +
                b_G1.RegimeMatched * G1 * M + b_G1.RegimeUnmatched * G1 * U)
  return(exp(out))
}


### Plot overall results for now
jpeg(file = "Figure 1.jpeg",width = 6.0, height = 5.0, units='in', res=800)


Gen = 4:9

p.Control = sapply(1:length(Gen), function(i) p_link(post= post, Generations = Gen[i], Matched = 0, Unmatched = 0))
p_means = apply(p.Control, 2, mean)
p_ci = (apply(p.Control, 2, HPDI, prob =.95))

length(p_means)

plot(p_means ~ Gen, pch = "", ylim = c(50, 170), ylab = "Number of Larvae", xlab = "Generations")
lines(p_means ~ Gen, col = 'orange', lwd = 2)
shade(p_ci, Gen, col = col.alpha('orange', 0.2))


p.Match = sapply(1:length(Gen), function(i) p_link(post= post, Generations = Gen[i], Matched = 1, Unmatched = 0))
p_means = apply(p.Match, 2, mean)
p_ci = (apply(p.Match, 2, HPDI, prob =.95))

lines(p_means ~ Gen, col = 'black', lwd = 2)
shade(p_ci, Gen, col = col.alpha('black', 0.2))


p.UnMatch = sapply(1:length(Gen), function(i) p_link(post= post, Generations = Gen[i], Matched = 0, Unmatched = 1))
p_means = apply(p.UnMatch, 2, mean)
p_ci = (apply(p.UnMatch, 2, HPDI, prob =.95))

lines(p_means ~ Gen, col = 'blue', lwd = 2)
shade(p_ci, Gen, col = col.alpha('blue', 0.2))



legend('topleft',c('Control','Matched', "Unmatched"), lty = 1,  col=c('orange', 'black', 'blue'),bty='n',cex=1.25)
graphics.off()
