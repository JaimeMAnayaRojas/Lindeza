####  
# library(brms)
library(rethinking)
library(data.table)
library(dplyr)
library(readxl)
library(tseries)
library(ggpubr)

rm(list=ls(all=TRUE))

LOS <- function(a){
  
  b = 100 * length(which(a>0))/length(a)
  return(b)
}
###########################################################################################################
# First get the data
getwd()
#setwd("~/Dropbox/Jaime M/Projects_JM/Muenster/Lindeza/")
list.files()
data <- read_excel("data/Larvae density_Larvae and adults.xlsx")
names(data)[c(4,5)]
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
# 
# # Exploratory model, I am using the treatment as random effects
# 
# df$G1 = df$Generation - 4
# df$G2 = df$Generation^2 - 4^2 
# 
# # Exploratory model, I am using the treatment as random effects
# 
# levels(df$Treat)
# # 
# # f1 <- bf(larvae ~  G1*Regime + G2  +  A1 + (1|Treat), family = poisson())
# # f2 <- bf()
# # mod4 <- brm(f1 + f2,  data = df, 
# #             control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
# # summary(mod4)
# # 
# # 
# # mod5 <- brm(f2,  data = df, 
# #             control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
# # summary(mod5)
# # 
# 
# # mod6 <- brm(A2 ~  adults * Regime  + (1|Treat), family = poisson(),  data = df, 
# #             control = list(max_treedepth = 10, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
# # summary(mod6)
# # 
# # plot((A2) ~ (larvae), col = Regime, df)
# # 
# # names(df)
# # cor(df[, c(4,9:12)])
# # plot(mod6)
# # p_means = apply(posterior_predict(mod6), 2, mean)
# # p_ci = t(apply(posterior_predict(mod6), 2, HPDI, prob =.95))
# # plot(df$larvae ~ p_means, xlim = c(30, 200))
# # abline(0,1)
# # segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)
# # 
# # 
# 
# #


data$Gen = factor(paste("G-",data$Generation, sep=""))
levels(data$Gen) = c("G3", "G4", "G5", "G6", "G7", "G8", "G9")

levels(data$Treat)

mod0 <- brm(larvae ~ 0 +  Gen:Regime  + (1|Treat), family = poisson() ,data = data, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)

summary(mod0)

conditional_effects(mod0, effects = "Gen:Regime")


###


post0 = posterior_samples(mod0)
str(post0)






post = posterior_samples(mod0)
post = exp(post[,1:21])
str(post)

sumPred = as.data.frame(precis(post, prob = .95))
sumPred$L68 = precis(post, prob = .68)[,3]
sumPred$U68 = precis(post, prob = .68)[,4]
sumPred$Gen = factor(rownames(sumPred))
sumPred$Gen <- rep(c("3","4","5","6","7","8", "9"), 3)
sumPred$Regime <-  c(rep("Control", 7), rep("Matched", 7), rep("Unmatched", 7))

# Get the predictions for each 


df_pred = sumPred


##

names(df_pred)
names(df_pred)[c(1,3,4)] <- c("larvae","L95", "U95")


pd = position_dodge(0.75)
(mPlot = ggplot(df_pred, aes(x=Gen, y=larvae, colour=Regime, group =Regime )) + 
    geom_errorbar(aes(ymin=L95, ymax=U95), width=.1, position=pd, alpha =1, size = 0.4) +
    geom_errorbar(aes(ymin=L68, ymax=U68), width=.1, position=pd, alpha =1, size = 1.2) +
    geom_line(position=pd, size = 0.5) +
    geom_point(position=pd, size = 2.5) 
  
)

pointdata = data
pointdata$larvae
levels(pointdata$Gen) <- c("3","4","5","6","7","8", "9")

levels(pointdata$Regime)

plotA <- mPlot + geom_point(data = pointdata, 
                             mapping = aes(x=Gen, y=larvae, colour=Regime),
                             position = position_jitterdodge(dodge.width=0.75),
                             size = 0.75, alpha = 0.3) +
    theme_bw() + theme(panel.grid.major = element_blank(), legend.position = c(0.85, .8)) + 
    scale_color_manual(values=c("#74a3dbff", "#aa2818ff", "#cb7f4eff")) + 
    ylab("Larvae (N)") +
    xlab("Generation") +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -2, ymax = 255), fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -2, ymax = 255),  fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -2, ymax = 255), fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -2, ymax = 255),  fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) 

  
plotA

##
df_contrast = data.frame(G3M = post$`b_GenG3:RegimeMatched` / post$`b_GenG3:RegimeControl`,
                         G4M = post$`b_GenG4:RegimeMatched` / post$`b_GenG4:RegimeControl`,
                         G5M = post$`b_GenG5:RegimeMatched` / post$`b_GenG5:RegimeControl`,
                         G6M = post$`b_GenG6:RegimeMatched` / post$`b_GenG6:RegimeControl`,
                         G7M = post$`b_GenG7:RegimeMatched`/ post$`b_GenG7:RegimeControl`,
                         G8M = post$`b_GenG8:RegimeMatched` / post$`b_GenG8:RegimeControl`,
                         G9M = post$`b_GenG9:RegimeMatched` / post$`b_GenG9:RegimeControl`,
                         
                         G3U = post$`b_GenG3:RegimeUnmatched` / post$`b_GenG3:RegimeControl`,
                         G4U = post$`b_GenG4:RegimeUnmatched` / post$`b_GenG4:RegimeControl`,
                         G5U = post$`b_GenG5:RegimeUnmatched` / post$`b_GenG5:RegimeControl`,
                         G6U = post$`b_GenG6:RegimeUnmatched` / post$`b_GenG6:RegimeControl`,
                         G7U = post$`b_GenG7:RegimeUnmatched` / post$`b_GenG7:RegimeControl`,
                         G8U = post$`b_GenG8:RegimeUnmatched` / post$`b_GenG8:RegimeControl`,
                         G9U = post$`b_GenG9:RegimeUnmatched` / post$`b_GenG9:RegimeControl`
                         
)

df_contrast = log(df_contrast)


contrast_tab = cbind(precis(df_contrast, prob = 0.95), precis(df_contrast, prob = 0.68)[,3:4])
names(contrast_tab)[c(3,4,6,7)] <- c("L95", "U95", "L68", "U68" )
contrast_tab$Gen = rownames(contrast_tab)
contrast_tab$Gen = rep(c("3","4","5","6","7","8","9"),2)
contrast_tab$Contrast = c(rep("Matched-Control",7), rep("Unmatched-Control",7))
contrast_tab$LOS = apply(df_contrast,2, LOS)



##
df_los = data.frame(     G3= post$`b_GenG3:RegimeMatched` / post$`b_GenG3:RegimeUnmatched`,
                         G4 = post$`b_GenG4:RegimeMatched` / post$`b_GenG4:RegimeUnmatched`,
                         G5 = post$`b_GenG5:RegimeMatched` / post$`b_GenG5:RegimeUnmatched`,
                         G6 = post$`b_GenG6:RegimeMatched` / post$`b_GenG6:RegimeUnmatched`,
                         G7 = post$`b_GenG7:RegimeMatched`/ post$`b_GenG7:RegimeUnmatched`,
                         G8 = post$`b_GenG8:RegimeMatched` / post$`b_GenG8:RegimeUnmatched`,
                         G9 = post$`b_GenG9:RegimeMatched` / post$`b_GenG9:RegimeUnmatched`
)


df_los = data.frame(Gen = 3:9, LOS = apply(log(df_los), 2, LOS))

df_los$Gen = factor(df_los$Gen)
df_los$y = -2
##

(plotB <- ggplot(data=contrast_tab, aes(x=Gen, y=mean, ymin=L95, ymax=U95, colour = Contrast)) +
  geom_pointrange(size = 0.5, position = pd) + 
  geom_errorbar(aes(ymin=L68, ymax=U68), width=.1, position=pd, alpha =1, size = 1.2) +
  geom_hline(yintercept=0, lty=1, colour = "#74a3dbff") +  # add a dotted line at x=1 after flip
  xlab("Generation") + ylab("Effect size (LRR)") +
  theme_bw() + theme(panel.grid.major = element_blank(), legend.position = "none" ) + 
  scale_color_manual(values=c("#aa2818ff", "#cb7f4eff")) 
)

(plotB <- plotB + geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -1.5, ymax = 1.2), fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -1.5, ymax = 1.2), fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -1.5, ymax = 1.2),  fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) +
  geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -1.5, ymax = 1.2), fill= "cyan",  inherit.aes = FALSE, alpha = 0.005) 
)  

plotB

figure <- ggarrange(plotA, plotB,
                    labels = c("A) ", "B) "), font.label = list(size= 11),
                    label.x = -0.01, label.y = 1.05, common.legend = T,
                    ncol = 1, nrow = 2)
figure

ggsave(filename = "plots/Figure-Larvae.svg", plot = figure, device = "svg", width = 105, height = 174, units = "mm" )
