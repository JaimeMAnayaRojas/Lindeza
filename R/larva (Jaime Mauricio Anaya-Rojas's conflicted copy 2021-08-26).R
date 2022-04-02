
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
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
setwd("~/Dropbox/Projects_JM/Muenster/Lindeza/")
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


op<-par(mfrow=c(5,3), mar = c(4, 5, 2, 1), oma = c(1, 1, 1, 1), cex = 0.75)

for(i in 1:15){
  
  pacf(data$larvae[which(data$Treat == levels(data$Treat)[i]  )], 9,  ylim=range(-1,1))
  
}


for(i in 1:15){
  
  pacf(data$adults[which(data$Treat == levels(data$Treat)[i]  )], 9,  ylim=range(-1,1))
  
}


sum_data <- data %>% 
  group_by(Regime, Generation) %>%
  summarise(larvae_m = mean(larvae), larvae_sd = sd(larvae), adults_m = mean(adults), adults_sd = sd(adults))


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

f1 <- bf(larvae ~  G1*Regime + G2  +  A1 + (1|Treat), family = poisson())
f2 <- bf()
mod4 <- brm(f1 + f2,  data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
summary(mod4)


mod5 <- brm(f2,  data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
summary(mod5)


mod6 <- brm(A2 ~  adults * Regime  + (1|Treat), family = poisson(),  data = df, 
            control = list(max_treedepth = 10, adapt_delta = 0.92), iter = 4000, cores = 4, chains = 4)
summary(mod6)


plot(mod4)
p_means = apply(posterior_predict(mod4), 2, mean)
p_ci = t(apply(posterior_predict(mod4), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)



# 

mod4 <- brm(larvae ~  G1*Regime + G2  +  A1 + (1|Treat), family = poisson() ,data = df, 
            control = list(max_treedepth = 15, adapt_delta = 0.92), iter = 4000, cores = 6, chains = 4)
summary(mod4)

plot(mod4)
p_means = apply(posterior_predict(mod4), 2, mean)
p_ci = t(apply(posterior_predict(mod4), 2, HPDI, prob =.95))
plot(df$larvae ~ p_means, xlim = c(30, 200))
abline(0,1)
segments(x0 = p_ci[,1], x1 = p_ci[,2], y0 = df$larvae, y1 = df$larvae)



post = data.frame(posterior_samples(mod4))

precis(post, prob = .95, digits = 3, depth = 2 )


##

p_link <- function(post, Generations, Matched, Unmatched, Adults){
  
  G = Generations - 4
  G2 =  Generations^2 - 4^2
  U = Unmatched
  M = Matched
  out <- with(post, b_Intercept + b_G1*G + b_RegimeMatched* M + b_RegimeUnmatched * U + b_G2*G2 + b_A1 * Adults )
  return(exp(out))
}


## I also have to model the number of adults produce
My_AR_pred <- function(post= post, Generations, Matched = 0, Unmatched = 0, Initial_adults=50 ){
  
  df = array(0,c(dim(post)[1],length(G)))  
  
  for( i in 1:length(Generations)){
    
    p_link(post= post, Generations = Gen, Matched = Matched, Unmatched = Unmatched, Adults = Initial_adults)  
    
    
  }
  
  
  
  
}



### Plot overall results for now

Generations = 3:9

p.Control = sapply(1:length(Gen), function(i) p_link(post= post, Generations = 4, Matched = 0, Unmatched = 0, Larvae = 100, Adults = 100))
p_means = apply(p.Control, 2, mean)
p_ci = (apply(p.Control, 2, HPDI, prob =.95))

length(p_means)

plot(p_means ~ Gen, pch = "", ylim = c(0, 200))
lines(p_means ~ Gen, col = 'orange', lwd = 2)
shade(p_ci, Gen, col = col.alpha('orange', 0.2))


p.Match = sapply(1:length(Gen), function(i) p_link(post= post, Generations = 4, Matched = 1, Unmatched = 0, Larvae = 100, Adults = 100))
p_means = apply(p.Match, 2, mean)
p_ci = (apply(p.Match, 2, HPDI, prob =.95))

lines(p_means ~ Gen, col = 'black', lwd = 2)
shade(p_ci, Gen, col = col.alpha('black', 0.2))


p.UnMatch = sapply(1:length(Gen), function(i) p_link(post= post, Generations = 4, Matched = 0, Unmatched = 1, Larvae = 100, Adults = 100))
p_means = apply(p.UnMatch, 2, mean)
p_ci = (apply(p.UnMatch, 2, HPDI, prob =.95))

lines(p_means ~ Gen, col = 'blue', lwd = 2)
shade(p_ci, Gen, col = col.alpha('blue', 0.2))
