

#--- Central Tendency and spread --- 

x=rnorm(20,mean=6,sd=0.7)

mean(x)
median(x)

var(x)
sd(x)
IQR(x)


#get the probability of P(X =< -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2)

#get the probability of P(X > -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2,lower.tail = FALSE)

#get 5 random numbers from norm. dist. | mean=0 and sd=2
rnorm(5, mean=0 , sd=2)


#--- precision of estimates ----




#--- hypothesis testing ----



#---- regression ----






