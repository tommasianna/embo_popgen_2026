# exercise day 4

# exercise 2a
# Consider a population of size 2N=100 and a current allele frequency of f=10/100=0.1.
#a) What is the probability that the allele frequency remains unchanged after one generation of random mating?

dbinom(x=10, size=100, prob=0.1)

# x=10 is the number of successes (the allele frequency remains unchanged), 
# size = 100 is the total number of trials (the population size), 
# and prob=0.1 is the probability of success on each trial (the current allele frequency).


# exercise 2b
# What is the probability that the allele is lost after one generation of random mating?

dbinom(x=0, size=100, prob=0.1)


# exercise 2c
# c) What is the probability that the allele is fixed after one generation of random mating?

dbinom(x=100, size=100, prob=0.1)


# exercise 2d
# What is the probability that the allele frequency decreases after one generation of random mating?

dbinom(x=9, size=100, prob=0.1)

# x=9 because we want to calculate the probability that the allele frequency decreases, 
# which means that there are 9 copies of the allele in the next generation (instead of 10).


# exercise 2e
# What is the probability that the allele frequency increases after one generation of random mating?

1 - pbinom(10, size = 100, prob = 0.1)

# We can calculate the probability of increase analogously as the integral from 11 to 100, 
# or as 1 minus the integral from 0 to 100 (?)
# Note that the probabilities for an increase and decrease are not equal, i.e. this distribution is not symmetric.



# QUESTION 3
# Plot the probability distribution on the allele frequency for a current allele 
# frequency f=0.1 and 2N=10, 2N=100, 2N=1000 and 2N=10000. 
# What is the effect of population size?

par(mfrow=c(2,2)) # this allow me to plot 4 graphs in one window

for(twoN in c(10, 100, 1000, 10000)) {
  p <- 0:twoN # this creates a sequence of numbers from 0 to twoN, which represent the possible allele frequencies in the next generation
  d <- dbinom(x=p, size=twoN, prob=0.1) # this calculates the probability of each allele frequency in the next generation using the binomial distribution
  plot(p, d, type = 'b')
}



# QUESTION 4
# Simulating neutral allele frequency trajectories
# Write a function simulateWF() to simulate allele trajectories under the 
# Wright-Fisher model. Your function should take as input 
# i) the population size 2N, 
# ii) the initial allele frequency f and 
# iii) the number G of generations to simulate. 
# It should then return the allele frequency (between 0 and 1) for each generation as a vector.


simulateWF <- function(twoN, f, G){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN
  }
  return(p)
}



# QUESTION 5
#  Use your function simulateWF() to simulate 1000 trajectories with 
# 2N=100 and f=0.1 for G=1000 and plot them in one plot. 
# In how many cases was the allele lost? 
# Does this match the expectation? Repeat for different population sizes 
# and initial allele frequencies. What is the effect of the population size?

trajectories <- replicate(1000, simulateWF(twoN = 100, f = 0.1, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was lost in", sum(trajectories[1000,] == 0), "/", ncol(trajectories), "cases."))
## [1] "Allele was lost in 897 / 1000 cases."

# The allele is expected to go to fixation with the probability of its initial allele 
# frequency f, and hence to be lost with probability 1−f.
# If you use larger population sizes, you should observe 
# less fluctuation and eventually many sites will remain polymorphic 
# even after G=1000 generations. However, if run for enough generations, 
# the same fraction of alleles will be lost.



# QUESTION 6
# Use your function simulateWF() to study fixation probability 
# of a new mutation under different population sizes (ensure G is large enough). 
# Does the fixation probability depend on the population size? 
# How does this affect the substitution rate?

simulateWF <- function(twoN, f, G){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN
  }
  return(p)
}

trajectories <- replicate(1000, simulateWF(twoN = 100, f = 1/twoN, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))



# QUESTION 7
# Write a function simulateWFWithSelection() that simulates allele trajectories 
# under both genetic drift and viability selection. Similar to your function simulateWF() 
# you wrote above, it should take as input
# i) the population size 2N, 
# ii) the initial allele frequency f, 
# iii) the number G of generations to simulate and 
# iv) also a vector v of viabilities for the genotypes AA, Aa and aa. 
# It should then return the allele frequency (between 0 and 1) for each generation as a vector. 
# In each generation, your function should apply selection to alter the allele frequency, 
# and then use binomial sampling to simulate genetic drift in that modified allele frequency.

simulateWFWithSelection <- function(twoN, f, G, v){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    #selection
    fA <- p[i]
    fa <- 1-fA
    fPrime <- (v[1]*fA*fA + v[2]*fA*fa)/(v[1]*fA*fA + v[2]*2*fA*fa + v[3]*fa*fa)
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN
  }
  return(p)
}
