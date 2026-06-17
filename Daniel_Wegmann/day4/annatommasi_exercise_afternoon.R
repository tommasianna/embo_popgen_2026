# Anna Tommasi day 4 afternoon exercises

# QUESTION 1
# Given a constant population of size 2N=100,
# what is the probability that two individuals pick the same parent one generation before? 
# What that they pick the same parent two generations before?

twoN = 100

p1 <- 1/twoN

p1

# 2 generation ago

p2 <- (1-p1)*p1



# QUESTION 2
# A researcher sequences a DNA fragment of length 10kb of a diploid individual 
# and detects 21 heterozygous sites. Assume that μ=10−9 per site for this segment 
# and that the assumptions of the standard coalescent and the infinites sites model hold.


# 2A
# What is the effective size of this population?

H <- 21
mi <- 10^-9

Ne <- H/(4*mi)

Ne


# 2B
# What is the expected time to the most recent common ancestor of the two sequences 
# in this individual?

ETMRCA <- 2*Ne # it's 2Ne because we are considering two sequences

ETMRCA



# QUESTION 3
# Consider a sample of two diploid individuals from a population for which 
# the assumptions of the standard coalescence model holds true. 
# What is the probability that the two sequences of the first individual 
# are more closely related to each other than they are to any of the sequences 
# of the other individual?

NO IDEA



# QUESTION 4A
# The following DNA sequences were obtained from a population:

nseq <- 3 # number of sequences
seqlength <- 20 # length of sequences

1:ATCGTGCACA ACTTGCAACA
2:ATCGTGGACC ACTTGCAACT
3:AGCGTGGACC ACTTGCAACT


# number of pairwise differences between sequences
d12 <- sum(strsplit("ATCGTGCACA ACTTGCAACA", "")[[1]] != strsplit("ATCGTGGACC ACTTGCAACT", "")[[1]])
d23 <- sum(strsplit("ATCGTGGACC ACTTGCAACT", "")[[1]] != strsplit("AGCGTGGACC ACTTGCAACT", "")[[1]])
d13 <- sum(strsplit("ATCGTGCACA ACTTGCAACA", "")[[1]] != strsplit("AGCGTGGACC ACTTGCAACT", "")[[1]])


# estimation of Tajima's theta
thetaT <- ((d12 + d23 + d13)/((nseq*(nseq-1))/2))/seqlength


# question 4B
# Calculate the Watterson estimator θw per site from this data.
# To do this I need to calculate the number of segregating sites (S) and the harmonic number (a1).

num_segregating_sites <- 5

thetaW <- num_segregating_sites * (1/sum()
                                   
                                   
                                   
                                   
                                   
                                   
# question 4D
# Plot a histogram of the site frequency specturm (SFS) for this data. 
# Do the data contain more, fewer or the same number of singletons as expected under 
# the standard neutral coalescent model?

barplot(c(3,1,4), xlab = "Derived allele frequency", ylab = "Counts")
