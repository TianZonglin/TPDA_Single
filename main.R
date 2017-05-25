
## Step.0
## These packages should be initialized at first.
## In addition , they only need to be initialized during the first execution of this program.
## If these packages do not exist,please copy these files to R's library manually.

# install.packages('Rgraphviz', dependencies = TRUE) 
# install.packages('graph',depend=TRUE) 
# install.packages('arules', depend=TRUE)
# install.packages('bnlearn', depend=TRUE)
# install.packages('igraph', depend=TRUE)
# install.packages('BiocGenerics', depend=TRUE)
# install.packages('grid', depend=TRUE)
 

## Set the current working directory.

setwd(getwd())
TEST_DATA <- "IN/Alarm1_s500_v1.txt"
COMP_DATA <- "IN/Alarm1_graph.txt"

## Step.1 
## In this step, we completed the following works:
## > sampling from source data,
## > pretreatment for the aprori, 
## > get association rules,
## > re-process the results.
source("Get_Rules.R")


## Step.2
## Implementation of TPDA algorithm, the algorithm input in two cases:
## > One is the result of direct sampling of source data.
## > The other is the result sorted by the association rule after sampling.
## Then make calculations.
source("TPDA_Function.R")


## Step.3
## Use nine algorithms to get the results.



## Step.4[1]
## Call the function :TPDA_Algorithm
## For the specified experimental data, 
## here we use different threshold, dispersion, discrete method to get a set of experimental results of TPDA.
## The result is showed in the console by default.
TPDA_Algorithm(0.04,0.06,0.06,TEST_DATA,COMP_DATA)


## Step.4[2]
## Call the function :NineAlgorithmTest.
## The nine algorithms are used to experiment with the data, 
## including gs,hc,iamb,mmpc,rsmax2,tabu,fastiamb,interiamb,mmhc.
## The result is stored in : OUT/CmpResult.txt by default.
source("Cmp_Algorithms.R")
CompareAlgorithm(TEST_DATA,COMP_DATA,'interval',5,"OUT/CmpResult.txt")




