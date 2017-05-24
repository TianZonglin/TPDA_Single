## Step.1 
## In this step, we completed the following works:
## > sampling from source data,
## > pretreatment for the aprori, 
## > get association rules,
## > re-process the results.



library("arules")



## This part is the pretreatment of association rules.
AllG <- read.table(paste(getwd(),TEST_DATA,sep = '/'))
jRow <- length(rownames(AllG))
jFin <- paste(getwd(),TEST_DATA,sep = '/')
jFout <- paste(getwd(),"EX/rules/Mult_data_abc.txt",sep = '/')
system("javac preTreatment.java")
system(paste("java preTreatment",jRow,jFin,jFout,sep = ' '))




##  This part use the encapsulated function 'apriori' to complete the analysis of association rules.

tr <-  read.transactions("EX/rules/Mult_data_abc.txt", format="basket", sep = ",")
rules = apriori(tr,parameter = list(support = 0.3,confidence = 0.6,minlen = 2,maxlen = 2)) #pre 0.3  0.6
summary(rules)     
outdata <- as(rules,"data.frame")   
write.table(outdata,"EX/rules/Result_abc.csv",sep = ",")   


# This part is the further processing of the results after association rules.
jRules <- length(rules)
jFin <- paste(getwd(),"EX/rules/Result_abc.csv",sep = '/')
jFout <- paste(getwd(),"EX/rules/sum_data_abc.csv",sep = '/')
system("javac furTreatment.java")
system(paste("java furTreatment",jRules,jFin,jFout,sep = ' '))


 
