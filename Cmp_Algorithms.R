
## Step.4
## Use nine algorithms to get the results.
## The nine algorithms are used to experiment with the data, 
## including gs,hc,iamb,mmpc,rsmax2,tabu,fastiamb,interiamb,mmhc.
 
library(bnlearn) 
library(igraph)
library(Rgraphviz) 

## Function : CompareAlgorithm
## methodd : discrete method.
## dimsize : discrete value.
## sourcee : input data / Source data.
## cmpdata : the right/standard dataset.
## outfile : the path of output file.
CompareAlgorithm <- function(sourcee,cmpdata,methodd,dimsize,outfile){  
  
    mydata <-  read.table(sourcee, header = FALSE)
    Alarm1_graph <- read.table(cmpdata)
    store_text <- outfile
    t <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
    rownames(t) <- rownames(mydata)
    colnames(t) <- colnames(mydata)
    mydata <- t
    
    ## discretize.
    mydata <- discretize(mydata,method = 'interval',breaks = 7) #ÀëÉ¢»¯
    
    rownames(Alarm1_graph) <- colnames(Alarm1_graph)
    Alarm1_graph <- data.frame(t(Alarm1_graph))
    Alarm1_graph <- as.matrix(Alarm1_graph)
    g1 <- graph_from_adjacency_matrix(Alarm1_graph)
    
    ## Standardization of standard datasets.
    result2 <- as_edgelist(g1, names = TRUE)
    paste_12 <- function(result){
        string1 <- c()
        for(i in 1:length(result[,1])){
            t <- sort(c(result[i,1],result[i,2]))
            t <- paste(t[1],t[2])
            string1 <- c(string1,t)
        }
        string1 <- union(string1,string1)
    }
    string1 <- paste_12(result2)
    
    ## Execute algorithm 'gs', save time and compare with standard set.
    timestart4 <- Sys.time()
    gs_2 <- gs(mydata)
    timeend4 <- Sys.time()
    time_gs_2 <- c(timeend4-timestart4)
    string_gs_2 <- paste_12(gs_2$arcs)
    length(intersect(string1, string_gs_2))
    
    ## Execute algorithm 'hc', save time and compare with standard set.
    timestart1 <- Sys.time()
    hc_2 <- hc(mydata)
    timeend1 <- Sys.time()
    time_hc_2 <- c(timeend1-timestart1)
    string_hc_2 <- paste_12(hc_2$arcs)
    length(intersect(string1, string_hc_2))
    
    ## Execute algorithm 'iamb', save time and compare with standard set.
    timestart5 <- Sys.time()
    iamb_2 <- iamb(mydata)
    timeend5 <- Sys.time()
    time_iamb_2 <- c(timeend5-timestart5)
    string_iamb_2 <- paste_12(iamb_2$arcs)
    length(intersect(string1, string_iamb_2))
    
    ## Execute algorithm 'mmpc', save time and compare with standard set.
    timestart8 <- Sys.time()
    mmpc_2 <- mmpc(mydata)
    timeend8 <- Sys.time()
    time_mmpc_2 <- c(timeend8-timestart8)
    string_mmpc_2 <- paste_12(mmpc_2$arcs)
    length(intersect(string1, string_mmpc_2))
    
    ## Execute algorithm 'rsmax2', save time and compare with standard set.
    timestart9 <- Sys.time()
    rsmax2_2 <- rsmax2(mydata)
    timeend9 <- Sys.time()
    time_rsmax2_2 <- c(timeend9-timestart9)
    string_rsmax2_2 <- paste_12(rsmax2_2$arcs)
    length(intersect(string1, string_rsmax2_2))
    
    ## Execute algorithm 'tabu', save time and compare with standard set.
    timestart2 <- Sys.time()
    tabu_2 <- tabu(mydata)
    timeend2 <- Sys.time()
    time_tabu_2 <- c(timeend2-timestart2)
    string_tabu_2 <- paste_12(tabu_2$arcs)
    length(intersect(string1, string_tabu_2))
    
    ## Execute algorithm 'fastiamb', save time and compare with standard set.
    timestart6 <- Sys.time()
    fastiamb_2 <- fast.iamb(mydata)
    timeend6 <- Sys.time()
    time_fastiamb_2 <- c(timeend6-timestart6)
    string_fastiamb_2 <- paste_12(fastiamb_2$arcs)
    length(intersect(string1, string_fastiamb_2))
    
    ## Execute algorithm 'interiamb', save time and compare with standard set.
    timestart7 <- Sys.time()
    interiamb_2 <- inter.iamb(mydata)
    timeend7 <- Sys.time()
    time_interiamb_2 <- c(timeend7-timestart7)
    string_interiamb_2 <- paste_12(interiamb_2$arcs)
    length(intersect(string1, string_interiamb_2))
    
    ## Execute algorithm 'mmhc', save time and compare with standard set.
    timestart3 <- Sys.time()
    mmhc_2 <- mmhc(mydata)
    timeend3 <- Sys.time()
    time_mmhc_2 <- c(timeend3-timestart3)
    string_mmhc_2 <- paste_12(mmhc_2$arcs)
    length(intersect(string1, string_mmhc_2))
    
    ## Store the comparison results.
    write.table(TEST_DATA,append = TRUE,file = store_text ,row.names = F, quote = F)
    d <- data.frame(
        gs = c( length(string_gs_2) ,length(intersect(string1, string_gs_2))), 
        hc = c( length(string_hc_2) ,length(intersect(string1, string_hc_2))),
        iamb = c( length(string_iamb_2) ,length(intersect(string1, string_iamb_2))),
        mmpc = c( length(string_mmpc_2) ,length(intersect(string1, string_mmpc_2))),
        rsmax = c( length(string_rsmax2_2) ,length(intersect(string1, string_rsmax2_2))),
        tabu = c( length(string_tabu_2) ,length(intersect(string1, string_tabu_2))),
        fastiamb = c( length(string_fastiamb_2) ,length(intersect(string1, string_fastiamb_2))),
        interiamb = c( length(string_interiamb_2) ,length(intersect(string1, string_interiamb_2))),
        mmhc = c( length(string_mmhc_2) ,length(intersect(string1, string_mmhc_2)) )
    )
    write.table(d, append = TRUE,file = store_text, row.names = F, quote = F,sep = "\t") 
  
}




