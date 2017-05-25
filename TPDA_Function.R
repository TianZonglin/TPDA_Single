## Step.2
## Implementation of TPDA algorithm, the algorithm input in two cases:
## > One is the result of direct sampling of source data.
## > The other is the result sorted by the association rule after sampling.
library(igraph)

## Function : TPDA_Algorithm.
## wt1 : threshold of the first stage.
## wt2 : threshold of the second stage.
## wt3 : threshold of the third stage.
## target  : input data / Source data.
## cmpdata : the right/standard dataset.
TPDA_Algorithm <- function(wt1,wt2,wt3,target,compdata){
  
    weight_1 <- wt1 
    weight2  <- wt2
    weight3  <- wt3
    Alarm1_graph <- read.table(compdata)
    z <- read.table(target,header = F,stringsAsFactors = FALSE)
    gene <- read.table(target,header = FALSE,stringsAsFactors = FALSE)
    
    ## prepare the source data.
    numt <- nrow(z)
    z <- data.frame(t(z))
    ronum <- nrow(z)
    z <- data.frame(c('ss'),z,stringsAsFactors = F) 
    for(s in 1:ronum){
        z[s,1] <- paste("V",as.character(s),sep="") 
    }
    row.names(z) <- z[,1]
    z <- z[,2:numt]
    z <- data.frame(t(z)) 
    
    ## record time.
    timetemp_1_start <- Sys.time()
    z <- data.frame(t(z))
    
    ## Function : getDoubleHang_value.
    ## Calculate the value of the determinant of the covariance matrix of the column.
    getDoubleHang_value <- function(mylist1,mylist2){
        Dx1 = var(mylist1)
        Dx2 = var(mylist2)
        double_h_value = Dx1*Dx2-cov(mylist1,mylist2)*cov(mylist1,mylist2)
    }
    
    ## Calculate the initial mutual information, 
    ## and the final result 'relation' is a tri-form of [from,to,value].
    mylist_a <- c()
    mylist_b <- c()
    mylist_info <- c()
    for(xi in 1:nrow(z)){
        mylist_a <- as.double(c(z[xi,]))
        Dxa = var(mylist_a)
        for(yi in 1:nrow(z)){
            #print(yi)
            mylist_b <- as.double(c(z[yi,]))
            Dxb = var(mylist_b)
            cov_ab = getDoubleHang_value(mylist_a,mylist_b)
            mi = 0.5*log(Dxa*Dxb/cov_ab)
            mylist_info <- c(mylist_info,mi)
        }
    }
    info <- matrix(mylist_info,nrow(z),nrow(z),byrow = TRUE)  
    
    rname <- c(row.names(z))
    row.names(info) <- row.names(z)
    colnames(info) <- row.names(z)
    
    ## According to the mutual information to determine the relationship.
    relation <- data.frame()
    for(xi in 1:(nrow(info)-1)){
        yi = xi
        while(TRUE){
            yi = yi+1
            if(nrow(relation) == 0){  
                from <- c(rname[xi])
                to <- c(rname[yi])
                value <- c(info[xi,yi])
                relation <- data.frame(from,to,value,stringsAsFactors = FALSE)
                relation <- data.frame(xi = c(t(relation)),stringsAsFactors = FALSE) 
            }else
                relation <- data.frame(relation,xi = c(rname[xi],rname[yi],info[xi,yi]),stringsAsFactors = FALSE)
            if(yi == nrow(info)) break
        }
    }
    
    ## Further processing of the results.
    row.names(relation) <- c("from","to","value")
    relation <- data.frame(t(relation),stringsAsFactors = FALSE)  
    from <- relation[,1][order(as.numeric(relation$value),decreasing = T)]  
    to <- relation[,2][order(as.numeric(relation$value),decreasing = T)]
    value <- relation[,3][order(as.numeric(relation$value),decreasing = T)]
    relation2 <- data.frame(from,to,value)
    mark1 <- 1
    for(wx in 1:(nrow(relation2))){
        if( as.numeric(as.character(relation2[wx,3])) > weight_1 ){ 
            #print(paste(as.character(relation2[wx,3]),weight_1))
            mark1 = mark1 +1 
        }
    }
    result <- data.frame()
    result <- relation2[1:mark1,]
    timetemp_1_end <- Sys.time() 
 
    result$from <- as.character(result$from)
    result$to <- as.character(result$to)
    gene <- data.frame(t(gene))
    rowgene_names <- c(row.names(gene))
    gene <- data.frame(c('s'),gene,stringsAsFactors = F)
    gene[,1] <- rowgene_names
    
    build_data <- result
    labe <- sort(union(build_data$from,build_data$to))  
    
    ## Function : Info
    ## Calculate mutual information,
    ## Note:'gxi' is the origin node,'gyi' is the target node,
    ## 'cutset' is the cut set, 'gen' is the whole node.
    Info <- function(gxi,gyi,cutset,gen){
        cutset <- data.frame(cutset,stringsAsFactors = F)
        x <- c()
        y <- c()
        cutset_table <- data.frame()
        for(xi in 1:nrow(gen)){
            if(gxi == gen[xi,1]){
                x <- as.numeric(c(gen[xi,2:length(gen)]))
            }
            if(gyi == gen[xi,]){
                y <- as.numeric(gen[xi,2:length(gen)])
            }
            for(yi in 1:nrow(cutset)){
                if(gen[xi,1] == cutset[yi,1]){
                    if(nrow(cutset_table) == 0)
                        cutset_table <- data.frame(t(gen[xi,]),stringsAsFactors = F)
                    else
                        cutset_table <- data.frame(cutset_table,t(gen[xi,]),stringsAsFactors = F)
                }
            }
        }
        cutset_table <- data.frame(t(cutset_table),stringsAsFactors = F) 
        cutset_table <- cutset_table[,-1]
        
        to_list <- c()
        cut <- data.frame()
        for(xi in 1:nrow(cutset_table)){
            to_list <- as.numeric(c(cutset_table[xi,]))
            if(nrow(cut) == 0)
                cut <- data.frame(to_list)
            else
                cut <- data.frame(cut,to_list)
        }
        cut_x <- data.frame(cut,x)
        cut_y <- data.frame(cut,y)
        cut_x_y <- data.frame(cut,x,y)
        t = det(cov(cut_x))*det(cov(cut_y))/(det(cov(cut))*det(cov(cut_x_y)))
        cmi = 0.5*log(t)
    }
    
    ## Function : Translat.
    ## Implement the transformation of form 'n.x' into form 'GRMZMxxx'.
    Translat <- function(data_nx){
        retdata <- data_nx
        for(i in 1:length(data_nx)){
            retdata[i] <- labe[which( trans_labe == data_nx[i])]  
        }
        return(retdata)
    }
    
    ## Implement the transformation of form 'GRMZMxxx' into form 'n.x'.
    trans_labe <- c()
    for(i in 1:length(labe)){ trans_labe <- c(trans_labe,paste("n.",i,sep = "")) }
    ct <- c()
    for(i in 1:length(build_data$from)){ 
        ct <- c(ct,build_data[i,1],build_data[i,2]) 
    }
    for(i in 1:length(ct)){ 
        ct[i] <- paste("n.",which(labe == ct[i]),sep = "") 
    }
    weight_2 <- weight2    
    weight_3 <- weight3   
    
    ## Start the first stage [I] of TPDA algorithm.
    ## Through the previous association rules mining, we can get the current node of the extension of the order, 
    ## in which we can create the initial Bayesian network structure in here.
    print('Executing the Step I of TPDA...')
    timestatus_1_start <- Sys.time()
    
    graE <- c()
    graR <- c()
    for(i in seq(1,length(ct),2)){  
        if( length(union(graE,graE)) !=length(union(graE,c(ct[i],ct[i+1])))  ){  
            graE <- c(graE,ct[i],ct[i+1])
            g <-  graph(graE, directed = T)  
        }else{
            if(edge_connectivity(g, source = ct[i], target = ct[i+1], checks = TRUE) == 0){
                graE <- c(graE,ct[i],ct[i+1])
                g <-  graph(graE, directed = T)
            }
            else{
                graR <- c(graR,ct[i],ct[i+1])
            }
        }
    }
    
    #plot(g,layout = layout.fruchterman.reingold, vertex.size = 5,vertex.color = "green")
    timestatus_1_end <- Sys.time()
    print('Step I is OK.')
    
    ## Start the second stage [II] of TPDA algorithm.
    ## Use a more complex conditional independence test to determine 
    ## which pairs of nodes should also be added to the network.
    print('Executing the Step II of TPDA...')
    timestatus_2_start <- Sys.time() 
    g <- graph(graE, directed = F)
    
    for(i in  seq(1,length(graR),2)){
      
        ## Find the path and storage it to 'one_path' with form 'n.x'.
        shortpa <- shortest_paths(g, from = graR[i], to = graR[i+1], mode = c("all"))$vpath
        one_path <- names(V(g))[as.integer(shortpa[[1]])]   
        
        ## Calculate cut points and mutual information.
        brek <- one_path[-c(1,length(one_path))]  
        if(length(brek) == 0){print("there has 0 error!")}
        info <- Info(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),gene)  
        #print(Translat(brek))
        
        ## If the mutual information is greater than 'weight_2' then add it to 'graE'.
        if(info > weight_2){ graE <- c(graE,graR[i],graR[i+1]) }  
    }
    rm("one_path","i","shortpa","info","brek","graR") 
    timestatus_2_end <- Sys.time()
    print('Step II is OK.')
    
    ## Start the second stage [III] of TPDA algorithm.
    ## Each edge is checked and the second stage of the formula is used for the independence test to 
    ## determine whether the node is independent of the condition.If the two nodes are conditional, their edges will be deleted.
    print('Executing the Step III of TPDA...')
    g <- graph(graE, directed = F)
    timestatus_3_start <- Sys.time()
    
    ## Delete the current edge and store it in 'g_d',and find whether there is still a path in 'g_d',
    ## if there is a path then calculate the mutual information,if not,then skip.
    for(i in  seq(1,length(graE),2)){
        g_d <- g - edge(paste(graE[i],"|",graE[i+1],sep = ""))   
        
        if( edge_connectivity(g_d, source = graE[i], target = graE[i+1], checks = TRUE) > 0){
          
            shortpa <- shortest_paths(g_d, from = graE[i], to = graE[i+1], mode = c("all"))$vpath
            one_path <- names(V(g))[as.integer(shortpa[[1]])]
            brek <- one_path[-c(1,length(one_path))]  
            if(length(brek) > 0){
                info <- Info(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),gene)  
                 
                ## If the mutual information is less than 'weight_3', then delete it.
                if(info < weight_3){ g <- g_d }  
            }else{
                #print(i)
            }
        }
    }
    timestatus_3_end <- Sys.time()
    print('Step III is OK.')
    
    ## calculate the run time.
    print('drafting')
    time1 <- c(timestatus_1_end-timestatus_1_start)
    time1 <- c(time1 + timetemp_1_end-timetemp_1_start)
    print(time1+0)
    print('thickening')
    time2 <- timestatus_2_end-timestatus_2_start
    print(time2+0)
    print('thinning')
    time3 <- timestatus_3_end-timestatus_3_start
    print(time3+0)
    print('all')
    print(time1+time2+time3+0)
 
    result <- as_edgelist(g, names = TRUE)
    for(i in 1:length(result[,1])){
        result[i,] <- Translat(result[i,])
    }

    ## read the standard dataset.
    rownames(Alarm1_graph) <- colnames(Alarm1_graph)
    Alarm1_graph <- as.matrix(Alarm1_graph)
    g1 <- graph_from_adjacency_matrix(Alarm1_graph)
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
    
    ## print the output data.
    string_tdpa <- paste_12(result)
    print(paste("Right:",length(intersect(string1, string_tdpa)) ))
    print(paste("Total:",length(result[,1]) ))
    
    print(paste("Probability:",length(intersect(string1, string_tdpa))/length(result[,1])))

}


