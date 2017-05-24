## Step.2
## To confirm whether the source data need to sort by the association rules, 
## if not, skip to the third step.



fin_1  <- TEST_DATA
fin_2  <- "EX/rules/sum_data_abc.csv"
fout_1 <- "EX/rules/genNode_abc.csv"
fout_2 <- "EX/rules/express_after_abc.csv"
 
ruleBack <- read.table(fin_1,header = TRUE)
ruleBack <- data.frame(ruleBack)
ruleBack <- t(ruleBack)
namess <- c()
for(s in 1:nrow(gen)){
  namess <- c(namess,paste('V',s,sep = ''))
}
rownames(ruleBack) <- namess

gen <- t(ruleBack)
gen <- data.frame(t(gen))
rownum <- nrow(gen)
gen <- data.frame(c('s'),gen,stringsAsFactors = F)
gen <- data.frame(gen[,1:2])
#gen[,2] <- namess

## give the name
for(s in 1:rownum){
    gen[s,1] <- paste(as.character(rownames(gen)[s]),sep = "")
}
gen[,1] <- as.character(gen[,1])
order <- c(nrow(gen))
gen <- data.frame(order,gen)
gen <- data.frame(gen[,1:2])


## Read the result of the Aprori.
hc <- read.csv(fin_2,header = T,stringsAsFactors = F)
hc <- hc[,-2] 
hc <- hc[,-3] 
hc <- hc[,-3]   

order <- c(nrow(hc))
from = hc[,1]
to = hc[,2]
row_names <- c(row.names(hc))
hc40ft <- data.frame(row_names,from,to,order,stringsAsFactors = F)


tt = 0;
count = 0;
gencount = 0;
mylist <- c()
mycount <- c();
mygencount <- c();
hc40ft <- data.frame(hc40ft,tag = c(0))

for(ix in 1:nrow(hc40ft)){
    for(iy in 1:nrow(hc40ft)){
        if( ix == 79 && iy == 98 ){
            print(hc40ft[79,2])
            print(hc40ft[98,3])
            if(hc40ft[ix,2] == hc40ft[iy,3]){
                print("mark...")
                print(hc40ft[ix,5])
            }
          
        }
        if(hc40ft[ix,2] == hc40ft[iy,3]&&hc40ft[ix,5] == 0){
            hc40ft[ix,5] = 1
            mylist <- c(mylist,hc40ft[ix,1])
          
        }
    }
}
mylist2 <- c()
for(xi in 1:nrow(hc40ft)){
    t = 0
    for(yi in 1:length(mylist)){
        if(hc40ft[xi,1] == mylist[yi]){
            t = 1
            break
        }
    }
    if(t == 0){
        mylist2 <- c(mylist2,hc40ft[xi,1])
    }
}
if(length(mylist2) == 0){
    print("----- There is no head node but ring -----")
    stop();
}


## Number for the head node of hc40ft_1
hc40ft_1 <- hc40ft[mylist2,] 
for(xi in 1:nrow(hc40ft_1)){
    for(yi in 1:nrow(hc40ft)){
        if(hc40ft_1[xi,1] == hc40ft[yi,1]&&hc40ft[yi,4] > count){
            count = count+1;
            hc40ft_1[xi,4] = count;
            hc40ft[yi,4] = count;
            break;
        }
    }
}
 

## Number for the head node of gen
for(xi in 1:nrow(hc40ft_1)){
    for(yi in 1:nrow(gen)){
        if(hc40ft_1[xi,2] == gen[yi,2]&&gen[yi,1]>gencount){
            gencount = gencount+1;
            gen[yi,1] = gencount;
            break;
        }
    }
}


mycount <- c(mycount,count);
mygencount <- c(mygencount,gencount);

## Number for the next node of the head node.(di er ci)
for(zz in 1:1000){
    #print(zz)
    for(xi in 1:nrow(hc40ft_1)){
          for(yi in 1:nrow(gen)){
            if(hc40ft_1[xi,3] == gen[yi,2]&&gen[yi,1]>gencount){
                gencount = gencount+1;
                gen[yi,1] = gencount;
                break;
            }
        }
    }
  
    ## Select the next node of the head node
    mygencount <- c(mygencount,gencount);
    tt = count;
    mylist3 <- c()
    mylist4 <- c()
    for(xi in 1:nrow(hc40ft_1)){
        for(yi in 1:nrow(hc40ft)){
            if(hc40ft_1[xi,3] == hc40ft[yi,2]&&hc40ft[yi,4]>count){
              mylist3 <- c(mylist3,hc40ft_1[xi,4])
              mylist4 <- c(mylist4,hc40ft[yi,1])
            }
        }
    }
    if(length(mylist3) == 0&&length(mylist4) == 0){ break;}
    mydata_f_t <- data.frame(mylist3,mylist4)
    
    ## hc40ft, excute the second number 
    for(xi in 1:nrow(mydata_f_t)){
        for(yi in 1:nrow(hc40ft)){
            if(hc40ft[yi,1] == mydata_f_t[xi,2]&&hc40ft[yi,4]>count){
                count = count+1;
                hc40ft[yi,4] = count;
                break;
            }
        }
    }
    
    ## remove duplicates
    p = 0
    mycount <- c(mycount,count);
    if(tt == count||gencount == (nrow(gen)-1)){break;}
    hc40ft_2 <- hc40ft[mylist4,]
    for(xi in 1:nrow(hc40ft_2)){
        q = 0
        ll <- c(strsplit(row.names(hc40ft_2)[xi-p],""))          
        for(yi in 1:length(ll[[1]])){     
            if(ll[[1]][yi] == "."){q = 1}
        }
        if(q == 1){
            hc40ft_2 <- hc40ft_2[-(xi-p),]
            p = p+1
        }
        if((xi-p) >= nrow(hc40ft_2)){break;}
    }
    hc40ft_1 <- data.frame();
    mydata_f_t <- data.frame();
    hc40ft_1 <- hc40ft_2;

}


## Now the gene table structure for the two columns, 
## left is the node name, the right is the node order, and then write the file.
cc = 0
for(ix in 1:nrow(gen)){
    if(gen[ix,1] == nrow(gen))
        #print(paste(gen[ix,1]," ",nrow(gen)))
    cc = cc+1;
}

names(gen) <- c('order','gene')
gen <- cbind(gen$gene, gen$order)

write.csv(gen,fout_1,row.names = F)
if(length(mylist3) == 0&&length(mylist4) == 0){
  print("----- There has noob gene -----")
}


## The last sort
## 2016-6-24 changed:
## Now,the output file is the result of the sequence by adjusting the source data. 
sx <- read.csv(fout_1,header = T) 
zero <- t(ruleBack)
zero <- data.frame(t(zero))
one <- data.frame(sx,zero)
two <- one[order(one[,2],decreasing = T),] #Order [ descending 100...1 ]
two <- two[,-2]
two <- t(two)
write.csv(two,fout_2,row.names = F)


print("----- Sort by aprori completed! -----")



