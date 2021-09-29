# Function to get modules from network adjacency matrix
findModules.speakeasy <- function(adj, iter = 10, timesteps = 20, module.pval = 0.05, module.min.size=30,doPar = TRUE){
  
  layers=1
  GeneNames <- colnames(ADJ)

  # Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




##This is the suggestedtt fuction through which to interact with SpeakEasy
#clustering.  The loops in this funtion are related to calling SpeakEasy on
#subsets of data.  This maybe useful when there are genuine hierarchies of clusters in the
#data, although you cannot always force SpeakEasy to split clusters (unlike hierarchical clustering).
#
#Description of inputs:
#"layers" is number of times to subcluster: should be integer >=1 (1 is no subclustering)  If you want to avoid subclustering clusters smaller than a given size, enter that as the 2nd entry i.e. [3,15] does three rounds of clustering and won't touch a cluster unless it has at least 15 members
#"iter" is number of replicate runs - depends on how finely you want to estimate cluster confidence, but 1 is min value and 10+ are useful especially when overlapping output is requested
#"ADJ" is an double, sparse or full adjacency matrix (weighted, unweighted; negative weights are fine too)
#"timesteps" is number of times to repeat the label selection - 30 tends to handle even huge matrices, as the number of labels collapses rapidly in real networks after the first few timesteps
#varargin - optional input is "multi_community_switch" in the range of [0 1], which if <1 allows overlapping clusters - closer to 0 is more lenient definition of overlapping nodes.
#We suggest a value of 1/(max# of overlapping communities per node desired).
#For instance if you want nodes to be in no more than 3 clusters, multi_community_switch=.33333
#
#Description of outputs:
#"partition_codes_overlapping" is a two column matrix with node ID's in the first column and numeric label ID's in the 2nd
#"cell_partition_overlapping" is a sequence of cells, each containing a list of nodes in a given cluster.  The clusters are ordered in size.
#convenient_node ordering{i} is imply vertcat(partition_n_idx_overlapping_cell{i}{:})
#even though the output variables contains the word "overlapping", if disjoint output is requested, that's what these will contain
#
#to check that your ADJ was clustered reasonably, try running: imagesc(ADJ(convenient_node_ordering{1},convenient_node_ordering{1})


virtual_cooccurrence <- function(ADJ,partitions,partitionID, main_iter, accept_multi){
    library(pracma)
    print(dim(partitions))
    print(dim(as.matrix(partitionID)))
    print(length(partitionID[[1]]))
#    partitions <- t(as.matrix(unlist(partitions)))
    
    adjustedrand <- function(partitionA,partitionB){
        
        ng1 <- max(partitionA)
        ng2 <- max(partitionB)
        #ctabmat=full(sparse(partitionA,partitionB,1,ng1,ng2)); slightly faster, but memory intense for large matrices
        ctabmat <- sparseMatrix(i=as.vector(partitionA),j=as.vector(t(partitionB)),x=1L,
                                symmetric = FALSE, repr='T',index1 = FALSE)#same format as lxn
       
        #dims = c(ng1,ng2),
         #SpeakEasy passes in sequentially labeled clusters, so this shouldn't be an issue
        
        n <- sum(colSums(ctabmat))
        nis <- sum(rowSums(ctabmat)^2)#sum of squares of sums of rows
        njs <- sum(colSums(ctabmat)^2)#sum of squares of sums of columns
        
        t1 <- combs(n,2)#total number of pairs of entities
        t2 <- sum(colSums(ctabmat^2))#sum over rows & columnns of nij^2
        t3 <- 0.5*(nis+njs)
        
        #Expected index (for adjustment)
        nc <- (n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1))
        
        A <- (t1+t2-t3)#no. agreements
        D <-   (-t2+t3)#no. disagreements
        
        if (is.na(nc)){
            ARI <- 0
        } else if(t1==nc){
            ARI <- 0 #avoid division by zero; if k=1, define Rand = 0
        } else {
            ARI <- (A-nc)/(t1-nc) #adjusted Rand - Hubert & Arabie 1985
        }
        return(ARI)
    }
    
    fast_ind2sub <- function(matsize,idx){
        
        nrs <- nrow(matsize)
        ncs <- ncol(matsize)
        r <- mod((idx-1),nrs) + 1
        c <- ((idx-r)/nrs) + 1
        res <- c(r,c)
        return(res)
        
    }
    
    for (i in 2:size(partitions,2)){           #make sure partition codes are distinct across partitions
        partitions[,i] <- partitions[,i]+ max(partitions[,i-1])
    }
    
    adjustedrand_pairwise <- zeros(size(partitions,2))#holds all possible adjusted rand index comparisons of partitions
    for (i in 1:size(partitions,2)){
        # Dev
         print(i)
        for (j in 1:size(partitions,2)){
            
            if (j<i){#metric is symetric, so save time by computing only half
                adjustedrand_pairwise[i,j] <- adjustedrand(partitionID[[as.character(i)]],partitionID[[as.character(j)]])  
                
                # Dev
                # print(j)
            }
        }
        print('done')
    }
    
    adjustedrand_pairwise <- adjustedrand_pairwise+ t(adjustedrand_pairwise)
    most_similar_var <-  max(colSums(adjustedrand_pairwise))
    most_similar_idx <- which.max(colSums(adjustedrand_pairwise))#select representative partition
    winning_partition <- partitions[,most_similar_idx]
    winning_members_unq <- unique(partitions[,most_similar_idx])
    print("winning_members_unq")
    print(length(winning_members_unq))
    
    cell_partition <- list()#get the indices of nodes which end up in the same cluster
    for (i in 1:length(winning_members_unq)){
        cell_partition[[i]]= c(which(winning_partition==winning_members_unq[i]))#seems like same contest as partitionID(most_similar_idx) but reordered
    }
    
    if (accept_multi<1){
        cat('getting overlapping clusters from virtual cooccurrence matrix',main_iter)
        all_partitions_unq <- unique(partitions)
        
        #markertitle is a list of clusterID's from all partitions
        #markerfound contains the positions of all the clusterID's in the ADJ
        #also, if clusterID==3, all the positions of nodes taked with clusterID==3 are in markerfound{3}
        sorted_labels <- sortrows(rbind(partitions,t(as.matrix(1:numel(partitions)))))
        sorted_labels <- as.matrix(sorted_labels)
        transitions=matrix(c(0, which(sorted_labels[2:nrow(sorted_labels),1]-sorted_labels[1:(nrow(sorted_labels)-1),1]!=0),size(sorted_labels,1)))#find sets of locations of the end of the previous label
        markertitle <- zeros(max(dim(sorted_labels)),1)#for all clusters, get their locations in the partitions and their numeric identifier
        markerfound = list()
        for (i in 1:(size(transitions,1)-1)){#puts the ith label in the ith cell of markerfound - remember that labels are unique across partitions
            markerfound[[i]] <- sorted_labels[c((transitions[i]+1):transitions[i+1]),1]#this actually gets very big too
            markertitle[i,1] <- sorted_labels[transitions[i]+1,1]
        }
        
        nodes_and_partition_identifiers_overlapping <- list()
        multi_store <- list()
        mutlicom_node_cell <- list() 
        cell_partition_overlapping <- list() 
        ban_module_size <- 3
        nodes_x_partitions_coocvals_holder <- list()
        
        banned_modules <- which(lengths(cell_partition) < ban_module_size)#sometimes a node may be half in it's "own cluster" (size one) and in some other, so not really multi-comm
        for (i in 1:length(winning_members_unq)){#for each cluster in winning partition
            
            if (!(i %in% (c(0, banned_modules)))){
                winning_partition_idxs=which(winning_partition==winning_members_unq[i],arr.ind = T)
                winning_partition_chunk <- partitions[winning_partition_idxs,]#all labels ever used for a cluster in all other partitions
                all_labels_for_winning_cluster_unq <- unique(winning_partition_chunk)
                countseach <- pracma::histc(as.matrix(winning_partition_chunk),all_labels_for_winning_cluster_unq)
                countseach <- countseach$cnt
                
                fraction_counts_labels <- list()
                for (j in 1:length(all_labels_for_winning_cluster_unq)){ #for all labels that have ever characterized any node in a winning cluster
                    
                    linear_idx = markerfound[[all_labels_for_winning_cluster_unq[j]]]
                    tt <- as.matrix(countseach[j,])/size(winning_partition_chunk,1)
                    tt <- tt*(ones(length(linear_idx) ,1))
                    fraction_counts_labels[[j]] <- tt
                }
                
                for (k in 1:length(fraction_counts_labels)){
                    if(k == 1){
                        fraction_counts_labels_db <- data.frame("fraction_counts_labels" = fraction_counts_labels[[1]])
                    }else {
                        fraction_counts_labels_db <- rbind(fraction_counts_labels_db, data.frame("fraction_counts_labels" = fraction_counts_labels[[k]]))
                    }
                }
                
                
                #fraction_counts_labels <- vertcat(fraction_counts_labels{:})
                
                for (k in all_labels_for_winning_cluster_unq){
                    if(k == all_labels_for_winning_cluster_unq[1]){
                        verted_markerfound <- data.frame("mark" = markerfound[[k]])
                    }else {
                        verted_markerfound <- rbind(verted_markerfound, data.frame("mark" = markerfound[[k]]))
                    }
                }
                
                fast_res <- fast_ind2sub(size(partitions),verted_markerfound)
                row_fs <- fast_res[1]
                col_fs <- fast_res[2]
                nodes_x_partitions_coocvals <- sparseMatrix(i=as.vector(row_fs),j=as.vector(t(col_fs)),x=fraction_counts_labels,
                                                            symmetric = FALSE, repr='T',dims = c(size(partitions,1), size(partitions,2)))
                # in nodes_x_partitions_coocvals rows are nodes and columns are partitions, values are mean of co-occurence with all members of a (winning) cluster using data from a non-winning partition
                nodes_x_partitions_coocvals[winning_partition_idxs,] <- 0
                nodes_x_partitions_coocvals_holder[[i]] <- nodes_x_partitions_coocvals
                
                mutlicom_node_cell[[i]]       <- which(rowMeans(nodes_x_partitions_coocvals)>accept_multi)#it's a bit slow to take these means sequentially, but I've tried doing them all at once and it ends up taking up too much ram for networks of 100K+
                cell_partition_overlapping[[i]] <- rbind(as.matrix(mutlicom_node_cell[[i]]),as.matrix(cell_partition[[i]]))#i is the target cluster and the source is whatever cluster the node was originally in
                nodes_and_partition_identifiers_overlapping[c(rbind(as.matrix(mutlicom_node_cell[[i]]),as.matrix(cell_partition[[i]])))] <- i
                
            } else {#for small/banned modules
                mutlicom_node_cell[[i]] <- c()
                cell_partition_overlapping[[i]] <- cell_partition[[i]]
                
            }
        }
        
        for (k in 1:length(mutlicom_node_cell)){
            if(k == 1){
                verted_mutlicom_node_cell <- data.frame("vert" = mutlicom_node_cell[[1]])
            }else {
                verted_mutlicom_node_cell <- rbind(verted_mutlicom_node_cell, data.frame("vert" =mutlicom_node_cell[[k]]))
            }
        }
        
        saveRDS(multicom_nodes_all,"multi_com_list.RDS")
    } else {
        cell_partition_overlapping <- cell_partition
        multicom_nodes_all <- c()
        saveRDS(multicom_nodes_all,"multi_com_list.RDS") #just so not confused with old result
    }
    
    if (accept_multi<1){
        cat('overlapping vs discrete length: ', as.character(sum(unlist(lapply(cell_partition_overlapping,length)))), ' vs ' ,str(sum(unlist(lapply(cell_partition,length)))), str(main_iter))
    }
    
    
    #from here on it's just arranging the output
    idx_large_partition <- c()
    for (t in 1:length(cell_partition)){
        idx <- pracma::size(cell_partition[[t]],1)
        idx_large_partition <- c(idx_large_partition,idx)
    }
    idx_large_partition <- order(idx_large_partition, decreasing = TRUE)
    cell_partition <- cell_partition[idx_large_partition]# reorder partitions by size, with no change to contents
    
    dx_large_partition <- c()
    for (t in 1:length(cell_partition_overlapping)){
        idx <- pracma::size(cell_partition_overlapping[[t]],1)
        dx_large_partition <- c(dx_large_partition,idx)
    }
    dx_large_partition <- c(order(dx_large_partition, decreasing = TRUE))
    cell_partition_overlapping <- cell_partition_overlapping[dx_large_partition]
    
    
    
    cluster_density <- list()#sort order of nodes within each cluster for display purposes
    for (i in 1:length(cell_partition)){
        cluster_density[[i]] <- colMeans(as.matrix(ADJ[unlist(cell_partition[[i]]),unlist(cell_partition[[i]])]))
        ind <- order(unlist(cluster_density[[i]]), decreasing=TRUE)
        cell_partition[[i]] <- cell_partition[[i]][ind]
    }
    #best_nodeorder_hard=vertcat(cell_partition{:});
    
    
    cluster_density_overlapping <- list() #sort order of nodes within each cluster for display purposes
    for (i in 1:length(cell_partition_overlapping)){
        cluster_density_overlapping[[i]] <- colMeans(as.matrix(ADJ[unlist(cell_partition_overlapping[[i]]),unlist(cell_partition_overlapping[[i]])]))
        ind <- order(unlist(cluster_density_overlapping[[i]]), decreasing=TRUE)
        cell_partition_overlapping[[i]] <- cell_partition_overlapping[[i]][ind]
    }
    
    #best_nodeorder_hard=vertcat(cell_partition_overlapping{:});
    
    
    partition_marker_sorted_hard <- list()
    partition_marker_sorted_overlapping <- list()
    for (i in 1:length(cell_partition)){
        partition_marker_sorted_hard <- c(partition_marker_sorted_hard,rep(i,length(cell_partition[[i]])))
        partition_marker_sorted_overlapping <- c(partition_marker_sorted_overlapping,rep(i,length(cell_partition_overlapping[[i]])))   
    }
    
    for (k in 1:length(cell_partition)){
        if(k == 1){
            cell_partition_db <-data.frame("cell_part" = cell_partition[[1]])
        }else {
            cell_partition_db <- rbind(cell_partition_db, data.frame("cell_part" = cell_partition[[k]]))
        }
    }
    
    
    print("Cell Partition Length")
    print(dim(cell_partition_db))
    print("partition_marker_sorted_hard Length")
    print(length(partition_marker_sorted_hard))
    
    print(dim(data.matrix(data.frame("V1" = cell_partition_db$cell_part,"V2" = unlist(partition_marker_sorted_hard)))))
    nodes_and_partition_identifiers_hard <- sortrows(data.matrix(data.frame("V1" = cell_partition_db$cell_part,"V2" = unlist(partition_marker_sorted_hard))))
    print("nodes_and_partition_identifiers_hard Length")
    
    print(dim(nodes_and_partition_identifiers_hard))
    
    partition_marker_sorted_hard <- list()
    for (i in 1:length(cell_partition_overlapping)){
        partition_marker_sorted_hard <- c(partition_marker_sorted_hard,rep(i,length(cell_partition_overlapping[[i]])))
    }
    for (k in 1:length(cell_partition_overlapping)){
        if(k == 1){
            print(length(cell_partition_overlapping[[1]]))
            cell_partition_db_2 <- data.frame("cell_part" =cell_partition_overlapping[[1]])
        }else {
            cell_partition_db_2 <- rbind(cell_partition_db_2, data.frame("cell_part" = cell_partition_overlapping[[k]]))
        }
    }
    print("Cell PartitionDB")
    print(dim(cell_partition_db_2))
    nodes_and_partition_identifiers_overlapping <- sortrows(data.matrix(data.frame("V1" = cell_partition_db_2$cell_part,"V2"= unlist(partition_marker_sorted_hard))))
    
    
    record_stuff <- 0#cluster cood density stats
    # if (record_stuff==1){
    # 
    #     partition_rows <- zeros(length(cell_partition),length(cooc))
    #     
    #     cooc_temp <- cooc
    #     for (i in 1:length(cell_partition)){
    # 
    #         partition_rows[i,:] <- sum(cooc_temp(cell_partition[[i]],:))./(length(cell_partition{i})-1)
    # 
    #     }
    # 
    #     [cluster row]=find(partition_rows!=0)
    #     value <- partition_rows(find(partition_rows>0))
    # 
    #     if (isempty(accept_multi)){
    #         csvwrite('SpeakEasy_cluster_assignment.csv',[row cluster value])
    #     }
    print('new res')
    results <- list("nodes_and_partition_identifiers_hard"=nodes_and_partition_identifiers_hard,"nodes_and_partition_identifiers_overlapping"=nodes_and_partition_identifiers_overlapping,
                    "cell_partition"=cell_partition, "cell_partition_overlapping"=cell_partition_overlapping, "multicom_nodes_all"=multicom_nodes_all)
    return(results)
}

SpeakEasycore <- function(ADJ,total_time,IC_store_relevant,nback,force_efficient){
    library(pracma)
    kin <- as.matrix(colSums(ADJ))
    
    
    
    aggregate_labels <- function(ADJ,labels,nback){
        library(pracma)
        labels <- data.matrix(as.vector(labels))
        labels_unq <- unique(labels)
        
        #sorted_labels=sortrows([labels (1:length(labels))']);  #two columns of labels and label locations
        indices <- sort(labels, index.return=TRUE)$ix#sort faster than sortrows
        labs <- c(1:length(labels))
        temp <- data.matrix(rbind(as.vector(labels),as.vector(labs)))
        temp <- t(temp)
        sorted_labels <- temp[as.vector(indices),]
        transitions=data.matrix(c(0, pracma::finds(unlist(sorted_labels[(2:nrow(sorted_labels)),1])-unlist(sorted_labels[(1:(nrow(sorted_labels)-1)),1]) != 0),dim(sorted_labels)[1]))#find sets of locations of the end of the previous label
        # node_identifiers=zeros(length(labels),1);
        # for i=1:size(transitions,1)-1 #relabel the... labels with sequential integers, starting with 1, which becomes a time-saver later
        #      node_identifiers(sorted_labels(transitions(i)+1:transitions(i+1),2))=i;  #maybe replace with cumsum on sprasemat based on transitions, if that's fast
        # end
        
        future_markers <- zeros(size(sorted_labels,1),1)
        future_markers[1] <- 1
        for (k in (transitions[2:length(transitions)-1])){
            if (k != 0){
                future_markers[k] = 1
            }
        }
        future_markers <- cumsum(future_markers)
        node_identifiers <- zeros(size(sorted_labels,1),1)
        node_identifiers[as.vector(sorted_labels[,2]),1] = future_markers
        
        
        #idea is to consider the labels from different time-steps to be different
        #(even though they are not) then add the row of temp created with these pseudo-different labels, to add up the rows that are infact related to the #same core label
        #the obvious way to do this is:
        #temp=zeros(length(labels_unq), length(node_identifiers));
        #temp=bsxfun(@eq, sparse(1:length(labels_unq)).', node_identifiers');
        #but a faster way is:
        j = t(as.matrix(1:length(node_identifiers)))
        nodes_by_labels_all_times <- sparseMatrix(i= as.vector(node_identifiers), j=t(c(1:length(node_identifiers))), 
                                                  x = as.vector(t(ones(length(node_identifiers),1))),
                                                  symmetric = FALSE,repr="T")
        #we can do that becaue node identifiers are numbers sequentially starting at 1
        #in each row of "nodes_by_labels_all_times" we tick off positions (using a 1) where that label occurs in the full list of labels
        
        section_length <- dim(nodes_by_labels_all_times)[2]/nback
        #we've arranged the labels from several timesteps chunks arranged horizontally, and we want to sum each of those chunks independently and add them all back up (because the rows in each chunch are in fact synchronized)
        nodes_by_labels_all_times = as.matrix(nodes_by_labels_all_times)
        for (i in 1:nback){
            if (i==1){
                running_sum <- nodes_by_labels_all_times[,1:section_length]
            } else {
                start = pracma::ceil((section_length)*(i-1))
                end = floor(start+section_length-1)
                running_sum <- running_sum+nodes_by_labels_all_times[,start:end]
            }
        }
        
        lxn <- running_sum%*%ADJ
        return(lxn)
        #size(runing_sum,2)==length(ADJ), lnx will be sparse if ADJ is sparse
        #lxn has counts of each label (each row of lxn represents a different label) for each node (columns of lxn)
    }
    
    
    
    #initialize matrix to store chosen labels
    listener_history <- matrix(0,nback+total_time,max(dim(ADJ)))#each column is history for a single node, starting at row 1
    listener_history[1:(nback+1),1:max(dim(ADJ))] <- IC_store_relevant
    #nback <- 
    for (i in (2+nback):(nback+total_time)){
        print(i)
        current_listener_history <- t(listener_history[(i-nback):(i-1),])
        current_listener_history = t(as.matrix(as.vector(current_listener_history)))
        temp_agr <- aggregate_labels(ADJ,current_listener_history,nback)#actual_counts is sparse for sparse input
        actual_counts <- temp_agr 
        active_labels <- unique(as.vector(current_listener_history))
        counts_normk <- colSums(actual_counts)  #',2)' needed for case of size-1 clusters
        print('line 395 comp!')
        count_normk_norm1 <- as.matrix(counts_normk/sum(counts_normk))#proportions of various labels normalized to 1
        
        #if matrix is very sparse, or too large to store a full ADJ, we only care about generating expected counts of labels if a node actually receives some of that label
        if (force_efficient == 1){
            x_y <- which(actual_counts != 0, arr.ind = T)#x will be labels and y will be nodeID
            x <- x_y[,1]
            y <- x_y[,2]
            #two lines below are easier to understand but slightly slower separately
            #scaled_kin=nback*([full(count_normk_norm1(x))]'.* kin(y)); #scales normalized counts by total input (some nodes have more inputs and thus you would expect more of all labels)
            #expected=sparse(x,y,scaled_kin, size(actual_counts,1),size(actual_counts,2));  #same format as lxn
            expected <- sparseMatrix(i=as.vector(x),j=as.vector(t(y)),x=as.vector(nback*(as.matrix(count_normk_norm1[x])* kin[y])),
                                     symmetric = FALSE, repr='T',dims = size(actual_counts))#same format as lxn
            print('line 408 comp!')
            
        } else {
            expected <- nback*(as.matrix(count_normk_norm1)*as.matrix(kin)) # %a bit slower for ADJ's with less than 3% density compared to option above (and more mem usage) but 10x faster on sparse
        }
        #,dims = c(length(unique(x)), dim(actual_counts)[2])
        # } else {
        #     # sum(expected)==kin*nback
        #     expected <- nback*t(as.matrix(count_normk_norm1))*(as.matrix(kin))#a bit slower for ADJ's with less than 3# density compared to option above (and more mem usage) but 10x faster on sparse
        # 
        #diagnostic
        #                 [x y]=find(actual_counts);   #x will be labels and y will be nodeID
        #                 expectedsparse=sparse(x,y,nback*([full(count_normk_norm1(x))]'.* kin(y)), size(actual_counts,1),size(actual_counts,2));  #same format as lxn
        #          full(actual_counts)
        #          full(expectedsparse)
        #          full(expected)
        #          full( expectedsparse-actual_counts)
        #          full( expected-actual_counts)
        #         length(find(min(( expectedsparse-actual_counts))>0))
        print(dim(actual_counts))
        print(dim(expected))
        tem <- actual_counts-expected
        
        tem <- t(tem)
        max_vals <- apply(tem,1,max)
        max_idx <- apply(tem,1,which.max)
        # length(find(maxvals==0))   #might worry that for force_efficient==1 case max of some column will be zero, indicating a non-elegible label, but that never happens... could add 100 to each non-zero entry of actual_counts if worried, but really not necessary
        listener_history[i,] <-  active_labels[max_idx]
        print('line 433 comp!')
        
    }#time-step loop
    
    #identify nodes in same clusters
    sorted_labels <- t(sortrows(rbind(t(listener_history[(dim(listener_history)[1]-2),]),t(1:length(listener_history[(dim(listener_history)[1]-2),])))))
    transitions=matrix(c(0, pracma::finds(sorted_labels[(2:dim(sorted_labels)[1]),1]-sorted_labels[(1:dim(sorted_labels)[1]-1),1]!=0), size(sorted_labels,1)))#find sets of locations of the end of the previous label
    partitionID <- zeros(max(dim(ADJ)),1)#for all clusters, get their locations in the partitions and their numeric identifier
    label_assignment <-list()
    for (i in 1:(size(transitions,1)-1)){
        ids <- c((transitions[i]+1):transitions[i+1])
        partitionID[sorted_labels[ids,2]] <- i
        
        #this actually gets very big too
        label_assignment[[i]] <- sorted_labels[((transitions[i]+1):transitions[i+1]),2]
    }
    final_res <- list("label_assignment" = label_assignment, "partitionID" = partitionID, "listener_history" = listener_history)
    return(final_res)
}


bootstrap_SpeakEasy <- function(iter,ADJ,timesteps,nback,is_ADJ_weighted, force_efficient, main_iter, multi_community_switch,varargin){     
    
    cat('starting to setup ICs for all runs',main_iter, '\n')
    
    
    genIC_full <- function(ADJ,how_many_runs,nback,varargin){
        
        IC_store <- matrix(0,how_many_runs*(nback+1),max(dim(ADJ)))
        
        for (m in 1:max(dim(ADJ))){#setup initial listener histories with randomly selected neighbor labels
            
            contacts=as.list(which(ADJ[,m]!=0))
            contacts = unlist(unname(contacts))
            if (length(contacts)==0){#particularly in very small modules, we might have no connections... usually diag(ADJ) is set to all ones, which also solves this, but in case it is not, this avoids an error message
                IC_store[,m] <- m
                
                
            } else {
                
                contacts_values <- ADJ[contacts,m]
                
                ## Made an edit here to change the contacts_values  to 0 rahter than 1 ##
                
                if (length(which(contacts_values==1))!=length(contacts_values)){#this should generally be the case for weighted networks
                    contacts_weights <- contacts_values-min(contacts_values)
                    # Commented the below sections since it was popping up NaN values -- Pradeep
                    contacts_weights <- contacts_weights/max(contacts_weights)#rescale weights for to [0 1]... not totally necessary, but seems to help in some cases
                    #should note that this is just for the selection of IC's - the
                    #links with negative weights in ADJ should continue to have
                    #different (not just rescaled) behavior from positive links -
                    #practically speaking tests bear this out as nodes with negative
                    #correlations to some other group of nodes end up outside of that
                    #given cluster
                    
                } else {
                    contacts_weights <- matrix(1,length(contacts),1)
                    
                }
                
                bins <- c(0,cumsum(contacts_weights/sum(contacts_weights)),1)
                bins <- as.matrix(bins)
                bins <- as.matrix(apply(bins, 1, FUN = min))
                
                library(pracma)
                bins <- sort(bins, decreasing = FALSE)
                ret_hist <- histc(rand(how_many_runs*(nback+1),1),bins)
                idxs <- c()
                for (i in 1:length(ret_hist$bin)){
                    temp <- ret_hist$bin[i]
                    idxs <- c(idxs,contacts[temp])
                }
                IC_store[,m] <- idxs
                #using cellfun (@randsample or bsxfn isn't any faster
                
            }
            
        }
        return(IC_store)
    }
    
    
    genIC_sparse_unweighted <- function(ADJ,how_many_runs,nback,varargin){
        library(pracma)
        k <- colSums(ADJ)
        kcumsum <- cumsum(k)
        basicrand <- rand((nback+1)*how_many_runs,max(dim(ADJ)))
        indices <- which(ADJ != 0, arr.ind = TRUE)
        row_ind <- indices[,1]
        basicrand_scaled <- ceil(basicrand*drop(repmat(k,size(basicrand,1),1)))
        basicrand_scaled_lifted <- basicrand_scaled+repmat(c(0,kcumsum[-length(kcumsum)]),size(basicrand_scaled,1),1)
        IC_store <- Reshape(row(basicrand_scaled_lifted),(nback+1)*how_many_runs,max(dim(ADJ)))
        return(IC_store)
    }
    
    if (is_ADJ_weighted!=0){#use the appropriate routine to generate the initial conditions for the run (we can do this more quickly for sparse unweighted networks
        IC_store <- genIC_full(ADJ,iter,nback)
    } else {
        IC_store <- genIC_sparse_unweighted(ADJ,iter,nback)
    }
    
    partition_columns <- matrix(nrow=max(dim(IC_store)),ncol=iter)
    for (i in 1:iter){
        
        cat(paste0('iteration #',as.character(i), ' of ',str(iter),main_iter,'\n'))
        
        
        IC_store_relevant <- IC_store[1:(nback+1),]
        IC_store <- IC_store[(nback+1):nrow(IC_store),] #set in SSLPA
        speakeasy_res <- SpeakEasycore(ADJ,timesteps,IC_store_relevant,nback,force_efficient)
        listener_history = speakeasy_res$listener_history
        partition_columns[,i] = speakeasy_res$partitionID  #i.e.check sparse and suppress graphics
    }
    
    partitionID_lst <- list()
    for( i in 1:ncol(partition_columns)){
        partitionID_lst[[as.character(i)]] <- partition_columns[,i]
    }
    
    ### Take a look into with @Jake
    save_results <- 0
    if (main_iter==1){
        if (save_results==1){
            cat('saving partitions')
            saveRDS(list(ADJ,partition_columns,partitionID,multi_community_switch),'record_of_all_partitions.RDS')#sometimes useful for very large networks if you want to save results before consensus clustering, especially if you want to try several multi-community cutoffs
        }
    }
    
    cat(paste0('started consensus clustering : ', as.character(main_iter),'\n'))
    
    virt_res <- virtual_cooccurrence(ADJ,partition_columns,partitionID_lst, main_iter, multi_community_switch)
    
    boot_res <- list("partition_codes"= virt_res["nodes_and_partition_identifiers_hard"],
                     "partition_codes_overlapping" = virt_res["nodes_and_partition_identifiers_overlapping"],
                     "cell_partition" =virt_res["cell_partition"],
                     "cell_partition_overlapping" =virt_res["cell_partition_overlapping"])
    return(boot_res)
    
}


layer_SpeakEasy <- function(layers,iter,ADJ,timesteps,varargin=NULL){ ##ok<NCOMMA>
## Get ADJ properties to optimize runtime
# values in this section do not change the SpeakEasy method, just help to
# run it efficiently and determine how it is applied.  You likely don't need to
# adjust anything here, ever.

    nback <- 5#min recommended value is 5 - if you increase communities evolve more slowly... never seen a need to change this
    force_efficient <- 0#recommended to leave==0, which adaptively selects different forms of computation based on ADJ density, setting to 1 forces slower alternative
    max_ADJ_size <- 50000
    
    #set "max_ADJ_size" to the length of the largest ADJ your machine can hold in available
    #RAM 2x.  You can of course cluster sparse matrices that are larger than
    #this.  The reason for this value is that we use a different (~60# faster)
    #way of estimated expected label frequncy for relatively dense matrices if their dimension is less than this value.  The values computed are identical in either case; it's only a matter of efficiency.
    #On my laptop with 16GB, I can use a 10Kx10K full ADJ, so I'd set this value to 10000.
    #for my 128gb machine, I set "max_ADJ_size" to x with no trouble.
    #if you machine has 8GB of ram, I recommend max_ADJ_size=8000, as a comfortable value if you're not running a bunch of other stuff on the machine
    
    #deal with optional inputs for overlapping clustering and
    if (is.null(varargin)){
        multi_community_switch <- 1#default is disjoint clusters
    } else {
        if (length(varargin)>=1){
            multi_community_switch <- varargin[1]
            if (is.null(multi_community_switch)){
                multi_community_switch <- 1
            }
    
            if (multi_community_switch>1){
                cat('error - multi_community_switch should be in range of [0 1]')
                stop()
            }
    
            cat('optional input multi_community_switch set to ',str(multi_community_switch))
            if (multi_community_switch==1){
                cat('indicating disjoint output')
            }
        }
    }
    
    
    if (length(layers)==2){
        subclustersize <- layers[2]
        cat('optional input subclustersize set to ', str(subclustersize))
    } else {
        subclustersize <- 5#default...clusters smaller than this will be left alone (NOT subclustered)
    }
    
    noinputs=which(colSums(abs(ADJ))==0)#because if a node has no inputs it can't receive labels, we set it to receive it's own label, i.e. remain isolated
    ADJ[(noinputs-1)*(max(dim(ADJ)))+noinputs] <- 1
    
    ADJ_characteristics <- function(ADJ,max_ADJ_size, force_efficient){
        
        if (dim(ADJ)[1] != dim(ADJ)[2]){
            cat('your ADJ is not square, please fix \n')
            stop()
        }
        
        if ((max(ADJ))>1 | min(ADJ)<(-1)){
            cat('your connection strength is outside [-1 1] please fix \n')
            stop()
        }
        
        is_ADJ_weighted <- c()
        fraction_to_be_full <- 0.3#should be around .2 or less - if the matrix is more dense than this it will be treated as full
        
        sample_of_links = ADJ
        ADJ_density=length(which(sample_of_links!=0))/length(sample_of_links)
        cat('approximate edge density is ',str(ADJ_density))
        
        if (length(which(sample_of_links!=0))/length(sample_of_links)>fraction_to_be_full){
            ADJ_is_full <- 1
            
            if (length(which(sample_of_links>.9999))+length(which(abs(sample_of_links)<.0001))!=length(sample_of_links)){#in case there are rounding errors
                cat('weighted full')#ADJ with links of only +1 and -1 will still be considered weighted
                is_ADJ_weighted <- 1
            } else {
                cat('unweighted full')
                is_ADJ_weighted <- 0
            }
            
        } else {
            ADJ_is_full <- 0
            
            #in this case we can afford to test all links
            c <- ADJ[which(ADJ != 0)]
            if (length(which(c>.9999))!=length(c)){#in case there are rounding errors
                cat('weighted sparse \n')
                is_ADJ_weighted <- 1
            } else {
                cat('unweighted_sparse \n')
                is_ADJ_weighted <- 0
            }
            
        }
        
        
        #force_efficient--1 if the ADJ is huge or if it is small&dense
        if (force_efficient != 1){
            
            if (dim(ADJ)[1]<max_ADJ_size){#if matrix is small enough that we could use full multiplicationn
                
                if (ADJ_density>0.03){             #do not change this value; this value was selected by observing the two multiplication strategies it switches between require equal time at .03 density){
                    
                    force_efficient = 0}#ADJ dense enough that full mult is faster than pairwise (and small enough it will fit in memory)
                else {
                    
                    force_efficient = 1#'ultrasparse and small ADJ, so pairwise multi is faster'
                }
                
            } else {
                if (ADJ_density<=.03){#you still might want ADJ sparse if it's huge and RAM usage is an issue, but hoping the user necessarily does that
                    
                    force_efficient = 1}
                
                else {
                    force_efficient = 0}#if you have a huge matrix is still may be wise to put it as sparse even if density >.03.  .03 is the point at which it becomes faster to compute with a full, but you may need to convert to sparse simply for memory reasons
                
            }
            
        }
        values <- c(is_ADJ_weighted, force_efficient)
        return(values)
    }
    
    values <- ADJ_characteristics(ADJ, max_ADJ_size, force_efficient)#this is used because there are different routines for setting up the initial conditions for weighted vs unweighted networks
    is_ADJ_weighted <- values[1] 
    force_efficient <- values[2]
    
    ##  Call SpeakEasy on primary clusters and again on each cluster (if layers>1)
    cluster_stats_store <- list()
    partition_codes <- list()
    partition_codes_overlapping <- list()
    cell_partition <- list()
    cell_partition_overlapping <- list()
    convenient_node_ordering <- list()
    for (i in 1:layers){
        main_iter <- i#for clarity when passing
        if (i==1){
            cat('calling main routine at level 1 \n')
            boot_res <- bootstrap_SpeakEasy(iter,ADJ,timesteps,nback,is_ADJ_weighted, force_efficient, main_iter,multi_community_switch)
            partition_codes[[i]] <- boot_res["partition_codes"]
            partition_codes_overlapping[[i]] <- boot_res["partition_codes_overlapping"]
            cell_partition[[i]] <- boot_res["cell_partition"]
            cell_partition_overlapping[[i]] <- boot_res["cell_partition_overlapping"]
    
        } else {
            cat('doing subclustering at level ', str(i))
            cat('terminal output reduced', str(main_iter-1))
            partition_codes_temp <- list()
            partition_codes_overlapping_temp <- list()
            cell_partition <- list()
            cell_partition_overlapping_temp <- list()
            for (j in 1:length(cell_partition_overlapping[[i-1]])){
                current_nodes <- cell_partition_overlapping[[i-1]][j] #needed becase first node in sub cluster may be 99th overal etc
                if (length(current_nodes)>subclustersize){  #don't even try to subcluster communites that are smaller than a certain size
                    boot_res2 <- bootstrap_SpeakEasy(iter,ADJ[c(current_nodes),c(current_nodes)],timesteps,nback,is_ADJ_weighted,force_efficient,main_iter,multi_community_switch)
                    partition_codes_temp[[j]] <- boot_res2["partition_codes"]
                    partition_codes_overlapping_temp[[j]]  <- boot_res2["partition_codes_overlapping"]
                    cell_partition[[j]] <- boot_res2["cell_partition"]
                    cell_partition_overlapping_temp[[j]]  <- boot_res2["cell_partition_overlapping"]
    
                    partition_codes_overlapping_temp[[j]][,1] <- current_nodes[partition_codes_overlapping_temp[[j]][,1]]
                    for (k in 1:length(cell_partition_overlapping_temp[[j]])){ #update to main node ID's
                        cell_partition_overlapping_temp[[j]][k] <- current_nodes[cell_partition_overlapping_temp[[j]][k]]
                    }
                } else {#if module is too small for subclustering
                    cell_partition_overlapping_temp[[j]] <- c(current_nodes)
                }
            }
            
            for (k in 1:length(cell_partition_overlapping_temp)){
                if(k == 1){
                    cell_partition_overlapping_temp_db <- as.data.frame(cell_partition_overlapping_temp[[1]])
                }else {
                    cell_partition_overlapping_temp_db <- rbind(cell_partition_overlapping_temp_db, as.data.frame(cell_partition_overlapping_temp[[k]]))
                }
            }
            cell_partition_overlapping[[i]]<- cell_partition_overlapping_temp_db
            cell_partition_overlapping_temp <- list()
    
            temp <- list()
            for (m in 1:length(cell_partition_overlapping[[i]])){
                
                temp[[m]] <- c(cell_partition_overlapping[[i]][m], kronecker(matrix(length(cell_partition_overlapping[[i]][m]),1),m))
            }
            for (k in 1:length(temp)){
                if(k == 1){
                    temp_db <- as.data.frame(temp[[1]])
                }else {
                    temp_db <- rbind(temp_db, as.data.frame(temp[[k]]))
                }
            }
            partition_codes_overlapping[[i]] <- temp_db
    
        }
        
        #Since its always going to be 1 layer
        for (k in 1:length(cell_partition_overlapping)){
            if(k == 1){
                print(length(cell_partition_overlapping[[1]]))
                cell_partition_db_3 <- data.frame("cell_part" =cell_partition_overlapping[[1]])
            }else {
                cell_partition_db_3 <- rbind(cell_partition_db_3, data.frame("cell_part" = cell_partition_overlapping[[k]]))
            }
        }
        
        convenient_node_ordering<- cell_partition_overlapping
        
        res_layer <- c("partition_codes_overlapping"=partition_codes_overlapping,
                       "cell_partition_overlapping"=cell_partition_overlapping,
                       "convenient_node_ordering"=convenient_node_ordering)
        return(res_layer)
    }
}

  result_temp <- layer_SpeakEasy(layers=1,iter=iter,ADJ=adj,timesteps=timesteps)

  gene_Modules <- result_temp$partition_codes_overlapping$partition_codes_overlapping$nodes_and_partition_identifiers_overlapping
  gene_Modules <- as.data.frame(gene_Modules)
  colnames(gene_Modules) <- c('Gene.ID','moduleNumber')
  if(!(is.null(GeneNames))){
    gene_Modules[Gene.ID] = GeneNames
  }
  gene_Modules['moduleSize'] <- 0
  count_mtx <- table(gene_Modules$moduleNumber)

  for (i in 1:nrow(gene_Modules)){
    temp = gene_Modules$moduleNumber[i]
    gene_Modules$moduleSize[i] = as.integer(count_mtx[temp])
  }

  gene_Modules = gene_Modules %>%
    group_by(Gene.ID) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                  moduleNumber = as.numeric(moduleNumber)) %>% filter(moduleSize > 30)


  # Change cluster number to color labels
  gene_Modules$moduleNumber = as.numeric(factor(gene_Modules$moduleNumber))
  gene_Modules$moduleLabel = WGCNA::labels2colors(gene_Modules$moduleNumber)
  mod = gene_Modules[c('Gene.ID','moduleNumber','moduleLabel')]

  return(mod)
}