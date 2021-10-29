####Generalized infrastructure for multiple family comparisons
#### This code was originally written for Jackson et al. Science, 2019
#### Only slightly adapted for Avery et al.
#### Most of the R libraries are easily available, except for the clues library, which isn't strictly needed for our analysis

library(igraph)
library(clues)
library(Rcpp)
library(RcppArmadillo)
library(gtools)
library(aricode)
sourceCpp("RandomWalkCode.cpp")

#Read in CLICS3 full network
graph<-as.undirected(read.graph("./network-1-families.gml",format="gml"))

#Remove Merry-Sad colex
#Merge Happy/Happiness
nodes = 1:2821
nodes[1916] = 2520
graph <- delete.edges(graph,E(graph)[get.edge.ids(graph, c(1667, 2051))])
graph <- contract(graph, mapping = nodes, vertex.attr.comb = toString)
graph = set_vertex_attr(graph, "Gloss",V(graph)[2520], "HAPPY")
graph = delete.vertices(graph, V(graph)[1916])

###This function performs all analyses on the colex network, only for the newest CLICS database.
#    walk_steps is the number of steps the random walk takes
#    concept_comm_threshold is the largest number of total unconnected components a concept is in.
#    for example, if for 13 concepts, there are 13 connected components, that means each concept will by definition
#    be an isolate, as each concept is in a separate component of the graph.
#   
#    WARNING: permute = T will perform strength and edge rewiring permutations on all language families. This takes 3+ days to complete with 1000 iterations

multi_family_comparison <- function(graph, concept_table, output_dir = "./Multi_fam_output",
                                    language_file = "languages_v2.csv", walk_steps = 5,
                                    concept_comm_threshold = 7, use_comm_threshold = F,
                                    permute = T,  permuteIter = 10, savePermute = T){
  require(igraph)
  require(clues)
  require(expm)
  d<-vertex.attributes(graph)
  d<-as.data.frame(d)
  dir.create(output_dir)
  
  e <- as.data.frame(edge.attributes(graph))
  
  
  targetNodes <- d$id[which(d$Gloss %in% concept_table[-1,2])]
  
  print(paste("Found ", length(targetNodes), " Concept Nodes"))
  
  targetIndices <<- which(d$Gloss %in% concept_table[-1,2])
  emoLabels = levels(d$Gloss)[as.numeric(d$Gloss[which(d$Gloss %in% concept_table[-1,2])])]
  
  
  ###Language Group Analysis Script
  langFams <- read.csv(language_file, stringsAsFactors = F)
  langMap = langFams[,c(1, 9)]
  langMap[,1] = paste(langFams$dataset_ID, langFams$ID, sep = "-")
  langs <-strsplit(x = levels(e$families), ";", fixed = T)
  
  langs <- do.call("c", langs)
  
  families <- levels(as.factor(langs))
  
  netFams <- list()
  
  compFunc <- function(vec,fam){
    targetLanguages <- langMap[which(langMap[,2] == fam),1]
    return(length(which(vec %in% targetLanguages)))
    
    
  }
  
  
  targetFams = names(table(langMap[,2])[table(langMap[,2])>2])
  
  netFams = list()
  counter = 1
  LangsCount = c()
  
  langs_by_fam = list()
  for(i in 1:length(targetFams)){
   
    netFams[[i]] <- subgraph.edges(graph, eids = which(grepl(pattern = targetFams[i],E(graph)$families, fixed = T)), delete.vertices = F)
    subLangs <- strsplit(E(netFams[[i]])$languages, ";", fixed = T)
    unique_langs = unique(do.call("c", subLangs))
    family_languages = length(langMap[which(langMap[,2] == targetFams[i]),1])
    used_languages = length(na.omit(match(langMap[which(langMap[,2] == targetFams[i]),1], unique_langs)))
    LangsCount = c(LangsCount,subLangs)
    E(netFams[[i]])$famLangWeight <-sapply(X = subLangs, FUN = compFunc, fam = targetFams[i])
    counter = counter+1
    langs_by_fam[[i]] = data.frame(family_languages, used_languages)
    
  }
  names(netFams) <- targetFams
  langs_by_fam <- do.call("rbind", langs_by_fam)
  
  print(paste("Found the following languages:", targetFams))
  
  #Compute information about Family Networks
  
  net_stat_func <- function(network){
    
    density = edge_density(network)
    n_con_components = components(network)
    closeness_centrality= centr_clo(network)
    concept_comms = length(unique(n_con_components$membership[targetIndices]))
    max_concepts = max(table(n_con_components$membership[targetIndices]))
    max_dist = distances(network, targetIndices, targetIndices, weights = NA)
    max_dist = unlist(max_dist[upper.tri(max_dist)])
    max_dist = max(max_dist[!is.infinite(max_dist)])
   # print(data.frame(density, n_con_components$no,concept_comms,max_concepts,max_dist, closeness_centrality$centralization))
    return(data.frame(density, n_con_components$no,concept_comms,max_concepts,max_dist, closeness_centrality$centralization))  
  }
  
  net_stats = lapply(netFams, net_stat_func)
  
  net_stats = do.call("rbind", net_stats)
  net_stats = cbind(langs_by_fam, net_stats)
  net_stats$family = names(netFams)
  
  LangsCount = do.call("c", LangsCount)
  if(use_comm_threshold){
  valid_fams = which(net_stats$concept_comms <= concept_comm_threshold & !(net_stats$family %in% c("", "Bookkeeping", "Kolopom", "Lakes Plain")))
  net_stats = net_stats[valid_fams,]
  }else{
    valid_fams = which(!(net_stats$family %in% c("", "Bookkeeping", "Kolopom", "Lakes Plain")))
    net_stats = net_stats[valid_fams,]
  }
  
  #net_stats = net_stats[-which(net_stats$family %in% c("Bookkeeping", "", "Kolopom")),]
  
  commExtract <- function(graph, target_attr = "famLangWeight"){
    require(expm)
    fullMat  <- as_adj(graph, attr = target_attr)
    # randWalkDist <- function(matrix, steps){
    #   probs <- list()
    #   sums <- list()
    #   probs[[1]] <- matrix
    #   sums[[1]] <- matrix
    #   for(i in 2:steps){
    #     print(i)
    #     probs[[i]] <- matrix%^%i
    #     sums[[i]] = probs[[i]]*(1-sums[[i-1]]) + sums[[i-1]]
    #   }
    #   return(sums)
    # }
    
    test =prob_walk_net(fullMat, walk_steps) 
    
    EmoDist <- test[targetIndices,targetIndices]
    
    EmoDist = (EmoDist +t(EmoDist))/2
    EmoDistGraph <- graph.adjacency(round(EmoDist*100,digits = 0), mode = "undirected", diag = F, weighted = T)
    comms <-cluster_optimal(EmoDistGraph, weights = E(EmoDistGraph)$weight)
    mod = modularity(comms)
    if(is.nan(mod)){
      mod = 0
    }
    return(list(EmoDist,membership(comms), mod))
  }
  
  universal_fam = commExtract(graph, "LanguageWeight")
  comms = universal_fam[[2]]
  net = universal_fam[[1]]
  net_mat = net + t(net)
  net = graph.adjacency(net_mat, "undirected", weighted = T, diag = F)
  V(net)$"emonames" = emoLabels
  V(net)$"comms" = comms
  V(net)$"strength" = rowSums(net_mat)
  write.graph(net, paste(output_dir, "/UniversalNet.gml",sep  =""), "gml")
  
  familyWiseComms <- list()
  print("Extracting Communities")
 
  netFams = netFams[valid_fams]
  net_stats$Mod = 0
  global_netfams <<- netFams
  for(i in 1:length(netFams)){
    print(names(netFams)[i])
    familyWiseComms[[i]] <- commExtract(netFams[[i]])
    net_stats[i,"Mod"] = familyWiseComms[[i]][[3]]
  }
  
  

  
  names(familyWiseComms) <- names(netFams)
  
  familyWiseComms2 <- familyWiseComms[which(sapply(X = familyWiseComms,FUN = function(x){length(table(x[[2]]))}) != nrow(concept_table))]
  net_stats <- net_stats[which(sapply(X = familyWiseComms,FUN = function(x){length(table(x[[2]]))}) != nrow(concept_table)),]
  
  
  EmoDistStatList = list()
  EmoDistStatList2 = list()
  emo_net_stats = function(fam, name){
    print(name)
    require(brainGraph)
    library(qgraph)
    temp_net <<- graph.adjacency(fam[[1]], mode = "undirected", diag = F, weighted = T)
    density = edge_density(temp_net)
    ave_degree = mean(degree(temp_net))
    ave_strength = mean(strength(temp_net), weights = E(temp_net)$weight)
    global_eff = efficiency(g = temp_net, type = "global")
    clust_coef = clustOnnela(x = temp_net)[,1]
    clust_coef[is.nan(clust_coef)] = 0
    clust_coef = mean(clust_coef)
    modularity = fam[[3]]
    name = name
    part1 = data.frame(name, density, ave_degree, ave_strength, global_eff, clust_coef, modularity)
    
    
    node_names = d$Gloss[targetIndices]
    degree = degree(temp_net)
    eg_cent = eigen_centrality(temp_net, weights = E(temp_net)$weight)$vector
    print(name)
   
    print(name)
    betw = betweenness(temp_net, weights = E(temp_net)$weights)
    part2 = data.frame(fam = name, node = node_names, degree =degree, eg_cent = eg_cent,  betw = betw)
    
    return(list(part1, part2))
  }
  
  for(i in 1:length(familyWiseComms2)){
    
    temp = emo_net_stats(familyWiseComms2[[i]], names(familyWiseComms2)[[i]])
    EmoDistStatList[[i]] = temp[[1]]
    EmoDistStatList2[[i]] = temp[[2]]
  }
  
  uniTemp =emo_net_stats(universal_fam, "Universal")
  
  EmoDistStatList[[length(EmoDistStatList) +1]] =uniTemp[[1]]
  EmoDistStatList2[[length(EmoDistStatList2) +1]] =uniTemp[[2]]
  EmoDistStatList = do.call("rbind", EmoDistStatList)
  EmoDistStatList2 = do.call("rbind", EmoDistStatList2)
  write.csv(EmoDistStatList,paste(output_dir, "/family_concept_network_statistics.csv", sep = ""))
  write.csv(EmoDistStatList2,paste(output_dir, "/family_node_concept_network_statistics.csv", sep = ""))
  glob_famcoms <<- familyWiseComms2
  write.csv(net_stats, paste(output_dir, "/network_statistics.csv", sep = ""))
  
  HARandMat <- matrix(NA, nrow = length(familyWiseComms2), ncol = length(familyWiseComms2))
  library(clues)
  for(i in 1:length(familyWiseComms2)){
    for(j in 1:length(familyWiseComms2)){
      intersect1 = names(table(familyWiseComms2[[i]][[2]]))[which(table(familyWiseComms2[[i]][[2]])>1)]
      nodes1 = which(familyWiseComms2[[i]][[2]] %in% intersect1)
      
      intersect2 = names(table(familyWiseComms2[[j]][[2]]))[which(table(familyWiseComms2[[j]][[2]])>1)]
      nodes2 = which(familyWiseComms2[[j]][[2]] %in% intersect2)
      intersect= intersect(nodes1, nodes2)
      intersect = 1:length(targetIndices)
      HARandMat[i,j] <- adjustedRand(as.numeric(familyWiseComms2[[i]][[2]][intersect]),
                                     as.numeric(familyWiseComms2[[j]][[2]][intersect]), "HA")     
    }
  }
  
  colnames(HARandMat) =names(familyWiseComms2)
  rownames(HARandMat) =names(familyWiseComms2)
  
  write.csv(HARandMat, paste(output_dir,"/LanguageComparisonARI.csv", sep = ""))
  HARandMat[which(HARandMat <0)] <- 0
  HARandMat[which(is.nan(HARandMat))] = 0
  HARandMatNet <- graph.adjacency(adjmatrix = HARandMat, weighted = T, mode = "undirected", diag = F)
  
  
  
  HARandMat <- matrix(NA, nrow = length(familyWiseComms2), ncol = length(familyWiseComms2))
  library(clues)
  for(i in 1:length(familyWiseComms2)){
    for(j in 1:length(familyWiseComms2)){
      intersect1 = names(table(familyWiseComms2[[i]][[2]]))[which(table(familyWiseComms2[[i]][[2]])>1)]
      nodes1 = which(familyWiseComms2[[i]][[2]] %in% intersect1)
      
      intersect2 = names(table(familyWiseComms2[[j]][[2]]))[which(table(familyWiseComms2[[j]][[2]])>1)]
      nodes2 = which(familyWiseComms2[[j]][[2]] %in% intersect2)
      intersect= intersect(nodes1, nodes2)
      intersect = 1:length(targetIndices)
      HARandMat[i,j] <- NMI(as.numeric(familyWiseComms2[[i]][[2]][intersect]),
                                     as.numeric(familyWiseComms2[[j]][[2]][intersect]))     
    }
  }
  
  colnames(HARandMat) =names(familyWiseComms2)
  rownames(HARandMat) =names(familyWiseComms2)
  
  write.csv(HARandMat, paste(output_dir,"/LanguageComparisonNMI.csv", sep = ""))
  
  
  HARandUni = vector()
  for(i in 1:length(familyWiseComms2)){
    HARandUni[i] = adjustedRand(as.numeric(familyWiseComms2[[i]][[2]]), as.numeric(comms), "HA")
  }
  
  
  V(HARandMatNet)$"Names" = names(familyWiseComms2)
  V(HARandMatNet)$"HARandUni" = HARandUni
  famComms <- walktrap.community(HARandMatNet, weights = E(HARandMatNet)$weight)
  V(HARandMatNet)$"comm" = membership(famComms)
  
  write.graph(HARandMatNet,paste(output_dir,"/FamilySimNet.gml", sep = ""), "gml")
  
  
  
  emotionCategories <- data.frame(emoLabels)

  GML_Net_Generator = function(name){
    net = familyWiseComms2[[name]][[1]]
    net_mat = net + t(net)
    net = graph.adjacency(net_mat, "undirected", weighted = T, diag = F)
    V(net)$"emonames" = emoLabels
    V(net)$"comms" = familyWiseComms2[[name]][[2]]
    V(net)$"strength" = rowSums(net_mat)
    write.graph(net, paste(output_dir, "/",name,".gml",sep  =""), "gml")
    emotionCategories[,name] <<- V(net)$comms
  }
  
  names(familyWiseComms2)
  global_comms <<- familyWiseComms2
  for(i in 1:length(familyWiseComms2)){
    GML_Net_Generator(names(familyWiseComms2)[i])  
  }
  write.csv(emotionCategories,paste(output_dir,"/AllFamConceptCategories.csv", sep = ""))
  
  print("Completed Data Analysis")
  
  #
  
  strength_rewire = function(graph, target_attr="famLangWeight"){
    
    out_strength = strength(graph,mode = "out", weights = get.edge.attribute(graph, target_attr))
    in_strength = strength(graph,mode = "out", weights = get.edge.attribute(graph, target_attr))
    
    out_deg = degree(graph, mode = "out")
    in_deg = degree(graph, mode = "in")
    in_out_seq = cbind(out_deg, in_deg)
    in_out_seq_u =unique(in_out_seq)
    
    for(i in 1:nrow(in_out_seq_u)){
      #print(in_out_seq_u[i,])
      k_ind = which((in_out_seq[,1] == in_out_seq_u[i,1]) & (in_out_seq[,2] == in_out_seq_u[i,2]))
      #k_ind = which((in_out_seq[,1] == 100) & (in_out_seq[,2] == 100))
      #print(k_ind)
      if(length(k_ind) >1){
        new_ind = gtools::permute(k_ind)
      }else{
        new_ind = k_ind
      }
      #print(new_ind)
      
      out_strength[new_ind]
      out_strength[k_ind] = out_strength[new_ind]
      
      if(any(out_strength[k_ind] == 0 &(in_out_seq_u[i, 1] !=0) )){
        
        print(in_out_seq[i, 1])
        break
      }
      
      in_strength[k_ind] = in_strength[new_ind]
      
    }
    
    mean_k = mean(in_out_seq)
    mean_s = mean(cbind(out_strength, in_strength))
    
    
    return(list(cbind(out_strength, in_strength), mean_k, mean_s, cbind(out_deg, in_deg)))
  }
  
  
  permute_graph_comm <- function(graph, iters = permuteIter, target_attr ="famLangWeight", tvec ){
    
    per_list = vector()
    for(i in 1:iters){
      print(i)
      strength_rewired = strength_rewire(graph, target_attr = target_attr)
      temp = rewire(graph, with = keeping_degseq(niter = vcount(graph)*500))
      temp_adj = as_adj(temp)
      temp = strength_preserve_rewire(temp_adj, strength_rewired[[1]], strength_rewired[[2]], strength_rewired[[3]], strength_rewired[[4]])
      temp_net = graph.adjacency(temp, mode = "undirected", diag = F, weighted = T)
      comms <- commExtract(temp_net, target_attr = "weight")[[1]]
      temp = graph.adjacency(round(comms*100,0),"undirected", diag = F, weighted = T)
      
      per_list[i] = modularity(temp, tvec, weights = E(temp)$weight)
      print(per_list[i])
    }
    return(per_list)
    
  }
  
  
  if(permute){
    mod_perm = vector()
    for(i in 1:length(familyWiseComms2)){
    print(names(netFams)[i])
    perm_list = permute_graph_comm(netFams[[i]], tvec = familyWiseComms2[[i]][[2]])
    
    distrib = ecdf(perm_list)
    mod_perm[i] = 1 - distrib(familyWiseComms2[[i]][[3]])
    print(perm_list)
    print(familyWiseComms2[[i]][[3]])
    print(mod_perm[i])
    net_stats$ModP[i] = mod_perm[i]
    write.csv(net_stats, paste(output_dir, "/network_statistics.csv", sep = ""))
    }
    
    pval_mat = matrix(0,nrow = length(familyWiseComms2),ncol = length(familyWiseComms2))
    for(i in 1:length(familyWiseComms2)){
      for(j in i:length(familyWiseComms2)){
        print(j)
        true_ARI = ARI(as.numeric(familyWiseComms2[[i]][[2]]), as.numeric(familyWiseComms2[[j]][[2]]))
        perm_ARI = vector()
        for(it in 1:permuteIter){
          perm_1 = permute(familyWiseComms2[[i]][[2]])
          perm_2 = permute(familyWiseComms2[[j]][[2]])
          
          perm_ARI[it] = ARI(as.numeric(perm_1), as.numeric(perm_2))+rnorm(1,0,.0000001)
          if(is.nan(perm_ARI[it])){
            perm_ARI[it] = 0+rnorm(1,0,.0000001)
          }
          
          
          
        }
        #print(true_ARI)
        #print(perm_ARI)
        perm_ecdf = ecdf(perm_ARI)
        pval_mat[i,j] =1 - perm_ecdf(true_ARI)
        
        
      }
    }
    
    rownames(pval_mat) = names(familyWiseComms2)
    colnames(pval_mat) = names(familyWiseComms2)
    
    write.csv(pval_mat, paste(output_dir, "/ARI_pvals.csv", sep = ""))
    
    }
  
  write.csv(net_stats, paste(output_dir, "/network_statistics.csv", sep = ""))
  
    ##L1 Matrix Norm calculations
    L1norm = matrix(0,length(familyWiseComms2), length(familyWiseComms2))
    for(i in 1:length(familyWiseComms2)){
      for(j in i:length(familyWiseComms2)){
        net = familyWiseComms2[[i]][[1]]
        net_mat1 = net + t(net)
        net = familyWiseComms2[[j]][[1]]
        net_mat2 = net + t(net)
        L1norm[i,j] = sum(abs(net_mat1-net_mat2))
      }
    }
    
    ##Additional Network Statistics
    av_str = vector()
    w_clust = vector()
    eff = vector()
    require(brainGraph)
    require(qgraph)
    for(i in 1:length(familyWiseComms2)){
      net = familyWiseComms2[[i]][[1]]
      net_mat = net + t(net)
      net = graph.adjacency(net_mat, "undirected", weighted = T, diag = F)
      eff[i] =  efficiency(net, "global",weights = E(net)$weights)
      w_clust[i] = mean(na.omit(clustOnnela(as.matrix(net_mat))[,1]))
      av_str[i] = mean(net_mat[upper.tri(net_mat)])
      
    }
    w_clust[which(is.nan(w_clust))] = 0
    net_stats$av_str = av_str
    net_stats$w_clust = w_clust
    net_stats$glob_eff = eff
    write.csv(net_stats, paste(output_dir, "/network_statistics.csv", sep = ""))
    
    
    ##Pairwise clustering coeff diff
    clust_diff = matrix(0,length(familyWiseComms2), length(familyWiseComms2))
    for(i in 1:length(familyWiseComms2)){
      for(j in i:length(familyWiseComms2)){
        clust_diff[i,j] = abs(w_clust[i]-w_clust[j])
      }
    }
    
    ##Pairwise average dist 
    dist_diff = matrix(0,length(familyWiseComms2), length(familyWiseComms2))
    for(i in 1:length(familyWiseComms2)){
      for(j in i:length(familyWiseComms2)){
        net1 = as.matrix(familyWiseComms2[[i]][[1]])
        net_mat1 = net1 + t(net1)
        net2 = as.matrix(familyWiseComms2[[j]][[1]])
        net_mat2 = net2 + t(net2)
        
        #Take the inverse
        z_loc = which(net_mat1 == 0, arr.ind = T)
        net_mat1 = 1/net_mat1
        net_mat1[z_loc] = 0
        
        z_loc = which(net_mat2 == 0, arr.ind = T)
        net_mat1 = 2/net_mat2
        net_mat1[z_loc] = 0
        
        net1 = graph.adjacency(net_mat1, "undirected", weighted = T, diag = F)
        net2 = graph.adjacency(net_mat2, "undirected", weighted = T, diag = F)
        
        dist1 = distances(net1, weights = E(net1)$weight)
        dist2 = distances(net2, weights = E(net2)$weight)
        dist1[which(is.infinite(dist1))] =0
        dist2[which(is.infinite(dist2))] =0
        
        dist_diff[i,j] = mean(abs(dist1-dist2))
      }
    }
    
    rownames(dist_diff) = names(familyWiseComms2)
    colnames(dist_diff) = names(familyWiseComms2)
    write.csv(dist_diff, paste(output_dir, "/DistDiff.csv", sep = ""))
    
    rownames(clust_diff) = names(familyWiseComms2)
    colnames(clust_diff) = names(familyWiseComms2)
    write.csv(clust_diff, paste(output_dir, "/ClustDiff.csv", sep = ""))
    
    rownames(L1norm) = names(familyWiseComms2)
    colnames(L1norm) = names(familyWiseComms2)
    write.csv(L1norm, paste(output_dir, "/L1norm.csv", sep = ""))
    
    
    }

#Perform Combined Emotion + Taste + Color Analysis
concept_table = read.table("concept_lists/ETC_concepts.txt", sep = "\t")
multi_family_comparison(graph, concept_table,walk_steps = 5, output_dir = "graph_files/ETC_concepts", use_comm_threshold = T,concept_comm_threshold = 49,  permute = F, permuteIter = 10)
