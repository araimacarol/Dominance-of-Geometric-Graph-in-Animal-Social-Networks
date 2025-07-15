# filezilla and cmd
#(1)  C:/Users/rcappaw/OneDrive - University of Tasmania/Desktop
#(2) ssh -i RaimaAppaw.pem raima@131.217.173.254
#(3) raima@crunchie:~/data/simulations/igraphEpiNew$
#(4) screen -S name
#(5) ctrl+A and D to detach
#(6) screen -r name to reattach

#setwd("C:/Users/rcappaw/OneDrive - University of Tasmania/Desktop/R/Workflow/igraphEpi-New")

library(hrbrthemes)
library(gdtools)
library(progress)
library(readr)
library(stringr)
library(tidyverse)
library(tidymodels)
library(vip)
library(flashlight)
library(doParallel)
library(janitor)
library(viridis)
library(GGally)
library(ggfortify)
library(future.apply)
library(parallel)
library(hrbrthemes)
library(iml)
library(plyr)
library(EpiModel)
library(igraph)
# install.packages("hrbrthemes")
# install.packages("gdtools")




generate_graphs <- function(network_size, num_replicates, EDGE_COUNT, Degseqsize, NEIGHBORS,
                            REWIRING_PROB, BARABASI_M, SPATIAL_R, SIZES_2_BLOCKS,
                            PROBABILITY_MATRIX_2_BLOCKS, SIZES_3_BLOCKS,
                            PROBABILITY_MATRIX_3_BLOCKS) {
  # ... (previous code remains the same)
  unweighted_graphs <- list()
  weighted_graphs <- list()
  unweighted_edgelists <- list()
  weighted_edgelists <- list()
  unweighted_adjacency_matrices <- list()
  weighted_adjacency_matrices <- list()
  unweighted_laplacian_matrices <- list()
  weighted_laplacian_matrices <- list()
  
  for (i in 1:num_replicates) {
    # Generate unweighted graphs
    repeat {
      out.deg <- sample(1:Degseqsize, network_size, replace = TRUE)
      if (sum(out.deg) %% 2 == 0) break
    }
    
    # Degree Sequence Graph
    ds_graph <- igraph::sample_degseq(out.deg, method = "simple")
    ds_graph <- igraph::simplify(ds_graph, remove.multiple = TRUE, remove.loops = TRUE)
    unweighted_graphs[[paste0("DS_graph_n", network_size, "_", i)]] <- ds_graph
    
    # Erdos-Renyi Graph
    p <- (2 * EDGE_COUNT) / (network_size * (network_size - 1))
    er_graph <- igraph::erdos.renyi.game(network_size, p, type = "gnp", directed = FALSE, loops = FALSE)
    unweighted_graphs[[paste0("ER_graph_n", network_size, "_", i)]] <- er_graph
    
    # Watts-Strogatz Graph
    sw_graph <- igraph::sample_smallworld(dim = 1, size = network_size, nei = NEIGHBORS, p = REWIRING_PROB, loops = FALSE, multiple = FALSE)
    unweighted_graphs[[paste0("SW_graph_n", network_size, "_", i)]] <- sw_graph
    
    # Barabasi-Albert Graph
    sf_graph <- igraph::sample_pa(n = network_size, m = BARABASI_M, directed = FALSE)
    unweighted_graphs[[paste0("SF_graph_n", network_size, "_", i)]] <- sf_graph
    
    # Spatial Graph
    sp_graph <- fastSpatialNetwork(n = network_size, r = SPATIAL_R, makeConnected = TRUE, keepCellsSeparate = FALSE)
    unweighted_graphs[[paste0("SP_graph_n", network_size, "_", i)]] <- sp_graph
    
    # Stochastic Block Model with 2 blocks
    sbm2_graph <- igraph::sample_sbm(n = sum(SIZES_2_BLOCKS), pref = PROBABILITY_MATRIX_2_BLOCKS, block.sizes = SIZES_2_BLOCKS, directed = FALSE, loops = FALSE)
    unweighted_graphs[[paste0("SBM2_graph_n", network_size, "_", i)]] <- sbm2_graph
    
    # Stochastic Block Model with 3 blocks
    sbm3_graph <- igraph::sample_sbm(n = sum(SIZES_3_BLOCKS), pref = PROBABILITY_MATRIX_3_BLOCKS, block.sizes = SIZES_3_BLOCKS, directed = FALSE, loops = FALSE)
    unweighted_graphs[[paste0("SBM3_graph_n", network_size, "_", i)]] <- sbm3_graph
    
    # Generate weighted graphs
    for (graph_name in names(unweighted_graphs)) {
      graph <- unweighted_graphs[[graph_name]]
      
      # Add weights based on topological features
      if (startsWith(graph_name, "DS_")) {
        weights <- sapply(E(graph), function(e) {
          igraph::degree(graph, V(graph)[ends(graph, e)[1]]) * 
            igraph::degree(graph, V(graph)[ends(graph, e)[2]])
        })
      } else if (startsWith(graph_name, "ER_")) {
        weights <- runif(ecount(graph), 1, 10)  # Random weights between 1 and 10
      } else if (startsWith(graph_name, "SW_")) {
        edge_list <- as_edgelist(graph, names = FALSE)
        num_edges <- nrow(edge_list)
        
        weights <- numeric(num_edges)
        for (i in 1:num_edges) {
          v1 <- edge_list[i, 1]
          v2 <- edge_list[i, 2]
          weights[i] <- 1 / (1 + abs(v1 - v2))  # Weight based on "distance" in the original ring
        }
        
        # Check for any remaining NAs and replace them
        na_indices <- which(is.na(weights))
        if (length(na_indices) > 0) {
          weights[na_indices] <- mean(weights, na.rm = TRUE)
        }
        
        # Ensure no zero weights
        weights[weights == 0] <- min(weights[weights > 0])
      } else if (startsWith(graph_name, "SF_")) {
        weights <- sapply(E(graph), function(e) {
          (igraph::degree(graph, V(graph)[ends(graph, e)[1]]) + 
             igraph::degree(graph, V(graph)[ends(graph, e)[2]]))
        })
      } else if (startsWith(graph_name, "SP_")) {
        coords <- layout_with_fr(graph)
        weights <- sapply(E(graph), function(e) {
          v1 <- ends(graph, e)[1]
          v2 <- ends(graph, e)[2]
          dist <- sqrt(sum((coords[v1,] - coords[v2,])^2))
          1 / (1 + dist)
        })
      } else if (startsWith(graph_name, "SBM")) {
        membership <- V(graph)$block
        if (is.null(membership)) {
          weights <- runif(ecount(graph), 1, 10)  # Random weights between 1 and 10
        } else {
          weights <- sapply(E(graph), function(e) {
            v1 <- ends(graph, e)[1]
            v2 <- ends(graph, e)[2]
            if (membership[v1] == membership[v2]) 10 else 1  # Stronger within-block connections
          })
        }
      }
      
      # Scale weights to 1-10 range
      if (length(weights) > 0) {
        weights <- 1 + 9 * (weights - min(weights)) / (max(weights) - min(weights))
      } else {
        weights <- numeric(0)
      }
      
      # Create weighted graph
      weighted_graph <- graph_from_edgelist(as_edgelist(graph), directed = FALSE)
      E(weighted_graph)$weight <- weights
      weighted_graphs[[graph_name]] <- weighted_graph
      
      # Create edgelists
      unweighted_edgelist <- as_edgelist(graph)
      weighted_edgelist <- cbind(as_edgelist(weighted_graph), weights)
      unweighted_edgelists[[graph_name]] <- unweighted_edgelist
      weighted_edgelists[[graph_name]] <- weighted_edgelist
      
      # Create adjacency and Laplacian matrices
      unweighted_adj_matrix <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
      weighted_adj_matrix <- igraph::as_adjacency_matrix(weighted_graph, attr = "weight", sparse = FALSE)
      
      unweighted_laplacian_matrix <- igraph::laplacian_matrix(graph, normalized = FALSE)
      weighted_laplacian_matrix <- igraph::laplacian_matrix(weighted_graph, normalized = FALSE)
      
      unweighted_adjacency_matrices[[graph_name]] <- unweighted_adj_matrix
      weighted_adjacency_matrices[[graph_name]] <- weighted_adj_matrix
      unweighted_laplacian_matrices[[graph_name]] <- unweighted_laplacian_matrix
      weighted_laplacian_matrices[[graph_name]] <- weighted_laplacian_matrix
    }
  }
  
  return(list(
    unweighted_graphs = unweighted_graphs,
    weighted_graphs = weighted_graphs,
    unweighted_edgelists = unweighted_edgelists,
    weighted_edgelists = weighted_edgelists,
    unweighted_adjacency_matrices = unweighted_adjacency_matrices,
    weighted_adjacency_matrices = weighted_adjacency_matrices,
    unweighted_laplacian_matrices = unweighted_laplacian_matrices,
    weighted_laplacian_matrices = weighted_laplacian_matrices
  ))
}


create_directories <- function(base_dir, network_sizes, network_types) {
  for (size in network_sizes) {
    # Create base directory for each network size
    size_dir <- file.path(base_dir, paste0("size_", size))
    dir.create(size_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Create separate directories for unweighted and weighted networks
    for (weight_type in c("unweighted", "weighted")) {
      weight_dir <- file.path(size_dir, weight_type)
      dir.create(weight_dir, showWarnings = FALSE)
      
      # Create subdirectories for each network type
      for (network_type in network_types) {
        dir.create(file.path(weight_dir, network_type), showWarnings = FALSE)
      }
    }
  }
}

# Function to save networks as edge lists
save_networks <- function(graphs, base_dir, size, weight_type, replicate) {
  for (name in names(graphs)) {
    # Get the model type (e.g., DS, SBM2, SBM3) from the name
    model_name <- sub("_graph_.*", "", name)
    
    # Create the directory path for this specific network type
    network_dir <- file.path(base_dir, paste0("size_", size), weight_type, model_name)
    
    if (weight_type == "unweighted") {
      # Save the unweighted graph with the appropriate naming convention, including replicate
      filename <- file.path(network_dir, paste0(model_name, "_unweighted_graph_n", size, "_", replicate, ".edges"))
      write.table(as_edgelist(graphs[[name]]), file = filename, row.names = FALSE, col.names = FALSE)
    } else {
      # Save the weighted graph with the appropriate naming convention, including replicate
      edgelist <- as_edgelist(graphs[[name]])
      weights <- E(graphs[[name]])$weight
      weighted_edgelist <- cbind(edgelist, weights)
      filename <- file.path(network_dir, paste0(model_name, "_weighted_graph_n", size, "_", replicate, ".edges"))
      write.table(weighted_edgelist, file = filename, row.names = FALSE, col.names = FALSE)
    }
  }
}





generate_and_save_networks <- function(network_sizes, num_replicates, EDGE_COUNT,
                                       Degseqsize, NEIGHBORS,
                                       REWIRING_PROB, BARABASI_M, SPATIAL_R,
                                       SIZES_2_BLOCKS,
                                       SIZES_3_BLOCKS,
                                       PROBABILITY_MATRIX_2_BLOCKS, 
                                       PROBABILITY_MATRIX_3_BLOCKS, 
                                       base_dir) {
  # Define network types
  network_types <- c("DS", "ER", "SBM2", "SBM3", "SF", "SP", "SW")
  
  # Create directories for all network sizes and types
  create_directories(base_dir, network_sizes, network_types)
  
  # Determine the number of cores to use (leave one core free)
  num_cores <- detectCores() - 1
  
  # Create a cluster
  cl <- makeCluster(num_cores)
  
  # Export necessary functions and variables to the cluster
  clusterExport(cl, c("generate_graphs", 
                      "save_networks", "create_directories", 
                      "fastSpatialNetwork","cellCoordsToIndex"))
  
  clusterEvalQ(cl, {
    library(igraph)
    #source("path/to/custom_functions.R")  # Replace with the actual path to your custom functions
  })
  
  # Function to process each size-replicate combination
  process_combination <- function(combination, params) {
    size <- combination$size
    replicate <- combination$replicate
    
    cat("Generating network of size", size, "replicate", replicate, "\n")
    
    # Generate the graphs
    result <- generate_graphs(size, 1, params$EDGE_COUNT, params$Degseqsize, params$NEIGHBORS,
                              params$REWIRING_PROB, params$BARABASI_M, params$SPATIAL_R, 
                              params$SIZES_2_BLOCKS, params$PROBABILITY_MATRIX_2_BLOCKS, 
                              params$SIZES_3_BLOCKS, params$PROBABILITY_MATRIX_3_BLOCKS)
    
    # Save the unweighted and weighted graphs
    save_networks(result$unweighted_graphs, params$base_dir, size, "unweighted", replicate)
    save_networks(result$weighted_graphs, params$base_dir, size, "weighted", replicate)
  }
  
  # Create a list of all size-replicate combinations
  combinations <- expand.grid(size = network_sizes, replicate = 1:num_replicates)
  
  # Create a list of parameters
  params <- list(
    EDGE_COUNT = EDGE_COUNT,
    Degseqsize = Degseqsize,
    NEIGHBORS = NEIGHBORS,
    REWIRING_PROB = REWIRING_PROB,
    BARABASI_M = BARABASI_M,
    SPATIAL_R = SPATIAL_R,
    SIZES_2_BLOCKS = SIZES_2_BLOCKS,
    SIZES_3_BLOCKS = SIZES_3_BLOCKS,
    PROBABILITY_MATRIX_2_BLOCKS = PROBABILITY_MATRIX_2_BLOCKS,
    PROBABILITY_MATRIX_3_BLOCKS = PROBABILITY_MATRIX_3_BLOCKS,
    base_dir = base_dir
  )
  
  # Use parLapply to process all combinations in parallel
  parLapply(cl, split(combinations, seq(nrow(combinations))), 
            function(comb) process_combination(comb, params))
  
  # Stop the cluster
  stopCluster(cl)
}


load_edgelists <- function(folder_path) {
  edgelists <- list()
  files <- list.files(folder_path, pattern = "\\.edges$", full.names = TRUE)
  
  print(paste("Found", length(files), ".edges files in the folder"))
  
  for (file in files) {
    print(paste("Processing file:", file))
    network_name <- tools::file_path_sans_ext(basename(file))
    edgelist <- tryCatch({
      df <- read_delim(file, delim = " ", col_names = FALSE, col_types = cols(.default = "c"))
      print(paste("File read successfully. Dimensions:", nrow(df), "x", ncol(df)))
      
      # Ensure we have at least two columns
      if (ncol(df) < 2) {
        stop("File does not have at least two columns")
      }
      
      # Convert first two columns to integer
      df[[1]] <- as.integer(df[[1]])
      df[[2]] <- as.integer(df[[2]])
      
      # If there's a third column, convert to numeric
      if (ncol(df) == 3) {
        df[[3]] <- as.numeric(df[[3]])
      }
      
      print(paste("Final dimensions:", nrow(df), "x", ncol(df)))
      as.matrix(df)
    }, error = function(e) {
      warning(paste("Error processing", file, ":", e$message))
      NULL
    })
    
    if (!is.null(edgelist)) {
      edgelists[[network_name]] <- edgelist
      print(paste("Added", network_name, "to edgelists. Dimensions:", nrow(edgelist), "x", ncol(edgelist)))
    }
  }
  
  print(paste("Total networks loaded:", length(edgelists)))
  return(edgelists)
}

getGraphFromFile <- function(file, simplify=TRUE, useBiggestComponent=TRUE, asUndirected=TRUE) {
  
  dat <- read.table(file) # just read static graph: ignoring third column
  G <- graph_from_data_frame(dat)
  if (asUndirected==TRUE) {
    G <- as.undirected(G, "collapse")
  }
  
  #  g_names <- gsub(".edges","",networks[i]) # one edge for each pair of connect vertices (not sure what this is for)
  
  if (useBiggestComponent==TRUE) {
    netclust <- components(G) #look for subgraphs
    gcc <- V(G)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    G <- induced.subgraph(G, gcc) #make it a igraph object.
  }
  if (simplify==TRUE) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  return(G)
}

###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}



calcGraphFeatures <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "edges",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "mean_degree",                # average degree of any vertex
    "minCut",                   # minimum cut weight of the graph (might take a while to compute)
    "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
    "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
    "closeness_centr",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is_connected))
  df$minCut = base::as.numeric(lapply(Graphs, min_cut))
  df$diameter = base::as.numeric(lapply(Graphs, diameter))
  df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs,igraph::degree )
  df$minDegree = base::as.numeric(lapply(degrees, min))
  df$maxDegree = base::as.numeric(lapply(degrees, max))
  df$mean_degree = base::as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  Fiedler.value=NULL
  norm.fiedler.value=NULL
  
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
    Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
    
    norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
    
    df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
    df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
    df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
    df$radius[i]<-radius(Graphs[[i]])
    df$mean_path_length[i]<-mean_distance(Graphs[[i]])
    #df$trace[i]<-sum(diag(Adj[[i]]))
    df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
    df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
    df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
    df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
    df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
    df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
    df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness_centr[i] = -1
    }
  }
  return (df)
}





countStates <- function(DF=NULL, states=NULL) {
  nr = dim(DF)[1]
  counts <- as.data.frame(matrix(0,ncol=length(states), nrow=nr))
  colnames(counts) = as.character(states)
  for (i in 1:nr) {
    col = 1
    for (state in states) {
      counts[i,col] = length(which(DF[i,]==state))
      col = col+1
    }
  }
  return(counts)
}


simPathE <- function(G, nTicks=100, beta=0.5, gamma=0.2, propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F) {
  nVertices = gorder(G)
  infectionState = as.data.frame(matrix(0, ncol = nVertices, nrow = nTicks+1))
  
  # Set up the initial state: all vertices are either exposed (state = 0) or infected (state = 1); none can be recovered yet (2).
  if (is.null(initialState)) {
    if (useProportion==T) {
      infectionState[1,] <- rbinom(nVertices, 1, propInfected) # set initial state of each nodes if known (note, no recovered node at initial state)
    } else {
      infected <- rep(1, nInfected) # just create a vector of the right number of 1s
      exposed <- rep(0, (nVertices - nInfected))
      infectionState[1,] <- sample(c(infected, exposed), nVertices, replace=FALSE)
    }
  } else {
    if (length(initialState) != nVertices) {
      return ("Initial state and order of Graph (number of vertices) are incompatible.  Check the sizes of your input.")
    }
    infectionState <- initialState # initial existing state.
  }
  
  adjacencyList <- as_adj_list(G) # this is a list of which vertices are adjacent to each vertex i
  
  # Now do the simulation through time:
  for (t in 1:nTicks) {
    # FIRST phase: transmission: S -> I
    for (i in which(infectionState[t,] == 0)) { # for all susceptible nodes (denoted 0) in previous time step
      infectionState[t+1,i] <- 0 # node remains as susceptible, since not all contact leads to an infection.
      for (j in adjacencyList[[i]]) { # for all neighbours of i
        if (infectionState[t,j][1] == 1) { # vertex j is infectious
          if ((runif(1)*G[i,j]) <= beta) { # ... and passes it on!
            infectionState[t+1,i] <- 1;
            break # assign node as infected if above condition is met, and break out of loop: we don't need
            # to check any more adjacent vertices.
          }
        }
      }
    }
    
    # SECOND phase: recovery: I -> R
    for (i in which(infectionState[t,] == 1)) { # for all infected nodes (denoted 1) in previous time step
      if (runif(1) <= gamma) { # compares a randomly generated uniform number to recovery rate
        infectionState[t+1,i] <- 2 # node is recovered
      } else {
        infectionState[t+1,i] <- 1 # node remains infected
      }
    }
    
    # THIRD phase: recovered stays recovered:
    for (i in which(infectionState[t,] == 2)) { # for all recovered nodes (denoted 2) in previous time step
      infectionState[t+1,i] <- 2 # node stays recovered
    }
  }
  rownames(infectionState) = 0:nTicks
  return(infectionState)
}


plan(multisession)
# Network_sim_pipeline <- function(folder_name, nsim = 2, nreps = 1, beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2), gamma_values = c(0.04, 0.2)) {
#   # List all .edges files in the directory
#   networks <- list.files(folder_name, pattern = "\\.edges$", full.names = TRUE)
#   
#   # Replicate each network nsim times
#   replicated_networks <- rep(networks, each = nreps)
#   
#   # Initialize an empty list to store results
#   Global_summary <- list()
#   
#   # Perform serialization processing (replace foreach with a loop)
#   for (i in 1:length(replicated_networks)) {
#     filename <- replicated_networks[i]
#     cat("Processing file:", filename, "\n")
#     
#     dat <- tryCatch({
#       read.table(filename, header = FALSE, stringsAsFactors = FALSE)
#     }, error = function(e) {
#       cat("Error reading file:", filename, "\n", e$message, "\n")
#       return(NULL)  # Return NULL if there's an error reading the file
#     })
#     
#     if (is.null(dat)) {
#       next  # Skip this iteration if dat is NULL
#     }
#     
#     # Determine if the network is weighted based on the number of columns
#     if (ncol(dat) == 2) {
#       weighted <- FALSE
#       g <- graph_from_data_frame(dat, directed = FALSE)
#     } else if (ncol(dat) == 3) {
#       weighted <- TRUE
#       g <- graph_from_data_frame(dat[, 1:2], directed = FALSE)
#       E(g)$weight <- dat[, 3]
#     } else {
#       cat("Invalid file format:", filename, "\n")
#       next
#     }
#     
#     # Simplify the graph and get the giant connected component
#     netclust <- components(g)
#     gcc <- V(g)[netclust$membership == which.max(netclust$csize)]
#     gccnet <- induced_subgraph(g, gcc)
#     gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
#     
#     network_size <- vcount(gcc_s)
#     network_edge <- ecount(gcc_s)
#     
#     # Adjacency matrix (with or without weights)
#     Adj <- as_adjacency_matrix(gcc_s, attr = if (weighted) "weight" else NULL)
#     
#     # Laplacian matrices (with or without weights)
#     L <- laplacian_matrix(gcc_s, normalized = FALSE, weights = if (weighted) E(gcc_s)$weight else NULL)
#     Norm_Lap <- laplacian_matrix(gcc_s, normalized = TRUE, weights = if (weighted) E(gcc_s)$weight else NULL)
#     
#     # Spectral radius (largest eigenvalue)
#     Adjacency_spectrum <- eigen(Adj)
#     Spectral_radius <- Adjacency_spectrum$values[1]
#     
#     spec <- eigen(L)
#     norm_spec <- eigen(Norm_Lap)
#     FiedlerValue <- spec$values[length(spec$values) - 1]
#     NormFiedlerValue <- norm_spec$values[length(norm_spec$values) - 1]
#     Eigen_central <- Adjacency_spectrum$vectors[, ncol(Adjacency_spectrum$vectors)]
#     
#     # Eigenvector centrality
#     if (weighted) {
#       eigen_cent <- eigen_centrality(gcc_s, directed = FALSE, scale = TRUE, weights = E(gcc_s)$weight)
#       eigen_centr <- (sum(max(eigen_cent$vector) - eigen_cent$vector)) / ((vcount(gcc_s) - 1) * (max(eigen_cent$vector) - min(eigen_cent$vector)))
#     } else {
#       eigen_centr <- centr_eigen(gcc_s, scale = T, directed = FALSE)$centralization
#     }
#     
#     Max_principal_eigenvector <- max(abs(Adjacency_spectrum$vectors[, 1]))
#     
#     # Mean degree
#     if (weighted) {
#       deg <- igraph::strength(gcc_s, weights = E(gcc_s)$weight)
#     } else {
#       deg <- igraph::degree(gcc_s)
#     }
#     
#     mean_deg <- mean(deg, normalized = TRUE)
#     Most_connected_node <- which.max(deg) # highest node degree
#     
#     # Spectral threshold (1/Spectral radius)
#     Rnot <- 1 / Spectral_radius
#     
#     g_names <- gsub(".edges", "", basename(filename))
#     connected <- as.numeric(is_connected(gcc_s))
#     most_infected_node <- match(max(Eigen_central), Eigen_central)  
#     cent <- eigen_centrality(gcc_s, directed = F, scale = T,
#                              weights = if (weighted) E(gcc_s)$weight else NULL)$value
#     
#     # Transitivity (global clustering coefficient)
#     trans <- transitivity(gcc_s, type = "global", weights = if (weighted) E(gcc_s)$weight else NULL)
#     
#     # Modularity and community detection
#     lc <- cluster_louvain(gcc_s, weights = if (weighted) E(gcc_s)$weight else NULL)
#     mod <- modularity(gcc_s, membership(lc), weights = if (weighted) E(gcc_s)$weight else NULL)
#     max_mod <- assortnet::assortment.discrete(Adj, types = as.factor(membership(lc)), weighted = TRUE)$r
#     Qrel <- mod / max_mod
#     
#     # Mean shortest path length
#     mean_dist <- mean_distance(gcc_s, directed = FALSE, weights = if (weighted) E(gcc_s)$weight else NULL)
#     
#     # Degree assortativity
#     if (weighted) {
#       strengths <- strength(gcc_s, mode = "all", weights = E(gcc_s)$weight)
#       deg_assort_coef <- assortativity(gcc_s, strengths, directed = FALSE)
#     } else {
#       deg_assort_coef <- assortativity_degree(gcc_s, directed = FALSE)
#     }
#     
#     # Degree centralization
#     if (weighted) {
#       strengths <- strength(gcc_s, mode = "all", weights = E(gcc_s)$weight)
#       deg_centr <- (sum(max(strengths) - strengths)) / ((vcount(gcc_s) - 1) * (vcount(gcc_s) - 2))
#     } else {
#       deg_centr <- centr_degree(gcc_s)$centralization
#     }
#     
#     # Eccentricity
#     if (weighted) {
#       dist_matrix <- distances(gcc_s, weights = E(gcc_s)$weight)
#       eccentricities <- apply(dist_matrix, 1, max)
#       mean_eccentr <- mean(eccentricities)
#     } else {
#       mean_eccentr <- mean(eccentricity(gcc_s))
#     }
#     
#     # Simulate SIR model over the network
#     simulation_results <- expand.grid(beta = beta_values, gamma = gamma_values)
#     simulation_results$avg_invasion_time <- NA
#     simulation_results$avg_prop_infected <- NA
#     
#     for (j in 1:nrow(simulation_results)) {
#       beta <- simulation_results$beta[j]
#       gamma <- simulation_results$gamma[j]
#       
#       inf_data <- list()
#       
#       for (k in 1:nsim) {
#         sm <- simPathE(gcc_s, nTicks = 100, beta = beta, gamma = gamma,
#                        propInfected = 0.05, initialState = NULL, nInfected = 1, useProportion = TRUE)
#         countEachState <- countStates(DF = sm, states = c(0:2))
#         
#         TwoPcInfected <- sum(sm[1,]) + ceiling((vcount(gcc_s) / 100) * 2)
#         invasion_time <- which(countEachState[, 2] >= TwoPcInfected)[1]
#         prop_infected <- max(countEachState[, 2]) / vcount(gcc_s)
#         
#         inf_data[[k]] <- c(invasion_time = invasion_time, prop_infected = prop_infected)
#       }
#       
#       complete_data <- as.data.frame(do.call(rbind, inf_data))
#       complete_data_noNA <- complete_data[complete.cases(complete_data),]
#       
#       simulation_results$avg_invasion_time[j] <- mean(complete_data_noNA$invasion_time)
#       simulation_results$avg_prop_infected[j] <- mean(complete_data_noNA$prop_infected)
#     }
#     
#     # Create the results data frame
#     results_df <- data.frame(
#       Network = g_names,
#       Network_size = network_size,
#       Edge = network_edge,
#       Connected = connected,
#       Fiedler = FiedlerValue,
#       Normalized_Fiedler = NormFiedlerValue,
#       Spectral_radius = Spectral_radius,
#       Eigen_centrality = eigen_centr,
#       Transitivity = trans,
#       Modularity = mod,
#       Relative_modularity = Qrel,
#       Mean_degree = mean_deg,
#       Mean_path_length = mean_dist,
#       Mean_eccentricity = mean_eccentr,
#       Degree_assortativity = deg_assort_coef,
#       Degree_centrality = deg_centr,
#       Max_principal_eigenvector = Max_principal_eigenvector,
#       Most_connected_node = Most_connected_node,
#       Rnot = Rnot,
#       Most_infected_node = most_infected_node
#     )
#     
#     # Add simulation results to the data frame
#     for (j in 1:nrow(simulation_results)) {
#       beta <- simulation_results$beta[j]
#       gamma <- simulation_results$gamma[j]
#       col_name_time <- paste0("avg_invasion_time_beta", beta, "_gamma", gamma)
#       col_name_prop <- paste0("avg_prop_infected_beta", beta, "_gamma", gamma)
#       results_df[[col_name_time]] <- simulation_results$avg_invasion_time[j]
#       results_df[[col_name_prop]] <- simulation_results$avg_prop_infected[j]
#     }
#     
#     # Append to the global summary
#     Global_summary[[i]] <- results_df
#   }
#   
#   # Combine all results into a single data frame
#   Global_summary <- do.call(rbind, Global_summary)
#   
#   return(Global_summary)
# }

# beta_values <- c(0.01, 0.025, 0.05, 0.1, 0.2)
# gamma_values <- c(0.04, 0.4)

Network_sim_pipeline <- function(folder_name, nsim = 2, nreps = 1, beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2), gamma_values = c(0.04, 0.2)) {
  # Create output directory for RDS files if it doesn't exist
  output_dir <- file.path(folder_name, "simulation_results")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # List all .edges files in the directory
  networks <- list.files(folder_name, pattern = "\\.edges$", full.names = TRUE)
  
  # Replicate each network nsim times
  replicated_networks <- rep(networks, each = nreps)
  
  # Initialize an empty list to store results
  Global_summary <- list()
  
  # Perform serialization processing
  for (i in 1:length(replicated_networks)) {
    filename <- replicated_networks[i]
    cat("Processing file:", filename, "\n")
    
    # Create unique identifier for this simulation
    sim_id <- paste0("sim_", basename(tools::file_path_sans_ext(filename)), "_", i)
    rds_path <- file.path(output_dir, paste0(sim_id, ".rds"))
    
    # Check if simulation already exists
    if (file.exists(rds_path)) {
      cat("Simulation", sim_id, "already exists, skipping...\n")
      # Load existing result into Global_summary
      Global_summary[[i]] <- readRDS(rds_path)
      next
    }
    
    dat <- tryCatch({
      read.table(filename, header = FALSE, stringsAsFactors = FALSE)
    }, error = function(e) {
      cat("Error reading file:", filename, "\n", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(dat)) {
      next
    }
    
    # Determine if the network is weighted
    if (ncol(dat) == 2) {
      weighted <- FALSE
      g <- graph_from_data_frame(dat, directed = FALSE)
    } else if (ncol(dat) == 3) {
      weighted <- TRUE
      g <- graph_from_data_frame(dat[, 1:2], directed = FALSE)
      E(g)$weight <- dat[, 3]
    } else {
      cat("Invalid file format:", filename, "\n")
      next
    }
    
    # Simplify the graph and get the giant connected component
    netclust <- components(g)
    gcc <- V(g)[netclust$membership == which.max(netclust$csize)]
    gccnet <- induced_subgraph(g, gcc)
    gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
    
    network_size <- vcount(gcc_s)
    network_edge <- ecount(gcc_s)
    
    # Adjacency matrix (with or without weights)
    Adj <- as_adjacency_matrix(gcc_s, attr = if (weighted) "weight" else NULL)
    
    # Laplacian matrices (with or without weights)
    L <- laplacian_matrix(gcc_s, normalized = FALSE, weights = if (weighted) E(gcc_s)$weight else NULL)
    Norm_Lap <- laplacian_matrix(gcc_s, normalized = TRUE, weights = if (weighted) E(gcc_s)$weight else NULL)
    
    # Spectral radius (largest eigenvalue)
    Adjacency_spectrum <- eigen(Adj)
    Spectral_radius <- Adjacency_spectrum$values[1]
    
    spec <- eigen(L)
    norm_spec <- eigen(Norm_Lap)
    FiedlerValue <- spec$values[length(spec$values) - 1]
    NormFiedlerValue <- norm_spec$values[length(norm_spec$values) - 1]
    Eigen_central <- Adjacency_spectrum$vectors[, ncol(Adjacency_spectrum$vectors)]
    # Eigenvector centrality
    if (weighted) {
      eigen_cent <- eigen_centrality(gcc_s, directed = FALSE, scale = TRUE, weights = E(gcc_s)$weight)
      eigen_centr <- (sum(max(eigen_cent$vector) - eigen_cent$vector)) / ((vcount(gcc_s) - 1) * (max(eigen_cent$vector) - min(eigen_cent$vector)))
    } else {
      eigen_centr <- centr_eigen(gcc_s, scale = T, directed = FALSE)$centralization
    }
    
    Max_principal_eigenvector <- max(abs(Adjacency_spectrum$vectors[, 1]))
    
    # Mean degree
    if (weighted) {
      deg <- igraph::strength(gcc_s, weights = E(gcc_s)$weight)
    } else {
      deg <- igraph::degree(gcc_s)
    }
    
    mean_deg <- mean(deg, normalized = TRUE)
    Most_connected_node <- which.max(deg)
    
    # Spectral threshold (1/Spectral radius)
    Rnot <- 1 / Spectral_radius
    
    g_names <- gsub(".edges", "", basename(filename))
    connected <- as.numeric(is_connected(gcc_s))
    most_infected_node <- match(max(Eigen_central), Eigen_central)  
    cent <- eigen_centrality(gcc_s, directed = F, scale = T,
                             weights = if (weighted) E(gcc_s)$weight else NULL)$value
    
    # Transitivity (global clustering coefficient)
    trans <- transitivity(gcc_s, type = "global", weights = if (weighted) E(gcc_s)$weight else NULL)
    
    # Modularity and community detection
    lc <- cluster_louvain(gcc_s, weights = if (weighted) E(gcc_s)$weight else NULL)
    mod <- modularity(gcc_s, membership(lc), weights = if (weighted) E(gcc_s)$weight else NULL)
    max_mod <- assortnet::assortment.discrete(Adj, types = as.factor(membership(lc)), weighted = TRUE)$r
    Qrel <- mod / max_mod
    
    # Mean shortest path length
    mean_dist <- mean_distance(gcc_s, directed = FALSE, weights = if (weighted) E(gcc_s)$weight else NULL)
    
    # Degree assortativity
    if (weighted) {
      strengths <- strength(gcc_s, mode = "all", weights = E(gcc_s)$weight)
      deg_assort_coef <- assortativity(gcc_s, strengths, directed = FALSE)
    } else {
      deg_assort_coef <- assortativity_degree(gcc_s, directed = FALSE)
    }
    
    # Degree centralization
    if (weighted) {
      strengths <- strength(gcc_s, mode = "all", weights = E(gcc_s)$weight)
      deg_centr <- (sum(max(strengths) - strengths)) / ((vcount(gcc_s) - 1) * (vcount(gcc_s) - 2))
    } else {
      deg_centr <- centr_degree(gcc_s)$centralization
    }
    
    # Eccentricity
    if (weighted) {
      dist_matrix <- distances(gcc_s, weights = E(gcc_s)$weight)
      eccentricities <- apply(dist_matrix, 1, max)
      mean_eccentr <- mean(eccentricities)
    } else {
      mean_eccentr <- mean(eccentricity(gcc_s))
    }
    # Simulate SIR model over the network
    simulation_results <- expand.grid(beta = beta_values, gamma = gamma_values)
    simulation_results$avg_invasion_time <- NA
    simulation_results$avg_prop_infected <- NA
    
    for (j in 1:nrow(simulation_results)) {
      beta <- simulation_results$beta[j]
      gamma <- simulation_results$gamma[j]
      
      inf_data <- list()
      
      for (k in 1:nsim) {
        sm <- simPathE(gcc_s, nTicks = 100, beta = beta, gamma = gamma,
                       propInfected = 0.05, initialState = NULL, nInfected = 1, useProportion = TRUE)
        countEachState <- countStates(DF = sm, states = c(0:2))
        
        TwoPcInfected <- sum(sm[1,]) + ceiling((vcount(gcc_s) / 100) * 2)
        invasion_time <- which(countEachState[, 2] >= TwoPcInfected)[1]
        prop_infected <- max(countEachState[, 2]) / vcount(gcc_s)
        
        inf_data[[k]] <- c(invasion_time = invasion_time, prop_infected = prop_infected)
      }
      
      complete_data <- as.data.frame(do.call(rbind, inf_data))
      complete_data_noNA <- complete_data[complete.cases(complete_data),]
      
      simulation_results$avg_invasion_time[j] <- mean(complete_data_noNA$invasion_time)
      simulation_results$avg_prop_infected[j] <- mean(complete_data_noNA$prop_infected)
    }
    
    # Create the results data frame
    results_df <- data.frame(
      Network = g_names,
      Network_size = network_size,
      Edge = network_edge,
      Connected = connected,
      Fiedler = FiedlerValue,
      Normalized_Fiedler = NormFiedlerValue,
      Spectral_radius = Spectral_radius,
      Eigen_centrality = eigen_centr,
      Transitivity = trans,
      Modularity = mod,
      Relative_modularity = Qrel,
      Mean_degree = mean_deg,
      Mean_path_length = mean_dist,
      Mean_eccentricity = mean_eccentr,
      Degree_assortativity = deg_assort_coef,
      Degree_centrality = deg_centr,
      Max_principal_eigenvector = Max_principal_eigenvector,
      Most_connected_node = Most_connected_node,
      Rnot = Rnot,
      Most_infected_node = most_infected_node
    )
    
    # Add simulation results to the data frame
    for (j in 1:nrow(simulation_results)) {
      beta <- simulation_results$beta[j]
      gamma <- simulation_results$gamma[j]
      col_name_time <- paste0("avg_invasion_time_beta", beta, "_gamma", gamma)
      col_name_prop <- paste0("avg_prop_infected_beta", beta, "_gamma", gamma)
      results_df[[col_name_time]] <- simulation_results$avg_invasion_time[j]
      results_df[[col_name_prop]] <- simulation_results$avg_prop_infected[j]
    }
    
    # Save individual simulation results
    saveRDS(results_df, rds_path)
    cat("Saved simulation results to:", rds_path, "\n")
    
    # Add to Global_summary
    Global_summary[[i]] <- results_df
  }
  
  # Combine all results into a single data frame
  final_summary <- do.call(rbind, Global_summary)
  
  
  # Save the combined results
  final_rds_path <- file.path(output_dir, "combined_results.rds")
  saveRDS(final_summary, final_rds_path)
  cat("Saved combined results to:", final_rds_path, "\n")
  
  return(final_summary)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Network size of 50
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##########################################################
# Erdos renyi
#########################################################
base_path <- "network_sim_50/ER_ALL_DATA_50"
folder_name=base_path

sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
#beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

er_sim_all_data_50=Network_sim_pipeline(base_path,
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(er_sim_all_data_50,"er_sim_all_data_50.rds")

###########################################################
# Spatial
###########################################################
base_path <- "network_sim_50/SP_ALL_DATA_50"
sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

sp_sim_all_data_50=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sp_sim_all_data_50,"sp_sim_all_data_50.rds")


###########################################################
# Small world
###########################################################
base_path <- "network_sim_50/SW_ALL_DATA_50"
sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

sw_sim_all_data_50=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sw_sim_all_data_50,"sw_sim_all_data_50.rds")

###########################################################
# SBM2
###########################################################
base_path <- "network_sim_50/SBM2_ALL_DATA_50"
sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

sbm2_sim_all_data_50=Network_sim_pipeline(base_path, 
                                       nsim = nsim,
                                       nreps = 1,
                                       beta_values = beta_values, 
                                       gamma_values = gamma_values)
 
saveRDS(sbm2_sim_all_data_50,"sbm2_sim_all_data_50.rds")

# ###########################################################
# # SBM3
# ###########################################################
# base_path <- "network_sim_50/SBM3_ALL_DATA_50"
# sizes <- c(50)  # Add more sizes if needed
# beta_values <- c(0.05,0.1, 0.25)
# gamma_values <- c(0.04, 0.4)
# nsim <- 50
# nreps <- 1
# folder_name=base_path
# 
# sbm3_sim_all_data_50=Network_sim_pipeline(base_path, 
#                                        nsim = nsim,
#                                        nreps = 1,
#                                        beta_values = beta_values, 
#                                        gamma_values = gamma_values)
# 
# saveRDS(sbm3_sim_all_data_50,"sbm3_sim_all_data_50.rds")
# 

###########################################################
# Scale-free
###########################################################
base_path <- "network_sim_50/SF_ALL_DATA_50"
sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

sf_sim_all_data_50=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sf_sim_all_data_50,"sf_sim_all_data_50.rds")

###########################################################
# Degree sequence
###########################################################
base_path <- "network_sim_50/DS_ALL_DATA_50"
sizes <- c(50)  # Add more sizes if needed
beta_values = c(0.01, 0.025, 0.05, 0.1, 0.2)
gamma_values <- c(0.04, 0.4)
nsim <- 100
nreps <- 1
folder_name=base_path

ds_sim_all_data_50=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(ds_sim_all_data_50,"ds_sim_all_data_50.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# END of Network size of 50
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#########################################################
# Erdos renyi
#########################################################
base_path <- "network_sim_100/ER_ALL_DATA_100"
folder_name=base_path

sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

er_sim_all_data_100=Network_sim_pipeline(base_path,
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(er_sim_all_data_100,"er_sim_all_data_100.rds")

###########################################################
# Spatial
###########################################################
base_path <- "network_sim_100/SP_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sp_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sp_sim_all_data_100,"sp_sim_all_data_100.rds")


###########################################################
# Small world
###########################################################
base_path <- "network_sim_100/SW_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sw_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sw_sim_all_data_100,"sw_sim_all_data_100.rds")

###########################################################
# SBM2
###########################################################
base_path <- "network_sim_100/SBM2_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm2_sim_all_data_100=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm2_sim_all_data_100,"sbm2_sim_all_data_100.rds")

###########################################################
# SBM3
###########################################################
base_path <- "network_sim_100/SBM3_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm3_sim_all_data_100=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm3_sim_all_data_100,"sbm3_sim_all_data_100.rds")


###########################################################
# Scale-free
###########################################################
base_path <- "network_sim_100/SF_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_pathss

sf_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sf_sim_all_data_100,"sf_sim_all_data_100.rds")

###########################################################
# Degree sequence
###########################################################
base_path <- "network_sim_100/DS_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

ds_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(ds_sim_all_data_100,"ds_sim_all_data_100.rds")





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Network size of 250
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#


##########################################################
# Erdos renyi
#########################################################
base_path <- "network_sim_250/ER-ALL-DATA"
folder_name=base_path

sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

er_sim_all_data=Network_sim_pipeline(base_path,
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(er_sim_all_data,"er_sim_all_data.rds")

###########################################################
# Spatial
###########################################################
base_path <- "network_sim_250/SP-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sp_sim_all_data=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sp_sim_all_data,"sp_sim_all_data.rds")


###########################################################
# Small world
###########################################################
base_path <- "network_sim_250/SW-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sw_sim_all_data=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sw_sim_all_data,"sw_sim_all_data.rds")

###########################################################
# SBM2
###########################################################
base_path <- "network_sim_250/SBM2-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm2_sim_all_data=Network_sim_pipeline(base_path, 
                                       nsim = nsim,
                                       nreps = 1,
                                       beta_values = beta_values, 
                                       gamma_values = gamma_values)

saveRDS(sbm2_sim_all_data,"sbm2_sim_all_data.rds")

###########################################################
# SBM3
###########################################################
base_path <- "network_sim_250/SBM3-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm3_sim_all_data=Network_sim_pipeline(base_path, 
                                       nsim = nsim,
                                       nreps = 1,
                                       beta_values = beta_values, 
                                       gamma_values = gamma_values)

saveRDS(sbm3_sim_all_data,"sbm3_sim_all_data.rds")


###########################################################
# Scale-free
###########################################################
base_path <- "network_sim_250/SF-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_pathss

sf_sim_all_data=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(sf_sim_all_data,"sf_sim_all_data.rds")

###########################################################
# Degree sequence
###########################################################
base_path <- "network_sim_250/DS-ALL-DATA"
sizes <- c(250)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

ds_sim_all_data=Network_sim_pipeline(base_path, 
                                     nsim = nsim,
                                     nreps = 1,
                                     beta_values = beta_values, 
                                     gamma_values = gamma_values)

saveRDS(ds_sim_all_data,"ds_sim_all_data.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Network size of 500
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##########################################################
# Erdos renyi
#########################################################
base_path <- "network_sim_500/ER_ALL_DATA_500"
folder_name=base_path

sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

er_sim_all_data_500=Network_sim_pipeline(base_path,
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(er_sim_all_data_500,"er_sim_all_data_500.rds")

###########################################################
# Spatial
###########################################################
base_path <- "network_sim_500/SP_ALL_DATA_500"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sp_sim_all_data_500=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sp_sim_all_data_500,"sp_sim_all_data_500.rds")


###########################################################
# Small world
###########################################################
base_path <- "network_sim_500/SW_ALL_DATA_500"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sw_sim_all_data_500=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sw_sim_all_data_500,"sw_sim_all_data_500.rds")

###########################################################
# SBM2
###########################################################
base_path <- "network_sim_500/SBM2_ALL_DATA_500"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm2_sim_all_data_500=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm2_sim_all_data_500,"sbm2_sim_all_data_500.rds")

###########################################################
# SBM3
###########################################################
base_path <- "network_sim_50/SBM3_ALL_DATA_50"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm3_sim_all_data_500=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm3_sim_all_data_500,"sbm3_sim_all_data_500.rds")


###########################################################
# Scale-free
###########################################################
base_path <- "network_sim_500/SF_ALL_DATA_500"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_pathss

sf_sim_all_data_500=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sf_sim_all_data_500,"sf_sim_all_data_500.rds")

###########################################################
# Degree sequence
###########################################################
base_path <- "network_sim_500/DS_ALL_DATA_500"
sizes <- c(500)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

ds_sim_all_data_500=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(ds_sim_all_data_500,"ds_sim_all_data_500.rds")



#########################################################
# Erdos renyi
#########################################################
base_path <- "network_sim_100/ER_ALL_DATA_100"
folder_name=base_path

sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

er_sim_all_data_100=Network_sim_pipeline(base_path,
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(er_sim_all_data_100,"er_sim_all_data_100.rds")

###########################################################
# Spatial
###########################################################
base_path <- "network_sim_100/SP_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sp_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sp_sim_all_data_100,"sp_sim_all_data_100.rds")


###########################################################
# Small world
###########################################################
base_path <- "network_sim_100/SW_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sw_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sw_sim_all_data_100,"sw_sim_all_data_100.rds")

###########################################################
# SBM2
###########################################################
base_path <- "network_sim_100/SBM2_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm2_sim_all_data_100=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm2_sim_all_data_100,"sbm2_sim_all_data_100.rds")

###########################################################
# SBM3
###########################################################
base_path <- "network_sim_100/SBM3_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

sbm3_sim_all_data_100=Network_sim_pipeline(base_path, 
                                           nsim = nsim,
                                           nreps = 1,
                                           beta_values = beta_values, 
                                           gamma_values = gamma_values)

saveRDS(sbm3_sim_all_data_100,"sbm3_sim_all_data_100.rds")


###########################################################
# Scale-free
###########################################################
base_path <- "network_sim_100/SF_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_pathss

sf_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(sf_sim_all_data_100,"sf_sim_all_data_100.rds")

###########################################################
# Degree sequence
###########################################################
base_path <- "network_sim_100/DS_ALL_DATA_100"
sizes <- c(100)  # Add more sizes if needed
beta_values <- c(0.05,0.1, 0.25)
gamma_values <- c(0.04, 0.4)
nsim <- 50
nreps <- 1
folder_name=base_path

ds_sim_all_data_100=Network_sim_pipeline(base_path, 
                                         nsim = nsim,
                                         nreps = 1,
                                         beta_values = beta_values, 
                                         gamma_values = gamma_values)

saveRDS(ds_sim_all_data_100,"ds_sim_all_data_100.rds")




####################################################################################
# Feature and Epi measure analysis
####################################################################################
x=readRDS("network_results_size_nsim100_50.rds")
colnames(x)
y= x%>%
dplyr::count(WeightType)


head(x)


###---Pearsons correlation----
cor(x$Spectral_radius, x$avg_invasion_time_beta0.025_gamma0.4, method = "pearson")


###---Multiple linear regression---
lm_multiple <- lm(avg_invasion_time_beta0.025_gamma0.4 ~ Spectral_radius + Fiedler + Modularity + Degree_assortativity, data = x)
summary(lm_multiple)

###---Stepwise selection----
stepwise_model <- step(lm(avg_invasion_time_beta0.025_gamma0.4 ~ 1, data = x), 
                       scope = list(lower = ~1, upper = ~ Spectral_radius + Fiedler + Modularity + Degree_assortativity),
                       direction = "forward")
summary(stepwise_model)

###---Partial correlation----
library(ppcor)
pcor_result <- pcor(cbind(x$Spectral_radius, x$avg_invasion_time_beta0.025_gamma0.4, x$Fiedler))
pcor_result$estimate

###---PCA----
pca_result <- prcomp(x[, c("Spectral_radius", "Fiedler", "Modularity", "Degree_assortativity")], scale = TRUE)
summary(pca_result)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Machine learning Disease Simulation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#from above Network_sim_pipeline function

 Simdata <- readRDS("network_results_size_nsim100_50.rds")
 
 Simdata[is.na(Simdata)] <- 0
 
 #note that - were changed to _ to avoid downstream issues in names. 
 colnames(Simdata)[1] <- c('Graph.Name')
 str(Simdata)  
 
 all_data <- Simdata
 
 str(all_data)
 

 
 mean(all_data$Network_size)
# 
# #most_infected_node removed as not a meaningful predictor in these models
# 
 Xworking <- all_data[-c(1:21)] #features
# 
# Ys <- all_data[c(2:21)] # labels all eta and gamma for prop infected and invasion time
# 
# Ys_prop <- Ys %>% select(contains("prop"))
# 
# Ys_time <- Ys %>% select(contains("time"))
# 
# Xfiltered <- subset(Xworking, select=-c(Source,  Edge_weight)) #these features aren't useful here
# 
# glimpse(Xfiltered)
# 
# dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Highest_degree', 'Adj_val','Diameter', 'Mean_degree' ))

##############################################################################################
#PCA 
##############################################################################################
# 
# pca_data <- Xworking[1:11]
# pca_res <- prcomp(pca_data, scale. = TRUE)
# 
# autoplot(pca_res)
# 
# #with labels to work out outliers
# 
# autoplot(pca_res, data = Xworking, colour = 'Class',
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 3,
#          label = TRUE)
# 
# #not labels
# 
# p <- autoplot(pca_res, data = Xworking, colour = 'Class',
#               loadings = TRUE, loadings.colour = 'black',
#               loadings.label = TRUE, loadings.label.size = 3)+
#   scale_colour_manual(values=c("azure3",'darkblue','green', 'deepskyblue','darkgray' ))
# 
# p+theme_bw()
# 
# # scale_fill_manual(colo)   
# 
# #Outlier (459) = complete voles network
# #Outlier (500) ~complete tortoise network.
# 
# #kmeans
# 
# autoplot(kmeans(pca_data, 3), data = pca_data, label = TRUE, label.size = 3)



##############################################################################################
#Correlations
#############################################################################################

# #species has too many levels so remove for the time being
# Xfiltered <- subset(Xworking, select=-c(Source, Species, Edge_weight, Class)) #these features aren't useful here
# ggpairs(Xfiltered) #lots of strong correlations
# 
# 
# dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Diameter','Modularity', 'Mean_pathLength', 'Mean_degree', 'Class', 'Highest_degree' ))
# ggpairs(dataNoCOr)
# #modularity and mean degree strongly correlated for example
# 
# glimpse(dataNoCOr)
# create_report(Xs)
# #beta 0.5 and beta 0.1 estimates strongly correlated. Choose 0.05 
# 
# 
# ##############################################################################################
# #test on a subset with the most common assocation type
# ##############################################################################################
# 
# #data_interactionType <- dataNoUnder10 %>% filter(Interaction_type=='Social_projection_bipartite')
# 
# data_interactionType <- dataNoUnder10 %>% filter(Interaction_type=='Physical_contact')
# 
# Xworking <- data_interactionType[-c(1:21)]
# 
# Ys <- data_interactionType[c(2:21)]
# 
# Ys_prop <- Ys %>% select(contains("prop"))
# 
# Ys_time <- Ys %>% select(contains("time"))
# 
# Xfiltered <- subset(Xworking, select=-c(Source,  Edge_weight)) #these features aren't useful here
# 
# glimpse(Xfiltered)
# 
# dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Highest_degree','Diameter', 'Mean_degree', 'Interaction_type',
#                                       'Modularity', 'Class', 'Mean_pathLength'))
# 
# 




##############################################################################################
#ML models
##############################################################################################

# model1 <- 
#   rand_forest(trees = 1000, mode = "regression", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity
#   set_engine("ranger", importance = "impurity") #random forest doesn't work well with dummy variables
# 
# model2 <- linear_reg() %>% 
#   set_engine("lm") %>% 
#   set_mode("regression")
# 
# library(kernlab)
# model3<-
#   svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
#   set_engine("kernlab") %>%
#   set_mode("regression")
# 
# 
# 
# #X <- dataNoUnder10 %>% select(avg_prop_infected_beta0.05)
# 
# # Define set the outcomes of interest
# 
# Xraw <- dataNoCOr #up to here  - need to add species again (no -s)
# 
# 
# #If needed - need to upadte MrIML
# library(fastDummies)
# 
# Xdummy <- dummy_cols(X, remove_first_dummy = TRUE, remove_selected_columns = TRUE)
# 
# #X %>% 
# # rename(`Captive_Semi-ranging` = Captive_Semi_ranging) #this isnt working hyphen
# 
# #colnames(Ydummy[6]) <- 'Captive_Semi_ranging' #r doesn't like '-' sometimes
# 
# #set up multicore
# cl <- parallel::makeCluster(5)
# plan(cluster, workers=cl)
# 
# --saveRDS(yhats_rf_prop , 'yhT_prop_infected_socialBipartit_July2023a')
# 
# #save(yhats, file='logreg_model')
# ModelPerf <- mrIMLperformance(yhats_rf_prop, model1, Y=Ys_prop) 
# ModelPerf 
# 
# #extract predictions
# predictionsR05 <- collect_predictions(yhats_rf_prop[[1]]$last_mod_fit, summarize = TRUE, , grid[1, ]) %>% arrange(.row)
# 
# training_pred <- 
#   predict(yhats_rf_prop[[1]]$mod1_k, as.data.frame(yhats_rf_prop[[1]]$data_train)) %>% 
#   bind_cols(predict(yhats_rf_prop[[1]]$mod1_k,yhats_rf_prop[[1]]$data_train)) %>% 
#   # Add the true outcome data back in
#   bind_cols(as.data.frame(yhats_rf_prop[[1]]$data_train) %>% 
#               select(class))
# 
# #predicts well on test data
# testing_pred <- 
#   predict(yhats_rf_prop[[1]]$mod1_k, as.data.frame(yhats_rf_prop[[1]]$data_testa)) %>% 
#   bind_cols(predict(yhats_rf_prop[[1]]$mod1_k,yhats_rf_prop[[1]]$data_testa)) %>% 
#   # Add the true outcome data back in
#   bind_cols(as.data.frame(yhats_rf_prop[[1]]$data_testa) %>% 
#               select(class))
# 
# #should now be fixed - work off mrvip_v2
# VI <- mrVip (yhats_rf_prop, X=Xraw)
# 
# p <- plot_vi(VI=VI,  X=Xraw,Y=Ys_prop, modelPerf=ModelPerf, cutoff= 0)+theme_bw()
# #checks which particular species etc
# 
# flashlightObj_rf <- mrFlashlight(yhats_rf_prop, X=Xraw,Y=Ys_prop, response = "multi", mode='regression')
# 
# #Interpretation
# 
# aledata_spec <-light_profile(flashlightObj_rf, v = "Spectral_radius", type = "ale",  n_bins =50)
# 
# mrProfileplot(aledata_spec , sdthresh = 0)+
#   theme_bw()+
#   geom_rug(data=Xraw, aes(x =Spectral_radius), inherit.aes = F)
# 
# 
# aledata_Fiedler <-light_profile(flashlightObj_rf, v = "Fiedler", type = "ale",  n_bins =50)
# 
# mrProfileplot(aledata_Fiedler , sdthresh = 0)+
#   theme_bw()+
#   geom_rug(data=Xraw, aes(x =Fiedler), inherit.aes = F)
# 
# aledata_Cent <-light_profile(flashlightObj_rf, v = "Centrality", type = "ale",  n_bins =50)
# 
# mrProfileplot(aledata_Cent , sdthresh = 0)+
#   theme_bw()+
#   geom_rug(data=X, aes(x =Centrality), inherit.aes = F)
# 
# aledata_Mod <-light_profile(flashlightObj_rf, v = "Qrel", type = "ale",  n_bins =50)
# 
# mrProfileplot(aledata_Mod , sdthresh = 0)+
#   theme_bw()+
#   geom_rug(data=X, aes(x =Modularity), inherit.aes = F)
# 
# #nteractions
# 
# #cl <- parallel::makeCluster(5)
# #plan(cluster, workers=cl)
# 
# #seemingly only work after running MrIML predicts
# interactions <-mrInteractions(yhats_rf_prop, Y=Ys_prop,X=X,  mode='regression') #this is computationally intensive so multicores are needed.
# 
# 
# mrPlot_interactions(interactions, Y=Ys_prop,X=X, top_ranking = 5, top_response=10)
# 
# mrIMLconverts_list <- MrIMLconverts(yhats_rf_prop,X=X, mode='regression')
# 
# featureA = 'Centrality' #doesn't do catergories
# featureB = 'Fiedler'
# 
# p2d <- mrProfile2D (mrIMLconverts_list, featureA,
#                     featureB,  mode='regression',
#                     grid.size=30, method = 'ale')
# plot(p2d) +theme_bw()
# 
# str(X)
# 
# #plot prediction scatter for all responses.
# 
# plot(light_profile(flashlightObj, v = "Modularity", type = "ale",  n_bins =25))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x =Modularity), inherit.aes = F)
# 
# plot(light_profile(flashlightObj, v = "network_size", type = "ale", n_bins =100))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x =network_size), inherit.aes = F)
# 
# plot(light_profile(flashlightObj, v = "Centrality", type = "ale",  n_bins =50))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x =Centrality), inherit.aes = F)
# 
# plot(light_profile(flashlightObj, v = "Fiedler", type = "ale", n_bins =100))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x=Fiedler), inherit.aes = F)
# 
# plot(light_ice(flashlightObj, v = "Fiedler", center = "first"))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x =Fiedler), inherit.aes = F)
# 
# plot(light_ice(flashlightObj, v = "Centrality", center = "first"))+
#   theme_bw()+
#   geom_rug(data=Y, aes(x = Centrality), inherit.aes = F)
# 
# plot(light_profile2d(flashlightObj, v = c("network_size", "Centrality")))
# 
# plot(light_profile2d(flashlightObj, v = c("network_size", "Mean_degree")))
# geom_rug(data=Y, aes(x = network_size, y=Mean_degree), inherit.aes = F)+
#   scale_fill_continuous(high = "#56B1F7", low = "#132B43")
# 
# plot(light_profile2d(flashlightObj, v = c("network_size", "Fiedler")))+
#   geom_rug(data=Y, aes(x = network_size, y=Fiedler), inherit.aes = F)+
#   #scale_fill_gradientn(colours = colorspace::heat_hcl(7))
#   scale_fill_continuous(high = "#132B43", low = "#56B1F7")
# 
# plot(light_global_surrogate(flashlightObj))
# 
# st <- light_interaction(flashlightObj, grid_size = 30, n_max = 50, seed = 42, pairwise=T)+
#   theme_bw()
# plot(st)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# End of Machine learning Disease Simulation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 # 
 # shap_plot_func <- function(trained_model = yhats_rf, 
 #                            feature = "Spectral_radius",
 #                            weighted = TRUE,
 #                            specific_type = "SW_unweighted_graph") {
 #   
 #   # Load the appropriate dataset based on `weight_type_data`
 #   if (weighted == T) {
 #     all_weighted_nets=readRDS("all_disease_simulation_data_weighted.rds")
 #     weighted_data=readRDS("weighted_data.rds")
 #     weight_type_data <- weighted_data
 #     network_type = factor(all_weighted_nets$Network)  # Reorder levels
 #     res <- weighted_data$Ys_prop$avg_prop_infected_beta0.05_gamma0.04
 #   }else{
 #     all_unweighted_nets=readRDS("all_disease_simulation_data_unweighted.rds")
 #     unweighted_data=readRDS("unweighted_data.rds")
 #     weight_type_data = unweighted_data
 #     network_type = factor(all_unweighted_nets$Network)  # Reorder levels
 #     res <- weighted_data$Ys_prop$avg_prop_infected_beta0.05_gamma0.04
 #   }
 #   
 #   Xfeat <- weight_type_data$Xraw
 #   network_type <-network_type  # Reorder levels
 #   res=res
 #   
 #   # Model selection
 #   model <- trained_model[1:6]  # proportioninfected responses only
 #   # model <- trained_model[7:12]  # Uncomment for invasion time responses only
 #   
 #   # Initialize storage for results
 #   shapley_results <- list()
 #   feature_effects <- list()
 #   
 #   # Loop through each network type
 #   for (net in levels(network_type)) {
 #     # Subset data for the current network type
 #     subset_indices <- which(network_type == net)
 #     Xfeat_subset <- Xfeat[subset_indices, , drop = FALSE]
 #     res_subset <- res[subset_indices]
 #     
 #     # Predictor setup
 #     if (inherits(model[[1]]$mod1_k$fit$fit$fit, "xgb.Booster")) {
 #       preds=Predictor$new(
 #         model = model[[1]]$mod1_k$fit$fit$fit,
 #         data = Xfeat_subset,
 #         predict.function = function(model, newdata) {
 #           newdata <- as.matrix(newdata)  # Ensure matrix format for xgboost
 #           predict(model, newdata)
 #         }
 #       )
 #     }else{
 #       preds=Predictor$new(
 #         model = model[[1]]$mod1_k$fit$fit$fit,
 #         data = Xfeat_subset,
 #         y = res_subset
 #       )
 #     }
 #     
 #     # Calculate Shapley values
 #     shapley <- Shapley$new(preds, x.interest = Xfeat_subset, sample.size = 50)
 #     shapley_results[[net]] <- shapley
 #     
 #     # Feature effects calculations
 #     cat("Calculating feature effects for network:", net, "\n")
 #     
 #     feature_effects[[net]] <- list(
 #       ale_feature = FeatureEffect$new(preds, 
 #                                       feature = feature, grid.size = 10),
 #       interaction_all = Interaction$new(preds, grid.size = 15),
 #       interaction_feature = Interaction$new(preds, 
 #                                             feature = feature, grid.size = 15),
 #       combined_effects = FeatureEffects$new(preds, grid.size = 10),
 #       tree = TreeSurrogate$new(preds, maxdepth = 2)
 #     )
 #   }
 #   
 #   # Test examples
 #   if (!is.null(shapley_results[[specific_type]])) {
 #     cat("Plotting results for specific network type:", specific_type, "\n")
 #     shapley_results[[specific_type]]$plot()
 #     feature_effects[[specific_type]]$ale_feature$plot()
 #     plot(feature_effects[[specific_type]]$combined_effects)
 #     plot(feature_effects[[specific_type]]$interaction_all)
 #     plot(feature_effects[[specific_type]]$interaction_feature)
 #     plot(feature_effects[[specific_type]]$tree)
 #   } else {
 #     warning("Specific network type not found in shapley_results.")
 #   }
 # }
 # 
 # 
 # xt=shap_plot_func(
 #   trained_model = yhats_rf,
 #   feature = "Spectral_radius",
 #   weighted = F,
 #   specific_type = "SP_unweighted_graph"
 # )
 # 
 # #easier to see with plots
 # plots <- mrPerformancePlot(ModelPerf1=ModelPerf_lm_weighted,
 #                            ModelPerf2 = ModelPerf_xgb_weighted,
 #                            mod_names=c('linear_reg',
 #                                        'boost_tree'),
 #                            mode='regression' ) 
 # 
 # plots[[1]]
 # 
 # plots[[2]]
 # 
 # ##---Variable Importance
 # VI=mrvip(yhats=yhats_xgb_weighted[1:6],
 #          X1=NULL,
 #          X=weighted_data$Xraw, 
 #          Y=weighted_data$Ys[1:6],
 #          ModelPerf=ModelPerf_xgb_weighted,
 #          mode="regression")
 # 
 # # VI_unweighted=mrvip(yhats=yhats_rf[1:6],
 # #                     X1=NULL,
 # #                     X=unweighted_data$Xraw, 
 # #                     Y=unweighted_data$Ys[1:6],
 # #                     ModelPerf=ModelPerf_rf,
 # #                     mode="regression")
 # # 
 # VI[[3]] #Importance plot
 # 
 # 
 # VI[[4]] #PCA 
 # 
 # 
 # flashlightObj <- mrFlashlight(yhats=yhats_xgb_weighted,
 #                               X=Xraw_weighted,
 #                               Y=Ys_prop_weighted,
 #                               response = "multi",
 #                               mode="regression")
 # 
 # profileData_pd <- light_profile(flashlightObj,
 #                                 v = "Spectral_radius") #partial dependencies
 # 
 # mrProfileplot(profileData_pd,
 #               sdthresh =0.001)
 # 
 # 
 # #>  Press [enter] to continue to the global summary plot
 # #> `geom_smooth()` using formula = 'y ~ x'
 # 
 # 
 # profileData_ale <- light_profile(flashlightObj,
 #                                  v = "Fiedler",
 #                                  type = "ale") #accumulated local effects
 # 
 # mrProfileplot(profileData_ale,
 #               sdthresh =0.01)
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # x=readRDS("network_results_size_nsim100_50.rds")
 # colnames(x)
 # y= x%>%
 #   dplyr::count(levels(as.factor(Network)))
 # 
 # 
 # head(x)
 # 
 # 
 # ###---Pearsons correlation----
 # cor(all_unweighted_nets$Spectral_radius,
 #     all_unweighted_nets$avg_invasion_time_beta0.05_gamma0.04,
 #     method = "pearson")
 # 
 # 
 # ###---Multiple linear regression---
 # lm_multiple <- lm(avg_invasion_time_beta0.05_gamma0.4 ~ Spectral_radius + Fiedler + Modularity + Degree_assortativity, 
 #                   data = all_unweighted_nets)
 # summary(lm_multiple)
 # 
 # ###---Stepwise selection----
 # stepwise_model <- step(lm(avg_invasion_time_beta0.05_gamma0.4 ~ 1, 
 #                           data = all_unweighted_nets), 
 #                        scope = list(lower = ~1, upper = ~ Spectral_radius + Fiedler + Modularity + Degree_assortativity),
 #                        direction = "forward")
 # summary(stepwise_model)
 # 
 # ###---Partial correlation----
 # library(ppcor)
 # pcor_result <- pcor(cbind(x$Spectral_radius, x$avg_invasion_time_beta0.025_gamma0.4, x$Fiedler))
 # pcor_result$estimate
 # 
 # ###---PCA----
 # pca_result <- prcomp(x[, c("Spectral_radius", "Fiedler", "Modularity", "Degree_assortativity")], scale = TRUE)
 # summary(pca_result)
 # 
 
 # 
 # st=yhats_xgb[1:6]#proportioninfected resonses only
 # class(st)
 # 
 # Xfeat=unweighted_data$Xraw
 # class(Xfeat)
 # res=unweighted_data$Ys$avg_prop_infected_beta0.05_gamma0.04
 # class(res)
 # 
 # preds<- Predictor$new(st[[1]]$mod1_k$fit$fit$fit, data = Xfeat,
 #                       y =unweighted_data$Ys_prop$avg_prop_infected_beta0.05_gamma0.04)
 # 
 # 
 # shapley <- Shapley$new(preds, 
 #                        x.interest = Xfeat,
 #                        sample.size = 50)
 # shapley$plot()
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # # Example data
 # Xfeat <- unweighted_data$Xraw
 # res <- unweighted_data$Ys$avg_prop_infected_beta0.05_gamma0.04
 # network_type <- factor(all_unweighted_nets$Network) # Assume this contains network types
 # 
 # # Check the structure
 # levels(network_type)
 # 
 # # Initialize a list to store Shapley results
 # shapley_results <- list()
 # 
 # # Loop through each network type
 # for (net in levels(network_type)) {
 #   # Subset data for the current network type
 #   subset_indices <- which(network_type == net)
 #   Xfeat_subset <- Xfeat[subset_indices, , drop = FALSE]
 #   res_subset <- res[subset_indices]
 #   
 #   # Train the model for the subset (modify the model as needed)
 #   # Here, we'll assume the model is pre-trained and directly usable
 #   preds <- Predictor$new(
 #     model = st[[1]]$mod1_k$fit$fit$fit, 
 #     data = Xfeat_subset,
 #     y = res_subset
 #   )
 #   
 #   # Calculate Shapley values
 #   shapley <- Shapley$new(preds, x.interest = Xfeat_subset[1, ], 
 #                          sample.size = 100)
 #   
 #   # Store results
 #   shapley_results[[net]] <- shapley
 # }
 # 
 # 
 # 
 # shapley_results$SW_unweighted_graph$plot()
 # 
 # 
 # 
 # # Example 2 data
 # Xfeat <- unweighted_data$Xraw
 # res <- unweighted_data$Ys_prop$avg_prop_infected_beta0.25_gamma0.04
 # network_type <- factor(all_unweighted_nets$Network) # Assume this contains network types
 # 
 # 
 # # Correct the order of levels in network_type to match the desired order
 # network_order <- c("ER_unweighted_graph", "SP_unweighted_graph", "SBM_unweighted_graph", 
 #                    "SW_unweighted_graph", "SF_unweighted_graph", "DS_unweighted_graph")
 # 
 # network_type <- factor(all_unweighted_nets$Network, levels = network_order)  # Reorder levels
 # 
 # # Check the structure of the updated levels
 # print(levels(network_type))
 # 
 # # Initialize a list to store Shapley results
 # shapley_results <- list()
 # 
 # # Loop through each network type
 # for (net in levels(network_type)) {
 #   # Subset data for the current network type
 #   subset_indices <- which(network_type == net)
 #   Xfeat_subset <- Xfeat[subset_indices, , drop = FALSE]
 #   res_subset <- res[subset_indices]
 #   
 #   # Train the model for the subset (modify the model as needed)
 #   # Here, we'll assume the model is pre-trained and directly usable
 #   preds <- Predictor$new(
 #     model = st[[1]]$mod1_k$fit$fit$fit, 
 #     data = Xfeat_subset,
 #     y = res_subset
 #   )
 #   
 #   # Calculate Shapley values using all data for the subsetted network
 #   shapley <- Shapley$new(preds, 
 #                          x.interest = Xfeat_subset[1, ], 
 #                          sample.size = 100)
 #   
 #   # Store results
 #   shapley_results[[net]] <- shapley
 # }
 # 
 # 
 # shapley_results$DS_unweighted_graph$plot()
 # 
 # 
 # 
 # ##---ALE PLOT
 # 
 # ale <- FeatureEffect$new(preds, 
 #                          feature = "Spectral_radius", 
 #                          grid.size = 10)
 # ale$plot() 
 # 
 # 
 # ale$set.feature("Modularity")
 # ale$plot() 
 # 
 # 
 # 
 # interact <- Interaction$new(preds, grid.size = 15) 
 # plot(interact)
 # 
 # interact <- Interaction$new(preds, 
 #                             feature = "Inverse_of_spectral_radius", grid.size = 15)
 # plot(interact)
 # 
 # 
 # effs <- FeatureEffects$new(preds, grid.size = 10)
 # plot(effs)
 # 
 # 
 # tree <- TreeSurrogate$new(preds, maxdepth = 2)
 # 
 # plot(tree)
 # 
 # #We can use the tree to make predictions:
 # 
 # head(tree$predict(Boston))
 # 
 # 
 # 
 # interactions <-mrInteractions(yhats_rf, unweighted_data$Xraw,
 #                               unweighted_data$Ys_prop,
 #                               feature = "Spectral_radius",
 #                               top.int = 10) #this is computationally intensive so multicores are needed.
 # 
 # mrPlot_interactions(interactions, unweighted_data$Xraw,
 #                     unweighted_data$Ys_prop, top_ranking = 2,
 #                     top_response=2) #can increase the number of interactions/SNPs ('responses') shown  
 # 
 # 
 # 
 # 
 # Xfeat <- unweighted_data$Xraw
 # res <- unweighted_data$Ys_prop$avg_prop_infected_beta0.25_gamma0.04
 # network_type <- factor(all_unweighted_nets$Network) 
 # 
 # st=yhats_rf[1:6]#proportioninfected responses only
 # 
 # # Correct the order of levels in network_type to match the desired order
 # network_order <- c("ER_unweighted_graph", "SP_unweighted_graph", "SBM_unweighted_graph", 
 #                    "SW_unweighted_graph", "SF_unweighted_graph", "DS_unweighted_graph")
 # 
 # network_type <- factor(all_unweighted_nets$Network, levels = network_order)  # Reorder levels
 # 
 # # Check the structure of the updated levels
 # print(levels(network_type))
 # 
 # # Initialize a list to store Shapley results and feature effects
 # shapley_results <- list()
 # feature_effects <- list()
 # 
 # # Loop through each network type
 # for (net in levels(network_type)) {
 #   # Subset data for the current network type
 #   subset_indices <- which(network_type == net)
 #   Xfeat_subset <- Xfeat[subset_indices, , drop = FALSE]
 #   res_subset <- res[subset_indices]
 #   
 #   # Train RF the model for the subset (modify the model as needed)
 #   preds <- Predictor$new(
 #     model = st[[1]]$mod1_k$fit$fit$fit,
 #     data = Xfeat_subset,
 #     y = res_subset
 #   )
 #   
 #   #----chnge when using xgboost 
 #   # preds <- Predictor$new(
 #   #   model = st[[1]]$mod1_k$fit$fit$fit,
 #   #   data = Xfeat_subset,
 #   #   predict.function = function(model, newdata) {
 #   #     newdata <- as.matrix(newdata)  # Convert to matrix for xgboost
 #   #     predict(model, newdata)
 #   #   }
 #   # )
 #   
 #   
 #   # Calculate Shapley values using all data for the subsetted network
 #   shapley <- Shapley$new(preds, x.interest = Xfeat_subset, sample.size = 50)
 #   
 #   # Store Shapley results
 #   shapley_results[[net]] <- shapley
 #   
 #   # Feature effects for a specific network type
 #   cat("Calculating feature effects for network:", net, "\n")
 #   
 #   # paste(cat("Calculating feature effects for',
 #   #           feature:,
 #   #           net, "\n"))
 #   
 #   # Example: Accumulated Local Effects (ALE) for a specific feature
 #   ale_feature <- FeatureEffect$new(preds, 
 #                                    feature = feature, grid.size = 10)
 #   
 #   # Feature effects for a specific network type
 #   cat("Calculating feature effects for all features:",
 #       net, "\n")
 #   
 #   feature_effects[[net]]$ale_feature <- ale_feature
 #   
 #   # # Example: Change the feature for ALE and plot
 #   #ale_modularity <- FeatureEffect$new(preds, feature = "Modularity", grid.size = 10)
 #   # feature_effects[[net]]$ale_modularity <- ale_modularity
 #   
 #   # Example: Interaction effects for all features
 #   interaction_all <- Interaction$new(preds, grid.size = 15)
 #   feature_effects[[net]]$interaction_all <- interaction_all
 #   
 #   # Example: Interaction effects for a specific feature
 #   interaction_feature <- 
 #     Interaction$new(preds, 
 #                     feature = feature, grid.size = 15)
 #   
 #   feature_effects[[net]]$interaction_feature <- interaction_feature
 #   # Example: Combined feature effects
 #   effs <- FeatureEffects$new(preds, grid.size = 10)
 #   feature_effects[[net]]$combined_effects <- effs
 #   
 #   # Example: Tree surrogate for understanding interactions
 #   tree <- TreeSurrogate$new(preds, maxdepth = 2)
 #   feature_effects[[net]]$tree <- tree
 # }
 # 
 # 
 # ## Test examples
 # feature="Spectral_radius"
 # # Example usage: Visualize results for a specific network type
 # specific_type <- "ER_unweighted_graph"  # Replace with the desired network type
 # 
 # # Shapley results plot
 # shapley_results[[specific_type]]$plot()
 # 
 # # ALE plots for specific features
 # feature_effects[[specific_type]]$ale_feature$plot()
 # 
 # # Combined feature effects
 # plot(feature_effects[[specific_type]]$combined_effects)
 # 
 # # Interaction plots
 # plot(feature_effects[[specific_type]]$interaction_all)
 # plot(feature_effects[[specific_type]]$interaction_feature)
 # # Tree surrogate plot
 # plot(feature_effects[[specific_type]]$tree)
 # 
 # # # Example: Use tree surrogate to make predictions (optional)
 # # tree_surrogate <- feature_effects[[specific_type]]$tree
 # # predictions <- tree_surrogate$predict(Xfeat_subset)
 # # print(predictions)