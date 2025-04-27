library(BoolNet)
library(gtools)
library(igraph)



# Generate adj matrix:
generate_matrices <- function(n) {
  values <- c(-1, 0, 1)
  all_combinations <- expand.grid(rep(list(values), n * n))
  matrices <- list()
  
  for (i in seq_len(nrow(all_combinations))) {
    mat <- matrix(as.numeric(all_combinations[i, ]), n, n)
    matrices <- c(matrices, list(mat))
  }
  
  cat("Found", length(matrices), "matrices.\n")
  return(matrices)
}

generate_matrices(2)

generate_unique_matrices <- function(n) {
  matrices <- generate_matrices(n)
  
  are_isomorphic <- function(A, B) {
    g1 <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
    g2 <- graph_from_adjacency_matrix(B, mode = "directed", weighted = TRUE)
    return(graph.isomorphic.vf2(g1, g2, edge.color1 = E(g1)$weight, edge.color2 = E(g2)$weight)[[1]] || all ((A == t(B))))
  }
  
  unique_graphs <- list()
  unique_graphs[[1]] <- matrices[[1]]
  
  for (i in 2:length(matrices)) {
    is_duplicate <- FALSE
    for (j in 1:length(unique_graphs)) {
      if (are_isomorphic(matrices[[i]], unique_graphs[[j]])) {
        is_duplicate <- TRUE
        break
      }
    }
    if (!is_duplicate) {
      unique_graphs[[length(unique_graphs) + 1]] <- matrices[[i]]
    }
  }
  
  cat("Found", length(unique_graphs), "unique adjacency matrices.\n")
  return(unique_graphs)
}

generate_unique_matrices(2)
#########################

# Generate all possible rules:
generate_all_boolean_rule_networks <- function(adj_matrix) {
  nodes <- LETTERS[1:nrow(adj_matrix)]
  
  # Build possible expressions for each node
  expr_options <- list()
  for (i in seq_along(nodes)) {
    node <- nodes[i]
    inputs_idx <- which(adj_matrix[, i] != 0)
    
    if (length(inputs_idx) == 0) {
      expr_options[[node]] <- "0"
      
    } else if (length(inputs_idx) == 1) {
      input_node <- nodes[inputs_idx]
      sign <- if (adj_matrix[inputs_idx, i] == -1) "!" else ""
      expr_options[[node]] <- paste0(sign, input_node)
      
    } else {
      n_inputs <- length(inputs_idx)
      n_combos <- 2^(n_inputs - 1)
      
      operators <- expand.grid(rep(list(c("&", "|")), n_inputs - 1))
      operators <- as.data.frame(operators, stringsAsFactors = FALSE)
      
      rules <- character(n_combos)
      
      for (j in 1:n_combos) {
        expr_parts <- c()
        for (k in 1:n_inputs) {
          sign <- if (adj_matrix[inputs_idx[k], i] == -1) "!" else ""
          expr_parts <- c(expr_parts, paste0(sign, nodes[inputs_idx[k]]))
          if (k < n_inputs) {
            expr_parts <- c(expr_parts, operators[j, k])
          }
        }
        rules[j] <- paste(expr_parts, collapse = " ")
      }
      
      expr_options[[node]] <- rules
    }
  }
  
  # Now, combine all possible networks
  network_combos <- expand.grid(expr_options, stringsAsFactors = FALSE)
  n_networks <- nrow(network_combos)
  
  # Format result
  results <- list()
  for (i in 1:n_networks) {
    network <- network_combos[i, ]
    lines <- character(length(nodes))
    for (j in seq_along(nodes)) {
      lines[j] <- paste0(nodes[j], ", ", network[[j]])
    }
    results[[paste0("network", i)]] <- lines
  }
  
  return(results)
}

#########################
# Remove duplicate rules:

remove_similar_networks <- function(networks) {
  unique_signatures <- list()
  unique_networks <- list()
  
  for (name in names(networks)) {
    rules <- sapply(networks[[name]], function(x) sub("^[A-Z], ", "", x))  # remove "A, " part
    signature <- paste(sort(rules), collapse = " | ")  # sorted canonical form
    
    if (!(signature %in% unique_signatures)) {
      unique_signatures <- c(unique_signatures, signature)
      unique_networks[[name]] <- networks[[name]]
    }
  }
  
  return(unique_networks)
}

# My Test:
adj_matrix = matrices[[33]]
gg = generate_all_boolean_rule_networks(adj_matrix)
gg
ll = remove_similar_networks(gg)
ll

#########################
## saving part:
save_networks_to_boolnet <- function(networks, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  
  for (i in seq_along(networks)) {
    net <- networks[[i]]
    file_name <- file.path(output_dir, paste0("network_", i, ".bn"))
    
    # BoolNet requires this header
    lines <- c("targets, factors")
    lines <- c(lines, net)
    
    writeLines(lines, file_name)
  }
  
  cat("Saved", length(networks), "networks to", output_dir, "\n")
} 

#########################
# Final Function:
generate_all_possible_rules <- function(n, output_dir) {
  matrices <- generate_unique_matrices(n)
  all_unique_networks <- list()
  network_counter <- 1
  
  for (mat in matrices) {
    colnames(mat) <- rownames(mat) <- LETTERS[1:n]
    
    all_networks <- generate_all_boolean_rule_networks(mat)
    unique_networks <- remove_similar_networks(all_networks)
    
    for (net in unique_networks) {
      all_unique_networks[[paste0("network_", network_counter)]] <- net
      network_counter <- network_counter + 1
    }
  }
  
  # Save to files
  save_networks_to_boolnet(all_unique_networks, output_dir)
  
  cat("Total unique BoolNet networks generated:", length(all_unique_networks), "\n")
  
  return(all_unique_networks)
}


result <- generate_all_possible_rules( n= 2, output_dir = "boolnet_networks")
