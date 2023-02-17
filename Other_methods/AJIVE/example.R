source("ajive.R")


blocks <- sample_toy_data(n=200, dx=100, dy=500)
scree_plot_blocks(blocks)
initial_signal_ranks <- c(2, 2)
ajive_output <- ajive(blocks, initial_signal_ranks,
                      n_wedin_samples = 100, n_rand_dir_samples = 100)



# Full matrix representation

k <- 1
J <- get_block_full(ajive_output, k, type='joint')
I <- get_block_full(ajive_output, k, type='individual')
E <- get_block_full(ajive_output, k, type='noise')
# check X = J + I + E
norm(blocks[[k]] - (J + I + E))




# Block Specific Scores


bss_1 <- get_block_scores(ajive_output, k, type='joint', normalized=FALSE)

normalized_bss_1 <- get_block_scores(ajive_output, k, type='joint', normalized=TRUE)





## Common Normalized Scores


cns <- get_common_normalized_scores(ajive_output)



## Individual Normalized Scores


ins_1 <- get_block_scores(ajive_output, k, type='individual', normalized=TRUE)


## Individual loadings

indiv_loadings_1 <- get_block_loadings(ajive_output, k, type='individual')

