set.seed(178)

library(dplyr)

# Simulate data according to chosen effect_size
sim_functions <- list(
  Exponential = function(effect_size) exp_sim(effect_size),
  Uniform = function(effect_size) unif_sim(effect_size),
  Beta = function(effect_size) beta_sim(effect_size)
)

exp_sim <- function(effect_size = 0){
  block <- as.factor(rep(1:5, each = 3)) # 5 blocks, which each get each treatment exactly once as per RCBD
  treatment <- rep(1:3, times = 5) # Each block gets each treatment once
  block_effects <- rnorm(5, mean = 0, sd = 0.25)
  block_effect <- block_effects[as.numeric(block)]
  obs <- rexp(15, rate = .8) + effect_size * (treatment - 1) + block_effect
  
  data.frame(block, treatment, obs)
}

unif_sim <- function(effect_size = 0, min = 0, max = 2){
  block <- as.factor(rep(1:5, each = 3)) # 5 blocks, which each get each treatment exactly once as per RCBD
  treatment <- rep(1:3, times = 5) # Each block gets each treatment once
  block_effects <- rnorm(5, mean = 0, sd = 0.25)
  block_effect <- block_effects[as.numeric(block)]
  obs <- runif(15, min = min, max = max) + effect_size * (treatment - 1) + block_effect

  data.frame(block, treatment, obs)
}

beta_sim <- function(effect_size = 0, shape1 = 1, shape2 = .5){
  block <- as.factor(rep(1:5, each = 3)) # 5 blocks, which each get each treatment exactly once as per RCBD
  treatment <- rep(1:3, times = 5) # Each block gets each treatment once
  block_effects <- rnorm(5, mean = 0, sd = 0.25)
  block_effect <- block_effects[as.numeric(block)]
  obs <- rbeta(15, shape1 = shape1, shape2 = shape2) + effect_size * (treatment - 1) + block_effect
  sd_blocks <- c(1, 5, 2.5, 1, 5) # Introduce unequal variances between observations
  obs <- obs * sd_blocks[as.numeric(block)]
  
  
  data.frame(block, treatment, obs)
}

anova_test <- function(data){
  summary(aov(obs ~ factor(treatment) + block, data = data))[[1]][1,5]
}

perm_test <- function(data, reps = 10000) {
  f_initial <- summary(aov(obs ~ factor(treatment) + block, data = data))[[1]][1, 4]
  
  perm_f <- replicate(reps, {
    perm_data <- data
    perm_data$obs <- ave(data$obs, data$block, FUN = sample)
    summary(aov(obs ~ factor(treatment) + block, data = perm_data))[[1]][1, 4]
  })
  mean(perm_f >= f_initial)
}

power_test <- function(iter = 1000, dist, effect_size = 0, test_type = 'anova'){
  anova_pvals <- numeric(iter)
  perm_pvals <- numeric(iter)
  
  for(i in 1:iter){
    data <- sim_functions[[dist]](effect_size)
    anova_pvals[i] <- anova_test(data)
    perm_pvals[i] <- perm_test(data, reps = 10000)
  }
  if(test_type == 'anova'){
    sum(anova_pvals < 0.05) / iter
  } else {
    sum(perm_pvals < 0.05) / iter
  }
}

results_list <- list()
counter <- 1
total_combinations <- length(c(0, .25, .5, .75)) * length(c('Exponential', 'Uniform', 'Beta')) * length(c('anova', 'permutation'))
save_interval <- ceiling(total_combinations * 0.01)

for(effect_size in c(0, .25, .5, .75)){
  for(dist in c('Exponential', 'Uniform', 'Beta')){
    for(test in c('anova', 'permutation')){
      power <- power_test(iter = 1000, dist, effect_size, test)
      
      results_list[[counter]] <- data.frame(
        distribution = dist,
        effect_size = effect_size, 
        test_type = test,
        power = power
      )
      counter <- counter + 1
      cat("Progress:", round(counter / total_combinations * 100, 2), "%\n")
      if (counter %% save_interval == 0) {
        saveRDS(results_list, file = "simulation_results_checkpoint.rds")
        cat("Checkpoint saved at iteration:", counter, "\n")
        }
    }
  }
}

results <- do.call(rbind, results_list)
save(results, file = "simulation_results.RData")
