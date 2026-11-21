library(future)
library(future.apply)

# Set up parallel execution
plan(multisession, workers = parallel::detectCores() - 1)

start <- Sys.time()
sims <- 200
years <- 1500  # Shared by simulation, error handling, and plotting

run_simulation <- function(sim) {
  tryCatch({
    
    # ---------------------------
    # Parameters, comment out for parameter testing or base value
    # ---------------------------
    dens   <- 0.9
    #   dens = runif(1, min = 0, max = 1) 
    MC     <- 0.3
    #   MC = runif(1, min = 0, max = 2)
    HC = 0.4
    #   HC = runif(1, min = 0, max = 2)
    #    thresh = -20    # For sedentary simulations
    # Only amendment: each family has an individual heritable migration threshold.
    thresh_min <- -5
    thresh_max <- 5
    thresh_mut_sd <- 0.05
    peak = 35   
    #   peak = runif(1, min = 25, max = 45)     
    RepT = 13.5  
    #   RepT = runif(1, min = 10, max = 14) 
    k = 4 
    #   k = runif(1, min = 0.5, max = 10)   
    age = 15    
    #   age = runif(1, min = 13, max = 18)  
    
    
    years <- 1500
    days  <- years * 365
    DC    <- 0.5
    f     <- 15
    mut   <- 0.005
    a     <- 2.7 * 10^-7
    b     <- 2.4 * 10^-4
    
    # -----------------------
    ## Yearly results vectors
    # -----------------------
    results_n <- rep(NA_real_, years)
    results_h <- rep(NA_real_, years)
    results_f <- rep(NA_real_, years)
    results_energy <- rep(NA_real_, years)
    results_group_mean <- rep(NA_real_, years)
    results_group_var <- rep(NA_real_, years)
    results_migrations <- rep(NA_real_, years)
    results_mean_thresh <- rep(NA_real_, years)
    results_sd_thresh <- rep(NA_real_, years)
    results_prop_reproducing <- rep(NA_real_, years)
    
    
    # ---------------------------
    # Probability of game (p)
    # ---------------------------
    num_groups <- 20
    game_prob <- matrix(0, nrow = days, ncol = num_groups)
    
    for (grp in 1:num_groups) {
      for (day_block in seq(1, days, by = 5)) {
        random_prob <- runif(1, min = 0.3, max = 1)
        game_prob[day_block:min(day_block + 4, days), grp] <- random_prob
      }
    }
    
    
    # ---------------------------
    # Agent state as vectors
    # ---------------------------
    id <- 1:100
    next_id <- 101L
    thresh_vec <- runif(100, thresh_min, thresh_max)
    initial_mean_thresh <- mean(thresh_vec)
    group <- sample(1:10, 100, replace = TRUE)
    energy <- runif(100, min = (0.5 * f), max = (0.7 * f))
    age_vec <- runif(100, min = 15, max = 50)
    gene <- rep("H", 100)
    w <- runif(100, min = 0.75, max = 1)
    last_reproduction <- rep(-365, 100)
    offspring_ages <- rep(list(numeric(0)), 100)
    offspring_energy_sum <- rep(list(numeric(0)), 100)
    offspring_energy_days <- rep(list(integer(0)), 100)
    last_migration <- rep(0, 100)
    migration_count <- rep(0, 100)
    migration_history <- vector("list", 100)
    
    
    ## Matrix for storing net energy for migration threshold ##
    # net_energy1, net_energy2, net_energy3
    net_energy <- matrix(0, nrow = 100, ncol = 3)
    
    # Add some random offspring ages for a subset of agents
    num_with_offspring <- 75
    agents_with_offspring <- sample(1:100, num_with_offspring)
    
    for (i in agents_with_offspring) {
      num_offspring <- sample(1:2, 1)
      offspring_ages[[i]] <- runif(num_offspring, min = 1, max = 14.8)
      offspring_energy_sum[[i]] <- rep(0, num_offspring)
      offspring_energy_days[[i]] <- rep(0L, num_offspring)
    }
    
    
    # Vectors for storing results and analysing outputs
    age_deaths_by_day <- rep(0, days)
    energy_deaths_by_day <- rep(0, days)
    births_by_day <- rep(0, days)
    maturity_by_day <- rep(0, days)
    population_by_day <- rep(NA_real_, days)
    energy_death_ages <- numeric(0)
    age_death_ages <- numeric(0)
    
    
    calculate_avg_individual_migration_interval <- function(migration_history) {
      individual_intervals <- c()
      
      for (i in seq_along(migration_history)) {
        if (length(migration_history[[i]]) > 1) {
          individual_intervals <- c(
            individual_intervals,
            diff(sort(migration_history[[i]]))
          )
        }
      }
      
      if (length(individual_intervals) > 0) {
        mean(individual_intervals)
      } else {
        NA_real_
      }
    }
    
    mean_occupied_group_size <- function(group_sizes) {
      occupied <- group_sizes[group_sizes > 0]
      if (length(occupied) > 0) {
        mean(occupied)
      } else {
        NA_real_
      }
    }
    
    # ---------------------------
    # Simulation loop
    # ---------------------------
    for (day in 1:days) {
      
      # Order agents by group
      ord <- order(group)
      id <- id[ord]
      thresh_vec <- thresh_vec[ord]
      group <- group[ord]
      energy <- energy[ord]
      age_vec <- age_vec[ord]
      gene <- gene[ord]
      last_reproduction <- last_reproduction[ord]
      offspring_ages <- offspring_ages[ord]
      offspring_energy_sum <- offspring_energy_sum[ord]
      offspring_energy_days <- offspring_energy_days[ord]
      last_migration <- last_migration[ord]
      net_energy <- net_energy[ord, , drop = FALSE]
      migration_count <- migration_count[ord]
      migration_history <- migration_history[ord]
      w <- w[ord]
      
      
      
      # Daily process for each group
      for (grp in unique(group)) {
        
        p <- game_prob[day, grp]
        
        group_rows <- which(group == grp)
        group_size <- length(group_rows)
        
        count_hunters <- 0
        success_hunts <- 0
        
        hunters_local <- which(runif(group_size) < w[group_rows])
        hunter_rows <- group_rows[hunters_local]
        hunters_id <- id[hunter_rows]
        
        # Gather hunters ages for peak hunting ability calculation
        hunters_age <- age_vec[group_rows][id[group_rows] %in% hunters_id]
        
        energy[hunter_rows] <- energy[hunter_rows] - HC
        count_hunters <- count_hunters + length(hunters_id)
        
        # Calculate successful hunts
        if (length(hunters_id) > 0) {
          for (hunter in hunters_id) {
            
            hunter_age <- age_vec[id == hunter]
            
            huntab <- 0.9 * exp(-(((hunter_age - peak) * 365)^2) / (2 * 7300^2))
            success_prob <- p * (dens ^ (count_hunters - 1)) * huntab
            
            if (runif(1) < success_prob) {
              success_hunts <- success_hunts + 1
            }
          }
        }
        
        # Meat sharing from all successful hunts
        meat_share <- (k * success_hunts) / group_size
        row_indexes <- group_rows
        energy[row_indexes] <- energy[row_indexes] + meat_share
        hunter_indexes <- which(id %in% hunters_id)
        
        base_energy <- meat_share - DC
        hunter_energy <- base_energy - HC
        
        
        ## Calculate average net energy income for last three days, for migration threshold ##
        
        slot <- day %% 3
        
        if (slot == 1) {
          net_energy[row_indexes, 1] <- base_energy
          net_energy[hunter_indexes, 1] <- hunter_energy
        } else if (slot == 2) {
          net_energy[row_indexes, 2] <- base_energy
          net_energy[hunter_indexes, 2] <- hunter_energy
        } else {
          net_energy[row_indexes, 3] <- base_energy
          net_energy[hunter_indexes, 3] <- hunter_energy
        }
      }
      
      # All agents age by 1 day
      offspring_ages <- lapply(offspring_ages, function(ages) ages + (1 / 365))
      age_vec <- age_vec + (1 / 365)
      
      # Daily energy cost
      energy <- energy - DC
      energy[energy > f] <- f
      
      #---------
      ## Deaths
      # --------
      ## All agents with 0 energy die
      remove_agents <- which(energy <= 0)
      if (length(remove_agents) >= 1) {
        
        energy_death_ages <- c(
          energy_death_ages,
          age_vec[remove_agents])
        
        energy_deaths_by_day[day] <- length(remove_agents)
        
        id <- id[-remove_agents]
        thresh_vec <- thresh_vec[-remove_agents]
        group <- group[-remove_agents]
        energy <- energy[-remove_agents]
        age_vec <- age_vec[-remove_agents]
        gene <- gene[-remove_agents]
        w <- w[-remove_agents]
        last_reproduction <- last_reproduction[-remove_agents]
        offspring_ages <- offspring_ages[-remove_agents]
        offspring_energy_sum <- offspring_energy_sum[-remove_agents]
        offspring_energy_days <- offspring_energy_days[-remove_agents]
        last_migration <- last_migration[-remove_agents]
        net_energy <- net_energy[-remove_agents, , drop = FALSE]
        migration_count <- migration_count[-remove_agents]
        migration_history <- migration_history[-remove_agents]
      }
      
      
      ## Stop if collapsed
      if (length(id) == 0) {
        total_migrations <- sum(migration_count)
        avg_days_between_migrations <- calculate_avg_individual_migration_interval(migration_history)
        
        
        mean_age_all_deaths <- if (
          length(c(energy_death_ages, age_death_ages)) > 0
        ) {
          mean(c(energy_death_ages, age_death_ages))
        } else {
          NA_real_
        }
        
        mean_age_energy_death <- if (length(energy_death_ages) > 0) {
          mean(energy_death_ages)
        } else {
          NA_real_
        }
        
        mean_age_natural_death <- if (length(age_death_ages) > 0) {
          mean(age_death_ages)
        } else {
          NA_real_
        }
        
        mean_prop_reproducing <- if (
          any(!is.na(results_prop_reproducing))
        ) {
          mean(results_prop_reproducing, na.rm = TRUE)
        } else {
          NA_real_
        }
        
        return(list(summary = c(
          dens, MC, HC, initial_mean_thresh, peak, RepT, k, age,
          day,
          0,
          NA_real_, NA_real_,
          total_migrations,
          avg_days_between_migrations,
          NA_real_,
          NA_real_,
          NA_real_,
          NA_real_,
          NA_real_,
          mean(results_n, na.rm = TRUE),
          mean(results_h, na.rm = TRUE),
          mean(results_f, na.rm = TRUE),
          mean(results_energy, na.rm = TRUE),
          mean(results_group_var, na.rm = TRUE),
          mean(results_group_mean, na.rm = TRUE),
          mean(results_migrations, na.rm = TRUE),
          NA_real_,
          NA_real_,
          mean_age_all_deaths,
          mean_age_energy_death,
          mean_age_natural_death,
          mean_prop_reproducing,
          NA_real_,
          NA_real_
        ), thresh_by_year = results_mean_thresh))
      }
      
      ## Deaths due to age
      age_days <- age_vec * 365
      age_death <- 1 - exp(-(a / b) * (exp(b * (age_days + 1)) - exp(b * (age_days - 1))))
      deaths <- which(runif(length(age_death)) < age_death)
      
      if (length(deaths) >= 1) {
        
        age_death_ages <- c(
          age_death_ages,
          age_vec[deaths]
        )
        
        id <- id[-deaths]
        thresh_vec <- thresh_vec[-deaths]
        group <- group[-deaths]
        energy <- energy[-deaths]
        age_vec <- age_vec[-deaths]
        gene <- gene[-deaths]
        w <- w[-deaths]
        last_reproduction <- last_reproduction[-deaths]
        offspring_ages <- offspring_ages[-deaths]
        offspring_energy_sum <- offspring_energy_sum[-deaths]
        offspring_energy_days <- offspring_energy_days[-deaths]
        last_migration <- last_migration[-deaths]
        net_energy <- net_energy[-deaths, , drop = FALSE]
        age_deaths_by_day[day] <- age_deaths_by_day[day] + length(deaths)
        migration_count <- migration_count[-deaths]
        migration_history <- migration_history[-deaths]
      }
      
      # Stop if civilisation has collapsed
      if (length(id) == 0) {
        total_migrations <- sum(migration_count)
        avg_days_between_migrations <- calculate_avg_individual_migration_interval(migration_history)
        
        
        valid_days <- !is.na(population_by_day) & population_by_day > 0
        
        
        mean_age_all_deaths <- if (
          length(c(energy_death_ages, age_death_ages)) > 0
        ) {
          mean(c(energy_death_ages, age_death_ages))
        } else {
          NA_real_
        }
        
        mean_age_energy_death <- if (length(energy_death_ages) > 0) {
          mean(energy_death_ages)
        } else {
          NA_real_
        }
        
        mean_age_natural_death <- if (length(age_death_ages) > 0) {
          mean(age_death_ages)
        } else {
          NA_real_
        }
        
        mean_prop_reproducing <- if (
          any(!is.na(results_prop_reproducing))
        ) {
          mean(results_prop_reproducing, na.rm = TRUE)
        } else {
          NA_real_
        }
        
        return(list(summary = c(
          dens, MC, HC, initial_mean_thresh, peak, RepT, k, age,
          day,
          0,
          NA_real_, NA_real_,
          total_migrations,
          avg_days_between_migrations,
          NA_real_,
          NA_real_,
          NA_real_,
          NA_real_,
          NA_real_,
          mean(results_n, na.rm = TRUE),
          mean(results_h, na.rm = TRUE),
          mean(results_f, na.rm = TRUE),
          mean(results_energy, na.rm = TRUE),
          mean(results_group_var, na.rm = TRUE),
          mean(results_group_mean, na.rm = TRUE),
          mean(results_migrations, na.rm = TRUE),
          NA_real_,
          NA_real_,
          mean_age_all_deaths,
          mean_age_energy_death,
          mean_age_natural_death,
          mean_prop_reproducing,
          NA_real_,
          NA_real_
        ), thresh_by_year = results_mean_thresh))
      }
      
      # Record parental energy for every dependent offspring, after mortality.
      # Each offspring accumulates its own history from birth to independence.
      for (j in seq_along(id)) {
        if (length(offspring_ages[[j]]) > 0) {
          offspring_energy_sum[[j]] <- offspring_energy_sum[[j]] + energy[j]
          offspring_energy_days[[j]] <- offspring_energy_days[[j]] + 1L
        }
      }
      
      # ---------------------------
      # Reproduction
      # ---------------------------
      can_reproduce <- which(energy > RepT & (day - last_reproduction) > 365)
      
      if (length(can_reproduce) > 0) {
        for (i in can_reproduce) {
          offspring_ages[[i]] <- c(offspring_ages[[i]], 0)
          offspring_energy_sum[[i]] <- c(offspring_energy_sum[[i]], 0)
          offspring_energy_days[[i]] <- c(offspring_energy_days[[i]], 0L)
          last_reproduction[i] <- day
        }
        births_by_day[day] <- length(can_reproduce)
      }
      
      # ---------------------------
      # Mature offspring: only become independent if their parent's
      # mean energy during the observed dependency period meets RepT.
      # ---------------------------
      n_parents <- length(id)  # Do not process newly created agents today
      for (i in seq_len(n_parents)) {
        mature <- which(offspring_ages[[i]] >= age)
        if (length(mature) == 0) next
        
        for (m in mature) {
          observed_days <- offspring_energy_days[[i]][m]
          mean_parent_energy <- if (observed_days > 0) {
            offspring_energy_sum[[i]][m] / observed_days
          } else {
            NA_real_
          }
          if (is.na(mean_parent_energy) || mean_parent_energy < RepT) next
          
          offspring_gene <- if (gene[i] == "H") {
            if (runif(1) < mut) "F" else "H"
          } else "F"
          offspring_w <- if (offspring_gene == "H") {
            runif(1, 0.75, 1)
          } else {
            runif(1, 0, 0.1)
          }
          
          offspring_thresh <- thresh_vec[i] + rnorm(1, 0, thresh_mut_sd)
          offspring_thresh <- min(max(offspring_thresh, thresh_min), thresh_max)
          id <- c(id, next_id)
          thresh_vec <- c(thresh_vec, offspring_thresh)
          next_id <- next_id + 1L
          group <- c(group, group[i])
          energy <- c(energy, energy[i])
          age_vec <- c(age_vec, age)
          gene <- c(gene, offspring_gene)
          w <- c(w, offspring_w)
          last_reproduction <- c(last_reproduction, day - 365)
          offspring_ages[[length(offspring_ages) + 1L]] <- numeric(0)
          offspring_energy_sum[[length(offspring_energy_sum) + 1L]] <- numeric(0)
          offspring_energy_days[[length(offspring_energy_days) + 1L]] <- integer(0)
          last_migration <- c(last_migration, day)
          net_energy <- rbind(net_energy, c(0, 0, 0))
          migration_count <- c(migration_count, 0)
          migration_history[[length(migration_history) + 1L]] <- integer(0)
          maturity_by_day[day] <- maturity_by_day[day] + 1
        }
        # Remove ALL mature offspring, including those failing the threshold.
        offspring_ages[[i]] <- offspring_ages[[i]][-mature]
        offspring_energy_sum[[i]] <- offspring_energy_sum[[i]][-mature]
        offspring_energy_days[[i]] <- offspring_energy_days[[i]][-mature]
      }
      
      # ---------------------------
      # Agents migrate if their net energy income over past 3 days is below the threshold
      # ---------------------------
      for (i in seq_along(id)) {
        agent_id <- id[i]
        
        if (day - last_migration[i] >= 3) {
          avg_meat <- mean(as.numeric(net_energy[i, ]))
          
          if (avg_meat < thresh_vec[i]) {
            possible_groups <- setdiff(1:num_groups, group[i])
            if (length(possible_groups) > 0) {
              new_group <- sample(possible_groups, 1)
              group[i] <- new_group
              energy[i] <- energy[i] - MC
              last_migration[i] <- day
              
              migration_count[i] <- migration_count[i] + 1
              migration_history[[i]] <- c(migration_history[[i]], day)
            }
          }
        }
      }
      
      population_by_day[day] <- length(id)
      
      
      # ---------------------------
      # Store results yearly
      # ---------------------------
      if (day %% 365 == 0) {
        year <- day / 365
        group_sizes <- as.numeric(table(factor(group, levels = 1:num_groups)))
        
        # Mean independent population during this year
        year_start <- ((year - 1) * 365) + 1
        year_end   <- year * 365
        
        mean_pop_this_year <- mean(
          population_by_day[year_start:year_end],
          na.rm = TRUE
        )
        
        # Total independent agents reproducing during this year
        births_this_year <- sum(
          births_by_day[year_start:year_end],
          na.rm = TRUE
        )
        
        # Proportion reproducing
        results_prop_reproducing[year] <-
          births_this_year / mean_pop_this_year
        
        results_n[year] <- length(id)
        results_h[year] <- mean(gene == "H")
        results_f[year] <- mean(gene == "F")
        results_energy[year] <- mean(energy)
        results_mean_thresh[year] <- mean(thresh_vec)
        results_sd_thresh[year] <- if (length(thresh_vec) > 1) sd(thresh_vec) else NA_real_
        results_group_var[year] <- var(group_sizes)
        results_group_mean[year] <- mean_occupied_group_size(group_sizes)
        results_migrations[year] <- sum(migration_count)
        
        
      }
    }
    
    # ---------------------------
    # End of simulation
    # ---------------------------
    total_migrations <- sum(migration_count)
    
    avg_days_between_migrations <- calculate_avg_individual_migration_interval(migration_history)
    final_group_sizes <- as.numeric(table(factor(group, levels = 1:num_groups)))
    final_mean_group_size <- mean_occupied_group_size(final_group_sizes)
    
    
    
    
    mean_prop_reproducing <- mean(
      results_prop_reproducing,
      na.rm = TRUE
    )
    
    mean_age_all_deaths <- mean(
      c(energy_death_ages, age_death_ages),
      na.rm = TRUE
    )
    
    mean_age_energy_death <- mean(
      energy_death_ages,
      na.rm = TRUE
    )
    
    mean_age_natural_death <- mean(
      age_death_ages,
      na.rm = TRUE
    )
    
    mean_w <- mean(w)
    sd_w <- if (length(w) > 1) sd(w) else NA_real_
    
    return(list(summary = c(
      dens, MC, HC, initial_mean_thresh, peak, RepT, k, age,
      days,
      length(id),
      mean(gene == "H"),
      mean(gene == "F"),
      total_migrations,
      avg_days_between_migrations,
      mean(energy),
      var(final_group_sizes),
      final_mean_group_size,
      mean(thresh_vec),
      if (length(thresh_vec) > 1) sd(thresh_vec) else NA_real_,
      mean(results_n, na.rm = TRUE),
      mean(results_h, na.rm = TRUE),
      mean(results_f, na.rm = TRUE),
      mean(results_energy, na.rm = TRUE),
      mean(results_group_var, na.rm = TRUE),
      mean(results_group_mean, na.rm = TRUE),
      mean(results_migrations, na.rm = TRUE),
      mean(results_mean_thresh, na.rm = TRUE),
      mean(results_sd_thresh, na.rm = TRUE),
      mean_age_all_deaths,
      mean_age_energy_death,
      mean_age_natural_death,
      mean_prop_reproducing,
      mean_w,
      sd_w
    ), thresh_by_year = results_mean_thresh))
    
  }, error = function(e) {
    print(e)
    return(list(summary = rep(NA_real_, 34), thresh_by_year = rep(NA_real_, years)))
  })
}


# Run simulations in parallel
results <- future_lapply(1:sims, run_simulation, future.seed = TRUE)

data <- do.call(rbind, lapply(results, `[[`, "summary"))
data <- as.data.frame(data)

colnames(data) <- c(
  "dens_val",
  "MC_val",
  "HC_val",
  "thresh_val",
  "peak_val",
  "RepT_val",
  "k_val",
  "age_val",
  "time_to_collapse",
  "final_population_size",
  "final_prop_hunters",
  "final_prop_free_riders",
  "total_migrations",
  "avg_days_between_migrations",
  "final_mean_energy",
  "final_group_size_variance",
  "final_mean_occupied_group_size",
  "final_mean_thresh",
  "final_sd_thresh",
  "mean_population_size",
  "mean_prop_hunters",
  "mean_prop_free_riders",
  "mean_energy",
  "mean_group_size_variance",
  "mean_occupied_group_size",
  "mean_total_migrations",
  "mean_thresh_over_time",
  "mean_sd_thresh_over_time",
  "mean_age_all_deaths",
  "mean_age_energy_death",
  "mean_age_natural_death",
  "mean_prop_reproducing",
  "mean_w",
  "sd_w"
)

# ------------------------------------------------------------
# Annual mean mobility threshold for every simulation
# Grey: individual simulations; black: mean across simulations
# ------------------------------------------------------------
threshold_plot_data <- do.call(rbind, lapply(seq_along(results), function(i) {
  x <- results[[i]]$thresh_by_year
  data.frame(sim = i, year = seq_along(x), mean_thresh = x)
}))

# Initial population thresholds are drawn from U(-5, 5).
# Record each simulation's actual initial mean, not an assumed zero.
threshold_initial <- data.frame(
  sim = seq_along(results),
  year = 0,
  mean_thresh = data$thresh_val
)
threshold_plot_data <- rbind(threshold_initial, threshold_plot_data)
threshold_plot_data <- threshold_plot_data[!is.na(threshold_plot_data$mean_thresh), ]

threshold_summary <- aggregate(mean_thresh ~ year, threshold_plot_data, mean)
n_sims_by_year <- aggregate(mean_thresh ~ year, threshold_plot_data, length)
names(n_sims_by_year)[2] <- "n_surviving_sims"
threshold_summary <- merge(threshold_summary, n_sims_by_year, by = "year")

library(ggplot2)
threshold_plot <- ggplot() +
  geom_line(
    data = threshold_plot_data,
    aes(x = year, y = mean_thresh, group = sim),
    colour = "grey75", linewidth = 0.3, alpha = 0.45
  ) +
  geom_line(
    data = threshold_summary,
    aes(x = year, y = mean_thresh),
    colour = "black", linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = 0.5, colour = "grey50", linewidth = 0.4) +
  annotate("text", x = years * 0.57, y = 0.58,
           label = "Lewis et al. default value", size = 3) +
  scale_x_continuous(limits = c(0, years), expand = expansion(mult = 0)) +
  labs(x = "Year", y = "Energy threshold at which migration occurs (thresh)") +
  theme_classic(base_size = 12)

print(threshold_plot)

end <- Sys.time()

print(end - start)

plan(sequential)


