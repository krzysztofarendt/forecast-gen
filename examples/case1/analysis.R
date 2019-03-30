# Set working directory -------------------------------------------------------
setwd("D:/code/forecast-gen/examples/case1")
rm(list=ls())

# Load libraries --------------------------------------------------------------
library(tidyverse)
library(lubridate)

# Read cases ------------------------------------------------------------------
wd <- file.path("results", "mpc", "r1c1_dymola_1e-11", "h6")  # Horizon might change

re0  <- file.path(wd, "re0")
re10 <- file.path(wd, "re10")
re20 <- file.path(wd, "re20")
re30 <- file.path(wd, "re30")

case_dirs <- list(re0, re10, re20, re30)
names(case_dirs) <- c("RE0", "RE10", "RE20", "RE30")

inp  <- list()
xctr <- list()
xemu <- list()
yemu <- list()
u    <- list()

for (n in names(case_dirs)) {

    case_path <- case_dirs[[n]]
    run_dirs  <- list.dirs(path = case_path, recursive = FALSE)
    
    inp[[n]]  <- list()
    xctr[[n]] <- list()
    xemu[[n]] <- list()
    yemu[[n]] <- list()
    u[[n]]    <- list()
    
    for (i in 1:length(run_dirs)) {
        inp[[n]][[i]]  <- read_csv(file.path(case_path, i, "inp_ctr.csv"))
        xctr[[n]][[i]] <- read_csv(file.path(case_path, i, "xctr.csv"))
        xemu[[n]][[i]] <- read_csv(file.path(case_path, i, "xemu.csv"))
        yemu[[n]][[i]] <- read_csv(file.path(case_path, i, "yemu.csv"))
        u[[n]][[i]]    <- read_csv(file.path(case_path, i, "u.csv"))
    }
}

# Read constraints
constr <- read_csv(file.path(case_dirs[[1]], 1, "constr.csv"))
constr$time <- constr$time / 3600
constr$Tmin <- constr$Tmin - 273.15
constr$Tmax <- constr$Tmax - 273.15

# Get variable (mean and sd) --------------------------------------------------
get_var <- function(v, l, offset = 0) {
    
    # Find RE cases
    re_cases <- names(l)
    
    # Find number of runs
    n_runs <- numeric(length(re_cases))
    names(n_runs) <- re_cases
    
    i <- 1
    for (case in re_cases) {
        n_runs[i] <- length(l[[case]])
        i <- i + 1
    }
    
    # Put all data into data.frame
    df <- tibble(time = l[[1]][[1]]$time / 3600)
    
    for (case in re_cases) {
        for (n in 1:n_runs[[case]]) {
            values <- l[[case]][[n]][[v]] + offset
            df[[paste(v, case, n, sep = "_")]] <- values
        }
    }
    
    # Collapse to means and sd
    m <- tibble(time = df$time)
    s <- tibble(time = df$time)
    
    for (case in re_cases) {
        cols <- paste(v, case, 1:n_runs[[case]], sep = "_")    
        m[[case]] <- rowMeans(df[cols])
        s[[case]] <- apply(df[cols], 1, sd)
    }
    
    # Return means and sd
    res <- list(m, s)
    names(res) <- c("mean", "sd")
    
    return(res)
}

# Get variable (all runs) -----------------------------------------------------
get_var_all <- function(v, l, offset = 0) {
    
    # Find RE cases
    re_cases <- names(l)
    
    # Find number of runs
    n_runs <- numeric(length(re_cases))
    names(n_runs) <- re_cases
    
    i <- 1
    for (case in re_cases) {
        n_runs[i] <- length(l[[case]])
        i <- i + 1
    }
    
    # Put all data into data.frame
    df <- tibble(time = l[[1]][[1]]$time / 3600)
    
    for (case in re_cases) {
        for (n in 1:n_runs[[case]]) {
            values <- l[[case]][[n]][[v]] + offset
            df[[paste(v, case, n, sep = "_")]] <- values
        }
    }
    
    return(df)
}


# Plotting function -----------------------------------------------------------
plot_time_vs_var <- function(v, l, offset = 0) {
    data <- get_var(v, l, offset)    
    
    cols <- c("0%" = 1, "10%" = 2, "20%" = 3, "30%" = 4)
    sd_cols <- c("10% +/- sd" = 2, "20% +/- sd" = 3, "30% +/- sd" = 4)
    
    p <- ggplot() +
        geom_line(data$mean, mapping = aes(time, RE0, color = "0%"), size = 2) +
        geom_line(data$mean, mapping = aes(time, RE10, color = "10%"), size = 1) +
        geom_ribbon(data$sd, mapping = aes(x = time, ymin = data$mean$RE10 - RE10,
                                           ymax = data$mean$RE10 + RE10,
                                           fill = names(sd_cols)[1]),
                    alpha = 0.2) +
        geom_line(data$mean, mapping = aes(time, RE20, color = "20%"), size = 1) +
        geom_ribbon(data$sd, mapping = aes(x = time, ymin = data$mean$RE20 - RE20,
                                           ymax = data$mean$RE20 + RE20,
                                           fill = names(sd_cols)[2]),
                    alpha = 0.2) +
        geom_line(data$mean, mapping = aes(time, RE30, color = "30%"), size = 1) +
        geom_ribbon(data$sd, mapping = aes(x = time, ymin = data$mean$RE30 - RE30,
                                           ymax = data$mean$RE30 + RE30,
                                           fill = names(sd_cols)[3]),
                    alpha = 0.2) +
        scale_color_manual(name = "Mean", values = cols) +
        scale_fill_manual(name = "Range", values = sd_cols) +
        guides(color = guide_legend(order = 1),
               fill = guide_legend(order = 2)) +
        ylab(v) +
        theme_bw() +
        scale_x_continuous(name = "time [h]", breaks = c(0, 24, 48, 72))
    
    return(p)
}

# Tout ------------------------------------------------------------------------
plot_time_vs_var("Tout", inp) +
    labs(title = "Oudoor temperature",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/outdoor_T.png")

# solrad ----------------------------------------------------------------------
plot_time_vs_var("solrad", inp) + coord_cartesian(ylim = c(45, 1000)) +
    labs(title = "Solar radiation",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/solar_rad.png")

# occ -------------------------------------------------------------------------
plot_time_vs_var("occ", inp) + coord_cartesian(ylim = c(2, 60)) +
    labs(title = "Occupancy schedules",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/occupancy.png")

# Tin -------------------------------------------------------------------------
p <- plot_time_vs_var("cair.T", xemu, -273.15)
p + ylab('Tin') +
    geom_line(constr, mapping = aes(time, Tmin), col = 1, linetype = "dashed") +
    geom_line(constr, mapping = aes(time, Tmax), col = 1, linetype = "dashed") +
    labs(title = "Indoor temperature vs. constraints",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/indoor_T.png")

# u ---------------------------------------------------------------------------
u_lim <- tibble(time = constr$time)
u_lim$u_max <- rep(100, times = nrow(u_lim))
u_lim$u_min <- rep(-100, times = nrow(u_lim))
    
p <- plot_time_vs_var("vpos", u)
p + ylab('u [%]') + 
    geom_line(u_lim, mapping = aes(time, u_min), col = 1, linetype = "dashed") +
    geom_line(u_lim, mapping = aes(time, u_max), col = 1, linetype = "dashed") +
    labs(title = "Heating/cooling signal",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/u_signal.png")

# Calculate total energy consumption Q [kWh] ----------------------------------
n_runs = length(u[[length(u)]])  # Last one, because RE0 doesn't need multiple runs

Q <- tibble(RE0  = numeric(n_runs),
            RE10 = numeric(n_runs),
            RE20 = numeric(n_runs),
            RE30 = numeric(n_runs))

get_Q <- function(x) {
    # x: data frame (columns: time, vpos)
    # return: scalar (Qtot)
    
    Q <- sum(abs(x$vpos)) / 1000
    
    return(Q)
}

for (i in 1:length(Q)) {
    Q[[i]] <- sapply(u[[i]], get_Q)
}

df <- Q %>% gather("RE", "Q") %>%
    mutate(RE = as.factor(paste0(substr(RE, 3, 5), "%")))

ggplot(df, aes(x=RE, y=Q)) + geom_boxplot(fill = 1:4, alpha = 0.3) +
    xlab("Relative Error [%]") +
    ylab("Q [kWh]") + theme_bw() +
    labs(title = "Energy demand Q vs. forecast relative error",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)")

ggsave("figs/energy_demand.png")

# Calculate total discomfort DC [Kh] ------------------------------------------
Tin <- get_var_all("cair.T", xemu, -273.15)

names(Tin) <- gsub("cair.T_", "", names(Tin))

# Temperature violation: lower constraint
Tvlo <- constr$Tmin - select(Tin, -time)
Tvlo[Tvlo < 0] <- 0
Tvlo <- colSums(Tvlo)
Tvlo

# Temperature violation: upper constraint
Tvhi <- select(Tin, -time) - constr$Tmax
Tvhi[Tvhi < 0] <- 0
Tvhi <- colSums(Tvhi)
Tvhi

# Total violation: lower + upper
Tv <- Tvlo + Tvhi

# Prepare for bar plot
Tv <- as.data.frame(Tv)
re_run <- strsplit(row.names(Tv), "_")
re <- sapply(re_n, function(x) x[1])
run <- sapply(re_n, function(x) x[2])
Tv$RE <- re
Tv$run <- run
row.names(Tv) <- 1:nrow(Tv)

Tv <- Tv %>% 
    mutate(DC = Tv) %>% 
    select(-Tv) %>%
    mutate(RE = as.factor(paste0(substr(RE, 3, 5), "%")),
           run = as.factor(run))

# Plot
ggplot(Tv, aes(RE, DC)) + geom_boxplot(fill = 1:4, alpha = 0.3) +
    xlab("Relative Error [%]") +
    ylab("Discomfort D [Kh]") +
    labs(title = "Discomfort vs. forecast relative error",
         subtitle = "Based on 5 MPC runs per error level (random error in each run)") +
    theme_bw()

ggsave("figs/discomfort.png")

# Actual Tout profiles (all runs) ---------------------------------------------
Tout <- tibble(
    time = inp$RE0[[1]]$time / 3600,
    ideal = inp$RE0[[1]]$Tout,
    run1 = inp$RE30[[1]]$Tout,
    run2 = inp$RE30[[2]]$Tout,
    run3 = inp$RE30[[3]]$Tout,
    run4 = inp$RE30[[4]]$Tout,
    run5 = inp$RE30[[5]]$Tout
)

Tout <- Tout %>% gather("Case", "Tout", -time)

cols <- c("ideal" = 'black', "run1" = 'gray', "run2" = 'gray',
          "run3" = 'gray', "run4" = 'gray', "run5" = 'gray')

ggplot(Tout, aes(time, Tout, color = Case)) +
    geom_line() +
    scale_color_manual(values = cols) +
    ggtitle("Relative Error: 30%") +
    theme_bw() +
    ylim(-13, 21) +
    scale_x_continuous(name = "time [h]", breaks = c(0, 24, 48, 72))

ggsave("figs/Tout_RE30.png")



