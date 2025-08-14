#packages
library(coda)
library(bdskytools)
library(beastio)
library(lubridate)

#paths

#Replication log-discarded
#Log_file <- "C:/Users/filip/Desktop/MCMC replication run 30 mil/Beauti_29_seq 30 mil MCMC.log"

#Second replication-discarded 
#logfile2 <- "C:/Users/filip/Desktop/Second replication/Beauti_29_seq RE 4.log"

#Third replication-used in the analysis
logfile4 <- "C:/Users/filip/Desktop/MCMC 108 REPLICATION/Beauti_29_seq 108.log"

#My analysis-discarded
#logfile <- "C:/Users/filip/Desktop/MCMC original analysis/Beauti_29_seq ORIGINAL RUN.log"

#My analysis 2- file used in the analysis

logfile3 <- "C:/Users/filip/Desktop/2 RE, 1 SAMPL/Beauti_29_seq 2 RE.log"

#read log with 10% burnin - original-discarded
#bdsky_trace <- beastio::readLog(Log_file, burnin = 0.10)
#all_names    <- colnames(bdsky_trace)

#  read log with 10% burnin - second replication-discarded
#bdsky_trace <- beastio::readLog(logfile2, burnin = 0.10)
#all_names    <- colnames(bdsky_trace)

# read log with 10% burnin - third replication-used in the analysis
bdsky_trace <- beastio::readLog(logfile4, burnin = 0.10)
all_names    <- colnames(bdsky_trace)

# - read log with 10% burnin - my analysis-discarded
#bdsky_trace <- beastio::readLog(logfile, burnin = 0.10)
#all_names    <- colnames(bdsky_trace)

#  read log with 10% burnin - my analysis2-used in the analysis
bdsky_trace <- beastio::readLog(logfile3, burnin = 0.10)
all_names    <- colnames(bdsky_trace)

# --- helper: numeric-sort columns that end .1, .2, _1, _2, etc. ---
order_by_interval_suffix <- function(nms){
  idx <- suppressWarnings(as.integer(sub(".*(?:\\.|_|\\[)(\\d+)\\]?\\s*$", "\\1", nms)))
  ord <- order(ifelse(is.na(idx), Inf, idx), nms)
  nms[ord]
}

# --- STEP 1: grab Re interval samples (iterations × nIntervals) ---
re_cols <- grep("reproductiveNumber", all_names, value = TRUE, ignore.case = TRUE)
if (!length(re_cols)) re_cols <- grep("Re(\\b|_)", all_names, value = TRUE, ignore.case = TRUE)
stopifnot(length(re_cols) > 0)

re_cols <- order_by_interval_suffix(re_cols)
Re_mat  <- as.matrix(bdsky_trace[, re_cols, drop = FALSE])  # iterations × nIntervals

# --- STEP 2: choose a regular time grid (years before present) ---
# Find a tmrca/origin parameter to set the *oldest* time (in years ago)
tmrca_cols <- grep("tmrca|treeheight|origin", all_names, value = TRUE, ignore.case = TRUE)
stopifnot(length(tmrca_cols) > 0)
tmrca_median <- median(as.numeric(bdsky_trace[, tmrca_cols[1]]))  # years before present (BEAST convention)

gridsize <- max(200, ncol(Re_mat) * 5)  # denser than intervals but not crazy
timegrid <- seq(from = tmrca_median, to = 0, length.out = gridsize)  # years before present


# Make sure tmrca_median is a single number
stopifnot(length(tmrca_median) == 1, !is.na(tmrca_median))

# 1) Ensure vectors go from oldest to present (in years before present)
# Example: tmrca_median = 6.5 → present = 0
interval_edges <- seq(from = 0, to = tmrca_median, length.out = ncol(Re_mat) + 1)
timegrid_sorted <- seq(from = 0, to = tmrca_median, length.out = gridsize)

# 2) Map each gridpoint to an interval
idx <- findInterval(timegrid_sorted, interval_edges,
                    rightmost.closed = TRUE, all.inside = TRUE)

# 3) Fill the grid with the Re values from the matching intervals
Re_grid <- Re_mat[, idx, drop = FALSE]   # iterations × gridsize



# --- STEP 4: HPD + median at every gridpoint -> shape 3 × gridsize (lower, median, upper) ---
hpd_3xN <- function(M){             # M = iterations × points
  stopifnot(is.matrix(M), nrow(M) > 0, ncol(M) > 0)
  out <- t(vapply(seq_len(ncol(M)), function(j){
    v <- M[, j]
    h <- HPDinterval(mcmc(v), prob = 0.95)
    c(lower = h[1], median = median(v), upper = h[2])
  }, numeric(3)))
  t(out)                            # 3 × N
}
Re_grid_hpd <- hpd_3xN(Re_grid)

# --- STEP 5: plot “years before present” (tutorial-style smooth skyline) ---
bdskytools::plotSkyline(timegrid, Re_grid_hpd, type = "step",
                        xlab = "Years before most recent sample",
                        ylab = "Effective reproductive number (Re)")

# --- STEP 6: set latest sample date and plot on calendar time ---

# 1) Set your most recent sampling time (decimal year):
most_recent_decimal <- 2013.379

# 2) Convert "years before present" -> calendar decimal years -> Date
decimal_to_date <- function(x){
  y  <- floor(x)
  d0 <- lubridate::ymd(sprintf("%d-01-01", y))
  d1 <- lubridate::ymd(sprintf("%d-01-01", y + 1))
  d0 + (x - y) * as.numeric(d1 - d0)
}

# NOTE: use the same vector you used to build Re_grid_hpd.
# If you followed the manual-grid fix, that's `timegrid_sorted` (0 -> tmrca).
calendar_decimal <- most_recent_decimal - timegrid_sorted
calendar_dates   <- decimal_to_date(calendar_decimal)

# 3) Make x increase left→right (older → newer)
ord <- order(calendar_dates)
cal_x <- calendar_dates[ord]
Re_med <- Re_grid_hpd[2, ord]
Re_lo  <- Re_grid_hpd[1, ord]
Re_hi  <- Re_grid_hpd[3, ord]

# 4) Plot like Zinsstag et al.: median solid, HPD dashed, threshold at Re=1

plot(cal_x, Re_med, type = "l", lwd = 2,
     xlab = "Calendar year",
     ylab = "Effective reproductive number (Re)",
     main = "Smoothed Re skyline (HPD on regular time grid)")
lines(cal_x, Re_lo, lty = 2)
lines(cal_x, Re_hi, lty = 2)
abline(h = 1, lty = 3)

#grep("reproductiveNumber|\\bRe\\b", nms, value = TRUE, ignore.case = TRUE)
#dim(Re_mat)



