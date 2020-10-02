args = commandArgs(trailingOnly = T)
sample_dir = args[1]
tune_sample_1 = as.numeric(args[2])
tune_sample_2 = as.numeric(args[3])
n.chr = args[4]
resolution = as.numeric(args[5])

file_list = list.files(sample_dir)
sample_files = grep("normalized", file_list, value = TRUE)
anchor_files = grep("anchors", file_list, value = TRUE)

n.samples = length(sample_files)

# Sanity check

sample_end = unlist(gregexpr("normalized", sample_files))
anchor_end = unlist(gregexpr("anchors", anchor_files))


samples = anchors = rep(NA, n.samples)
for (i in 1:n.samples) {
	samples[i] = substr(sample_files[i], 1, sample_end[i] - 2)
	anchors[i] = substr(anchor_files[i], 1, anchor_end[i] - 2)
}

if (!all(sort(samples) == sort(anchors))) {
	print("Each sample must have corresponding anchor list!")
	quit("no")
}


library(plyr)
library(data.table)


n.bins = 1e6/resolution * 2


# Perform HPRep procedure

get_anchor_pair_counts <- function(anchor, pairs, resolution, n_bins) {
  pair_vec <- c(seq(anchor - n_bins/2 * resolution, anchor - resolution, resolution),
                seq(anchor + resolution, anchor + n_bins/2 * resolution, resolution))
  res_vec <- rep(0, n_bins)
  for (i in 1:n_bins) {
    count_location <- pairs[, 1] %in% pair_vec[i]
    if (sum(count_location == 1)) {
      res_vec[i] <- pairs[count_location, 2]
    }
  }
  return(res_vec)
}

generate_contact_matrix <- function(sample, anchors, resolution, n_bins) {
  n_anchors <- length(anchors)
  contact_mat <- matrix(0, nrow = n_anchors, ncol = n_bins)
  for (i in 1:n_anchors) {
    c1 <- sample[sample[, 1] %in% anchors[i], c(2, 5)]
    c2 <- sample[sample[, 2] %in% anchors[i], c(1, 5)]
    colnames(c1) <- colnames(c2) <- c("Bin", "Normalized_count")
    pairs = rbind(c1, c2)
    contact_mat[i, ] <- get_anchor_pair_counts(anchors[i], pairs, resolution, n_bins)
  }
  rownames(contact_mat) <- anchors
  return(contact_mat)
}

matrix_smooth = function(sample, h) {
  j_range <- dim(sample)[2]
  sample_smooth <- sample
  if (h > 0) {
    for (i in 1:j_range) {
      j_min = max(1, i-h)
      j_max = min(j_range, i+h)
      sample_smooth[, i] = apply(sample[ ,j_min:j_max], 1, mean)
    }
  }
  return(sample_smooth)
}

hicrep <- function(sample1, sample2, resolution) {
  n_bins <- 1e6/resolution * 2
  #anchors <- data.table::fread(anchor_file, data.table = F)
  #anchors <- sort(anchors[anchors$chr == paste0("chr", chr), 2])

  #cmat1 <- generate_contact_matrix(sample1, anchors, resolution, n_bins)
  #cmat2 <- generate_contact_matrix(sample2, anchors, resolution, n_bins)

  cmat_smooth1 <- sample1
  cmat_smooth2 <- sample2

  rho_k <- r2_k <- N_k <- rep(NA, n_bins/2)
  for (i in 1:(n_bins/2)) {
    x <- c(cmat_smooth1[, n_bins/2-i+1], cmat_smooth1[, n_bins/2+i])
    y <- c(cmat_smooth2[, n_bins/2-i+1], cmat_smooth2[, n_bins/2+i])
    x1 <- x[x != 0 | y != 0]
    y1 <- y[x != 0 | y != 0]
    n <- length(x1)
    if (n != 0) {
      rho_k[i] <- cor(x1, y1)
      r2_k[i] <- sqrt( (sum(x1^2)/n - mean(x1)^2) * (sum(y1^2)/n - mean(y1)^2) )
      N_k[i] <- n
    }
  }
  weights <- r2_k * N_k / sum(r2_k * N_k, na.rm = T)
  rho_s <- sum(weights * rho_k, na.rm = T)

  return(rho_s)
}


matrix_flatten <- function(sample, resolution, n_bins) {
  n_row <- dim(sample)[1]
  n_col <- dim(sample)[2]
  sample_flat <- matrix(0, nrow = n_row * n_col, ncol = 3)
  sample_flat[, 3] <- as.vector(sample)
  sample_flat[, 1] <- sample(as.integer(row.names(samples)), n_col)
  pass <- 1
  for (i in c(-(n_bins/2):-1, 1:(n_bins/2))) {
    sample_flat[((pass-1) * n_row + 1):(pass * n_row), 2] <- as.integer(row.names(sample)) + i * resolution
    pass <- pass + 1
  }
  bin1 <- pmin(sample_flat[,1], sample_flat[,2])
  bin2 <- pmax(sample_flat[,1], sample_flat[,2])
  sample_flat[,1] <- bin1
  sample_flat[,2] <- bin2
  colnames(sample_flat) = c("Bin1", "Bin2", "Normalized_contact")
  return(sample_flat)
}


hicrep_tune <- function(sample1, sample2, anchors, resolution, seed) {
  set.seed(seed)
  n_bins <- 1e6/resolution * 2
  #anchors <- data.table::fread(anchor_file, data.table = F)
  #anchors <- sort(anchors[anchors$chr == paste0("chr", chr), 2])

  cmat1 <- generate_contact_matrix(sample1, anchors, resolution, n_bins)
  cmat2 <- generate_contact_matrix(sample2, anchors, resolution, n_bins)

  n_row <- dim(cmat1)[1]
  cvec1 <- as.vector(cmat1)
  idx <- which(cvec1 != 0)
  cvec2 <- as.vector(cmat2)

  res = rep(NA, 21)
  for (h in 0:20) {
    res_temp = rep(NA, 10)
    for (j in 1:10) {
      idx_sample <- sample(idx, floor(0.25 * length(idx)), replace = F)
      cvec1_new = rep(0, length(cvec1))
      cvec1_new[idx_sample] = cvec1[idx_sample]
      cmat1_sampled <- matrix(cvec1_new, nrow = n_row)
      rm(cvec1_new)

      cvec2_new = rep(0, length(cvec2))
      cvec2_new[idx_sample] = cvec2[idx_sample]
      cmat2_sampled <- matrix(cvec2_new, nrow = n_row)
      rm(cvec2_new)

      cmat_smooth1 <- matrix_smooth(cmat1_sampled, h)
      cmat_smooth2 <- matrix_smooth(cmat2_sampled, h)

      rho_k <- r2_k <- N_k <- rep(NA, n_bins/2)
      for (i in 1:(n_bins/2)) {
        x <- c(cmat_smooth1[, n_bins/2-i+1], cmat_smooth1[, n_bins/2+i])
        y <- c(cmat_smooth2[, n_bins/2-i+1], cmat_smooth2[, n_bins/2+i])
        x1 <- x[x != 0 | y != 0]
        y1 <- y[x != 0 | y != 0]
        n <- length(x1)
        if (n != 0) {
          rho_k[i] <- cor(x1, y1)
          r2_k[i] <- sqrt( (sum(x1^2)/n - mean(x1)^2) * (sum(y1^2)/n - mean(y1)^2) )
          N_k[i] <- n
        }
      }
      weights <- r2_k * N_k / sum(r2_k * N_k, na.rm = T)
      rho_s <- sum(weights * rho_k, na.rm = T)

      res_temp[j] = rho_s
    }
    res[h+1] = mean(res_temp)
    print(paste0("Completed training for h = ", h))
    if (h > 0) {
      if (res[h+1] / res[h] < 1.01) {
        break
      }
    }
  }
  print(paste0("Optimized h is ", h-1))
  return(h-1)
}


#hicrep_tune(get(samples[tune_sample_1]), get(samples[tune_sample_2]), anchors, resolution, 1)


h = 6  # Yin's data and hichip
#h = 7  # mouse
runs = combn(n.samples, 2)
n.runs = dim(runs)[2]
res = matrix(NA, nrow = n.runs, ncol = 2 + as.numeric(n.chr))
#s1_name = s2_name = SCC = rep(NA, n.runs)


for (chr in 1:as.numeric(n.chr)) {

# Read in anchors
  anchors = data.frame()
  for (i in 1:n.samples) {
        a =  read.table(anchor_files[i], header = T, stringsAsFactors = F)
        anchors = rbind(anchors, a[a$chr == paste0("chr", chr), ])
  } 


  anchors = sort(unique(anchors)[, 2])

  for (i in 1:n.samples) {
        a = fread(sample_files[i], header = T, data.table = F)
        a = a[a$Chr == paste0("chr", chr), 2:6]
        assign(samples[i], a)
  }

  for (i in 1:n.runs) {
	res[i, 1] = samples[runs[1, i]]
        res[i, 2] = samples[runs[2, i]]
  }

  n.bins = 1e6/resolution * 2

  for (i in 1:n.samples) {
	assign(paste0("contact", i), generate_contact_matrix(get(samples[i]), anchors, resolution, n.bins))
	assign(paste0("smoothed", i), matrix_smooth(get(paste0("contact", i)), h))
	print(date())
  }

  for (i in 1:n.runs) {
	res[i, 2 + chr] = hicrep(get(paste0("smoothed", runs[1, i])), get(paste0("smoothed", runs[2, i])), resolution)
  }
  print(paste0("Completed chromosome ", chr, " ", date()))
}
write.table(res, file = paste0("SCC.txt"), row.names = F, col.names = F, sep = '\t', quote = F)




