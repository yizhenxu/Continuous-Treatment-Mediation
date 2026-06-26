############################################################
## Revised exact RKHS/CME conditional density estimator
## with internal 3-fold cross-validation for lambda.
##
## Drop-in replacement for fit_cme_conditional_density().
##
## Supports either:
##   fit_cme_conditional_density(df, A_col = "A", X_cols = ...)
## or:
##   fit_cme_conditional_density(df, y_col = "A", x_cols = ...)
##
## CV criterion:
##   mean held-out negative log conditional density
##   -mean(log f_hat(Y_i | X_i))
############################################################

.gaussian_kernel_matrix <- function(X1, X2 = NULL, sigma) {
  X1 <- as.matrix(X1)
  if (is.null(X2)) X2 <- X1
  X2 <- as.matrix(X2)

  if (ncol(X1) != ncol(X2)) {
    stop("X1 and X2 must have the same number of columns.")
  }
  if (!is.finite(sigma) || sigma <= 0) {
    stop("sigma must be positive and finite.")
  }

  sq1 <- rowSums(X1^2)
  sq2 <- rowSums(X2^2)

  D2 <- outer(sq1, sq2, "+") - 2 * tcrossprod(X1, X2)
  D2[D2 < 0] <- 0

  exp(-D2 / (2 * sigma^2))
}

.gaussian_density_kernel_A <- function(a_eval, A_train, h_A) {
  outer(as.numeric(a_eval), as.numeric(A_train), function(a, A) {
    dnorm(a, mean = A, sd = h_A)
  })
}

.cme_trapz_normalize <- function(Fhat, A_grid, eps = 1e-12) {
  if (length(A_grid) < 2) stop("A_grid must have length at least 2.")

  dx <- diff(A_grid)
  area <- rowSums(
    (Fhat[, -1, drop = FALSE] + Fhat[, -length(A_grid), drop = FALSE]) *
      matrix(dx, nrow = nrow(Fhat), ncol = length(dx), byrow = TRUE) / 2
  )

  Fhat / pmax(area, eps)
}

fit_cme_conditional_density <- function(df,
                                        A_col = NULL,
                                        X_cols = NULL,
                                        y_col = NULL,
                                        x_cols = NULL,
                                        sigma_X = NULL,
                                        h_A = NULL,
                                        lambda = NULL,
                                        lambda_grid = c(1e-1, 3e-2, 1e-2,
                                                        3e-3, 1e-3, 3e-4,
                                                        1e-4),
                                        cv_lambda = TRUE,
                                        cv_folds = 3,
                                        cv_seed = 1,
                                        standardize_X = TRUE,
                                        jitter = 1e-8,
                                        eps = 1e-12,
                                        verbose = TRUE) {
  ##########################################################
  ## Column-name compatibility
  ##########################################################

  if (is.null(y_col)) y_col <- A_col
  if (is.null(x_cols)) x_cols <- X_cols

  if (is.null(y_col)) stop("Specify y_col or A_col.")
  if (is.null(x_cols)) stop("Specify x_cols or X_cols.")

  if (!y_col %in% names(df)) stop("y_col/A_col not found in df.")
  if (!all(x_cols %in% names(df))) stop("Some x_cols/X_cols are not found in df.")

  Y_raw <- as.numeric(df[[y_col]])
  X_raw <- as.matrix(df[, x_cols, drop = FALSE])

  keep <- stats::complete.cases(Y_raw, X_raw)
  if (!all(keep)) {
    if (verbose) message("Dropping ", sum(!keep), " rows with missing outcome or conditioning variables.")
    Y_raw <- Y_raw[keep]
    X_raw <- X_raw[keep, , drop = FALSE]
  }

  n <- nrow(X_raw)
  p <- ncol(X_raw)

  if (n < 5) stop("Need at least 5 complete observations.")
  if (cv_folds < 2) stop("cv_folds must be at least 2.")
  if (cv_folds > n) stop("cv_folds cannot exceed n.")
  if (any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    stop("lambda_grid must contain positive finite values.")
  }

  ##########################################################
  ## Standardize X once on full training data
  ## For causal cross-fitting, this whole function is called inside
  ## each nuisance-training fold, so this does not leak information
  ## from the inference fold.
  ##########################################################

  if (standardize_X) {
    X_center <- colMeans(X_raw)
    X_scale <- apply(X_raw, 2, stats::sd)
    X_scale[!is.finite(X_scale) | X_scale == 0] <- 1
    Xs <- scale(X_raw, center = X_center, scale = X_scale)
  } else {
    X_center <- rep(0, p)
    X_scale <- rep(1, p)
    Xs <- X_raw
  }

  ##########################################################
  ## Bandwidths
  ##########################################################

  if (is.null(sigma_X)) {
    id <- sample(seq_len(n), min(n, 1000))
    Dmat <- as.matrix(stats::dist(Xs[id, , drop = FALSE]))
    sigma_X <- stats::median(Dmat[Dmat > 0], na.rm = TRUE)
    if (!is.finite(sigma_X) || sigma_X <= 0) sigma_X <- 1
  }

  if (is.null(h_A)) {
    h_A <- 1.06 * stats::sd(Y_raw) * n^(-1 / 5)
    if (!is.finite(h_A) || h_A <= 0) h_A <- 1
  }

  ##########################################################
  ## Internal low-level fitter for a fixed lambda
  ##########################################################

  fit_fixed_lambda <- function(Y_train, X_train_s, lambda_value) {
    n_train <- nrow(X_train_s)
    K_X <- .gaussian_kernel_matrix(X_train_s, sigma = sigma_X)
    Reg <- K_X + (n_train * lambda_value + jitter) * diag(n_train)

    chol_Reg <- tryCatch(
      chol(Reg),
      error = function(e) {
        chol(Reg + 10 * jitter * diag(n_train))
      }
    )

    solve_alpha <- function(k_x_mat) {
      backsolve(chol_Reg, forwardsolve(t(chol_Reg), k_x_mat))
    }

    list(
      Y_train = Y_train,
      X_train_s = X_train_s,
      K_X = K_X,
      chol_Reg = chol_Reg,
      solve_alpha = solve_alpha,
      lambda = lambda_value
    )
  }

  predict_fixed_at_y <- function(fixed_fit,
                                 X_new_s,
                                 y_eval,
                                 A_grid = NULL,
                                 n_grid = 300,
                                 grid_probs = c(0.001, 0.999),
                                 pred_chunk_size = 1000,
                                 normalize = TRUE,
                                 nonnegative = TRUE) {
    Y_train <- fixed_fit$Y_train
    X_train_s <- fixed_fit$X_train_s

    if (is.null(A_grid)) {
      q <- stats::quantile(Y_train, probs = grid_probs, na.rm = TRUE, names = FALSE)
      margin <- 3 * h_A
      A_grid <- seq(q[1] - margin, q[2] + margin, length.out = n_grid)
    } else {
      A_grid <- as.numeric(A_grid)
    }

    if (length(A_grid) < 2) stop("A_grid must have length at least 2.")

    Psi_A_grid <- .gaussian_density_kernel_A(A_grid, Y_train, h_A)

    n_new <- nrow(X_new_s)
    Fhat <- matrix(NA_real_, nrow = n_new, ncol = length(A_grid))

    starts <- seq(1, n_new, by = pred_chunk_size)

    for (s in starts) {
      e <- min(s + pred_chunk_size - 1, n_new)
      idx <- s:e

      K_train_new <- .gaussian_kernel_matrix(
        X_train_s,
        X_new_s[idx, , drop = FALSE],
        sigma = sigma_X
      )

      Alpha <- fixed_fit$solve_alpha(K_train_new)
      Fhat[idx, ] <- t(Psi_A_grid %*% Alpha)
    }

    if (nonnegative) Fhat <- pmax(Fhat, eps)
    if (normalize) Fhat <- .cme_trapz_normalize(Fhat, A_grid, eps = eps)

    out <- numeric(n_new)
    for (i in seq_len(n_new)) {
      val <- stats::approx(
        x = A_grid,
        y = Fhat[i, ],
        xout = y_eval[i],
        rule = 1
      )$y
      out[i] <- ifelse(is.na(val), eps, val)
    }

    pmax(out, eps)
  }

  ##########################################################
  ## 3-fold CV for lambda
  ##########################################################

  cv_table <- NULL

  if (isTRUE(cv_lambda)) {
    if (is.null(lambda)) {
      lambda_candidates <- lambda_grid
    } else {
      ## Include user-supplied lambda in the candidate set.
      lambda_candidates <- sort(unique(c(lambda, lambda_grid)), decreasing = TRUE)
    }

    set.seed(cv_seed)
    folds <- sample(rep(seq_len(cv_folds), length.out = n))

    cv_loss <- rep(NA_real_, length(lambda_candidates))

    if (verbose) {
      message("Selecting CME lambda by ", cv_folds, "-fold CV")
      message("  candidates: ", paste(signif(lambda_candidates, 4), collapse = ", "))
    }

    for (ell in seq_along(lambda_candidates)) {
      lam <- lambda_candidates[ell]
      fold_loss <- rep(NA_real_, cv_folds)

      for (k in seq_len(cv_folds)) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)

        fixed_k <- fit_fixed_lambda(
          Y_train = Y_raw[train_idx],
          X_train_s = Xs[train_idx, , drop = FALSE],
          lambda_value = lam
        )

        f_test <- predict_fixed_at_y(
          fixed_fit = fixed_k,
          X_new_s = Xs[test_idx, , drop = FALSE],
          y_eval = Y_raw[test_idx],
          pred_chunk_size = 1000,
          normalize = TRUE,
          nonnegative = TRUE
        )

        fold_loss[k] <- mean(-log(pmax(f_test, eps)), na.rm = TRUE)
      }

      cv_loss[ell] <- mean(fold_loss, na.rm = TRUE)

      if (verbose) {
        message("  lambda = ", signif(lam, 4),
                ", CV neg-log-density = ", signif(cv_loss[ell], 5))
      }
    }

    best_id <- which.min(cv_loss)
    lambda_selected <- lambda_candidates[best_id]

    cv_table <- data.frame(
      lambda = lambda_candidates,
      cv_neg_log_density = cv_loss
    )
  } else {
    if (is.null(lambda)) lambda <- 1e-3
    lambda_selected <- lambda
  }

  if (verbose) {
    message("Fitting exact CME conditional density estimator")
    message("  n = ", n, ", p = ", p)
    message("  sigma_X = ", signif(sigma_X, 4))
    message("  h_A     = ", signif(h_A, 4))
    message("  lambda  = ", signif(lambda_selected, 4))
  }

  ##########################################################
  ## Final fit on all data using selected lambda
  ##########################################################

  final_fit <- fit_fixed_lambda(
    Y_train = Y_raw,
    X_train_s = Xs,
    lambda_value = lambda_selected
  )

  predict_density_grid <- function(new_df,
                                   A_grid = NULL,
                                   n_grid = 200,
                                   grid_probs = c(0.005, 0.995),
                                   pred_chunk_size = 1000,
                                   normalize = TRUE,
                                   nonnegative = TRUE) {
    X_new <- as.matrix(new_df[, x_cols, drop = FALSE])
    X_new_s <- scale(X_new, center = X_center, scale = X_scale)

    n_new <- nrow(X_new_s)

    if (is.null(A_grid)) {
      q <- stats::quantile(Y_raw, probs = grid_probs, na.rm = TRUE, names = FALSE)
      margin <- 3 * h_A
      A_grid <- seq(q[1] - margin, q[2] + margin, length.out = n_grid)
    } else {
      A_grid <- as.numeric(A_grid)
    }

    if (length(A_grid) < 2) stop("A_grid must have length at least 2.")

    Psi_A_grid <- .gaussian_density_kernel_A(A_grid, Y_raw, h_A)

    Fhat <- matrix(NA_real_, nrow = n_new, ncol = length(A_grid))
    starts <- seq(1, n_new, by = pred_chunk_size)

    for (s in starts) {
      e <- min(s + pred_chunk_size - 1, n_new)
      idx <- s:e

      K_train_new <- .gaussian_kernel_matrix(
        Xs,
        X_new_s[idx, , drop = FALSE],
        sigma = sigma_X
      )

      Alpha <- final_fit$solve_alpha(K_train_new)
      Fhat[idx, ] <- t(Psi_A_grid %*% Alpha)
    }

    if (nonnegative) Fhat <- pmax(Fhat, eps)
    if (normalize) Fhat <- .cme_trapz_normalize(Fhat, A_grid, eps = eps)

    colnames(Fhat) <- paste0("Y=", signif(A_grid, 5))

    list(
      A_grid = A_grid,
      y_grid = A_grid,
      density = Fhat
    )
  }

  predict_at_A <- function(new_df,
                           A_eval = NULL,
                           y_eval = NULL,
                           A_grid = NULL,
                           n_grid = 300,
                           grid_probs = c(0.001, 0.999),
                           pred_chunk_size = 1000,
                           normalize = TRUE,
                           nonnegative = TRUE,
                           eps_tail = eps) {
    if (is.null(y_eval)) y_eval <- A_eval

    if (is.null(y_eval)) {
      if (!y_col %in% names(new_df)) {
        stop("A_eval/y_eval is NULL and y_col/A_col is not found in new_df.")
      }
      y_eval <- as.numeric(new_df[[y_col]])
    } else {
      y_eval <- as.numeric(y_eval)
    }

    if (length(y_eval) != nrow(new_df)) {
      stop("length(A_eval/y_eval) must equal nrow(new_df).")
    }

    grid_fit <- predict_density_grid(
      new_df = new_df,
      A_grid = A_grid,
      n_grid = n_grid,
      grid_probs = grid_probs,
      pred_chunk_size = pred_chunk_size,
      normalize = normalize,
      nonnegative = nonnegative
    )

    A_grid <- grid_fit$A_grid
    Fhat <- grid_fit$density

    out <- numeric(nrow(new_df))
    for (i in seq_len(nrow(new_df))) {
      val <- stats::approx(
        x = A_grid,
        y = Fhat[i, ],
        xout = y_eval[i],
        rule = 1
      )$y
      out[i] <- ifelse(is.na(val), eps_tail, val)
    }

    pmax(out, eps_tail)
  }

  ## Alias with clearer generic name.
  predict_density_at <- function(new_df,
                                 y_eval,
                                 A_grid = NULL,
                                 n_grid = 300,
                                 grid_probs = c(0.001, 0.999),
                                 pred_chunk_size = 1000,
                                 normalize = TRUE,
                                 nonnegative = TRUE,
                                 eps_tail = eps) {
    predict_at_A(
      new_df = new_df,
      y_eval = y_eval,
      A_grid = A_grid,
      n_grid = n_grid,
      grid_probs = grid_probs,
      pred_chunk_size = pred_chunk_size,
      normalize = normalize,
      nonnegative = nonnegative,
      eps_tail = eps_tail
    )
  }

  list(
    A_train = Y_raw,
    y_train = Y_raw,
    X_train_scaled = Xs,
    X_center = X_center,
    X_scale = X_scale,
    X_cols = x_cols,
    x_cols = x_cols,
    A_col = y_col,
    y_col = y_col,
    sigma_X = sigma_X,
    h_A = h_A,
    lambda = lambda_selected,
    lambda_selected = lambda_selected,
    lambda_grid = lambda_grid,
    cv_lambda = cv_lambda,
    cv_folds = cv_folds,
    cv_table = cv_table,
    eps = eps,
    K_X = final_fit$K_X,
    chol_Reg = final_fit$chol_Reg,
    solve_alpha = final_fit$solve_alpha,
    predict_density_grid = predict_density_grid,
    predict_at_A = predict_at_A,
    predict_density_at = predict_density_at
  )
}
