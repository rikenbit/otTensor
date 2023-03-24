OTT <- function(X, Y, f, ps=NULL, qs=NULL,
    loss=.absolute_error, num.sample=1000,
    num.iter=200, epsilon=1e-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkOTT(X, Y, f, ps, qs, loss, num.sample, num.iter, epsilon, verbose)
    ######################################
    # Initialization
    ######################################
    int <- .initOTT(X, Y, f, ps, qs)
    ps <- int$ps
    qs <- int$qs
    Ts <- int$Ts
    A <- int$A
    D <- int$D
    Is <- int$Is
    Ks <- int$Ks
    ######################################
    # Iteration
    ######################################
    iter <- 1
    for(s in 0:(num.iter-1)){
        for(a in seq_len(A)){
            # The size of nabla_T_a_epsilon is the same as T_a.
            nabla_T_a_epsilon <- matrix(0, length(ps[[a]]), length(qs[[a]]))
            # d_prime_list is a list of d' (dimension) that satisfies f(d') == a
            d_prime_list <- which(f == a)
            for(d_prime in d_prime_list){
                # probability_list is a list of scalars that take values between 0 and 1.
                # This is the probability of each sample.
                probability_list <- numeric(num.sample)
                # A list of matrices that correspond to the left side of L.
                c_d_prime_list <- vector("list", num.sample)
                for (i_sample in seq_len(num.sample)) {
                    # Sample i1, i2, ..., iD and
                    # k1, k2, ..., kD and
                    # id and kd are used in the right side.
                    # the rest are used in the left side.
                    # i_list[d] means id.
                    # k_list[d] means kd.
                    # The size of i_list is D
                    # The size of k_list is D
                    i_list <- integer(D)
                    k_list <- integer(D)
                    for (d in 1:D) {
                        i_list[d] <- sample(1:Is[d], 1)
                        k_list[d] <- sample(1:Ks[d], 1)
                    }
                    # Right hand side of the equation.
                    # products is a scalar.
                    products <- 1.0
                    for (d in 1:D) {
                        if (d == d_prime) next
                        i_d <- i_list[d]
                        k_d <- k_list[d]
                        products <- products * Ts[[f[d]]][i_d, k_d]
                    }
                    rhs <- products
                    # Compute the left side of the equation.
                    Xs <- .extract_vector_from_array(X, i_list, d_prime)
                    Ys <- .extract_vector_from_array(Y, k_list, d_prime)
                    L <- outer(Xs, Ys, loss)
                    # Store the probability and the corresponding matrix.
                    probability_list[i_sample] <- rhs
                    c_d_prime_list[[i_sample]] <- L
                }
                # normalized_probability_list is a list of scalars that take values between 0 and 1.
                normalized_probability_list <- probability_list / sum(probability_list)
                # Generate E(C_d_prime) by summing matrices with normalized probabilities.
                hat_expectation_c_d_prime <- NULL
                for (i in 1:num.sample) {
                    if (is.null(hat_expectation_c_d_prime)) {
                        hat_expectation_c_d_prime <- normalized_probability_list[i] * c_d_prime_list[[i]]
                    } else {
                        hat_expectation_c_d_prime <- hat_expectation_c_d_prime + normalized_probability_list[i] * c_d_prime_list[[i]]
                    }

                }
                nabla_T_a_epsilon <- nabla_T_a_epsilon + hat_expectation_c_d_prime
            }
            # Sinkhorn algorithm
            Ts[[a]] <- .compute_transport_plan(p = rowSums(Ts[[a]]), q = colSums(Ts[[a]]), M = nabla_T_a_epsilon)
        }
        iter <- iter + 1
    }
    # Output
    list(Ts=Ts)
}

# absolute error
.absolute_error <- function(x, y){
    abs(x - y)
}

.checkOTT <- function(X, Y, f, ps, qs, loss, num.sample,
    num.iter, epsilon, verbose){
    stopifnot(is(X)[1] == "Tensor")
    stopifnot(is(Y)[1] == "Tensor")
    stopifnot(length(dim(X)) == length(dim(Y)))
    stopifnot(is.vector(f))
    stopifnot(length(unique(f)) <= length(dim(X)))
    stopifnot(length(unique(f)) <= length(dim(Y)))
    if(!is.null(ps)){
        stopifnot(is.list(ps))
        stopifnot(length(ps) == length(unique(f)))
    }
    if(!is.null(qs)){
        stopifnot(is.list(qs))
        stopifnot(length(qs) == length(unique(f)))
    }
    stopifnot(is(loss)[1] == "function")
    stopifnot(is.numeric(num.sample))
    stopifnot(num.sample >= 1)
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter >= 1)
    stopifnot(is.numeric(epsilon))
    stopifnot(epsilon >= 0)
    stopifnot(is.logical(verbose))
}

.initOTT <- function(X, Y, f, ps, qs){
    D <- length(dim(X))
    A <- length(unique(f))
    Is <- dim(X)
    Ks <- dim(Y)
    if (is.null(ps)){
        ps <- list()
        for (a in seq_len(A)) {
            ds <- which(f == a)
            d <- ds[1]
            length_of_p_a <- dim(X)[d]
            ps[[a]] <- rep(1, length_of_p_a)
            ps[[a]] <- ps[[a]] / sum(ps[[a]])
        }
    }
    if (is.null(qs)){
        qs <- list()
        for (a in 1:A) {
            ds <- which(f == a)
            d <- ds[1]
            length_of_q_a <- dim(Y)[d]
            qs[[a]] <- rep(1, length_of_q_a)
            qs[[a]] <- qs[[a]] / sum(qs[[a]])
        }
    }
    stopifnot(length(ps) == A)
    stopifnot(length(qs) == A)
    Ts <- list()
    for (a in 1:A) {
        Ts[[a]] <- outer(ps[[a]], qs[[a]])
    }
    list(ps=ps, qs=qs, Ts=Ts, A=A, D=D, Is=Is, Ks=Ks)
}

.extract_vector_from_array <- function(target_array, fixed_indices, variable_dimension){
    indices <- as.list(fixed_indices)
    indices[[variable_dimension]] <- 1:dim(target_array@data)[variable_dimension]
    indices_as_matrix <- do.call(cbind, indices)
    as.vector(target_array@data[indices_as_matrix])
}

# Sinkhorn algorithm
.compute_transport_plan <- function(p, q, M, lam = 10, max_iter = 100){
    # vector
    p_full <- p
    p_full <- p_full / sum(p_full)
    p <- p_full
    # vector, size = (N)
    q <- q / sum(q)
    N <- length(q)
    # vector size = Q
    large_i <- which(p > 0)
    # vector size = Q
    p <- p[large_i]
    # matrix size = (Q, N)
    M <- M[large_i, ]
    # matrix size = (Q, N)
    K <- exp(- lam * M)
    # matrix size = (Q, N)
    u <- matrix(1, length(p), N) / length(p)
    # matrix size = (Q, N)
    K_tilde <- diag(1 / p) %*% K
    for(i_iter in 1:max_iter){
        #         (Q,N)       (N) / (N,Q) %*% (Q, N)
        # size = (Q, N)
        u <- 1 / (K_tilde %*% (q / (t(K) %*% u)))
    }
    #   (N)   (N, Q) %*% (Q, N)
    # size = (N, N)
    v <- q / (t(K) %*% u)
    # Extract only the first column.
    P_lam <- diag(u[, 1]) %*% K %*% diag(v[, 1])
    P_lam_full <- matrix(0, length(p_full), N)
    P_lam_full[large_i, ] <- P_lam
    P_lam_full
}
