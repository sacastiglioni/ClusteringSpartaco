main <- function(x, coordinates,
                 K,
                 column.labels,
                 Delta.constr = 10,
                 max.iter = 10^3,
                 estimate.iterations = 10,
                 conv.criterion = NULL,
                 input.values = NULL,
                 save.options = NULL,
                 verbose = F,
                 verbose.parallel.label = NULL,
                 lambda.ridge = 0,
                 lambda.lasso = 0){
  Dist <- as.matrix(stats::dist(coordinates))
  Ds <- column.labels
  R <- length(unique(column.labels))
  if(is.null(input.values)){
    cur.Cs <- best.Cs <- sample(1:K, size = nrow(x), replace = T)
    cur.phi <- best.phi <- runif(R, 1, 5)
    cur.mu <- best.mu <- matrix(runif(K*R, 1, 10), K, R)
    cur.tau <- best.tau <- matrix(runif(K * R, 1e-7, Delta.constr), K, R)
    cur.xi <- best.xi <- Delta.constr - best.tau
    cur.alpha <- best.alpha <- matrix(runif(K*R, 1, 3), K, R)
    cur.beta <- best.beta <- matrix(runif(K*R, 1, 3), K, R)} else {
      cur.Cs <- best.Cs <- input.values$Cs
      cur.phi <- best.phi <- input.values$phi
      cur.mu <- best.mu <- input.values$mu
      cur.tau <- best.tau <- input.values$tau
      cur.xi <- best.xi <- Delta.constr - best.tau
      cur.alpha <- best.alpha <- input.values$alpha
      cur.beta <- best.beta <- input.values$beta}
  Cs <- matrix(0, nrow(x), max.iter)
  Cs[,1] <- cur.Cs

  Uglob <- list()
  Dglob <- numeric(ncol(x))
  sapply(1:R, function(r){
    eigK <- eigen(exp(-Dist[Ds == r, Ds == r]/cur.phi[r]))
    Uglob[[r]] <<- eigK$vec
    Dglob[Ds == r] <<- eigK$val
  })
  ll <- rep(-1e+40, max.iter)
  logL.values <- matrix(0, K, R)
  i <- 1
  if(!is.null(conv.criterion)) counter.conv <- 0
  if(verbose == "progress"){
    progressr::handlers(global = T)
    P <- progressr::progressor(along = 1:max.iter)}
  while(T){
    if(i == max.iter) break
    i <- i + 1
    #if(!verbose) svMisc::progress(i/max.iter*100)
    #P(message = sprintf("Added %g", i))
    #if(verbose == "parallel" & i %% 10 == 0) cat(paste(verbose.parallel.label,": Iteration ",i," of ",max.iter,"\n",sep=""))
    if(verbose == "progress") P()
    if(verbose == "full") cat(paste("---Iteration",i,"\n"))

    # ---M Step
    if(verbose == "full") cat("M Step/")
    goodK <- sort(unique(cur.Cs))
    goodR <- sort(unique(Ds))
    sapply(goodR, function(r){
      traceDelta_r <- Delta.constr * sum(Ds == r)
      sapply(goodK, function(k){
        #cat(k, r, sum(cur.Cs == k), sum(Ds == r), "\n")
        #cat(cur.mu[k,r], "\n")
        estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, Ds == r],
                                                                                         traceDelta = traceDelta_r,
                                                                                         U = Uglob[[r]],
                                                                                         d = Dglob[Ds == r],
                                                                                         mu0 = cur.mu[k,r],
                                                                                         alpha0 = cur.alpha[k,r],
                                                                                         beta0 = cur.beta[k,r],
                                                                                         tau0 = cur.tau[k,r],
                                                                                         maxit = estimate.iterations,
                                                                                         lambda.ridge = lambda.ridge,
                                                                                         lambda.lasso = lambda.lasso
        )
        cur.mu[k,r] <<- estimation.parameters$mu
        cur.tau[k,r] <<- estimation.parameters$tau
        cur.xi[k,r] <<- estimation.parameters$xi
        cur.alpha[k,r] <<- estimation.parameters$alpha
        cur.beta[k,r] <<- estimation.parameters$beta
      })
      cur.phi[r] <<- updatePhi_r_marginal(x = x[,Ds == r],
                                          Cs = cur.Cs,
                                          Dist = Dist[Ds == r, Ds == r],
                                          Mu = cur.mu[,r],
                                          Tau = cur.tau[,r],
                                          Xi = cur.xi[,r],
                                          Alpha = cur.alpha[,r],
                                          Beta = cur.beta[,r],
                                          phi.old = cur.phi[r])
      EigenK <- eigen(exp(-Dist[Ds == r, Ds == r]/cur.phi[r]))
      Uglob[[r]] <<- EigenK$vec
      Dglob[Ds == r] <<- EigenK$val
    })

    # ---CE Step
    if(verbose == "full") cat("CE Step/")
    cur.cs <- RowClustering(x = x, Ds = Ds, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi, Uglob = Uglob, Dglob = Dglob)
    cur.Cs <- cur.cs$allocation
    goodK <- sort(unique(cur.Cs))
    goodR <- sort(unique(Ds))
    sapply(goodR, function(r){
      sapply(goodK, function(k){
        logL.values[k,r] <<- logL.Cocluster(x = x[cur.Cs == k, Ds == r],
                                            Mu = cur.mu[k,r],
                                            Tau = cur.tau[k,r],
                                            Xi = cur.xi[k,r],
                                            Alpha = cur.alpha[k,r],
                                            Beta = cur.beta[k,r],
                                            U = Uglob[[r]],
                                            d = Dglob[Ds == r],
                                            lambda.ridge = 0,
                                            lambda.lasso = 0)
      })
    })
    ll[i] <- sum(logL.values)
    Cs[,i] <- cur.Cs

    if(ll[i] == max(ll)){
      best.phi <- cur.phi
      best.mu <- cur.mu
      best.tau <- cur.tau
      best.xi <- cur.xi
      best.alpha <- cur.alpha
      best.beta <- cur.beta
      best.Cs <- cur.Cs
    }

    if(verbose == "full"){
      cat(paste("diff(logL) =",round(diff(ll)[i-1],5),"\n"))
      cat(paste("Row clusters size =", paste(table(cur.Cs), collapse = ", "),"\n"))
      cat(paste("Column clusters size =", paste(table(Ds), collapse = ", "),"\n"))
    }

    if(!is.null(conv.criterion)){
      if(ll[i] >= ll[i-1] & (ll[i] - ll[i-1] < conv.criterion$epsilon)){
        counter.conv <- counter.conv + 1
        if(counter.conv == conv.criterion$iterations){
          cat("Converged\n")
          break}
      } else {
        counter.conv <- 0
      }
    }

    # save the result in the given location
    if(!is.null(save.options)){
      if(i %% save.options$after == 0){
        ICL <- max(ll) - nrow(x)*K - ncol(x)*R - .5*(4*K*R+R)*log(nrow(x) * ncol(x))
        results <- list(
          mu = best.mu,
          tau = best.tau,
          xi = best.xi,
          alpha = best.alpha,
          beta = best.beta,
          phi = best.phi,
          Cs = best.Cs,
          Ds = Ds,
          logL = ll[c(2:i)],
          ICL = ICL,
          x = x,
          coordinates = coordinates
        )
        class(results) <- "spartaco"
        save(results, file = save.options$file.name)
      }
    }
  }

  ICL <- max(ll) - nrow(x)*log(K) - ncol(x)*log(R) - .5*(4*K*R+R)*log(nrow(x) * ncol(x))
  # save the result in the given location
  if(!is.null(save.options)){
    results <- list(
      mu = best.mu,
      tau = best.tau,
      xi = best.xi,
      alpha = best.alpha,
      beta = best.beta,
      phi = best.phi,
      Cs = best.Cs,
      Ds = best.Ds,
      logL = ll[c(2:i)],
      ICL = ICL,
      x = x,
      coordinates = coordinates)
    class(results) <- "spartaco"
    save(results, file = save.options$file.name)
  }
  results <- list(
    mu = best.mu,
    tau = best.tau,
    xi = best.xi,
    stn.ratio = best.tau/best.xi,
    alpha = best.alpha,
    beta = best.beta,
    phi = best.phi,
    Cs = best.Cs,
    Ds = Ds,
    logL = ll[c(2:i)],
    ICL = ICL,
    x = x,
    coordinates = coordinates
  )
  class(results) <- "spartaco"
  return(results)
}

