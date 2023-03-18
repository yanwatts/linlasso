# Linear Lasso function
# Authors : Yan Watts and Mylene Bedard
# Update : 2023-03-18

stand <- function(data){

  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  ###################### Standardizing data ########################

  # Center data (mean=0)
  mean_col <- apply(data, 2, mean)
  data_cent <- data - matrix(mean_col,n,p,byrow=T)

  # Standardize data (sd=1); we use the MLE estimator of the variance (i.e. the biased version of the variance)
  sd_col <- sqrt(apply(data_cent^2, 2, sum)/n)
  data_stand <- data_cent/matrix(sd_col, n, p, byrow=T)

  return(list(mean=mean_col, sd=sd_col, cent=data_cent, stand=data_stand))
}


correl <- function(resp, pred) {

  n <- nrow(pred)
  p <- ncol(pred)

  ###################### Computing correlations ########################

  # Compute the correlations between y and each of the x's
  #cj <- (t(pred) %*% resp)/n
  cj <- crossprod(pred, resp)/n

  # Compute the correlations between the x's
  #Cjj <- (t(pred) %*% pred)/n
  Cjj <- crossprod(pred)/n

  # Transform vectors so that they lie in the positive half space
  pred_pos <- matrix((cj >= 0) - (cj < 0), n, p, byrow=T) * pred

  # Compute the positive correlations between y and each of the x's
  cj_pos <- crossprod(pred_pos, resp)/n

  # Compute the correlations between the x's (in positive half space)
  # This inverts the sign of some of the correlations (but not all of them)
  #Cjj_pos <- (t(pred_pos) %*% pred_pos)/n
  Cjj_pos <- crossprod(pred_pos)/n

  # Note: the LS estimates associated to inverted x's also invert their sign
  return(list(cj = cj, Cjj=Cjj,  modif= (cj >= 0) - (cj < 0), pred_pos = pred_pos, cj_pos = cj_pos, Cjj_pos = Cjj_pos))

}


xlasso_b2 <- function (resp, pred, cj, Cjj, m) {

  #resp=y.stand ; pred=x.stand ; m=sum(data_cor$cj_pos <= gamma) ; cj = data_cor$cj_pos ; Cjj= data_cor$Cjj_pos

  n <- nrow(pred)
  p <- ncol(pred)

  # 	m <- m * (m < p-1) + (p-1) * ( m >= (p-1))
  m <- m + (p - 1 - m) * (m >= (p-1))

  # Compute the variance term for the complete model
  var_term <- t(cj) %*% solve(Cjj) %*% cj
  #	print(var_term)

  # LS estimate full model
  beta_5 <- solve(t(pred) %*% pred) %*% (t(pred) %*% resp)
  #	beta_5 <- solve(Cjj*n) %*% (cj*n)
  LS_list <- list(beta_5)


  # We start with a vector of 0's. The value 1 will represent the first variable to be discarded from the model (for instance, the value 1 in the 4th position indicates that the 4th variable is the first to exit the model), etc.
  ord <- rep(0, p)

  # Enumerate the explanatory variables in the order in which they exit the model
  pos <- c()

  # Vector of 1's initially; then, when a variable is discarded from the model, a 0 appears at the corresponding position in the vector
  #	delt <- rep(1, p)

  # Enumerate the positions that are still in the model
  upd <- 1:p

  if(m>0){
    for(j in 1:m){
      # Discards an explanatory variable according to the cj's ordering
      cond <- as.vector(cj)  == sort(as.vector(cj))[j]
      #			delt <- delt - delt * cond
      ord[cond] <- j
      pos <- c(pos, (1:p)[cond])
      upd <- (1:p)[-pos]

      # LS estimate of reduced model
      beta_i <- solve(t(pred[,upd]) %*% pred[,upd]) %*% (t(pred[,upd]) %*% resp)
      #beta_i <- solve(Cjj[upd, upd]*n) %*% (cj[upd]*n)

      LS_list <- c(LS_list, list(beta_i))

    }
  }

  # Compute the variance term for the complete model
  # beta_i <- solve(Cjj[upd, upd]*n, cj[upd]*n)
  #LS_list <- list(beta_i)
  var_term <- crossprod(cj[upd], solve(Cjj[upd, upd], cj[upd]))

  if(m<(p-1)){
    for(j in (m+1):(p-1)){

      var_term_i <- rep(-Inf, p)

      for(i in 1:(p-j+1)){

        # Compute the variance term for the reduced model
        var_term_i[upd[i]] <- crossprod(cj[-c(upd[i], pos)], solve(Cjj[-c(upd[i], pos), -c(upd[i], pos)], cj[-c(upd[i], pos)]))

      }

      # Discards an explanatory variable based on variance criterion (we discard the variable that minimizes the drop in the variance)
      cond <- (as.vector(var_term) - var_term_i) == min(as.vector(var_term) - var_term_i)
      #			delt <- delt - delt * cond
      ord[cond] <- j
      pos <- c(pos, (1:p)[cond])
      upd <- (1:p)[-pos]


      # LS estimate reduced model
      beta_i <- solve(crossprod(pred[,upd]), crossprod(pred[,upd], resp))
      #beta_i <- solve(t(pred[,upd]) %*% pred[,upd]) %*% (t(pred[,upd]) %*% resp)
      LS_list <- c(LS_list, list(beta_i))

      var_term <- var_term_i[pos[j]]

    }
  }

  # Includes the last variable to exit the model
  ord[upd] <- p
  pos <- c(pos, upd)

  return(list(Pos=pos, "LS_estimates"=LS_list))

}

pred_err <- function (L = 10, K = 10, y, pred, m){

  #y=y.stand ; pred=x.stand ; m=sum(data_cor$cj_pos <= gamma)

  n <- nrow(pred)
  p <- ncol(pred)

  # 	m <- m * (m < p-1) + (p-1) * ( m >= (p-1))
  m <- m + (p - 1 - m) * (m >= (p-1))

  CVMSE_rep = matrix(nrow = L, ncol = p)
  CVVAR_rep = matrix(nrow = L, ncol = p)
  #CVSE_rep = matrix(nrow = L, ncol = p - m)


  # Complete a total of ? it ? cross-validation cycles
  for(l in 1:L){

    CVMSE <- matrix(nrow=K, ncol = p)
    #if (l %% 1 == 0) {
    #  cat('Iteration:', l, "\p")
    #}
    # Shuffle the observations so as to eventually divide them into K groups of size z
    shuf <- sample(1:n, n, replace=F)
    #		shuf <- 1:n

    #CVMSE <- rep(0, p-m)
    #CVMSE2 <- rep(0, p-m)

    z <- n/K
    z <- round(z)

    # One complete cross-validation cycle (i.e. each one of the K groups is  considered as the test set)
    for(i in 1:K){

      # If unequal groups
      idx = (i-1)*z + 1:z
      if (i == K){idx = ((K-1)*z+1):n ; z = n - (K-1)*z}

      # Remove one of the K groups of size z, and define the rest as training set
      y_train <- y[shuf[-idx]]
      pred_train <- pred[shuf[-idx], ]
      # Compute correlations
      cor_i <- correl(y_train, pred_train)

      # Fit model based on remaining data
      fit_i <- xlasso_b2(y_train, cor_i$pred_pos, cor_i$cj_pos, cor_i$Cjj_pos, m)

      # Compute the residuals sum of squares for each of the p models (of sizes p to 1)
      for(j in 1:p){

        # Select the appropriate explanatory variables (according to fit_i)
        # Change the sign of some explanatory vectors to ensure positive correlations with the response
        pred_test <- matrix(cor_i$modif[sort(fit_i$Pos[j:p])], z, (p+1) - j, byrow=T) * pred[ shuf[idx], sort(fit_i$Pos[j:p])]

        # Compute the z predictions
        y_test <- apply(matrix(fit_i$LS_estimates[[j]], z, (p+1)- j , byrow=T) * pred_test, 1, sum)
        CVMSE[i,j] <- mean((y[shuf[idx]] - y_test)^2)
        #CVMSE2[j] <- CVMSE2[j] + sum((y[shuf[(i-1)*z + 1:z]] - y_test)^4)

      }
    }

    # For each complete cross-validation cycle, we end up with a MSE vector of length p; these ? it ? vectors are put end to end and will eventually be transformed into a matrix.
    # SD of predictions within one complete cross-validation cycle
    #SD_vec <- c(SD_vec, CVMSE2/n - (CVMSE/n)^2)
    CVMSE_rep[l,] <- apply(CVMSE, 2, mean)
    CVVAR_rep[l,] <- apply(CVMSE, 2, var)
    #CVSE_rep[l,] <- sqrt(CVVAR_rep[l,])/sqrt(K)

  }

  MSE = apply(CVMSE_rep, 2, mean)
  SD = sqrt(apply(CVVAR_rep, 2, mean))
  SE = SD/sqrt(K)

  # SD of predictions between the ? it ? cross-validation cycles
  #SD_MSE <- sqrt(apply(matrix(CVMSE_vec, l, p-m, byrow=T), 2, var))

  return(list(MSE=MSE, SD=SD, SE = SE))

}

beta_full <- function(p, var.left, beta.hat, var.names.x){
  beta.hat.matrix = matrix(nrow = p, ncol = 1)
  rownames(beta.hat.matrix) <- var.names.x
  count = 1
  for(i in 1:p){
    if(i %in% var.left){
      beta.hat.matrix[i,1] = beta.hat[count]
      count = count + 1
    }
    else beta.hat.matrix[i,1] = 0
  }
  return(beta.hat.matrix)
}

graph.one.by.one <- function(MSE, gamma, index.1se, table.MSE, K, L, french){


  #par(mfrow=c(1,1))

  minimum <- which.min(MSE)
  col <- c()
  for(k in 1:length(MSE)){
    if(k == minimum) col<- c(col, 2)#red
    else if(k == index.1se) col<- c(col, 3)#green
    else{col<- c(col, "black")}
  }

  if(french){
    main = ""
  }else{
    main = "Mean squared error of rejected variables for the One-by-one procedure
       (Red line indicates final model)"
  }
  plot(MSE, pch = 16, cex = 0.5,
       xlab = paste0("Order in which variables are rejected (gamma = ", gamma,")"),
       ylab = "MSE",
       main = main,
       cex.main = 0.8,
       col = col,
       ylim = c(min(MSE), max(MSE)+0.1),
       cex.lab = 0.8,
       xaxt = 'n')
  axis(1, at=1:length(MSE), labels = seq(1,length(MSE)))
  lines(seq(minimum,length(MSE)), MSE[minimum:length(MSE)], col = "red")
  max.char = max(nchar(row.names(table.MSE)))
  if(minimum == index.1se){
    legend("topleft", legend = "Minimum MSE et 1se MSE",
           col = 2, pch = 16, cex = 0.8)
  }else{
    legend("topleft", legend = c("Minimum MSE", "1se MSE"),
           pch = c(16,16), col = c(2,3), cex = 0.8)
  }
  if(max.char < 4){
    text(MSE, row.names(table.MSE), pos=3, cex = 0.7)
  }else{
    text(MSE + 0.002, row.names(table.MSE), pos=4, cex = 0.7, offset = 0, srt=90)
  }
}


LL <- function(y, x, gamma = 0.2, K = 10, L = 10, plot = F, french = F, max.cor = F){

  #y = as.matrix(diabetes[,11]) ;  x <- model.matrix( ~ .-1, diabetes[,-11]) ; L = 50 ; gamma = 0.2 ; K = 13 ; cor.only = F

  #y = trim32[,1]; x = trim32[,-1]; L = 10; plot = T ; gamma = 0.2 ; K = 10

  x <- as.data.frame(x)
  x <- model.matrix(~., x) ;  x <- x[,-1]
  x <- apply(x, 2,function(x) as.numeric(as.character(x)))

  # Error message if x is not a matrix
  if (is.null(n <- nrow(x))) stop("'x' must be a matrix")

  # Error message if K is not numeric or L is not numeric
  if (K%%1 != 0 | K <= 1) stop("'K' must be an integer greater than 1")
  if (L%%1 != 0) stop("'L' must be an integer")

  n <- nrow(x)
  p <- ncol(x)

  #Standardize the data
  data_stand <- stand(data = data.frame(y,x))
  y.stand <- data_stand$stand[,1] ; x.stand <- data_stand$stand[,2:(p+1)]

  #Calculate the correlations of C (between predictors) and c (between predictor and response)
  data_cor <- correl(y.stand, x.stand)

  x.names = colnames(x)

  if(p >= n){
    c <- data_cor$cj_pos
    #c.sort <- sort(data_cor$cj_pos, decreasing = T)
    c.sort <- c[order(c[,1],decreasing=T),]
    d = round(n/log(n),0)
    #d = round(0.5*n,0)
    print(paste0("The chosen d is ", d,", because of high dimensionnality. This leaves us with ", d," variables for the Linear Lasso"))
    gamma.temp = round(c.sort[d],4)

    x.names = colnames(x)
    idx = which(c > gamma.temp)
    x.stand <- x.stand[,idx]
    x.mod <- x[,idx]
    data_cor <- correl(y.stand, x.stand)

    if(gamma.temp > gamma){
      gamma = round(gamma.temp,2)
    }
  }

  # Treat categorical data by adding small increment so they are not the same variables
  #if(TRUE %in% duplicated(round(data_cor$cj_pos, 8))){
  #  c <- data_cor$cj_pos
  #  var.name <- rownames(c)[(1:p)[duplicated(c)]]
  #  var.name.sub <- substring(var.name,1,nchar(var.name) - 1)
  #  stop(paste(c("Model.matrix was not performed right. Identical variables in design matrix. Variable(s) in question :", var.name.sub), collapse=" "))
  #}


  # Error message if the gamma chosen by the user is too high
  #if (gamma > max(data_cor$cj_pos)) stop("The chosen gamma is too high. All coefficients are zero.")

  # The m that is chosen by the user with the gamma by correlations

  m=sum(data_cor$cj_pos <= gamma)

  # If the gamma chosen by the user is too high, automatically 0 coefficients
  #if (gamma > max(data_cor$cj_pos)){
  #  beta.min = matrix(rep(0,p), nrow = p, ncol = 1)
  #  c = c(data_cor$cj_pos)
  #  names(c) = colnames(x)
  #  rownames(beta.min) = colnames(x)
  #  paste("The chosen gamma is too high. All coefficients are zero.")
  #  return(list(beta.min = beta.min, c.pos = sort(c, decreasing = T)))
  #}

  # Cross-validation to find the optimal variables
  data_CV <- pred_err(L = L, K = K, y=y.stand, pred=x.stand, m=m)

  # The fit on the whole model with the specified m
  fit <- xlasso_b2(y.stand, data_cor$pred_pos, data_cor$cj_pos, data_cor$Cjj_pos, m=m)

  if(p >= n) p = ncol(x.stand)

  # Results
  MSE = data_CV$MSE
  SD = data_CV$SD
  SE = data_CV$SE
  table.MSE <- matrix(nrow = length(MSE), ncol = 4)
  colnames(table.MSE) <- c("Length", "MSE.CV", "SD.CV", "SE.CV")
  rownames(table.MSE) <- colnames(x)[fit$Pos[1:p]]

  #rownames(table.MSE)<-vec_names
  table.MSE[,1] <- seq(nrow(table.MSE),1,-1)
  table.MSE[,2] <- round(as.numeric(MSE),4)
  table.MSE[,3] <- round(as.numeric(SD), 4)
  table.MSE[,4] <- round(as.numeric(SE), 4)

  names(dimnames(table.MSE)) <- c("Rejected variables in order", "")
  #row.names(table.MSE) <- colnames(x)[fit$Pos[(m+1):p]]

  index.1se <- max(which(MSE<=MSE[which.min(MSE)] + SE[which.min(MSE)]))

  # The number of variables left after cutoff
  var.left.cutoff <- sort(fit$Pos[(m+1):p])

  # The number of variables left with MSE
  var.left <- sort(fit$Pos[which.min(MSE):p])

  # The number of variables left with MSE (1se)
  var.left.1se <- sort(fit$Pos[index.1se:p])

  if(ncol(x) >= n){

    x.mod = as.matrix(x.mod)
    beta.hat.min = solve(t(x.mod[,var.left])%*%x.mod[,var.left])%*%t(x.mod[,var.left])%*%y
    beta.min <- beta_full(ncol(x), var.left, beta.hat.min, x.names)
    beta.hat.1se = solve(t(x.mod[,var.left.1se])%*%x.mod[,var.left.1se])%*%t(x.mod[,var.left.1se])%*%y
    beta.1se <- beta_full(ncol(x), var.left.1se, beta.hat.1se, x.names)

    c = c(data_cor$cj_pos)
    names(c) = colnames(x.mod)
  }else{
    x = as.matrix(x)
    beta.hat.min = solve(t(x[,var.left])%*%x[,var.left])%*%t(x[,var.left])%*%y
    beta.min <- beta_full(p, var.left, beta.hat.min, x.names)
    beta.hat.1se = solve(t(x[,var.left.1se])%*%x[,var.left.1se])%*%t(x[,var.left.1se])%*%y
    beta.1se <- beta_full(p, var.left.1se, beta.hat.1se, x.names)

    c = c(data_cor$cj_pos)
    names(c) = x.names
  }


  if(plot) graph.one.by.one(MSE, gamma, index.1se, table.MSE, K, L, french)

  #print(paste("Variables with inferior correlation to gamma =", gamma, ":", paste(colnames(x)[var.left.cutoff], collapse = " + ")))
  #print(paste(c("Variables left after cutoff =", gamma, "are", colnames(x)[var.left.cutoff]), collapse=" ", sep = "\n"))
  return(list(m=m, gamma=gamma,
              c.pos = sort(c, decreasing = T),
              table.MSE = table.MSE,
              "Variables with minimum MSE" = colnames(x)[var.left],
              beta.min = beta.min,
              "Variables with 1se of minimum MSE" = colnames(x)[var.left.1se],
              beta.1se = beta.1se))

}

#crime_data = read.csv( "C:\\Users\\waya7500\\Desktop\\Memoire\\Code R\\Exemples\\Datasets\\crime_data.txt",sep = "")
#crime_data = crime_data[,-2]
#y = crime_data[,1] ; x = crime_data[,2:6] ; L = 10 ; gamma = 0.2 ; K = 3 ; cor.only = F
#LL(y,x, K = 10, L = 50)
