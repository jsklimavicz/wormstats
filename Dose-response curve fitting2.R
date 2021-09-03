library(tidyverse)
library(ggplot2)
library(foreach)
library(doParallel)
library(plyr)
library(MASS)

curve_fit <- function(alive_counts, dead_counts, conc, plateID, param_count=2, n_trials = 1000, ctrl_ave=.95){
   #fitting functions
   ll_worms2 <- function(b, probs, conc, sigma = 1000){
      b0 <- b[1]
      b1 <- b[2]
      p_sum <- sum(probs)
      p_conc_sum <- sum(probs*conc)
      ll <- -(b0^2 + b1^2)/(2*sigma^2) + 
         b0*p_sum + 
         b1*p_conc_sum -
         sum(log(1 + exp(b0+b1*conc)))
      return(ll)
   }
   ll_worms2_grad <- function(b, probs, conc, sigma = 1000){
      b0 <- b[1]
      b1 <- b[2]
      p_sum <- sum(probs)
      p_conc_sum <- sum(probs*conc)
      xi = exp(b0+b1*conc)
      l = xi/(xi+1)
      g1 = -b0/sigma^2 + sum(probs) - sum(l)
      g2 = -b1/sigma^2 + sum(conc*probs) - sum(conc*l)
      return(c(g1,g2))
   }
   ll_worms3 <- function(b, probs, conc, sigma = 1000, gamma_param = c(1.5, 0.5)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      if (b2 <= 0) {return(-1e10)}
      alpha <- 1+exp(b0+b1*conc)
      ll <-  -(b0^2 + b1^2)/(2*sigma^2) +(gamma_param[1] -1)* log(b2) - 
         gamma_param[2]*b2 +
         sum(probs*log(alpha^b2-1)) - 
         b2 * sum(log(alpha))
      return(ll)
   }
   ll_worms3_grad <- function(b, probs, conc, sigma = 1000, gamma_param = c(1.5, 0.5)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      # if (b2 <= 0) {return(-1e308)}
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      lambda <- log(alpha)
      mu = (1+xi)^(b2-1)
      nu = mu*(1+xi)
      m = probs*xi * mu / (nu - 1)
      l = xi/alpha
      
      g1 = b0/sigma^2 + b2 * sum(m) - b2*sum(l)
      g2 = b1/sigma^2 + b2 * sum(conc*m) - b2*sum(conc*l)
      g3 = (gamma_param[1]-1)/b2 - gamma_param[2] + sum(probs*nu*lambda/(nu -1)) - sum(lambda)
      
      return(c(g1,g2,g3))
   }
   ll_worms4 <- function(b, probs, conc, sigma = 1000, gamma_param = c(1.5, 0.5), beta_param = c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      b3 <- b[4]
      if (b2 <= 0) {return(-1e10)}
      if (b3 <= 0 | b3 >= 1) {return(-1e10)}
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      a2 = alpha^2
      lambda = log(alpha)
      if (min(a2)-b3 <=0) {return(-1e10)}
      
      ll <- -(b0^2 + b1^2)/(2*sigma^2) + #MVN prior
         (gamma_param[1] -1)* log(b2) - gamma_param[2]*b2 + #gamma prior
         (beta_param[1] -1)* log(b3) + (beta_param[2] -1)*log(1-b3)  + #\beta prior
         sum(probs*log(a2-b3)) + sum((1-probs)*log(b3)) - 
         b2 * sum(lambda)
      
      
      return(ll)
   }
   ll_worms4_grad <- function(b, probs, conc,  sigma = 1000, gamma_param = c(1.5, 0.5), beta_param = c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      b3 <- b[4]
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      lambda <- log(alpha)
      mu = (1+xi)^(b2-1)
      nu = mu*(1+xi)
      d = (nu - b3)
      m = probs*xi * mu / d
      l = xi/alpha
      
      g1 = b0/sigma^2 + b2 * sum(m) - b2*sum(l)
      g2 = b1/sigma^2 + b2 * sum(conc*m) - b2*sum(conc*l)
      g3 = (gamma_param[1]-1)/b2 - gamma_param[2] + sum(probs*nu*lambda/d) - sum(lambda)
      g4 = (beta_param[1]-1)/b3 -(beta_param[2]-1)/(1-b3) - sum(probs/d) + sum((1-probs)/b3)
      
      return(c(g1,g2,g3,g4))
   }
   ll_worms3_wb <- function(b, probs, conc, sigma = 1000,  weibull_param=c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b3 <- b[3]
      if (b3 <= 0 | b3 > 1.001) {return(-1e10)}
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      lambda = log(alpha)
      if (min(alpha)-b3 <=0) {return(-1e10)}
      wk = weibull_param[1]
      wl = weibull_param[2]*(wk/(wk-1))^(1/wk)
      
      ll <- -(b0^2 + b1^2)/(2*sigma^2) - #MVN prior
         (b3/wl)^wk + (wk-1)*log(b3) + #weibull prior
         sum(probs*log(alpha-b3)) + sum((1-probs)*log(b3)) - sum(lambda)
      
      return(ll)
   }
   ll_worms3_grad_wb <- function(b, probs, conc,  sigma = 1000,  weibull_param=c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b2<-1
      b3 <- b[3]
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      lambda <- log(alpha)
      d = (alpha - b3)
      m = probs*xi / d
      l = xi/alpha
      wk = weibull_param[1]
      wl = weibull_param[2]*(wk/(wk-1))^(1/wk)
      
      g1 = -b0/sigma^2 + b2 * sum(m) - b2*sum(l)
      g2 = -b1/sigma^2 + b2 * sum(conc*m) - b2*sum(conc*l)
      g4 = (wk-1)/b3 - (wk/b3)*(b3/wl)^wk - sum(probs/d) + sum((1-probs)/b3)
      
      return(c(g1,g2,g4))
   }
   ll_worms4_wb <- function(b, probs, conc, sigma = 1000, gamma_param = c(1.5, 0.5),  weibull_param=c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      b3 <- b[4]
      if (b2 <= 0) {return(-1e10)}
      if (b3 <= 0 | b3 > 1.05) {return(-1e10)}
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      a2 = alpha^2
      lambda = log(alpha)
      if (min(a2)-b3 <=0) {return(-1e10)}
      wk = weibull_param[1]
      wl = weibull_param[2]*(wk/(wk-1))^(1/wk)
      
      ll <- -(b0^2 + b1^2)/(2*sigma^2) + #MVN prior
         (gamma_param[1] -1)* log(b2) - gamma_param[2]*b2 - #gamma prior
         (b3/wl)^wk + (wk-1)*log(b3) + #weibull prior
         sum(probs*log(a2-b3)) + sum((1-probs)*log(b3)) - 
         b2 * sum(lambda)
      
      return(ll)
   }
   ll_worms4_grad_wb <- function(b, probs, conc,  sigma = 1000, gamma_param = c(1.5, 0.5),  weibull_param=c(2,1)){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      b3 <- b[4]
      xi = exp(b0+b1*conc)
      alpha <- 1+xi
      lambda <- log(alpha)
      mu = (1+xi)^(b2-1)
      nu = mu*(1+xi)
      d = (nu - b3)
      m = probs*xi * mu / d
      l = xi/alpha
      wk = weibull_param[1]
      wl = weibull_param[2]*(wk/(wk-1))^(1/wk)
      
      g1 = -b0/sigma^2 + b2 * sum(m) - b2*sum(l)
      g2 = -b1/sigma^2 + b2 * sum(conc*m) - b2*sum(conc*l)
      g3 = (gamma_param[1]-1)/b2 - gamma_param[2] + sum(probs*nu*lambda/d) - sum(lambda)
      g4 = (wk-1)/b3 - wk/b3*(b3/wl)^wk - sum(probs/d) + sum((1-probs)/b3)
      
      return(c(g1,g2,g3,g4))
   }
   #fitting goodness
   curve_ll <- function(b, probs, conc, sigma = 1000, gamma_param = c(1.5, 0.5), weibull_param = c(2,1)) {
      v<- c(0,0,1,1)
      v[1:length(b)]<-b
      return (ll_worms4_wb(b, probs, conc, sigma = sigma, gamma_param = gamma_param,  weibull_param=weibull_param))
      
   }
   curve_R2 <- function(b, probs, conc) {
      SS_total <- sum((probs-.5)^2)
      v<- c(0,0,1,1)
      v[1:length(b)]<-b
      fitted_p <- 1-v[4]/(1+exp(v[1]+v[2]*conc))^v[3]
      ss_fit <- sum((probs-fitted_p)^2)
      return(1-ss_fit/SS_total)
   }

   #detect cores for parallel
   numCores <- detectCores()
   registerDoParallel(numCores)
   #calculate "correct" number of trials close to and >= requested number that allows for blocking
   max_blocks = 31
   block_alloc = 8
   if ((n_trials/numCores)/block_alloc < max_blocks) {
      blocks <- floor((n_trials/numCores)/block_alloc)
   } else {
      blocks <- max_blocks
   }
   if (blocks == 0) {
      bs=0
      last_bs <- ceiling(n_trials/numCores)* numCores
      n_trials2 <- last_bs
   } else {
      bs <- floor(floor(n_trials/blocks))
      last_bs <- n_trials - (bs*blocks)
      last_bs <- ceiling(last_bs/numCores)* numCores
      n_trials2 <- blocks* bs + last_bs 
   }
   if (n_trials2 != n_trials) {
      print(paste("Using n_trials =", n_trials2, "instead of the requested", n_trials, "to optimize CPU usage."))
   }
   n_trials <- n_trials2
   
   #get log(conc)
   log_conc <- log(conc)
   
   #Make correlated beta variables
   rho <- 0.25
   covar_mat <- matrix(0, nrow = length(plate), ncol = length(plate))
   for (p in unique(plate)){
      curr_plate <- which(plate %in% p)
      for (i in curr_plate) {
         for (j in curr_plate){
            if (i == j) {
               covar_mat[i,i] = 1
            } else {
               covar_mat[i,j] = rho/2^(abs(i-j)-1)
               # covar_mat[i,j] = rho
            }
         }
      }
   }
   probs <- as.data.frame(mvrnorm(n = n_trials, mu = rep(0, length(plate)), Sigma = covar_mat, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
   for (i in 1:length(alive_counts)){
      probs[,i] <- pnorm(probs[,i])
      probs[,i] <- qbeta(probs[,i], dead_counts[i]+.1, alive_counts[i]+.1)
   }
   
   all_params = data.frame(matrix(ncol = 3, nrow = 0))
   for (block in 1:(blocks+1)){
      print(paste("On block", block, "of", blocks))
      bs_curr = ifelse(block == (blocks + 1), last_bs,bs)
      params <- foreach (i=1:bs_curr, .combine = rbind) %dopar% {
         iter = (block-1)*bs + i
         if (param_count<=2){
            p1 <- optim(c(0,1), fn = ll_worms2, gr = ll_worms2_grad, probs=probs[iter,],
                        conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
            return(p1$par)
         } else if (param_count<=3){
            p1 <- optim(c(0,1), fn = ll_worms2, gr = ll_worms2_grad, probs=probs[iter,],
                        conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
            p2 <- optim(c(p1$par,ctrl_ave), fn = ll_worms3_wb, gr = ll_worms3_grad_wb, probs=probs[iter,], 
                        conc =  log_conc, weibull_param=c(2,ctrl_ave), 
                        method = 'BFGS', control = list(fnscale=-1))
            return(c(p1$par, p2$par))
         } else {
            p1 <- optim(c(0,1), fn = ll_worms2, gr = ll_worms2_grad, probs=probs[iter,],
                        conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
            p2 <- optim(c(p1$par,ctrl_ave), fn = ll_worms3_wb, gr = ll_worms3_grad_wb, probs=probs[iter,], 
                        conc =  log_conc, weibull_param=c(2,ctrl_ave), 
                        method = 'BFGS', control = list(fnscale=-1))
            out <- optim(c(p1$par,1,ctrl_ave), fn = ll_worms4_wb,  gr = ll_worms4_grad_wb,
                         probs=probs[iter,], conc =  log_conc, weibull_param=c(2,ctrl_ave),
                         method = 'BFGS', control = list(fnscale=-1,maxit=2000))
            return(c(p1$par, p2$par,out$par))
         }
      }
      if (empty(all_params)){
         all_params = params
      } else {
         all_params = rbind(all_params, params)
      }
   }
   params2 <- as.data.frame(all_params)
   if  (param_count>=4){
      return(list(p2=params2[,1:2], p3=params2[,3:5], p4=params2[,6:9]))
   } else if (param_count>=3){
      return(list(p2=params2[,1:2], p3=params2[,3:5]))
   } else {
      return(list(p2=params2))
   }
}



repcount = 4
conc = rev(rep(c(128,64, 32,16,8,4,2,1,.5, .25),repcount))
log_conc <- log(conc)
bp_low_mort = 1/(1+exp(1.5-1.5*log_conc)) #baseline_mortality
high_mort_baseline <- 0.3#30% baseline mortality

worm_counts = rbinom(n = 10*repcount, size = 30, prob = .5)
d_low = rbinom(n = 10*repcount, size = worm_counts, prob = bp_low_mort)
a_low = worm_counts-d_low

#now kill 30% of surviving worms... :(
d_high = rbinom(n = 10*repcount, size = a_low, prob = high_mort_baseline)
d_high = d_high + d_low
a_high = worm_counts-d_high
plate = rep(seq(1,repcount),each=10)

ctrl_ct = rbinom(n = 8, size = 30, prob = .5)
ctrl_mort = rbinom(n = 8, size = ctrl_ct, prob = high_mort_baseline)
ctrl_ave = mean((ctrl_ct - ctrl_mort)/ctrl_ct)

low <- curve_fit(a_low, d_low, conc, plate, param_count=3, n_trials = 100, ctrl_ave=.95)
high <- curve_fit(a_high, d_high, conc, plate, param_count=3, n_trials = 100, ctrl_ave=ctrl_ave)

low <- list(p2=data.frame(matrix(ncol = 2, nrow = 0)), p3=data.frame(matrix(ncol = 3, nrow = 0)))
high = list(p2=data.frame(matrix(ncol = 2, nrow = 0)), p3=data.frame(matrix(ncol = 3, nrow = 0)))
for (i in 1:20){
   print(paste("Rep",i))
   worm_counts = rbinom(n = 10*repcount, size = 30, prob = .5)
   d_low = rbinom(n = 10*repcount, size = worm_counts, prob = bp_low_mort)
   a_low = worm_counts-d_low
   
   #now kill 30% of surviving worms... :(
   d_high = rbinom(n = 10*repcount, size = a_low, prob = high_mort_baseline)
   d_high = d_high + d_low
   a_high = worm_counts-d_high
   plate = rep(seq(1,repcount),each=10)
   
   ctrl_ct = rbinom(n = 8, size = 30, prob = .5)
   ctrl_mort = rbinom(n = 8, size = ctrl_ct, prob = high_mort_baseline)
   ctrl_ave = mean((ctrl_ct - ctrl_mort)/ctrl_ct)
   low_curr <- curve_fit(a_low, d_low, conc, plate, param_count=3, n_trials = 500, ctrl_ave=.95)
   high_curr <- curve_fit(a_high, d_high, conc, plate, param_count=3, n_trials = 500, ctrl_ave=ctrl_ave)
   if (empty(high$p2)) {
      high = high_curr
      low = low_curr
   } else {
      low$p2 = rbind(low$p2, low_curr$p2)
      low$p3 = rbind(low$p3, low_curr$p3)
      high$p2 = rbind(high$p2, high_curr$p2)
      high$p3 = rbind(high$p3, high_curr$p3)
   }
}


minx = log(min(conc)) - 2
maxx = log(max(conc)) + 2
xs = seq(minx, maxx, length = 257)

low_lines <- foreach (i = 1:dim(low$p2)[1], .combine=cbind) %dopar% {
   b = low$p2[i,]
   b = unlist(b)
   num =  1- 1/(1+exp(b[1]+b[2]*xs))
}
high_lines <- foreach (i = 1:dim(high$p3)[1], .combine=cbind) %dopar% {
   b = high$p3[i,]
   b = unlist(b)
   num =  1- b[3]/(1+exp(b[1]+b[2]*xs))
}
l <- adply(low_lines, 1, function(d) quantile(d, c(0.025, .5, .975)))
h <- adply(high_lines, 1, function(d) quantile(d, c(0.025, .5, .975)))
l$X1 <- xs
h$X1 <- xs

med = median(high$p3[,length(high$p3)])


livedead_low = data.frame(x=log_conc, y = d_low/(a_low+d_low))
livedead_high = data.frame(x=log_conc, y = d_high/(a_high+d_high))
ggplot() +
   geom_ribbon(data = l, aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'steelblue1', alpha=0.3) +
   geom_line(data = l,aes(x = X1, y = `50%`), color = 'blue',alpha=0.8)+ 
   geom_jitter(data=livedead_low, aes(x=x,y=y), color = 'blue', width=0.1, height = 0,alpha=0.7) +
   geom_ribbon(data = h, aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'tomato1', alpha=0.3) +
   geom_line(data = h,aes(x = X1, y = `50%`), color = 'red',alpha=0.8)+ 
   geom_jitter(data=livedead_high, aes(x=x,y=y), color = 'red', width=0.1, height = 0,alpha=0.7) +
   ylab("Mortality") +
   scale_x_continuous(breaks = log(c(256,128,64, 32,16,8,4,2,1,.5, .25,.125)), 
                    labels = c(256,128,64, 32,16,8,4,2,1,.5, .25,.125))+
   xlab("Concentration (ppm)")+ 
   theme(
      panel.grid.minor.x = element_blank() 
   )


calc_EC <- function(param_df, cc = c(0.1,0.5,0.9), conf = 0.95){
   conf_low = 0.5 - conf/2
   conf_high = 0.5 + conf/2
   c = matrix(cc, nrow=1)
   if (dim(param_df)[2]==4){
      v1=matrix(rep(param_df[,1],length(cc)),ncol=length(cc))
      v2=matrix(rep(param_df[,2],length(cc)),ncol=length(cc))
      v3=matrix(rep(param_df[,4],length(cc)),ncol=length(cc))
      c2=matrix(rep(1-c,dim(param_df)[1]),ncol=length(cc), byrow=TRUE)
      vals <- (log(c2^(-1/v3)-1)-v1)/v2
   }
   if (dim(param_df)[2] %in% c(2,3)){
      v1=matrix(rep(param_df[,1],length(cc)),ncol=length(cc))
      v2=matrix(rep(param_df[,2],length(cc)),ncol=length(cc))
      c2=matrix(rep(1-c,dim(param_df)[1]),ncol=length(cc), byrow=TRUE)
      vals <- (log(1/c2-1)-v1)/v2
   }
   m <- adply(vals, 2, function(d) quantile(d, c(0.025, .5, .975)))
   m$X1 <- cc*100
   names(m)[names(m) == "X1"] <- "EC"
   vals=data.frame(vals)
   names(vals)<- paste("EC",cc*100,sep="")
   return(list(EC=m, vals=vals))
}

j = calc_EC(low$p2)
k = calc_EC(high$p3)

ggplot() + 
   geom_density(data = j$vals, aes(x=EC50, y = ..scaled..), color="blue") + 
   geom_density(data = k$vals, aes(x=EC50, y = ..scaled..), color="red") +
   scale_x_continuous(breaks = log(c(8,6,5,4,3,2.5,2,1.5,1.25,1,.5, .25)), 
                       labels = c(8,6,5,4,3,2.5,2,1.5,1.25,1,.5, .25)) +
   ylab("Scaled Density") + 
   xlab("LC50 (ppm)") + 
   theme(
      panel.grid.minor.x = element_blank() 
   )

lines2<-data.frame(low_lines)[seq(1,10000,25)]
lines2$X <- xs
lines2 = lines2%>% pivot_longer(cols = starts_with("result"),names_to = "result", values_to = "Y") 
lines3<-data.frame(high_lines)[seq(1,10000,25)]
lines3$X <- xs
lines3 = lines3%>% pivot_longer(cols = starts_with("result"),names_to = "result", values_to = "Y") 

# 
# lines2 %>% ggplot(aes(x=X, y=Y)) + 
#    geom_density_2d_filled(bins=15, contour_var = "ndensity") + 
#    theme(legend.position = "none") + 
#    geom_line(data = l, aes(x = X1, y = `50%`))+ 
#    geom_point(data=livedead, aes(x=x,y=y)) +
#    ylab("Mortality") +
#    xlab("ln(Concentration)")

ggplot() + 
   geom_line(data = lines2, aes(x=X, y=Y,group=result), size=0.25, alpha = 0.2, color = 'blue') + 
   geom_line(data = lines3, aes(x=X, y=Y,group=result), size=0.25, alpha = 0.2, color = 'red') + 
   theme(legend.position = "none") + 
   # geom_line(data = l, aes(x = X1, y = `50%`))+ 
   # geom_point(data=livedead_low, aes(x=x,y=y)) +
   ylab("Mortality") +
   scale_x_continuous(breaks = log(c(256,128,64, 32,16,8,4,2,1,.5, .25,.125)), 
                      labels = c(256,128,64, 32,16,8,4,2,1,.5, .25,.125))+
   xlab("Concentration (ppm)")

ny=1024
mat = matrix(data = NA, nrow = ny, ncol = dim(lines)[1])
for (xval in 1:dim(lines)[1]){
   k = density(unlist(lines[xval,]), from=0, to=1, n = ny)
   mat[xval,] =round(k$y/max(k$y),4)
}
matdf = data.frame(mat) 
names(matdf) <- seq(0,1,length=ny)
matdf$x <- xs
matdf = matdf%>% pivot_longer(cols = starts_with(c("0","1")),names_to = "y", values_to = "z") 
matdf$y <- as.numeric(matdf$y)
matdf$z <- as.numeric(matdf$z)
ggplot() + 
   geom_tile(data = matdf, aes(x, y, fill= 1-z)) + 
   scale_fill_viridis(discrete=FALSE, option = "magma")+
   theme(legend.position = "none") +
   geom_line(data = l, aes(x = X1, y = `50%`), color = 'green')+ 
   geom_point(data=livedead, aes(x=x,y=y), color = 'green') +
   ylab("Mortality") +
   xlab("ln(Concentration)")



minx = log(min(conc)) - 2
maxx = log(max(conc)) + 2
xs = seq(minx, maxx, length = 141)
df = data.frame(x = xs)
df <- df %>% mutate(y = b[4]/((1+exp(b[1]+b[2]*x))^b[3]))


ggplot() +
   geom_line(data = df, aes(x=x, y=y)) + 
   geom_point(data=livedead, aes(x=x,y=y))

# livedead = data.frame(x=log_conc, y = d/(a+d))
livedead = data.frame(x=log_conc, y = probs)
minx = log(min(conc)) - 2
maxx = log(max(conc)) + 2
xs = seq(minx, maxx, length = 141)
df = data.frame(x = xs)

p1 <- optim(c(3,-1.6), fn = ll_worms2, gr = ll_worms2_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
p2 <- optim(c(p1$par,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
p2b <- optim(c(p1$par,ctrl_ave), fn = ll_worms3_wb, gr = ll_worms3_grad_wb, probs=probs[i,], 
             conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
p3 <- optim(c(p2$par,ctrl_ave), fn = ll_worms4,  gr = ll_worms4_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
p4 = optim(c(1,-1,1,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #bad
p5 = optim(c(p1$par,1,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #bad
p6 =optim(c(p1$par,1,.999), fn = ll_worms4,  gr = ll_worms4_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #good
p7 =optim(c(p2$par,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #bad
p8 =optim(c(p2$par,.999), fn = ll_worms4,  gr = ll_worms4_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #good
p9 = optim(c(1,-1,1,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #poor
p10 = optim(c(p1$par,1,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #good
p11 =optim(c(p1$par,1,1), fn = ll_worms4_wb,  gr = ll_worms4_grad_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #bad
p12 =optim(c(p2$par,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #bad
p13 =optim(c(p2$par,1), fn = ll_worms4_wb,  gr = ll_worms4_grad_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)) #good


b<- c(0,0,1,1)
v <- p2b
b[1] = v$par[1]
b[2] = v$par[2]
b[4] = v$par[3]
b[1:length(v$par)]<-v$par
df <- df %>% mutate(y = 1- b[4]/(1+exp(b[1]+b[2]*x))^b[3])
df <- df %>% mutate(y = b[4]/(1+exp(b[1]+b[2]*x))^b[3])
curve_R2(b,probs=-probs,conc=log_conc)
ggplot() +
   geom_line(data = df, aes(x=x, y=y)) + 
   geom_point(data=livedead, aes(x=x,y=y))





library(microbenchmark)
j <- microbenchmark(
   # "ll3" = optim(c(p1$par,1), fn = ll_worms3, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll3g" =optim(c(p1$par,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   # "L-ll3" = optim(c(p1$par,1), fn = ll_worms3, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),
   "L-ll3g" =optim(c(p1$par,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),
   "S-ll3g" =optim(c(1,-1,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "S-L-ll3g" =optim(c(1,-1,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),
   # "ll2" = optim(c(3,-2), fn = ll_worms2, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll2g" = optim(c(1,-1), fn = ll_worms2,gr = ll_worms2_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   # "L-ll2" = optim(c(3,-2), fn = ll_worms2, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),
   "L-ll2g" = optim(c(1,-1), fn = ll_worms2,gr = ll_worms2_grad, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),
   "S-ll4" = optim(c(1,-1,1,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll4" = optim(c(p1$par,1,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll4-g" =optim(c(p1$par,1,.999), fn = ll_worms4,  gr = ll_worms4_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "H-ll4" =optim(c(p2$par,1), fn = ll_worms4, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "H-ll4-g" =optim(c(p2$par,.999), fn = ll_worms4,  gr = ll_worms4_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "S-ll4-WB" = optim(c(1,-1,1,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll4-WB" = optim(c(p1$par,1,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "ll4-g-WB" =optim(c(p1$par,1,1), fn = ll_worms4_wb,  gr = ll_worms4_grad_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "H-ll4-WB" =optim(c(p2$par,1), fn = ll_worms4_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   "H-ll4-g-WB" =optim(c(p2$par,1), fn = ll_worms4_wb,  gr = ll_worms4_grad_wb, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
   times=100,
   control = list(warmup=10),
   setup =  list(probs <- rbeta(length(a), a+.1, d+.1), 
                 p1 <-  optim(c(3,-2), fn = ll_worms2,gr = ll_worms2_grad, probs=probs, conc =  log_conc, method = 'BFGS', control = list(fnscale=-1)),
                 p2 <- optim(c(p1$par,1), fn = ll_worms3, gr = ll_worms3_grad, probs=probs, conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)))
)









params <- foreach (i=1:n_trials, .combine = rbind) %dopar% {
   out <- optim(c(-1,1), ll_worms2, probs=probs[i,], conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1))
   return(out$par)
}
params <- as.data.frame(params)

minx = log(min(conc)) - 2
maxx = log(max(conc)) + 2
xs = seq(minx, maxx, length = 141)

lines <- foreach (i = 1:n_trials, .combine=cbind) %dopar% {
   num = exp(params[i,1] + params[i,2]*xs)
   1/(1+num)
}
l <- adply(lines, 1, function(d) quantile(d, c(0.025, .5, .975)))
l$X1 <- xs
lines <- as.data.frame(lines)

l %>% ggplot() +
   geom_ribbon(aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'gray70') +
   geom_line(aes(x = X1, y = `50%`))


# params2 <- data.frame(matrix(NA, nrow = n_trials, ncol = 3))
# for (i in 1:n_trials) {
#   out <- try(optim(c(1,1,1), ll_worms3, probs=probs[i,], conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),silent = T)
#   if(inherits(out, "try-error")) {
#     print(i)
#     out <- optim(c(-1,1,1), ll_worms3, probs=probs[i,], conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
#   } 
#   params2[i,] <- out$par
# }




params2 <- foreach (i=1:n_trials, .combine = rbind) %dopar% {
   out <- try(optim(c(-1,1,1), ll_worms3, probs=probs[i,], conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),silent = T)
   if(inherits(out, "try-error")) {
      out <- optim(c(-1,1,1), ll_worms3, probs=probs[i,], conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
   } 
   return(out$par)
} 
params2 <- as.data.frame(params2)

lines2 <- foreach (i = 1:n_trials, .combine=cbind) %dopar% {
   num = exp(params2[i,1] + params2[i,2]*xs)
   1/(1+num)^params2[i,3]
}

l2 <- adply(lines2, 1, function(d) quantile(d, c(0.025, .5, .975)))
lines2 <- as.data.frame(lines2)
l2$X1 <- xs



params3 <- foreach (i=1:n_trials, .combine = rbind) %dopar% {
   out <- try(optim(c(-1,1,1,1), ll_worms4, probs=probs[i,], conc =  log_conc, method = 'L-BFGS-B', control = list(fnscale=-1)),silent = T)
   if(inherits(out, "try-error")) {
      out <- optim(c(-1,1,1,1), ll_worms4, probs=probs[i,], conc =  log_conc, method = 'BFGS', control = list(fnscale=-1))
   } 
   return(out$par)
} 
params3 <- as.data.frame(params3)

lines3 <- foreach (i = 1:n_trials, .combine=cbind) %dopar% {
   num = exp(params3[i,1] + params3[i,2]*xs)
   1/(1+num)^params3[i,3]
}

l3 <- adply(lines3, 1, function(d) quantile(d, c(0.025, .5, .975)))
lines3 <- as.data.frame(lines3)
l3$X1 <- xs



# LL3_LC50 <- params2 %>% mutate(lc50 = -(log((1-.5^V3)/(.5^V3)) + V1)/V2)
# # LL3_LC50 <- params2 %>% mutate(lc50 = -(log((1-.5^X3)/(.5^X3)) + X1)/X2)
# 
# # params2 %>% ggplot() + 
# #   geom_density_2d_filled(aes(x = -V1/V2, y =V3), n = 100, bins = 25) + 
# #   theme(legend.position = "none")
# 
# ggplot()+
#   geom_density(data = LL3_LC50, aes(x=lc50), color = 'red') + 
#   geom_density(data = params, aes(x=-V1/V2), color = "black")

ggplot() + 
   geom_ribbon(data = l,aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'springgreen', alpha = 0.5) +
   geom_ribbon(data = l2,aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'lightskyblue', alpha = 0.5) +
   geom_ribbon(data = l3, aes(x=X1, ymin = `2.5%`, ymax = `97.5%`), fill = 'thistle', alpha = 0.5) +
   geom_line(data = l,aes(x = X1, y = `50%`), color = "springgreen4", alpha = 0.5) +
   geom_line(data = l2,aes(x = X1, y = `50%`), color = "lightskyblue4", alpha = 0.5)+
   geom_line(data = l3,aes(x = X1, y = `50%`), color = "thistle4", alpha = 0.5)

