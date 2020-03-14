#model: Z_i ~ N(0,exp(h_i/2)), i=1, ...,T
#h_1 ~ N(mu,sigma^2/(1-xi^2))
#h_i|h_i-1,mu,xi,sigma ~ N(mu+xi*(h_i-1-mu), sigma^2), i=2, ...,T
#parameters of sv model: mean mu, persistance xi, standard deviation sigma, latent variables h    
#parameter transfromations: xi_cont=f_z(xi) (fisher's z transform), sigma_cont=log(sigma)
#let h_std ~ N(0,I) (T dimensional)
#h=numeric(T)
#h[1]=h_std[1]*sigma/sqrt(1-xi^2)+mu
#for(i in 2:T){
#  h[i]=mu+h_std[i]*sigma+xi*(h[i-1]-mu)
#}
# then h follows AR(1) process as described above, to go from h to h_std multply h with cov(h)^(-1/2) (cov(h) is the covariance matrix defined by the AR(1) process, for the sampler it is not necessary to calculate it)

fz=function(x) 0.5*log((1-x)^(-1)*(1+x)) #fisher's z transform
fz_inv=function(x) (exp(2*x)-1)* (exp(2*x)+1)^(-1)#inverse fisher's z transform

#prior densities and its derivatives, priors chosen as in stochvol
log_prior_mu=function(x) dnorm(x,sd=100, log=TRUE) 
gradient_log_prior_mu=function(x) -x/100^2

a_xi=5
b_xi=1.5
log_prior_xi_cont=function(x) (a_xi-1)*log(fz_inv(x)+1)+(b_xi-1)*log(1-fz_inv(x))+log(1-fz_inv(x)^2) 
gradient_log_prior_xi_cont=function(x) (a_xi-1)*1/(fz_inv(x)+1)*(1-fz_inv(x)^2)+(b_xi-1)*1/(1-fz_inv(x))*(-1)*(1-fz_inv(x)^2)+1/(1-fz_inv(x)^2)*2*(fz_inv(x))*(-1)*(1-fz_inv(x)^2)

log_prior_sigma_cont=function(x) -x-exp(2*x)/2+log(2)+2*x -1/2*log(2) -log(gamma(1/2))
gradient_log_prior_sigma_cont=function(x) 1-exp(2*x)

#log posterior density
log_post_sv=function(mu,xi_cont,sigma_cont,h_std,Z){
  xi=fz_inv(xi_cont)
  sigma=exp(sigma_cont)
  T=length(h_std)
  h=numeric(T)
  h[1]=h_std[1]*sigma/sqrt(1-xi^2)+mu
  for(i in 2:T){
    h[i]=mu+h_std[i]*sigma+xi*(h[i-1]-mu)
  }
  sum(dnorm(Z,mean=0,sd=exp(h/2), log=TRUE))+sum(dnorm(h_std,log=TRUE))+log_prior_sigma_cont(sigma_cont)+log_prior_xi_cont(xi_cont)+log_prior_mu(mu)
}

#include h_std, mu, xi_cont, sigma_cont in parameter vector par
log_post_dens=function(par,Z)
{
  mu=par[1]
  xi_cont=par[2]
  sigma_cont=par[3]
  h_std=par[4:length(par)]
  log_post_sv(mu,xi_cont,sigma_cont,h_std,Z)
}

#derivatives of the log posterior wrt all parameters
gradient_par=function(par,Z){
  mu=par[1]
  xi_cont=par[2]
  sigma_cont=par[3]
  h_std=par[4:length(par)]
  T=length(h_std)
  xi=fz_inv(xi_cont)
  sigma=exp(sigma_cont)
  gradient=numeric(T+3)
  h=numeric(T)
  h[1]=h_std[1]*sigma/sqrt(1-xi^2)+mu
  for(i in 2:T){
    h[i]=mu+h_std[i]*sigma+xi*(h[i-1]-mu)
  }
  gradient[1]=sum(((-1/exp(h/2)+Z^2/(exp(h/2))^3)*exp(h/2)*1/2))+gradient_log_prior_mu(mu)  
  x=numeric(T) #contains derivatives of h wrt. xi
  x[1]=h_std[1]*sigma*xi*(1-xi^2)^(-3/2)
  for(i in 2:T)
  {
    x[i]=h[i-1]-mu+xi*x[i-1]
  }
  gradient[2]=sum(((-1/exp(h/2)+Z^2/(exp(h/2))^3)*exp(h/2)*1/2)*x*(1-fz_inv(xi_cont)^2))+gradient_log_prior_xi_cont(xi_cont) 
  
  x=numeric(T) #contains derivatives of h wrt. sigma
  x[1]=h_std[1]*1/sqrt(1-xi^2)
  for(i in 2:T)
  {
    x[i]=h_std[i]+xi*x[i-1]
  }
  gradient[3]=sum(((-1/exp(h/2)+Z^2/(exp(h/2))^3)*exp(h/2)*1/2)*x)*exp(sigma_cont)+gradient_log_prior_sigma_cont(sigma_cont)  
  
  m=matrix(0,nrow=T,ncol=T) #contains derivaties of h wrt h_std
  m[,1]=xi^(0:(T-1))*sigma/sqrt(1-xi^2)
  for(j in 2:T)
  {
    m[j:T,j]=xi^(0:(T-j))*sigma
  }
  gradient[4:(T+3)]=(0.5*(Z^2*exp(-h)-1))%*%m - h_std #derivatives wrt. h_std
  return(gradient)
}  

#check if gradient is correct
#library(numDeriv)
#par=rnorm(23)
#Z=rnorm(20)
#gradient_par(par=par,Z=Z)-grad(log_post_dens,x=par,Z=Z)


#one iteration of hmc with stepsize epsilon, number of steps L, mass vector M (elements of diagonal mass matrix)
hmc_iteration = function(par,Z,epsilon,L,M)
{
  nan_count=0
  M_inv = 1/M
  d_par = length(par)
  phi = rnorm (d_par, 0, sqrt(M))
  par_old=par
  log_post_dens_old = log_post_dens(par,Z) - 0.5*sum(M_inv*phi^2)
  phi = phi + 0.5*epsilon*gradient_par(par, Z)
  for (l in 1:L){
    if((any(is.nan(c(par,phi)))))
    {
      p_jump=0
      par_new=par_old
      return (list (par=par_new, p_jump=p_jump,nan_count=1))
    }
    par = par + epsilon*M_inv*phi
    phi = phi + (if (l==L) 0.5 else 1)*epsilon*gradient_par(par,Z)
  }
  phi = -phi
  log_post_dens_new = log_post_dens(par,Z) - 0.5*sum(M_inv*phi^2)
  r = exp (log_post_dens_new - log_post_dens_old)
  if (is.nan(r)){
    r = 0
    nan_count=1
  }
  p_jump = min(r,1)
  par_new = if (runif(1) < p_jump) par else par_old
  return (list (par=par_new, p_jump=p_jump,nan_count=nan_count,log_post_dens_new=log_post_dens_new))
}

#hmc where stepsize epsilon is chosen uniformly between epsilon_1, epsilon_2, number of steps L is chosen uniformly between L_1,L_2
hmc_sv= function (Z,starting_values,M, iter=1000, epsilon_1=0,epsilon_2=0.1, L_1=0,L_2=100) {
  nan_count=0
  d_par=length(starting_values)
  sims = array (NA, c(iter, d_par))
  p_jump = numeric(iter)
  log_post = numeric(iter)
    par_temp = starting_values
    for (t in 1:iter){
      epsilon = epsilon_1+runif(1)*(epsilon_2-epsilon_1)
      L = ceiling(L_1+runif(1)*(L_2-L_1))
      temp = hmc_iteration(par=par_temp,Z, epsilon, L, M)
      p_jump[t] = temp$p_jump
      sims[t,] = temp$par
      nan_count=nan_count+temp$nan_count
      par_temp = temp$par
      if(t %% 100==0) {
        cat(paste0("iteration: ", t, "\n"))
      }
    }
    #go back to standard parametrization
    T=length(Z)
    mu=sims[,1]
    xi=fz_inv(sims[,2])
    sigma=exp(sims[,3])
    h_std=sims[,4:d_par]
    h=matrix(nrow=iter,ncol=T)
    h[,1]=h_std[,1]*sigma/sqrt(1-xi^2)+mu
    for(i in 2:T){
      h[,i]=mu+h_std[,i]*sigma+xi*(h[,i-1]-mu)
    }
  return (list (sims_o=cbind(mu,xi,sigma,h),sims=sims, p_jump=p_jump,nan_count=nan_count,log_post=log_post))
}

#function to simulate from sv model
stochvol_sim=function(T,mu,xi,sigma)
{
  h=numeric(T)
  eps=rnorm(T)
  delta=rnorm(T)
  h[1]=rnorm(1,mean=mu,sd=sigma/(1-xi^2))
  for(t in 2:T)
  {
    h[t]=mu+xi*(h[t-1]-mu)+delta[t-1]*sigma
  }
  Z=eps*exp(h/2)
  return(list(Z=Z,T=T,h=h))
}


#example:

#simulate data
sim=stochvol_sim(T=300,mu=-2,xi=0.9,sigma=0.3)

#run hmc
run=hmc_sv(Z=sim$Z,starting_values=rep(0,sim$T+3),M=rep(1,sim$T+3))
 
#avg acceptance prob
mean(run$p_jump)

#some traceplots 
plot(run$sims_o[,1], main="traceplot: mu", type="l")
plot(run$sims_o[,2], main="traceplot: xi", type="l")
plot(run$sims_o[,3], main="traceplot: sigma", type="l")
plot(run$sims_o[,20], main="traceplot: h", type="l")
plot(run$sims_o[,100], main="traceplot: h", type="l")
