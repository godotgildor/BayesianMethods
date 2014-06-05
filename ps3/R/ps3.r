########################################################################################################
# Question 2
sat_prep_data <- read.table("../data/sat_prep.txt",header=T)
attach(sat_prep_data)

############
# Sample directly
# This is the marginal posterior distribution of p(mu, tau^2 | yi, sigma^2)  It will actually
# just allow us to evaluate this at given points so that we can do grid sampling.
marginal_posterior_dist <- function(treatment_effects, treatment_effects_errors, mu, tau)
{
	prob <- 0

	num_treatments <- length(treatment_effects)

	for(j in 1:num_treatments)
	{
		curr_sd = sqrt( treatment_effects_errors[j]^2 + tau^2)
		prob <- prob +  dnorm(treatment_effects[j], mean=mu, sd=curr_sd, log = TRUE)
	}

	prob <- prob - log(tau)

	prob
}

# This will sample theta_i given mu, tau, y, and sigma
sample_theta <- function(treatment_effect, treatment_effect_error, mu, tau)
{
	Vi <- 1/(1/treatment_effect_error^2 + 1/tau^2)
	theta_hat <- (treatment_effect / treatment_effect_error^2 + mu/tau^2)* Vi

	out <- rnorm(1,mean=theta_hat, sd=sqrt(Vi))
	out
}

#######
# Sample our marginal posterior over a 1000 x 1000 grid
num_samples <- 1000
marg_post_dist <- array(NA,dim=c(num_samples, num_samples))
tau_sqd_grid = ppoints(num_samples)*10
mu_grid = ppoints(num_samples)*20 - 2.5
for(i in 1:num_samples)
{
	print (i)
	for(j in 1:num_samples)
	{
		marg_post_dist[i, j] <- marginal_posterior_dist(Estimate, StdErr, mu_grid[i], sqrt(tau_sqd_grid[j]))
	}
}
marg_post_dist <- exp(marg_post_dist - max(marg_post_dist))
marg_post_dist <- marg_post_dist / sum(marg_post_dist)

#######
# Calculate the marginal dist of tau
tau_sqd_marginal_dist <- array(NA,dim=num_samples)
for(i in 1:num_samples)
{
	tau_sqd_marginal_dist[i] <- sum(marg_post_dist[,i])
}
tau_sqd_marginal_dist <- tau_sqd_marginal_dist / sum(tau_sqd_marginal_dist)
# And now calculate the conditional dist of mu
mu_conditional_dist <- array(NA, dim=c(num_samples, num_samples))
for(i in 1:num_samples)
{
	for(j in 1:num_samples)
	{
		mu_conditional_dist[i, j] <- marg_post_dist[i,j] / tau_sqd_marginal_dist[j]
	}
}

#######
# Question 3 and 4
# Now take 1000 samples of mu and tau and calculate intervals.
num_samples <-1000
mu_samples <- array(NA, dim=num_samples)
tau_sqd_samples <- array(NA, dim=num_samples)
theta_samples <- array(NA, dim=c(length(Estimate),num_samples))
grid_size <- length(tau_sqd_marginal_dist)
for(i in 1:num_samples)
{
	a <- sample(1:grid_size, size=1, prob=tau_sqd_marginal_dist)
	b <- sample(1:grid_size, size=1, prob=mu_conditional_dist[,a])
	tau_sqd_samples[i] <- tau_sqd_grid[a]
	mu_samples[i] <- mu_grid[b]

	for (j in 1:length(Estimate) )
	{
		theta_samples[j, i] <- sample_theta(Estimate[j], StdErr[j], mu_samples[i], sqrt(tau_sqd_samples[i]))
	}
}
mean(mu_samples)
median(mu_samples)
quantile(mu_samples, 0.025)
quantile(mu_samples, 0.975)
mean(tau_sqd_samples)
median(tau_sqd_samples)
quantile(tau_sqd_samples, 0.025)
quantile(tau_sqd_samples, 0.975)
for(i in 1:length(Estimate))
{
	print(mean(theta_samples[i,]))
}

#############
# Second way to calculate mu and tau
# Note: not using this code for any answers on the PS.
# Sample tau first then mu
log_tau_sqd_posterior <- function(treatment_effects, treatment_effects_errors, tau)
{
	num_treatments <- length(treatment_effects)

	V_mu <- 0
	mu_hat <- 0
	for(i in 1:num_treatments)
	{
		temp <- 1/(treatment_effects_errors[i]^2 + tau^2)
		mu_hat <- mu_hat + temp * treatment_effects[i]		
		V_mu <- V_mu + temp		
	}
	V_mu <- 1/V_mu	
	mu_hat <- mu_hat * V_mu	

	out <- 0
	# Is it out - 0.5 log(V_mu) (notes) or out + 0.5 log(V_mu) (book)
	out <- out + 0.5 * log(V_mu)
	for(i in 1:num_treatments)
	{
		temp <- treatment_effects_errors[i]^2 + tau^2
		out <- out + -0.5 * log(temp)
		out <- out + -0.5 * (treatment_effects[i] - mu_hat)^2/temp
	}

	# What about p(tau^2) ?
	out <- out - log(tau)

	out
}

sample_mu_posterior <- function(treatment_effects, treatment_effects_errors, tau)
{
	num_treatments <- length(treatment_effects_errors)

	V_mu <- 0
	mu_hat <- 0 
	for(i in 1:num_treatments)
	{
		temp <- 1/(treatment_effects_errors[i]^2 + tau^2)
		mu_hat <- mu_hat + temp * treatment_effects[i]		
		V_mu <- V_mu + temp		
	}
	V_mu <- 1/V_mu
	mu_hat <- mu_hat * V_mu
	

	out <- rnorm(1, mean=mu_hat, sd=sqrt(V_mu))
	out
}

num_samples <- 1000
# Get a distribution for tau
tau_sqd_grid <- ppoints(num_samples)*10
tau_sqd_posterior <- array(NA,dim=num_samples)
for(i in 1:num_samples)
{
	curr_tau <- sqrt(tau_sqd_grid[i])
	tau_sqd_posterior[i] <- log_tau_sqd_posterior(Estimate, StdErr, curr_tau)
}
tau_sqd_posterior <- exp(tau_sqd_posterior - max(tau_sqd_posterior))
tau_sqd_posterior <- tau_sqd_posterior / sum(tau_sqd_posterior)
plot(tau_sqd_grid, tau_sqd_posterior, type="l")

tau_sqd_samples <- array(NA, dim=num_samples)
mu_samples <- array(NA, dim=num_samples)
theta_samples <- array(NA, dim=c(length(Estimate),num_samples))
for (i in 1:num_samples)
{
	# Sampling tau
	tau_sqd_samples[i] <- sample(tau_sqd_grid, size=1, prob=tau_sqd_posterior)
	# Sample mu
	mu_samples[i] <- sample_mu_posterior(Estimate, StdErr, sqrt(tau_sqd_samples[i]))

	for (j in 1:length(Estimate) )
	{
		theta_samples[j, i] <- sample_theta(Estimate[j], StdErr[j], mu_samples[i], sqrt(tau_sqd_samples[i]))
	}
}

mean(mu_samples)
median(mu_samples)
quantile(mu_samples, 0.025)
quantile(mu_samples, 0.975)

mean(tau_sqd_samples)
median(tau_sqd_samples)
quantile(tau_sqd_samples, 0.025)
quantile(tau_sqd_samples, 0.975)

for(i in 1:length(Estimate))
{
	print(mean(theta_samples[i,]))
}

# End of second way.

##############
# Questions 5 and 6
# Gibbs sampler portion
# Fix the value of tau^2.  This is the value of the median
# of tau^2 that I calculated above for one of my runs.
fixed_tau_sqd <- 2.345

num_iters <- 100500
thetas_from_gibbs <- array(NA, dim=c(length(Estimate),num_iters))
mus_from_gibbs <- array(NA, dim=num_iters)

# Initialize the first values
for(i in 1:length(Estimate))
{
	thetas_from_gibbs[i, 1] <- i*100
}
mus_from_gibbs[1] <- 300

for (i in 2:num_iters)
{
	# Sample the thetas
	for(j in 1:length(Estimate))
	{
		Vj <- 1 / (1/StdErr[j]^2 + 1/fixed_tau_sqd)
		theta_hat <- Estimate[j] / StdErr[j]^2 + mus_from_gibbs[i-1] / fixed_tau_sqd
		theta_hat <- theta_hat * Vj
		thetas_from_gibbs[j, i] <- rnorm(1,mean=theta_hat, sd=sqrt(Vj))
	}

	# Now sample the mu's
	avg_theta <- mean(thetas_from_gibbs[,i])
	mus_from_gibbs[i] <- rnorm(1,mean=avg_theta,sd=sqrt(fixed_tau_sqd/length(Estimate)))
}
par(mfrow=c(2,1))
y_min <- min(thetas_from_gibbs)
y_max <- max(thetas_from_gibbs)
plot(1:num_iters,thetas_from_gibbs[1,],type="l",main="thetas",ylim=c(y_min,y_max))
for(i in 2:length(Estimate))
{
	lines(1:num_iters,thetas_from_gibbs[i,],col=i)
}
y_min <- min(mus_from_gibbs)
y_max <- max(mus_from_gibbs)
plot(1:num_iters,mus_from_gibbs,type="l",main="mu",ylim=c(y_min,y_max))

# trim off the first 500 samples for burn-in
mus_from_gibbs <- mus_from_gibbs[501:num_iters]
thetas_from_gibbs <- thetas_from_gibbs[,501:num_iters]
acf(mus_from_gibbs)
acf(thetas_from_gibbs[1,])

# Sub-sample to get 1000 samples
temp <- length(mus_from_gibbs)/1000*c(1:1000)
mus_from_gibbs_trimmed <- mus_from_gibbs[temp]
thetas_from_gibbs_trimmed <- thetas_from_gibbs[,temp]
acf(mus_from_gibbs_trimmed)
acf(thetas_from_gibbs_trimmed[1,])

for(i in 1:length(Estimate))
{
	print(mean(thetas_from_gibbs_trimmed[i,]))
	print(quantile(thetas_from_gibbs_trimmed[i,],0.025))
	print(quantile(thetas_from_gibbs_trimmed[i,], 0.975))
}

# Question 7
# Let's see how many times program A had the better mean.
num_times_prog_a_best <- 0
num_samples <- length(thetas_from_gibbs_trimmed[1,])
num_schools <- length(thetas_from_gibbs_trimmed[,1])
sampled_yis = array(NA, dim=num_schools)
for(i in 1:num_samples)
{
	for(j in 1:num_schools)
	{
		sampled_yis[j] <- rnorm(1, mean=thetas_from_gibbs_trimmed[j,i],sd=StdErr[j])
	}
	if(sampled_yis[1] >= max(sampled_yis))
	{
		num_times_prog_a_best <- num_times_prog_a_best + 1
	}
}
prob_prog_a_best <- num_times_prog_a_best / num_samples
prob_prog_a_best

########################################################################################
# Questions 11 and 12
# Bicycle problems

bicycle_data <- read.table("../data/bicycles2.txt",header=T)
attach(bicycle_data)
# We are only interested in the first row of data.

num_iters <- 200500
thetas_from_gibbs <- array(NA, dim=c(length(bikes),num_iters))
alphas_from_gibbs <- array(NA, dim=length(bikes))
betas_from_gibbs <- array(NA, dim=length(bikes))

# Set some init values.
alphas_from_gibbs[1] <- 30
betas_from_gibbs[1] <- 23
for(i in 1:length(bikes))
{
	thetas_from_gibbs[i,1] <- i*33
}

for(i in 2:num_iters)
{
	print(i)
	for(j in 1:length(bikes))
	{
		thetas_from_gibbs[j, i] <- rgamma(1, alphas_from_gibbs[i-1] + bikes[j] + vehicles[j], rate=betas_from_gibbs[i-1]+1)
	}

	# Now do Metropolis Hastings for alpha and betas
	prev_alpha <- alphas_from_gibbs[i-1]
	prev_beta <- betas_from_gibbs[i-1]
	
	# Alpha
	#candidate_alpha <- rgamma(1,2,rate=(2/prev_alpha))
	candidate_alpha <- rnorm(1,mean=prev_alpha,sd=5)
	if(candidate_alpha <0)
	{
		candidate_alpha <- 0
	}
	a1 <- 0
	for(j in 1:length(bikes))
	{
		total_vehicles <- thetas_from_gibbs[j,i]
		a1 <- a1 + (candidate_alpha - prev_alpha)*log(prev_beta)-lgamma(candidate_alpha)+lgamma(prev_alpha)+(candidate_alpha - prev_alpha)*log(total_vehicles)
	}
	a1 <- exp(a1)
	#a2 <- dgamma(prev_alpha,2,rate=(2/candidate_alpha))/dgamma(candidate_alpha,2,rate=(2/prev_alpha))
	a2 <- 1
	a <- a1*a2
	temp <- runif(1,min=0, max=1)
	if(temp < a)
	{
		alphas_from_gibbs[i] <- candidate_alpha
	}
	else
	{
		alphas_from_gibbs[i] <- prev_alpha
	}
	
	# Beta
	curr_alpha <- alphas_from_gibbs[i]
	#candidate_beta <- rgamma(1,2,rate=(2/prev_beta))
	candidate_beta <- rnorm(1,mean=prev_beta,sd=5)
	if(candidate_beta <0)
	{
		candidate_beta <- 0
	}
	a1 <-0
	for(j in 1:length(bikes))
	{
		total_vehicles <- thetas_from_gibbs[j,i]
		a1 <- a1 + curr_alpha*(log(candidate_beta) - log(prev_beta)) + total_vehicles*(prev_beta-candidate_beta)
	}
	a1 <- exp(a1)
	#a2 <- dgamma(prev_beta,2,rate=(2/candidate_beta))/dgamma(candidate_beta,2,rate=(2/prev_beta))
	a2 <- 1
	a <- a1*a2
	temp <- runif(1,min=0, max=1)
	if(temp < a)
	{
		betas_from_gibbs[i] <- candidate_beta
	}
	else
	{
		betas_from_gibbs[i] <- prev_beta
	}
}

par(mfrow=c(3,1))
plot(1:length(alphas_from_gibbs),alphas_from_gibbs,type="l")
plot(1:length(betas_from_gibbs),betas_from_gibbs,type="l")
plot(1:length(thetas_from_gibbs[1,]), thetas_from_gibbs[1,],type="l")
for(i in 2:length(bikes))
{
	lines(1:length(thetas_from_gibbs[1,]), thetas_from_gibbs[1,],type="l",col=i)
}

# Trim down to get rid of burn in
alphas_from_gibbs <- alphas_from_gibbs[501:num_iters]
betas_from_gibbs <- betas_from_gibbs[501:num_iters]
thetas_from_gibbs <- thetas_from_gibbs[,501:num_iters]

acf(alphas_from_gibbs)
acf(betas_from_gibbs)
acf(thetas_from_gibbs[1,])

# Trim down to only 1000 samples
curr_size <- length(alphas_from_gibbs)
temp <- curr_size/1000*c(1:1000)
alphas_from_gibbs_trimmed <- alphas_from_gibbs[temp]
betas_from_gibbs_trimmed <- betas_from_gibbs[temp]
thetas_from_gibbs_trimmed <- thetas_from_gibbs[,temp]
acf(alphas_from_gibbs_trimmed)
acf(betas_from_gibbs_trimmed)
acf(thetas_from_gibbs_trimmed[1,])

mean(alphas_from_gibbs_trimmed)
quantile(alphas_from_gibbs_trimmed,0.025)
quantile(alphas_from_gibbs_trimmed,0.975)

mean(betas_from_gibbs_trimmed)
quantile(betas_from_gibbs_trimmed,0.025)
quantile(betas_from_gibbs_trimmed,0.975)

# Predict traffic on new residential street
num_iters <- 1000
total_vehicles <- array(NA, dim=num_iters)
for(i in 1:num_iters)
{
	curr_theta <- rgamma(1,alphas_from_gibbs_trimmed[i],rate=betas_from_gibbs_trimmed[i])
	total_vehicles[i] <- rpois(1,lambda=curr_theta)
}
mean(total_vehicles)
quantile(total_vehicles,0.025)
quantile(total_vehicles,0.975)
