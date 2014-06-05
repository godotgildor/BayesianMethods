#Problem 2
#Some initial setup
plane_data <- read.table("../data/planes.txt",header=T)
attach(plane_data)
num_replications <- 1000
num_years <- length(fatal)

fatal_obs <- array(NA, dim=c(1,num_years))
deaths_obs <- array(NA, dim=c(1,num_years))
for (curr_year in 1:num_years)
{
	fatal_obs[1, curr_year] <- fatal[curr_year]
	deaths_obs[1, curr_year] <- deaths[curr_year]
}

calc_t_stat_for_independence <- function(replicated_data)
{
	num_reps <- dim(replicated_data)[1]
	num_samples <- dim(replicated_data)[2]

	t_stat_ind <- array(NA,dim=num_reps)
	for (curr_rep in 1:num_reps)
	{
		curr_data <- replicated_data[curr_rep,]
		my_acf <- acf(curr_data,plot=FALSE)$acf[2]
		t_stat_ind[curr_rep] <- abs(my_acf)

# Previously I was correlating the even and odd samples.  Acf is probably better.
#		group1 <- (1:(num_samples/2))*2
#		group2 <- group1 - 1	
#		t_stat_ind[curr_rep] <- abs(cor(curr_data[group1], curr_data[group2]))
	}
	
	t_stat_ind
}

calc_t_stat_for_trend_with_time <- function(replicated_data, use_corr=FALSE)
{
	num_reps <- dim(replicated_data)[1]
	num_samples <- dim(replicated_data)[2]

	t_stat_trend_with_time <- array(0,dim=num_reps)
	# Let's get the slope of the best fit line with regards to years.
	for (curr_rep in 1:num_reps)
	{
		curr_data <- replicated_data[curr_rep,]
		sample_nums <- 1:num_samples
		
		if(use_corr)
		{
			t_stat_trend_with_time[curr_rep] <- abs(cor(curr_data, sample_nums))
		}
		else
		{
			# Calculate the linear fit of data.  Then take
			# slope and divide by std. err.  I found the equation for the
			# std err at: http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8467.html
			my_linear_model <- lm(curr_data~sample_nums)
			summary_my_lm <- summary(my_linear_model)
			X <- model.matrix(my_linear_model)
			var.betas <- solve(crossprod(X)) * summary_my_lm$sigma^2
			std_err <- sqrt(diag(var.betas))

			t_stat_trend_with_time[curr_rep] <- abs(my_linear_model$coeff[2] / std_err[2])
		}
	}

	t_stat_trend_with_time
}


#Part a
sum_fatal_accidents <- sum(fatal)

# Now loop over the number of replications and draw a new set of 
# fatalities.
fatal_accidents_replicated <- array(NA,dim=c(num_replications, num_years))
for (curr_replication in 1:num_replications)
{
	for (curr_year in 1:num_years)
	{
		sample_of_theta = rgamma(1, shape=(sum_fatal_accidents + 0.0001), rate = (num_years+.0001))
		fatal_accidents_replicated[curr_replication, curr_year] <- rpois(1,sample_of_theta)
	}
}

# Calculate the T stat for independence.
tstat_ind_rep <- calc_t_stat_for_independence(fatal_accidents_replicated)
tstat_ind_obs <- calc_t_stat_for_independence(fatal_obs)
p_val <- sum(tstat_ind_rep < tstat_ind_obs[1])/num_replications
print(p_val)

# Calculate the T stat for trend with time.
tstat_twt_slope_rep <- calc_t_stat_for_trend_with_time(fatal_accidents_replicated)
tstat_twt_slope_obs <- calc_t_stat_for_trend_with_time(fatal_obs)
mean(tstat_twt_slope_rep)
p_val <- sum(tstat_twt_slope_rep < tstat_twt_slope_obs[1])/num_replications
print(p_val)

tstat_twt_corr_rep <- calc_t_stat_for_trend_with_time(fatal_accidents_replicated, TRUE)
tstat_twt_corr_obs <- calc_t_stat_for_trend_with_time(fatal_obs, TRUE)
mean(tstat_twt_corr_rep)
p_val <- sum(tstat_twt_corr_rep < tstat_twt_corr_obs[1])/num_replications
print(p_val)

par(mfrow=c(3,1))
hist(tstat_ind_rep,main='Independence T Stat (Cor) for Part A')
abline(v=tstat_ind_obs,col=2)
hist(tstat_twt_slope_rep,main='Trend with Time T Stat (linear fit slope) for Part A')
abline(v=tstat_twt_slope_obs,col=2)
hist(tstat_twt_corr_rep,main='Trend with Time T Stat (correlation) for Part A')
abline(v=tstat_twt_corr_obs,col=2)

#Part b
# Get the passenger_miles (in 100 millions)
passenger_miles = deaths/rate
sum_passenger_miles = sum(passenger_miles)

# Now loop over the number of replications and draw a new set
# of fatalities.
fatal_accidents_replicated <- array(NA,dim=c(num_replications, num_years))
for (curr_replication in 1:num_replications)
{
	for (curr_year in 1:num_years)
	{
		sample_of_theta <- rgamma(1, shape=(sum_fatal_accidents + 0.0001), rate = (sum_passenger_miles+.0001))
		fatal_accidents_replicated[curr_replication, curr_year] <- rpois(1,passenger_miles[curr_year]*sample_of_theta)
	}
}

# Calculate the T stat for independence.
tstat_ind_rep <- calc_t_stat_for_independence(fatal_accidents_replicated)
tstat_ind_obs <- calc_t_stat_for_independence(fatal_obs)
p_val <- sum(tstat_ind_rep < tstat_ind_obs[1])/num_replications
print(p_val)

# Calculate the T stat for trend with time.
tstat_twt_slope_rep <- calc_t_stat_for_trend_with_time(fatal_accidents_replicated)
tstat_twt_slope_obs <- calc_t_stat_for_trend_with_time(fatal_obs)
mean(tstat_twt_slope_rep)
p_val <- sum(tstat_twt_slope_rep < tstat_twt_slope_obs[1])/num_replications
print(p_val)

tstat_twt_corr_rep <- calc_t_stat_for_trend_with_time(fatal_accidents_replicated, TRUE)
tstat_twt_corr_obs <- calc_t_stat_for_trend_with_time(fatal_obs, TRUE)
mean(tstat_twt_corr_rep)
p_val <- sum(tstat_twt_corr_rep < tstat_twt_corr_obs[1])/num_replications
print(p_val)

par(mfrow=c(3,1))
hist(tstat_ind_rep,main='Independence T Stat (Cor) for Part B')
abline(v=tstat_ind_obs,col=2)
hist(tstat_twt_slope_rep,main='Trend with Time T Stat (linear fit slope) for Part B')
abline(v=tstat_twt_slope_obs,col=2)
hist(tstat_twt_corr_rep,main='Trend with Time T Stat (correlation) for Part B')
abline(v=tstat_twt_corr_obs,col=2)

#Part c
sum_deaths <- sum(deaths)

pass_deaths_replicated <- array(NA,dim=c(num_replications, num_years))
for (curr_replication in 1:num_replications)
{
	for(curr_year in 1:num_years)
	{
		sample_of_theta <- rgamma(1, shape=(sum_deaths + 0.0001), rate = (num_years+.0001))
		pass_deaths_replicated[curr_replication, curr_year] <- rpois(1,sample_of_theta)
	}
}

# Calculate the T stat for independence.
tstat_ind_rep <- calc_t_stat_for_independence(pass_deaths_replicated)
tstat_ind_obs <- calc_t_stat_for_independence(deaths_obs)
p_val <- sum(tstat_ind_rep < tstat_ind_obs[1])/num_replications
print(p_val)

# Calculate the T stat for trend with time.
tstat_twt_slope_rep <- calc_t_stat_for_trend_with_time(pass_deaths_replicated)
tstat_twt_slope_obs <- calc_t_stat_for_trend_with_time(deaths_obs)
mean(tstat_twt_slope_rep)
p_val <- sum(tstat_twt_slope_rep < tstat_twt_slope_obs[1])/num_replications
print(p_val)

tstat_twt_corr_rep <- calc_t_stat_for_trend_with_time(pass_deaths_replicated, TRUE)
tstat_twt_corr_obs <- calc_t_stat_for_trend_with_time(deaths_obs, TRUE)
mean(tstat_twt_corr_rep)
p_val <- sum(tstat_twt_corr_rep < tstat_twt_corr_obs[1])/num_replications
print(p_val)

par(mfrow=c(3,1))
hist(tstat_ind_rep,main='Independence T Stat (Cor) for Part C')
abline(v=tstat_ind_obs,col=2)
hist(tstat_twt_slope_rep,main='Trend with Time T Stat (linear fit slope) for Part C')
abline(v=tstat_twt_slope_obs,col=2)
hist(tstat_twt_corr_rep,main='Trend with Time T Stat (correlation) for Part C')
abline(v=tstat_twt_corr_obs,col=2)

#Part d
millions_passenger_miles <- deaths/rate
sum_millions_passenger_miles <- sum(millions_passenger_miles)

pass_deaths_replicated <- array(NA,dim=c(num_replications, num_years))
for (curr_replication in 1:num_replications)
{
	for(curr_year in 1:num_years)
	{
		sample_of_theta <- rgamma(1, shape=(sum_deaths + 0.0001), rate = (sum_millions_passenger_miles+.0001))
		pass_deaths_replicated[curr_replication, curr_year] <- rpois(1,millions_passenger_miles[curr_year]*sample_of_theta)
	}
}

# Calculate the T stat for independence.
tstat_ind_rep <- calc_t_stat_for_independence(pass_deaths_replicated)
tstat_ind_obs <- calc_t_stat_for_independence(deaths_obs)
p_val <- sum(tstat_ind_rep < tstat_ind_obs[1])/num_replications
print(p_val)

# Calculate the T stat for trend with time.
tstat_twt_slope_rep <- calc_t_stat_for_trend_with_time(pass_deaths_replicated)
tstat_twt_slope_obs <- calc_t_stat_for_trend_with_time(deaths_obs)
mean(tstat_twt_slope_rep)
p_val <- sum(tstat_twt_slope_rep < tstat_twt_slope_obs[1])/num_replications
print(p_val)

tstat_twt_corr_rep <- calc_t_stat_for_trend_with_time(pass_deaths_replicated, TRUE)
tstat_twt_corr_obs <- calc_t_stat_for_trend_with_time(deaths_obs, TRUE)
mean(tstat_twt_corr_rep)
p_val <- sum(tstat_twt_corr_rep < tstat_twt_corr_obs[1])/num_replications
print(p_val)

par(mfrow=c(3,1))
hist(tstat_ind_rep,main='Independence T Stat (Cor) for Part D')
abline(v=tstat_ind_obs,col=2)
hist(tstat_twt_slope_rep,main='Trend with Time T Stat (linear fit slope) for Part D')
abline(v=tstat_twt_slope_obs,col=2)
hist(tstat_twt_corr_rep,main='Trend with Time T Stat (correlation) for Part D')
abline(v=tstat_twt_corr_obs,col=2)

# Problem 3
# Part A
# Let's try 10,000 replications this time.
num_replications <- 10000
num_samples_per_iter <- 100

replicated_samples <- array(NA,dim=c(num_replications,num_samples_per_iter))
for (curr_replication in 1:num_replications)
{
	curr_theta <- rnorm(1,mean=5.1, sd=sqrt(1/100))
	for (curr_sample in 1:num_samples_per_iter)
	{
		replicated_samples[curr_replication, curr_sample] <- rnorm(1,mean=curr_theta,sd=1)
	}
}

tstat_max <- array(NA,dim=num_replications)
for (curr_replication in 1:num_replications)
{
	tstat_max[curr_replication] <- max(abs(replicated_samples[curr_replication,]))
}
hist(tstat_max,main='Posterior Predictive Distribution of T(Y_rep)')
abline(v=8.1,col=2)
p_val <- sum(tstat_max > 8.1)/num_replications

# Part B
replicated_samples <- array(NA,dim=c(num_replications,num_samples_per_iter))
for (curr_replication in 1:num_replications)
{
	curr_theta <- runif(1,min=-10000,max=10000)
	for (curr_sample in 1:num_samples_per_iter)
	{
		replicated_samples[curr_replication, curr_sample] <- rnorm(1,mean=curr_theta,sd=1)
	}
}

tstat_max <- array(NA,dim=num_replications)
for (curr_replication in 1:num_replications)
{
	tstat_max[curr_replication] <- max(abs(replicated_samples[curr_replication,]))
}
hist(tstat_max,main='Prior Predictive Distribution of T(Y_rep)')
abline(v=8.1,col=2)
p_val <- sum(tstat_max > 8.1)/num_replications

# Problem 4
rat_tumor_data <- read.table("../data/rat_tumors.txt",header=T)
attach(rat_tumor_data)

transformed_prior_dist <- function(log_alpha_over_beta, log_alpha_plus_beta, num_rats, treatment_effect)
{
	alpha_plus_beta <- exp(log_alpha_plus_beta)
	alpha_over_beta <- exp(log_alpha_over_beta)
	beta <- alpha_plus_beta/(alpha_over_beta+1)
	alpha <- alpha_plus_beta - beta

	num_samples <- length(num_rats)

	out <- -5/2*log(alpha+beta)
	for(i in 1:num_samples)
	{
		out <- out + lgamma(alpha+beta) + lgamma(alpha+treatment_effect[i]) + lgamma(beta + num_rats[i] - treatment_effect[i])
		out <- out - lgamma(alpha) - lgamma(beta) - lgamma(alpha + beta + num_rats[i])
	}
	
	one_over_jac <- alpha*beta	
	out <- out + log(abs(one_over_jac))
	out
}

num_samples <- 1000

log_alpha_over_beta <- ppoints(num_samples) - 2.3
log_alpha_plus_beta <- ppoints(num_samples)*4 + 1

marg_post_dist <- array(NA, dim=c(num_samples, num_samples))
for( i in 1:num_samples)
{
	for( j in 1:num_samples)
	{
		marg_post_dist[i, j] <- transformed_prior_dist(log_alpha_over_beta[i], log_alpha_plus_beta[j], num_rats, treatment_effect)
	}
}

marg_post_dist <- marg_post_dist - max(marg_post_dist)
marg_post_dist <- exp(marg_post_dist)
marg_post_dist <- marg_post_dist / sum(marg_post_dist)
contour(x=log_alpha_over_beta, y=log_alpha_plus_beta, marg_post_dist)

log_alpha_over_beta_margin_dist <- array(NA,dim=num_samples)
log_alpha_plus_beta_cond_dist <- array(NA, dim=c(num_samples, num_samples))
for(i in 1:num_samples)
{
	log_alpha_over_beta_margin_dist[i] <- sum(marg_post_dist[i,])
}
for( i in 1:num_samples)
{
	for( j in 1:num_samples)
	{
		log_alpha_plus_beta_cond_dist[i,j] <- marg_post_dist[i,j] / log_alpha_over_beta_margin_dist[i]
	}
}

num_replicates <- 1000
treatment_effects_replicated <- array(NA, dim=c(num_replicates,length(num_rats)))
for( i in 1:num_replicates)
{
	log_alpha_over_beta_indx <- sample(1:num_samples, size=1, prob=log_alpha_over_beta_margin_dist)
	log_alpha_plus_beta_indx <- sample(1:num_samples, size=1, prob=log_alpha_plus_beta_cond_dist[log_alpha_over_beta_indx,])

	alpha_plus_beta <- exp(log_alpha_plus_beta[log_alpha_plus_beta_indx])
	alpha_over_beta <- exp(log_alpha_over_beta[log_alpha_over_beta_indx])
	beta <- alpha_plus_beta/(alpha_over_beta+1)
	alpha <- alpha_plus_beta - beta

	for( j in 1:length(num_rats))
	{
		theta <- rbeta(1, alpha + treatment_effect[j], beta + num_rats[j] - treatment_effect[j])
		treatment_effects_replicated[i, j] <- rbinom(1, num_rats[j], theta)
	}
}

tstat_num_zeros_obs <- sum(treatment_effect == 0)
tstat_max_obs <- max(treatment_effect)
tstat_mean_obs <- mean(treatment_effect)
tstat_sd_obs <- sqrt(var(treatment_effect))

tstat_num_zeros_rep <- array(NA, dim=num_replicates)
tstat_max_rep <- array(NA, dim=num_replicates)
tstat_mean_rep <- array(NA, dim=num_replicates)
tstat_sd_rep <- array(NA, dim=num_replicates)

for( i in 1:num_replicates)
{
	tstat_num_zeros_rep[i] <- sum(treatment_effects_replicated[i,] == 0)
	tstat_max_rep[i] <- max(treatment_effects_replicated[i,])
	tstat_mean_rep[i] <- mean(treatment_effects_replicated[i,])
	tstat_sd_rep[i] <- sqrt(var(treatment_effects_replicated[i,]))
}

p_val_num_zeros <- sum(tstat_num_zeros_rep > tstat_num_zeros_obs) / num_replicates
p_val_max <- sum(tstat_max_rep > tstat_max_obs) / num_replicates
p_val_mean <- sum(tstat_mean_rep > tstat_mean_obs) / num_replicates
p_val_sd <- sum(tstat_sd_rep > tstat_sd_obs) / num_replicates

par(mfrow=c(2,2))
hist(tstat_num_zeros_rep,main='Hist of Number of Rats with 0 Tumors')
abline(v=tstat_num_zeros_obs,col=2)
hist(tstat_max_rep,main='Hist of Max Number of Rats with Tumors')
abline(v=tstat_max_obs,col=2)
hist(tstat_mean_rep,main='Hist of Mean number of Rats with  Tumors')
abline(v=tstat_mean_obs,col=2)
hist(tstat_sd_rep,main='Hist of SD Tumors')
abline(v=tstat_sd_obs,col=2)

#################################
# Now look at the same stats but as ratio of tumors
tstat_max_obs <- max(treatment_effect/num_rats)
tstat_mean_obs <- mean(treatment_effect/num_rats)
tstat_sd_obs <- sqrt(var(treatment_effect/num_rats))

tstat_max_rep <- array(NA, dim=num_replicates)
tstat_mean_rep <- array(NA, dim=num_replicates)
tstat_sd_rep <- array(NA, dim=num_replicates)
for( i in 1:num_replicates)
{
	tstat_max_rep[i] <- max(treatment_effects_replicated[i,]/num_rats)
	tstat_mean_rep[i] <- mean(treatment_effects_replicated[i,]/num_rats)
	tstat_sd_rep[i] <- sqrt(var(treatment_effects_replicated[i,]/num_rats))
}
p_val_max <- sum(tstat_max_rep > tstat_max_obs) / num_replicates
p_val_mean <- sum(tstat_mean_rep > tstat_mean_obs) / num_replicates
p_val_sd <- sum(tstat_sd_rep > tstat_sd_obs) / num_replicates

par(mfrow=c(3,1))
hist(tstat_max_rep,main='Hist of Max Ratio of Rats with Tumors')
abline(v=tstat_max_obs,col=2)
hist(tstat_mean_rep,main='Hist of Mean Ratio of Rats with  Tumors')
abline(v=tstat_mean_obs,col=2)
hist(tstat_sd_rep,main='Hist of SD of Ratio of Rats with Tumors')
abline(v=tstat_sd_obs,col=2)


# Problem 5
football_data <- read.table("../data/football.txt",header=T)
attach(football_data)

num_replicates <- 1000
num_games <- length(favorite)

x_noise <- runif(num_games, min=-.1, max=.1)
y_noise <- runif(num_games, min=-.2, max=.2)

fav_score_minus_underdog_score_rep <- array(NA, dim=c(num_replicates, num_games))
for( i in 1:num_replicates)
{
	for(j in 1:num_games)
	{
		fav_score_minus_underdog_score_rep[i, j] <- rnorm(1, mean=spread[j], sd=14)
	}
}

stride <- floor(num_replicates/7)
par(mfrow=c(4,2))
plot(spread+x_noise,favorite-underdog+y_noise,main="Actual Data",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[1,]+y_noise,main="Replicated Data 1",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[stride+1,]+y_noise,main="Replicated Data 2",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[2*stride+1,]+y_noise,main="Replicated Data 3",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[3*stride+1,]+y_noise,main="Replicated Data 4",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[4*stride+1,]+y_noise,main="Replicated Data 5",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[5*stride+1,]+y_noise,main="Replicated Data 6",xlab="point spread", ylab="outcome")
plot(spread+x_noise, fav_score_minus_underdog_score_rep[6*stride+1,]+y_noise,main="Replicated Data 7",xlab="point spread", ylab="outcome")

stride <- floor(num_replicates/3)
par(mfrow=c(4,2))
plot(spread+x_noise,favorite-underdog-spread+y_noise,main="Actual Data",xlab="point spread", ylab="outcome-point spread")
hist(favorite-underdog-spread+y_noise,xlab="outcome - point spread",main="Actual Data", breaks=38)
plot(spread+x_noise, fav_score_minus_underdog_score_rep[1,]-spread+y_noise,main="Replicated Data 1",xlab="point spread", ylab="outcome-point spread")
hist(fav_score_minus_underdog_score_rep[1,]-spread+y_noise,xlab="outcome - point spread",main="Replicated Data 1", breaks=38)
plot(spread+x_noise, fav_score_minus_underdog_score_rep[stride+1,]-spread+y_noise,main="Replicated Data 2",xlab="point spread", ylab="outcome-spread")
hist(fav_score_minus_underdog_score_rep[stride+1,]-spread+y_noise,xlab="outcome - point spread",main="Replicated Data 2", breaks=38)
plot(spread+x_noise, fav_score_minus_underdog_score_rep[2*stride+1,]-spread+y_noise,main="Replicated Data 3",xlab="point spread", ylab="outcome-spread")
hist(fav_score_minus_underdog_score_rep[2*stride+1,]-spread+y_noise,xlab="outcome - point spread",main="Replicated Data 3", breaks=38)

calc_t_stat_var_simple <- function(replicated_data, point_spread, use_corr=FALSE)
{
	out <- 0

	point_spread_range <- (0:40)
	point_spread_range <- point_spread_range / 2

	temp <- array(NA, dim=length(point_spread_range) )
	for( i in 1:length(point_spread_range) )
	{
		games_of_interest <- which(point_spread == point_spread_range[i])

		if( length(games_of_interest) > 1)
		{
			temp[i] <- var(replicated_data[games_of_interest])
		}
		else
		{
			temp[i] <- NaN
		}
	}

	point_spread_range <- point_spread_range[is.finite(temp)]
	temp <- temp[is.finite(temp)]

	if(use_corr)
	{
		out <- cor(temp, point_spread_range)
	}
	else
	{		
		my_lm = lm(temp~point_spread_range)
		out <- my_lm$coeff[2]
	}

	out
}

calc_t_stat_var_equal_size_groups <- function(replicated_data, point_spread, use_corr=FALSE)
{
	out <- 0

	num_bins <- 40
	num_samples_per_bin = floor(length(replicated_data)/num_bins)

	sort_info <- sort(point_spread, index.return=TRUE)
	sort_indices <- sort_info$ix
	
	var_temp <- array(NA, dim=num_bins )
	ps_temp <- array(NA, dim=num_bins )
	
	for( i in 1:num_bins )
	{
		indices_of_interest <- sort_indices[((i-1)*num_samples_per_bin+1):(i*num_samples_per_bin)]
		var_temp[i] <- var(replicated_data[indices_of_interest])
		ps_temp[i] <- mean( point_spread[indices_of_interest] )
	}

	if(use_corr)
	{
		out <- cor(var_temp, ps_temp)
	}
	else
	{
		my_lm = lm(var_temp~ps_temp)
		out <- my_lm$coeff[2]
	}

	out
}

calc_t_stat_range <- function(replicated_data, point_spread, use_corr=FALSE)
{
	out <- 0

	point_spread_range <- (0:40)
	point_spread_range <- point_spread_range / 2

	temp <- array(NA, dim=length(point_spread_range) )
	for( i in 1:length(point_spread_range) )
	{
		games_of_interest <- which(point_spread == point_spread_range[i])

		if( length(games_of_interest) > 1)
		{
			temp[i] <- max(replicated_data[games_of_interest]) - min(replicated_data[games_of_interest])
		}
		else
		{
			temp[i] <- NaN
		}
	}

	point_spread_range <- point_spread_range[is.finite(temp)]
	temp <- temp[is.finite(temp)]
	
	if(use_corr)
	{
		out <- cor(temp, point_spread_range)
	}
	else
	{		
		my_lm = lm(temp~point_spread_range)
		out <- my_lm$coeff[2]
	}

	out
}

t_stat_var_simple_corr_obs <- calc_t_stat_var_simple(favorite-underdog,spread, TRUE)
t_stat_var_simple_lm_obs <- calc_t_stat_var_simple(favorite-underdog,spread, FALSE)
t_stat_equal_size_groups_corr_obs <- calc_t_stat_var_equal_size_groups(favorite-underdog,spread, TRUE)
t_stat_equal_size_groups_lm_obs <- calc_t_stat_var_equal_size_groups(favorite-underdog,spread, FALSE)
t_stat_range_corr_obs <- calc_t_stat_range(favorite-underdog,spread, TRUE)
t_stat_range_lm_obs <- calc_t_stat_range(favorite-underdog,spread, FALSE)

t_stat_var_simple_corr_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])
t_stat_var_simple_lm_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])
t_stat_equal_size_groups_corr_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])
t_stat_equal_size_groups_lm_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])
t_stat_range_corr_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])
t_stat_range_lm_rep <- array(NA,dim=dim(fav_score_minus_underdog_score_rep)[1])

for( i in 1:dim(fav_score_minus_underdog_score_rep)[1] )
{
	t_stat_var_simple_corr_rep[i] <- calc_t_stat_var_simple(fav_score_minus_underdog_score_rep[i,],spread, TRUE)
	t_stat_var_simple_lm_rep[i] <- calc_t_stat_var_simple(fav_score_minus_underdog_score_rep[i,],spread, FALSE)

	t_stat_equal_size_groups_corr_rep[i] <- calc_t_stat_var_equal_size_groups(fav_score_minus_underdog_score_rep[i,],spread, TRUE)
	t_stat_equal_size_groups_lm_rep[i] <- calc_t_stat_var_equal_size_groups(fav_score_minus_underdog_score_rep[i,],spread, FALSE)

	t_stat_range_corr_rep[i] <- calc_t_stat_range(fav_score_minus_underdog_score_rep[i,],spread, TRUE)
	t_stat_range_lm_rep[i] <- calc_t_stat_range(fav_score_minus_underdog_score_rep[i,],spread, FALSE)
}

p_val_var_simple_corr <- sum(t_stat_var_simple_corr_rep > t_stat_var_simple_corr_obs)/length(t_stat_var_simple_corr_rep)
p_val_var_simple_lm <- sum(t_stat_var_simple_lm_rep > t_stat_var_simple_lm_obs)/length(t_stat_var_simple_lm_rep)

p_val_equal_size_groups_corr <- sum(t_stat_equal_size_groups_corr_rep > t_stat_equal_size_groups_corr_obs)/length(t_stat_equal_size_groups_corr_rep)
p_val_equal_size_groups_lm <- sum(t_stat_equal_size_groups_lm_rep > t_stat_equal_size_groups_lm_obs)/length(t_stat_equal_size_groups_lm_rep)

p_val_range_corr <- sum(t_stat_range_corr_rep > t_stat_range_corr_obs)/length(t_stat_range_corr_rep)
p_val_range_lm <- sum(t_stat_range_lm_rep > t_stat_range_lm_obs)/length(t_stat_range_lm_rep)

par(mfrow=c(3,2))
hist(t_stat_var_simple_corr_rep,main='Correlation of Variance with Spread')
abline(v=t_stat_var_simple_corr_obs,col=2)
hist(t_stat_var_simple_lm_rep,main='Slope of Linear Fit of Variance with Spread')
abline(v=t_stat_var_simple_lm_obs,col=2)

hist(t_stat_equal_size_groups_corr_rep,main='Correlation of Variance with Spread (Equal Groups)')
abline(v=t_stat_equal_size_groups_corr_obs,col=2)
hist(t_stat_equal_size_groups_lm_rep,main='Slope of Linear Fit of Variance with Spread (Equal Groups)')
abline(v=t_stat_equal_size_groups_lm_obs,col=2)

hist(t_stat_range_corr_rep,main='Correlation of Range with Spread')
abline(v=t_stat_range_corr_obs,col=2)
hist(t_stat_range_lm_rep,main='Slope of Linear Fit of Range with Spread')
abline(v=t_stat_range_lm_obs,col=2)

print(p_val_var_simple_corr)
print(p_val_var_simple_lm)
print(p_val_equal_size_groups_corr)
print(p_val_equal_size_groups_lm)
print(p_val_range_corr)
print(p_val_range_lm)
