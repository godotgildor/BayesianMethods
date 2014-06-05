# Problem number 1
plane_data <- read.table("../data/planes.txt",header=T)
attach(plane_data)

# Subtract off the first year to make our math easier.
year <- year - year[1]

# Fit a linear model to our data.
my_linear_model = lm(fatal~year);
summary(my_linear_model)

# Plot the linear model and our actual data
plot(year, fatal, type="l", col=1, lwd=2, main="Exercise 1e: Linear Regression Fit",xlab="Year",ylab="Number of Fatalities")
lines(year, my_linear_model$fitted.values, type="l",col=2, lwd=2)

# Our crude estimate of alpha and beta
crude_alpha <- my_linear_model$coefficients[1];
crude_alpha_sd <- coef(summary(my_linear_model))[3]
crude_beta <- my_linear_model$coefficients[2];
crude_beta_sd <- coef(summary(my_linear_model))[4]

# Part b Sketch informative prior contours
num_grid_points <- 200
informative_prior_dist <- array(NA,dim=c(num_grid_points,num_grid_points))
my_alphas <- (ppoints(num_grid_points)-0.5)*crude_alpha_sd*4 + crude_alpha
my_betas <- (ppoints(num_grid_points)-0.5)*crude_beta_sd*4 + crude_beta 
for(i in 1:num_grid_points)
{
	for(j in 1:num_grid_points)
	{
		p_alpha <- dnorm(my_alphas[i], mean=crude_alpha, sd=crude_alpha_sd)
		p_beta <- dnorm(my_betas[j], mean=crude_beta,sd=crude_beta_sd)
		informative_prior_dist[i,j] = p_alpha*p_beta
	}
}

contour(my_alphas,my_betas,informative_prior_dist,main="Informative Prior for Alpha and Beta",xlab="alpha", ylab="Beta")


# We will use a non-informative prior, so our
# posterior will actually just be our likelihood.
posterior_prob <- function(alpha, beta)
{
	prob <- 0

	lambda <- alpha + beta * year;	

	if(min(lambda) < 0)
	{
		prob <- -Inf
	}
	else
	{
		for(i in 1:length(fatal))
		{
			prob <- prob + dpois(fatal[i], lambda[i], log=TRUE)
		}
	}

	prob
}

# How dense should our grid be?
num_grid_points <- 200
posterior_dist <- array(NA,dim=c(num_grid_points,num_grid_points))

# Make our sampling of alpha and beta values with each parameter
# centered at the crude estimate provided by our linear model.
alpha_ranges <- (ppoints(num_grid_points)-0.5)* 15 + crude_alpha
beta_ranges <- (ppoints(num_grid_points)-0.5) * 3 + crude_beta

# Create our grid of samples.
for(i in 1:num_grid_points)
{
	for(j in 1:num_grid_points)
	{
		posterior_dist[i,j] <- posterior_prob(alpha_ranges[i],beta_ranges[j])	
	}
}

# Our posterior distribution actually provided the log of the probabilities,
# not the actually probabilities, so take the inverse log now.
posterior_dist <- exp(posterior_dist - max(posterior_dist))
# Normalize the probabilities so they sum to 1.
posterior_dist <- posterior_dist / sum(posterior_dist)
# Plot the contour.
contour(alpha_ranges,beta_ranges,posterior_dist,main="Posterior Density using Non-informative Prior",xlab="alpha", ylab="Beta")

# Now calculate the marginal distribution of alpha
# by summing over all the betas for a given alpha
alpha_marginal <- rep(NA,num_grid_points)
for (i in 1:num_grid_points)
{
	alpha_marginal[i] <- sum(posterior_dist[i,])
}

# Now calculate the conditional distribution of beta
# which is just the posterior distribution divided by
# the alpha marginal distribution
# P(beta | alpha) = P(alpha, beta) / P(alpha)
beta_conditional <- matrix(NA,nrow=num_grid_points,ncol=num_grid_points)
for (i in 1:num_grid_points)
{
	for (j in 1:num_grid_points)
	{
		beta_conditional[i,j] <- posterior_dist[i,j]/alpha_marginal[i]
	}
}

# Part F
# Now we can sample
num_samples <- 1000
alpha_samp <- rep(NA,num_samples)
beta_samp <- rep(NA,num_samples)
for (m in 1:num_samples)
{
	curr_alpha_sample <- sample(1:num_grid_points,size=1,replace=F,prob=alpha_marginal)
	curr_beta_sample  <- sample(1:num_grid_points,size=1,replace=F,prob=beta_conditional[curr_alpha_sample,])
	alpha_samp[m] <- alpha_ranges[curr_alpha_sample]
	beta_samp[m] <- beta_ranges[curr_beta_sample]
}

# Now let's look at the samples that we drew
points(alpha_samp, beta_samp)
hist(alpha_samp)
hist(beta_samp)

# Get some stats about our distribution
mean(alpha_samp)
mean(beta_samp)

# Sort our alphas and betas so that we can
# calculate 95% intervals
alpha_samp_sorted <- sort(alpha_samp)
beta_samp_sorted <- sort(beta_samp)
lower_index <- round(num_samples * 0.025)
upper_index <- round(num_samples * 0.975)
alpha_samp_sorted[lower_index]
alpha_samp_sorted[upper_index]
beta_samp_sorted[lower_index]
beta_samp_sorted[upper_index]

# Part G
# Predict for next year
predicted_rate <- alpha_samp + beta_samp * (length(year) + 1)
hist(predicted_rate,main="Histogram of Predicted number of Fatal Accidents in 1986",xlab="Number of Fatal Accidents")

# Part H
num_accidents_next_year <- rep(NA,num_samples)
for (m in 1:num_samples)
{
	num_accidents_next_year[m] <- rpois(1, predicted_rate[m])
}
num_accidents_next_year_sorted <- sort(num_accidents_next_year)
num_accidents_next_year_sorted[lower_index]
num_accidents_next_year_sorted[upper_index]
detach(plane_data)
########################################################################################################
# Question 2
library(MASS)
radon_data <- read.table("../data/radon.csv",sep=",",header=TRUE)
# We could leave off the Goodhue data since Blue_Earth, Clay, and Goodhue are mutually exclusive.
# Thus, we only need to state 2 of the counties information. I'm not going to include it
# to avoid getting NA's in the various coefficient matrices etc.
first_floor <- 1 - radon_data$Basement
my_linear_model <- lm(log(radon_data$Radon_Level)~radon_data$Blue_Earth+radon_data$Clay+first_floor)
summary(my_linear_model)

qqnorm(my_linear_model$residuals)
qqline(my_linear_model$residuals)

# Now let's simulate sampling another house from
# Blue Earth county.  We'll do 2 samplings,
# one from the basement, one from the first floor.
new_Blue_Earth_basement <- t(c(1, 1, 0, 0))
new_Blue_Earth_first_floor <- t(c(1, 1, 0, 1))

# Get parameter values from our linear regression model.
beta_hat <- my_linear_model$coef
sigma <- summary(my_linear_model)$sigma
sigma_squared <- sigma^2
v_beta <- summary(my_linear_model)$cov.unscaled

# Basic parameters describing number of samples etc.
num_samples <- 5000
num_params <- length(beta_hat)
num_observations <- length(radon_data$Radon_Level)
num_degs_of_freedom <- num_observations - num_params

# Vars where we will store our sample beta and sigma vals.
sigma_squared_samples <- rep(NA,num_samples)
beta_samples <- matrix(NA,nrow=num_samples,ncol=num_params)

# Finally, the vars where we will store the predicted radon levels.
predicted_radon_level_Blue_Earth_basement <- rep(NA,num_samples)
predicted_radon_level_Blue_Earth_first_floor <- rep(NA,num_samples)

for (i in 1:num_samples)
{
	# Sample sigma squared, and then using this, sample beta.
	temp <- rchisq(1,num_degs_of_freedom)
	curr_sigma_squared <- num_degs_of_freedom * sigma_squared / temp
	curr_var_beta <- curr_sigma_squared * v_beta
	curr_beta <- mvrnorm(1,beta_hat,curr_var_beta)

	# Now use our sampled values of sigma squared and beta to 
	# calculate our new predicted radon values.
	new_sample_times_beta <- new_Blue_Earth_basement %*% t(t(curr_beta))
	predicted_radon_level_Blue_Earth_basement[i] <- rnorm(1,mean=new_sample_times_beta,sd=sqrt(curr_sigma_squared))

	# Now for first floor
	new_sample_times_beta <- new_Blue_Earth_first_floor %*% t(t(curr_beta))
	predicted_radon_level_Blue_Earth_first_floor[i] <- rnorm(1,mean=new_sample_times_beta,sd=sqrt(curr_sigma_squared))

	# Store away the sampled vals of sigma squared an beta, just for fun.
	sigma_squared_samples[i] <- curr_sigma_squared
	beta_samples[i,] <- curr_beta
}
# Exponentiate the values since we were actually predicting the log of the radon level.
predicted_radon_level_Blue_Earth_basement <- exp(predicted_radon_level_Blue_Earth_basement)
predicted_radon_level_Blue_Earth_first_floor<- exp(predicted_radon_level_Blue_Earth_first_floor)

# Calculate the 2.5% and 97.5% samples so that we can
# easily get our 95% interval.
lower <- round(0.025 * num_samples)
upper <- round(0.975 * num_samples);

sorted_vals <- sort(predicted_radon_level_Blue_Earth_basement)
sorted_vals[lower]
sorted_vals[upper]
sorted_vals <- sort(predicted_radon_level_Blue_Earth_first_floor)
sorted_vals[lower]
sorted_vals[upper]
########################################################################################################
# Baseball questions
# Question 3
baseball_data <- read.table("../data/hitters.post1975.txt",skip=1,sep=",",row.names=NULL)
at_bats <- baseball_data[,8]
# Take out entries that have less than 100 at bats and less than 1 home run.
baseball_data <- baseball_data[at_bats >=100,]
home_runs <- baseball_data[,13]
baseball_data <- baseball_data[home_runs >=1,]

# Extract fields of interest
years <- baseball_data[,3]
teams <- baseball_data[,5]
hits <- baseball_data[,10]
home_runs <- baseball_data[,13]
at_bats <- baseball_data[,8]
players <- baseball_data[,2]

batting_average <- hits / at_bats
hist(batting_average,main="Histogram of Batting Averages",xlab="Batting Average")
max(batting_average)
players[batting_average==max(batting_average)]
years[batting_average==max(batting_average)]
teams[batting_average==max(batting_average)]
min(batting_average)
players[batting_average==min(batting_average)]
years[batting_average==min(batting_average)]
teams[batting_average==min(batting_average)]
mean(batting_average)
avg_batting_average <- mean(batting_average)

# Question 4
num_grid_points <- 200
my_alphas = ppoints(num_grid_points)*500
my_betas = ppoints(num_grid_points)*1000
#num_grid_points <- 200
#my_alphas = (ppoints(num_grid_points)-0.5)*6+43
#my_betas = (ppoints(num_grid_points)-0.5)*12+120.5


my_log_likelihood <- array(NA,dim=c(num_grid_points,num_grid_points))
for(curr_alpha in 1:num_grid_points)
{
	for(curr_beta in 1:num_grid_points)
	{
		my_log_likelihood[curr_alpha,curr_beta] <- sum(dbeta(batting_average, my_alphas[curr_alpha], my_betas[curr_beta],log=TRUE))
	}
}
contour(my_alphas,my_betas,my_log_likelihood,main="Log Likelihood of Batting Average",xlab="Alpha", ylab="Beta")
# Question 5
max_indices <- which(my_log_likelihood==max(my_log_likelihood),arr.ind=TRUE)
max_alpha <- my_alphas[max_indices[1]]
max_beta <- my_betas[max_indices[2]]
expected_var <- max_alpha*max_beta/((max_alpha+max_beta+1)*(max_alpha+max_beta)^2)
var(batting_average)

# Question 6
beta_2d_nr <- function(old_alpha,old_beta,y)
{
	num_samples <- length(y)
	Jaa <- num_samples * trigamma(old_alpha + old_beta) - num_samples*trigamma(old_alpha)
	Jab <- num_samples * trigamma(old_alpha + old_beta)
	Jbb <- num_samples * trigamma(old_alpha + old_beta) - num_samples*trigamma(old_beta)
	J <- rbind(c(Jaa,Jab),c(Jab,Jbb))
	Jinv <- solve(J)
	ga <- num_samples * digamma(old_alpha + old_beta) - num_samples * digamma(old_alpha) + sum(log(y))
	gb <- num_samples * digamma(old_alpha + old_beta) - num_samples * digamma(old_beta) + sum(log(1-y))
	g_vector <- t(t(c(ga,gb)))
	old_vector <- t(t(c(old_alpha, old_beta)))
	new_vector <- old_vector - Jinv%*%g_vector
	new_vector
}

# What are our initial values?
init_alpha <- 2
init_beta <- 2
old_params <- c(init_alpha, init_beta)
iters_matrix <- NULL
iters_matrix <- rbind(iters_matrix, old_params)
new_params <- c(0,0)
diff <- sum(abs(new_params - old_params))
while (diff > 0.0000001)
{
	new_params <- beta_2d_nr(old_params[1],old_params[2],batting_average)
	iters_matrix <- rbind(iters_matrix,c(new_params[1],new_params[2]))
	diff <- sum(abs(new_params-old_params))
	old_params <- new_params
}

init_alpha <- 20
init_beta <- 50
old_params <- c(init_alpha, init_beta)
iters_matrix <- NULL
iters_matrix <- rbind(iters_matrix, old_params)
new_params <- c(0,0)
diff <- sum(abs(new_params - old_params))
while ( (diff > 0.0000001))
{
	new_params <- beta_2d_nr(old_params[1],old_params[2],batting_average)
	iters_matrix <- rbind(iters_matrix,c(new_params[1],new_params[2]))
	diff <- sum(abs(new_params-old_params))
	old_params <- new_params
}

# Question 7
home_run_average <- home_runs / at_bats
hist(home_run_average,main="Histogram of Home Run Averages",xlab="Home Run Average")
max(home_run_average)
players[home_run_average==max(home_run_average)]
years[home_run_average==max(home_run_average)]
teams[home_run_average==max(home_run_average)]
min(home_run_average)
players[home_run_average==min(home_run_average)]
years[home_run_average==min(home_run_average)]
teams[home_run_average==min(home_run_average)]
mean(home_run_average)
avg_home_run_average <- mean(home_run_average)


# Do grid sampling for log likelihood.
num_grid_points <- 100
my_alphas = (ppoints(num_grid_points)-0.5)*2 + 2
my_betas = (ppoints(num_grid_points)-0.5)*50+70

my_log_likelihood <- array(NA,dim=c(num_grid_points,num_grid_points))
for(curr_alpha in 1:num_grid_points)
{
	for(curr_beta in 1:num_grid_points)
	{
		my_log_likelihood[curr_alpha,curr_beta] <- sum(dbeta(home_run_average, my_alphas[curr_alpha], my_betas[curr_beta],log=TRUE))
	}
}
contour(my_alphas,my_betas,my_log_likelihood,main="Log Likelihood of Home Run Average",xlab="Alpha", ylab="Beta")

max_indices <- which(my_log_likelihood==max(my_log_likelihood),arr.ind=TRUE)
max_alpha <- my_alphas[max_indices[1]]
max_beta <- my_betas[max_indices[2]]
expected_var <- max_alpha*max_beta/((max_alpha+max_beta+1)*(max_alpha+max_beta)^2)
var(home_run_average)

init_alpha <- 1
init_beta <- 1
old_params <- c(init_alpha, init_beta)
iters_matrix <- NULL
iters_matrix <- rbind(iters_matrix, old_params)
new_params <- c(0,0)
diff <- sum(abs(new_params - old_params))
while (diff > 0.0000001)
{
	new_params <- beta_2d_nr(old_params[1],old_params[2],home_run_average)
	iters_matrix <- rbind(iters_matrix,c(new_params[1],new_params[2]))
	diff <- sum(abs(new_params-old_params))
	old_params <- new_params
}

init_alpha <- 1
init_beta <- 50
old_params <- c(init_alpha, init_beta)
iters_matrix <- NULL
iters_matrix <- rbind(iters_matrix, old_params)
new_params <- c(0,0)
diff <- sum(abs(new_params - old_params))
while ( (diff > 0.0000001))
{
	new_params <- beta_2d_nr(old_params[1],old_params[2],home_run_average)
	iters_matrix <- rbind(iters_matrix,c(new_params[1],new_params[2]))
	diff <- sum(abs(new_params-old_params))
	old_params <- new_params
}

# Question 8
library(msm)
mixture_log_likelihood <- function(alpha, mu0, sigma0, sigma1, x)
{
	
	elite_dist <- dnorm(x, mu0, sigma0)
	std_dist <- dtnorm(x, mean=0, sd=sigma1, lower=0)
	log_likelihood <- sum(log(alpha*elite_dist + (1-alpha)*std_dist))
	log_likelihood
}

# Expectation function
expectation_step <- function(alpha, mu0, sigma0, sigma1, x)
{
	num_samples <- length(x)
	indicator_vars <- rep(NA,num_samples)

	for (i in 1:num_samples)
	{
		prob0 <- alpha * dnorm(x[i], mean=mu0, sd=sigma0)
		prob1 <- (1 - alpha) * dtnorm(x[i], mean=0, sd=sigma1, lower=0)
		indicator_vars[i] <- prob0/(prob0+prob1)
	}
	indicator_vars
}

# Maximization Function
maximization_step <- function(indicator_vars, x)
{
	alpha <- mean(indicator_vars)
	mu0 <- sum(indicator_vars * x) / sum(indicator_vars)

	sigma0 <- sqrt(sum(indicator_vars*((x-mu0)^2))/sum(indicator_vars))
	sigma1 <- sqrt(sum((1-indicator_vars)*(x^2))/sum(1-indicator_vars))

	c(alpha, mu0, sigma0, sigma1)
}

## Starting values for EM algorithm:
# Both sets of the following init vals produce similar results.
#curr_alpha <- 0.5
#curr_mu0 <- .1
#curr_sigma0 <- .5
#curr_sigma1 <- .5
curr_alpha <- 0.95
curr_mu0 <- .05
curr_sigma0 <- .02
curr_sigma1 <- .02

iters_matrix <- NULL
iters_matrix <- c(curr_alpha, curr_mu0, curr_sigma0, curr_sigma1)
diff <- 1
num_iters <- 1
while (diff > 0.000001)
{
	num_iters <- num_iters+1
	curr_indicator_vars <- expectation_step(curr_alpha, curr_mu0, curr_sigma0, curr_sigma1, home_run_average)
	curr_params <- maximization_step(curr_indicator_vars, home_run_average)
	curr_alpha <- curr_params[1]
	curr_mu0 <- curr_params[2]
	curr_sigma0 <- curr_params[3]
	curr_sigma1 <- curr_params[4]
	iters_matrix <- rbind(iters_matrix, curr_params)
	diff <- sum(abs(iters_matrix[num_iters,]-iters_matrix[num_iters-1,])/abs(iters_matrix[num_iters,]))
	print (num_iters)
}

# Question 9
hist(home_run_average,prob=T,main="Histogram of Home Run Averages",xlab="Home Run Average")
my_points <- ppoints(1000)*max(home_run_average)
elite_density <- curr_alpha * dnorm(my_points, curr_mu0, curr_sigma0)
std_density <- (1 - curr_alpha) * dtnorm(my_points, mean=0, sd=curr_sigma1, lower=0)
lines(my_points, elite_density)
lines(my_points, std_density)
lines(my_points, elite_density+std_density)

# Question 10
player_indices <- which(players=="abreubo01",arr.ind=TRUE)
curr_indicator_vars[player_indices]

# Question 11

# Expectation function
expectation_step <- function(alpha, mu0, sigma, x)
{
	num_samples <- length(x)
	indicator_vars <- rep(NA,num_samples)

	for (i in 1:num_samples)
	{
		prob0 <- alpha * dnorm(x[i], mean=mu0, sd=sigma)
		prob1 <- (1 - alpha) * dtnorm(x[i], mean=0, sd=sigma, lower=0)
		indicator_vars[i] <- prob0/(prob0+prob1)
	}
	indicator_vars
}

# Maximization Function
maximization_step <- function(indicator_vars, x)
{
	alpha <- mean(indicator_vars)
	mu0 <- sum(indicator_vars * x) / sum(indicator_vars)

	sigma_squared <- (sum(indicator_vars * ((x-mu0)^2)) + sum(x^2*(1-indicator_vars))) / length(indicator_vars)
	sigma <- sqrt(sigma_squared)			   

	c(alpha, mu0, sigma)
}

## Starting values for EM algorithm:
curr_alpha <- 0.95
curr_mu0 <- .05
curr_sigma <- .02
# The following values converge to an alpha of near 0.
#curr_alpha <- 0.25
#curr_mu0 <- .1
#curr_sigma <- .02

iters_matrix <- NULL
iters_matrix <- c(curr_alpha, curr_mu0, curr_sigma)
diff <- 1
num_iters <- 1
while (diff > 0.000001)
{
	num_iters <- num_iters+1
	curr_indicator_vars <- expectation_step(curr_alpha, curr_mu0, curr_sigma, home_run_average)
	curr_params <- maximization_step(curr_indicator_vars, home_run_average)
	curr_alpha <- curr_params[1]
	curr_mu0 <- curr_params[2]
	curr_sigma <- curr_params[3]
	iters_matrix <- rbind(iters_matrix, curr_params)
	diff <- sum(abs(iters_matrix[num_iters,]-iters_matrix[num_iters-1,])/abs(iters_matrix[num_iters,]))
	print (num_iters)
}

hist(home_run_average,prob=T,main="Histogram of Home Run Averages",xlab="Home Run Average")
my_points <- ppoints(1000)*max(home_run_average)
elite_density <- curr_alpha * dnorm(my_points, curr_mu0, curr_sigma)
std_density <- (1 - curr_alpha) * dtnorm(my_points, mean=0, sd=curr_sigma, lower=0)
lines(my_points, elite_density)
lines(my_points, std_density)
lines(my_points, elite_density+std_density)

player_indices <- which(players=="abreubo01",arr.ind=TRUE)
curr_indicator_vars[player_indices]

# Question 12 
# Newton-Raphson for Gamma
gamma_NR <- function(old_a, old_b, indicator_vars, y)
{	
	num_samples <- length(y)
	one_minus_indicator_vars <- 1 - indicator_vars
	sum_one_minuse_indicator_vars <- sum(one_minus_indicator_vars)
	Jaa <- -sum_one_minuse_indicator_vars * trigamma(old_a)
	Jab <- sum_one_minuse_indicator_vars / old_b
	Jbb <- -sum_one_minuse_indicator_vars * old_a / (old_b^2)
	J <- rbind(c(Jaa,Jab),c(Jab,Jbb))
	Jinv <- solve(J)
	ga <- sum_one_minuse_indicator_vars*log(old_b) - sum_one_minuse_indicator_vars*digamma(old_a) + sum(one_minus_indicator_vars*log(y))
	gb <- sum_one_minuse_indicator_vars*old_a / old_b - sum(one_minus_indicator_vars*y)
	g_vector <- t(t(c(ga,gb)))
	old_vector <- t(t(c(old_a, old_b)))
	new_vector <- old_vector - Jinv%*%g_vector

	new_vector
}


# Expectation function
expectation_step <- function(alpha, mu, sigma, a, b, x)
{
	num_samples <- length(x)
	indicator_vars <- rep(NA,num_samples)

	for (i in 1:num_samples)
	{
		prob0 <- alpha * dnorm(x[i], mean=mu, sd=sigma)
		prob1 <- (1 - alpha) * dgamma(x[i], a, rate=b)
		indicator_vars[i] <- prob0/(prob0+prob1)
	}
	indicator_vars
}

# Maximization Function
maximization_step <- function(indicator_vars, old_a, old_b, x)
{
	alpha <- mean(indicator_vars)
	mu <- sum(indicator_vars * x) / sum(indicator_vars)
	sigma <- sqrt(sum(indicator_vars*((x-mu)^2))/sum(indicator_vars))

	#####################################################################
	# Run Newton-Raphson to find new estimates for gamma params	
	diff <- 1
	
	old_params <- c(old_a, old_b) # starting values
	while (diff > 0.0000001)
	{
		new_params <- gamma_NR(old_params[1], old_params[2], indicator_vars, x)
		diff <- sum(abs(new_params - old_params))
		old_params <- new_params
	}

	a <- new_params[1]
	b <- new_params[2]

	c(alpha, mu, sigma, a, b)
}

## Starting values for EM algorithm:
curr_alpha <- 0.95
curr_mu <- .05
curr_sigma <- .02
curr_a <- 1
curr_b <- 100
iters_matrix <- NULL
iters_matrix <- c(curr_alpha, curr_mu, curr_sigma, curr_a, curr_b)
diff <- 1
num_iters <- 1
while (diff > 0.000001)
{
	num_iters <- num_iters+1
	curr_indicator_vars <- expectation_step(curr_alpha, curr_mu, curr_sigma, curr_a, curr_b, home_run_average)
	curr_params <- maximization_step(curr_indicator_vars, curr_a, curr_b, home_run_average)
	curr_alpha <- curr_params[1]
	curr_mu0 <- curr_params[2]
	curr_sigma <- curr_params[3]
	curr_a <- curr_params[4]
	curr_b <- curr_params[5]
	iters_matrix <- rbind(iters_matrix, curr_params)
	diff <- sum(abs(iters_matrix[num_iters,]-iters_matrix[num_iters-1,])/abs(iters_matrix[num_iters,]))
	print (num_iters)
}

my_points <- ppoints(1000)*max(home_run_average)
elite_density <- curr_alpha * dnorm(my_points, curr_mu0, curr_sigma)
std_density <- (1 - curr_alpha) * dgamma(my_points, curr_a, rate=curr_b)
hist(home_run_average,prob=T,main="Histogram of Home Run Averages",xlab="Home Run Average",ylim=c(0, 30))
lines(my_points, std_density)
lines(my_points, elite_density+std_density)
lines(my_points, elite_density)


player_indices <- which(players=="abreubo01",arr.ind=TRUE)
curr_indicator_vars[player_indices]
