# Problem number 1a
theta <- (ppoints(20000)-.5)*25
theta.sampdist.mean1 <- dnorm(theta,mean=1,sd=2)
theta.sampdist.mean2 <- dnorm(theta, mean=2, sd=2)
theta.sampdist.combined <- (theta.sampdist.mean1 + theta.sampdist.mean2) /2
plot(theta,theta.sampdist.combined ,type="l",col=1,lwd=2,main="Excercise 1",xlab="y", ylab="Density")
lines(theta,theta.sampdist.mean1,type="l",col=2,lwd=2)
lines(theta,theta.sampdist.mean2,type="l",col=3,lwd=2)
legend(locator(1), c("Marginal probability density","Normal with mean = 1","Normal with mean = 2"),col=c(1,2,3),lwd=2)

# Problem 1c
y <- (ppoints(20000)-.5)*20
num <- 1
sigma <- 1
denom1 <- 1 + exp(-1/(2*sigma^2)*((y-2)^2-(y-1)^2))
z1 = num/denom1
sigma <- 2
denom2 <- 1 + exp(-1/(2*sigma^2)*((y-2)^2-(y-1)^2))
z2 = num/denom2
sigma <-3
denom3 <- 1 + exp(-1/(2*sigma^2)*((y-2)^2-(y-1)^2))
z3 = num/denom3
plot(y,z1,type="l",col=1,lwd=2,main="Excercise 1c",xlab="y", ylab="P(theta = 1 |y)")
lines(y, z2,type="l",col=2,lwd=2)
lines(y, z3,type="l",col=3,lwd=2)
legend(locator(1), c("sigma = 1","sigma = 2","sigma = 3"),col=c(1,2,3),lwd=2)

#Problem 5
# Set a variable for the time parameter for people coming into the clinic
time_parameter <- 10
# First calculate the number of minutes the clinic is open (From 9 until 4 currently)
minutes_open <- (16-9)*60

# Init statistics that we'll track
total_num_patients <- 0
num_patients_have_to_wait <- 0
avg_wait <- 0
office_close <- 0

num_sims <- 100

# Run simulations
for (curr_sim in 1:num_sims)
{
	# Now calculate the arrival times of each person from open to close.
	total_num_patients[curr_sim] <- 1
	#arrive_time <- rpois(1,time_parameter)
	arrive_time <- rexp(1,rate=1/time_parameter)
	while(arrive_time[total_num_patients[curr_sim]] < minutes_open)
	{
		total_num_patients[curr_sim] <- total_num_patients[curr_sim] + 1
		#arrive_time[total_num_patients[curr_sim]] <- arrive_time[total_num_patients[curr_sim]-1] + rpois(1,time_parameter)
		arrive_time[total_num_patients[curr_sim]] <- arrive_time[total_num_patients[curr_sim]-1] + rexp(1,rate=1/time_parameter)
	}
	# The last person exceeded our total time, so take them off
	total_num_patients[curr_sim] <- total_num_patients[curr_sim] - 1
	arrive_time <- arrive_time[1:total_num_patients[curr_sim]]

	# Now calculate the time each patient spends with the doctor.
	time_with_doctor <- runif(total_num_patients[curr_sim], min=5, max=20)
	# Keep a vector of when each doctor is free.
	dr <- c(0,0,0)
	# Keep a vector of wait times
	wait_times <- c(0)
	num_patients_have_to_wait[curr_sim] <- 0
	for (curr_patient in 1:total_num_patients[curr_sim])
	{
		# Get the earliest available doctor.
		curr_doctor <- which(dr==min(dr))
		# In case of a tie, always use the lower numbered doctor.
		curr_doctor <- curr_doctor[1]
		# Does this patient have to wait?
		if(dr[curr_doctor] > arrive_time[curr_patient])
		{
			num_patients_have_to_wait[curr_sim] <- num_patients_have_to_wait[curr_sim] + 1
			wait_times[num_patients_have_to_wait[curr_sim]] <- dr[curr_doctor] - arrive_time[curr_patient]
			# Book the time with this doctor.
			dr[curr_doctor] <- dr[curr_doctor] + time_with_doctor[curr_patient]

		}
		else
		{
			# Book the time with this doctor.
			dr[curr_doctor] <- arrive_time[curr_patient] + time_with_doctor[curr_patient]
		}
	
	}
	# What was the average wait time? (Just those that had to wait)
	if(sum(wait_times) > 0)
	{
		avg_wait[curr_sim] <- sum(wait_times)/num_patients_have_to_wait[curr_sim]
	}
	else
	{
		avg_wait[curr_sim] <- 0
	}
	# When did office close?  This time is given in minutes after 4:00pm
	office_close[curr_sim] <- max(max(dr) - minutes_open, 0)
}

# Print stats.
test <- sort(total_num_patients)
print(median(total_num_patients))
print(test[round(.25*num_sims)])
print(test[round(.75*num_sims)])
test <- sort(num_patients_have_to_wait)
print(median(num_patients_have_to_wait))
print(test[round(.25*num_sims)])
print(test[round(.75*num_sims)])
test <- sort(avg_wait)
print(median(avg_wait))
print(test[round(.25*num_sims)])
print(test[round(.75*num_sims)])
test <- sort(office_close)
print(median(office_close))
print(test[round(.25*num_sims)])
print(test[round(.75*num_sims)])

# Problem 8
#Part a
# Sketch a Beta distribution with alpha = 1 and beta = 2/3
x = ppoints(1000)
my_beta = dbeta(x, shape1=1, shape2 = 2/3)
plot(x,my_beta ,type="l",col=1,lwd=2,main="Excercise 8",xlab="theta", ylab="Density")

#Part b
# Sketch a Beta distribution with alpha = 651 and beta = 350+2/3
x = ppoints(1000)
my_beta = dbeta(x, shape1=651, shape2 = (350+2/3))
plot(x,my_beta ,type="l",col=1,lwd=2,main="Excercise 8 part b",xlab="theta", ylab="Density")

#Problem 9
#Some initial setup
plane_data <- read.table("../data/planes.txt",header=T)
attach(plane_data)
num_samples_of_theta = 10000

#Part a
num_years <- length(fatal)
sum_fatal_accidents <- sum(fatal)
theta = ppoints(1000)*50
posterior_dist = dgamma(theta, shape=(sum_fatal_accidents + 0.0001), rate = (num_years+.0001))
# Now that we have our posterior distribution, take samples of our theta (or lambda since its a poisson)
samples_of_theta = rgamma(num_samples_of_theta, shape=(sum_fatal_accidents + 0.0001), rate = (num_years+.0001))
num_accidents_in_1986 = 0;
# Now loop over all of our values of theta, and draw a value for the number of fatal accidents we expect.
for (curr_sample in 1:num_samples_of_theta)
{
	num_accidents_in_1986[curr_sample] = rpois(1,samples_of_theta[curr_sample])
}
quantile(num_accidents_in_1986,0.025)
quantile(num_accidents_in_1986,0.975)

#Part b
# Get the passenger_miles (in 100 millions)
passenger_miles = deaths/rate
sum_passenger_miles = sum(passenger_miles)

# Calculate the posterior dist and then sample to get values of theta (lambda)
posterior_dist = dgamma(theta, shape=(sum_fatal_accidents + 0.0001), rate = (sum_passenger_miles+.0001))
samples_of_theta = rgamma(num_samples_of_theta, shape=(sum_fatal_accidents + 0.0001), rate = (sum_passenger_miles+.0001))

# Now do a poisson draw for each of our value of theta and calculate the 95% interval for our posterior predictive.
num_accidents_in_1986 = 0;
for (curr_sample in 1:num_samples_of_theta)
{
	num_accidents_in_1986[curr_sample] = rpois(1,8000*samples_of_theta[curr_sample])
}
quantile(num_accidents_in_1986,0.025)
quantile(num_accidents_in_1986,0.975)

#Part c
sum_deaths <- sum(deaths)
theta <- ppoints(1000)*50
posterior_dist <- dgamma(theta, shape=(sum_deaths + 0.0001), rate = (num_years+.0001))
samples_of_theta <- rgamma(10000, shape=(sum_deaths + 0.0001), rate = (num_years+.0001))

num_accidents_in_1986 = 0;
for (curr_sample in 1:num_samples_of_theta)
{
	num_accidents_in_1986[curr_sample] = rpois(1,samples_of_theta[curr_sample])
}
quantile(num_accidents_in_1986,0.025)
quantile(num_accidents_in_1986,0.975)

#Part d
millions_passenger_miles <- deaths/rate
sum_millions_passenger_miles <- sum(millions_passenger_miles)
posterior_dist <- dgamma(theta, shape=(sum_deaths + 0.0001), rate = (sum_millions_passenger_miles+.0001))
samples_of_theta <- rgamma(10000, shape=(sum_deaths + 0.0001), rate = (sum_millions_passenger_miles+.0001))
num_accidents_in_1986 = 0;
for (curr_sample in 1:num_samples_of_theta)
{
	num_accidents_in_1986[curr_sample] = rpois(1,8000*samples_of_theta[curr_sample])
}
quantile(num_accidents_in_1986,0.025)
quantile(num_accidents_in_1986,0.975)
