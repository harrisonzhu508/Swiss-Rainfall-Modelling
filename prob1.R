#prepare the second-order Markov chain data
rain.yes <- (rain.ts > 0.05)
rain.no <- (rain.yes <= 0.05)

#model responses
rain.model.yes <- rep(0, 366)
rain.model.no <- rep(0, 366)

#yes response
rain.00.1 <- rep(0, 366)
rain.01.1 <- rep(0, 366)
rain.10.1 <- rep(0, 366)
rain.11.1 <- rep(0, 366)

#all responses
rain.00 <- rep(0, 366)
rain.01 <- rep(0, 366)
rain.10 <- rep(0, 366)
rain.11 <- rep(0, 366)

#simple rain or no rain
count <- 1

for (i in 1:13880){
  c <- 0
  
  if (count %% 365 == 0){
    count <- 1
  }
  
  if ( rain.yes[i] == TRUE) 
  {
    
    rain.model.yes[days.366[i]] = rain.model.yes[days.366[i]] + 1
  }
  else{
    
    rain.model.no[days.366[i]] = rain.model.no[days.366[i]] + 1
    
  }
  
  count <- count + 1
}

# (h, i) = (0, 0)
count <- 1

for (i in 3:13880){
  c <- 0
  
  if (count %% 365 == 0){
    count <- 1
  }
  
  if ( (rain.yes[i - 2] == FALSE) & (rain.yes[i - 1] == FALSE) ) 
  {
    rain.00[days.366[i]] = rain.00[days.366[i]] + 1
    if ( rain.yes[i] == TRUE){
      rain.00.1[days.366[i]] = rain.00.1[days.366[i]] + 1
    }
  }
  
  count <- count + 1
}


# (h, i) = (0, 1)
count <- 1

for (i in 3:13880){
  c <- 0
  
  if (count %% 365 == 0){
    count <- 1
  }
  
  if ( (rain.yes[i - 2] == FALSE) & (rain.yes[i - 1] == TRUE) ) 
  {
    rain.01[days.366[i]] = rain.01[days.366[i]] + 1
    if ( rain.yes[i] == TRUE){
      rain.01.1[days.366[i]] = rain.01.1[days.366[i]] + 1
    }
  }
  
  count <- count + 1
}

# (h, i) = (1, 0)
count <- 1

for (i in 3:13880){
  c <- 0
  
  if (count %% 365 == 0){
    count <- 1
  }
  
  if ( (rain.yes[i - 2] == TRUE) & (rain.yes[i - 1] == FALSE) ) 
  {
    rain.10[days.366[i]] = rain.10[days.366[i]] + 1
    if ( rain.yes[i] == TRUE){
      rain.10.1[days.366[i]] = rain.10.1[days.366[i]] + 1
    }
  }
  
  count <- count + 1
}
# (h, i) = (1, 1)
count <- 1

for (i in 3:13880){
  c <- 0
  
  if (count %% 365 == 0){
    count <- 1
  }
  
  if ( (rain.yes[i - 2] == TRUE) & (rain.yes[i - 1] == TRUE) ) 
  {
    rain.11[days.366[i]] = rain.11[days.366[i]] + 1
    if ( rain.yes[i] == TRUE){
      rain.11.1[days.366[i]] = rain.11.1[days.366[i]] + 1
    }
  }
  
  count <- count + 1
}


#pool the observations

rain.p11 <- rep(0, 73) 
rain.p01 <- rep(0, 73) 
rain.p10 <- rep(0, 73)
rain.p00 <- rep(0, 73)

rain.p11.1 <- rep(0, 73)
rain.p01.1 <- rep(0, 73)
rain.p10.1 <- rep(0, 73)
rain.p00.1 <- rep(0, 73)

for (i in 1:73){
  
  rain.p11[i] <- ceiling(sum(rain.11[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p01[i] <- ceiling(sum(rain.01[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p10[i] <- ceiling(sum(rain.10[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p00[i] <- ceiling(sum(rain.00[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  
  rain.p11.1[i] <- ceiling(sum(rain.11.1[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p01.1[i] <- ceiling(sum(rain.01.1[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p10.1[i] <- ceiling(sum(rain.10.1[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  rain.p00.1[i] <- ceiling(sum(rain.00.1[(5 * (i - 1) + 1): (5*(i-1) + 5)]) / 5)
  
}

#proportions r_{hi}
r.11 <- rain.11.1 / rain.11
r.10 <- rain.10.1 / rain.10
r.01 <- rain.01.1 / rain.01
r.00 <- rain.00.1 / rain.00

r.p11 <- rep(0, 73)
r.p10 <- rep(0, 73)
r.p01 <- rep(0, 73)
r.p00 <- rep(0, 73)

#proportions pooled over 5 days

for (i in 1:73){
  
  r.p11[i] <- sum(r.11[(5 * (i - 1)): (5*(i-1) + 4)]) / 5
  r.p10[i] <- sum(r.10[(5 * (i - 1)): (5*(i-1) + 4)]) / 5
  r.p01[i] <- sum(r.01[(5 * (i - 1)): (5*(i-1) + 4)]) / 5
  r.p00[i] <- sum(r.00[(5 * (i - 1)): (5*(i-1) + 4)]) / 5
  
}

rain.hist <- rep(0,366)
rain.hist.no <- rep(0,366)
#non-Markov chain frequency
for (i in 1:13880){
  rain.hist[days.366[i]] <- rain.hist[days.366[i]] + as.numeric(rain.yes[i])
  rain.hist.no[days.366[i]] <- rain.hist.no[days.366[i]] + (1 - as.numeric(rain.yes[i]))
}

