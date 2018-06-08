#prepare Markov chains for rainfall amount modelling
#calculate X'(t)
rain.yes <- (rain.ts > 0.05)
delta <- min(rain.ts[rain.ts > 0.05])
rain.ts <- rain.ts - delta


nj.t <- rep(0, 366)
X <- vector("list", 366)


for (i in 1:13880){
  
  if (rain.yes[i] == TRUE){
    
    nj.t[days.366[i]] <- nj.t[days.366[i]] + 1
    X[[days.366[i]]] <- c(X[[days.366[i]]], rain.ts[i])
    
  }
  
}

t.prime <- 2 * pi * c(1:366)/ 366
#we use the sufficient stats xbar

xbar <- rep(0, 366)

for (i in 1:366){
  
  xbar[i] <- sum(X[[i]]) / length(X[[i]])
  
}

#pooling with 5 days
xpool <- rep(0, 73)
npool <- rep(0, 73)

for (i in 0:72){
  
  times <- X[(i*5 + 1):(i*5 + 5)]
  somme <- 0
  nb <- 0
  for (j in 1:5){
    
    somme <- somme + sum(times[[j]])
    nb <- nb + length(times[[j]])
    
  }
  
  npool[i + 1] <- nb
  xpool[i + 1] <- somme / nb
  
}

t.pool <-  2 * pi * c(1:73)/ 73

#####Markov Chains

nmarkov11 <- rep(0, 366)
Xmarkov11 <- vector("list", 366)
nmarkov01 <- rep(0, 366)
Xmarkov01 <- vector("list", 366)
nmarkov10 <- rep(0, 366)
Xmarkov10 <- vector("list", 366)
nmarkov00 <- rep(0, 366)
Xmarkov00 <- vector("list", 366)
rain11 <- rep(0, 13880)
rain01 <- rep(0, 13880)
rain10 <- rep(0, 13880)
rain00 <- rep(0, 13880)

#(h,i)
for (i in 3:13880){
  
  if (rain.yes[(i-2)] == TRUE & rain.yes[(i-1)] == TRUE & rain.yes[i] == TRUE){
    rain11[i] = 1
  }
  if (rain.yes[(i-2)] == FALSE & rain.yes[(i-1)] == TRUE & rain.yes[i] == TRUE){
    rain01[i] = 1
  }
  if (rain.yes[(i-2)] == TRUE & rain.yes[(i-1)] == FALSE & rain.yes[i] == TRUE){
    rain10[i] = 1
  }
  if (rain.yes[(i-2)] == FALSE & rain.yes[(i-1)] == FALSE & rain.yes[i] == TRUE){
    rain00[i] = 1
  }
  
}

for (i in 3:13880){
  
  if (rain11[i] == 1){
    nmarkov11[days.366[i]] <- nmarkov11[days.366[i]] + 1
    Xmarkov11[[days.366[i]]] <- c(Xmarkov11[[days.366[i]]], rain.ts[i])
  }
  
  if (rain10[i] == 1){
    nmarkov10[days.366[i]] <- nmarkov10[days.366[i]] + 1
    Xmarkov10[[days.366[i]]] <- c(Xmarkov10[[days.366[i]]], rain.ts[i])
  }
  
  if (rain01[i] == 1){
    nmarkov01[days.366[i]] <- nmarkov01[days.366[i]] + 1
    Xmarkov01[[days.366[i]]] <- c(Xmarkov01[[days.366[i]]], rain.ts[i])
  }
  
  if (rain00[i] == 1){
    nmarkov00[days.366[i]] <- nmarkov00[days.366[i]] + 1
    Xmarkov00[[days.366[i]]] <- c(Xmarkov00[[days.366[i]]], rain.ts[i])
  }
  
}

#treat missing values
#missingval <- c(44, 60, 127, 156, 173, 180, 190, 222, 340, 356)

#for (i in missingval){
#  Xmarkov00[i] <- 0

#}

#we use the sufficient stats xbar

xbar11 <- rep(0, 366)
xbar01 <- rep(0, 366)
xbar10 <- rep(0, 366)
xbar00 <- rep(0, 366)

for (i in 1:366){
  
  xbar11[i] <- sum(Xmarkov11[[i]]) / length(Xmarkov11[[i]])
  xbar01[i] <- sum(Xmarkov01[[i]]) / length(Xmarkov01[[i]])
  xbar10[i] <- sum(Xmarkov10[[i]]) / length(Xmarkov10[[i]])
  xbar00[i] <- sum(Xmarkov00[[i]]) / length(Xmarkov00[[i]])
  
}

#pooling with 5 days
xpool11 <- rep(0, 73)
npool11 <- rep(0, 73)
xpool01 <- rep(0, 73)
npool01 <- rep(0, 73)
xpool10 <- rep(0, 73)
npool10 <- rep(0, 73)
xpool00 <- rep(0, 73)
npool00 <- rep(0, 73)

for (i in 0:72){
  
  times11 <- Xmarkov11[(i*5 + 1):(i*5 + 5)]
  somme11 <- 0
  nb11 <- 0
  
  times01 <- Xmarkov01[(i*5 + 1):(i*5 + 5)]
  somme01 <- 0
  nb01 <- 0
  
  times10 <- Xmarkov10[(i*5 + 1):(i*5 + 5)]
  somme10 <- 0
  nb10 <- 0
  
  times00 <- Xmarkov00[(i*5 + 1):(i*5 + 5)]
  somme00 <- 0
  nb00 <- 0
  
  
  for (j in 1:5){
    
    somme11 <- somme11 + sum(times11[[j]])
    nb11 <- nb11 + length(times11[[j]])
    somme01 <- somme01 + sum(times01[[j]])
    nb01 <- nb01 + length(times01[[j]])
    somme10 <- somme10 + sum(times10[[j]])
    nb10 <- nb10 + length(times10[[j]])
    somme00 <- somme00 + sum(times00[[j]])
    nb00 <- nb00 + length(times00[[j]])
    
  }
  
  npool11[i + 1] <- nb11
  xpool11[i + 1] <- somme11 / nb11
  npool01[i + 1] <- nb01
  xpool01[i + 1] <- somme01 / nb01
  npool10[i + 1] <- nb10
  xpool10[i + 1] <- somme10 / nb10
  npool00[i + 1] <- nb00
  xpool00[i + 1] <- somme00 / nb00
  
}