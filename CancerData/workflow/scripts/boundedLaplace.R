## Bounded Laplace

set.seed(8180)

## Acceptance Rejection Method for Bounded Laplce Density
boundlap <- function(x, mu, b, C) {
  if (x >0) {
    return(1/(2*C*b)*exp(-abs(x-mu)/b))
  } else {
    return(0)
  }
}

acc_rej <- function(n, mu, b, C) {
  accepted <- 1
  x <- numeric(n)
  while (accepted <= n) {
    y <- rlaplace(1, mu, b)
    u <- runif(1)
    if (u <= C*boundlap(y, mu, b, C)/(dlaplace(y, mu, b))) {
      x[accepted] <- y
      accepted = accepted + 1
    } 
  }
  return(x)
}

testl <- acc_rej(10000, 0, 1, 1/2)
hist(testl, breaks=seq(-20, 20, 0.25), xlab = '', ylab= '', main='Laplace (mu=0, sigma=1)')
hist(rlaplace(10000, 0, 1))
