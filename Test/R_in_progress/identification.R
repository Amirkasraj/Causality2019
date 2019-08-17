library(hdrcde)
library(KernSmooth)
library(Rlab)

generate.coditions <- function(Xs, k = 7) {
  # ToDo interval should be small ?! use hdr
  m = mean(Xs[[1]])
  s = sqrt(var(Xs[[1]]))
  seq(from = m-1*s, to = m+1*s, length.out = k)
}

generate.pdf <- function(xs,ys) {
  l=min(xs)
  r=max(xs)
  d = splinefun(xs,ys)
  i = integrate(d,lower = l ,upper = r)$value
  single <- function(x) {
    y = 0
    if (x<l) y = 0
    else if (x>r) y = 0
    else y = d(x)
    y/i
  }
  function(x) {
    vapply(FUN = single, X=x, FUN.VALUE = 1)
  }
}

w <-function(pdfs,x) {
  single <- function(x) {
    ans = NULL
    for (d in pdfs)
      ans = c(ans, d(x))
    ans
  }
  lapply(FUN = single,X = x)#,FUN.VALUE = double(length = length(pdfs)))
}

generate.C <- function(pdfs, l,r,power=10){
  #assume 2 dimensions
  groups = length(pdfs)
  len = r - l
  n = 2^power*100
  P = seq(from = l, to = r, length.out = n)
  Z = w(pdfs,P)
  ans = vector(mode = 'numeric',2^power)
  obs = vector(mode = 'numeric')
  for (z in Z){
    phi = z / sum(z)
    location = sqrt(phi[[1]]^2 + (1-phi[[2]])^2)/sqrt(2)
    index = ceiling(location*2^power)
    b = (mean(z))/2
    ans[index] = ans[index] + b
    #obs = c(obs, rnorm(n = floor(b*100),mean = location, sd = 10/(2^power)))
  }
  #print(length(obs))
  #points(density(obs),type = 'l')
  ans = ans/(sum(ans)/(2^power))
  plot(ans,type = 'h')
  ans
}

kl <- function(a,b) {
  a = a+1e-20
  b = b+1e-20
  mean(a*(log(a)-log(b)) + b*(log(b)-log(a)))
}

distance <- function(a,b) {
  n = length(a)
  I = n/10
  ans = 1e20
  for (i in -I:I) {
    l = max(1,i)
    r = min(n,n+i)
    a.shifted = a[l:r]
    null = vector('numeric',n-(r-l))
    if (i > 0) a.shifted = c(null,a.shifted)
    else a.shifted = c(a.shifted,null)
    ans = min(ans, kl(a.shifted,b))
  }
  ans
}

difference <- function(C) {
  K = length(C)
  ans = 0
  for (i in 1:(K-1))
    for (j in (i+1):K)
      ans = max(ans,min(distance(C[[i]],C[[j]]),distance(C[[j]],C[[i]])))
  ans
}

identification <- function(D) {
  groups = length(D)
  Xs = list(groups)
  for (g in 1:groups){
    Xs[[g]] = D[[g]]$x
    s = sqrt(var(D[[g]]$x))
    m = mean(D[[g]]$x)
    D[[g]] = D[[g]][abs(D[[g]]$x-m) < 1.5*s,]
  }
  conditions = generate.coditions(Xs)
  cdes = list(groups)
  for (g in 1:groups) {
    for (condition in conditions) {
      cdes[[g]] = cde(D[[g]]$x,D[[g]]$y,x.margin = conditions,rescale = T,nymargin = 1000)
    }
  }
  cdes.length = length(cdes[[1]]$z[1,])
  power = 10
  # :(
  C = data.frame(sample = 1:2^power)
  C = C[-1]
  K = length(conditions)
  #plot(0,0,pch = 16, cex = 0,xlim = c(0,1),ylim = c(0,5))
  ent = vector(mode = 'numeric',length = K)
  for (c in 1:K) {
    pdfs = list(groups)
    l = 1e6
    r = -1e6
    for (g in 1:groups) {
      x = cdes[[g]]$y
      y = cdes[[g]]$z[c,]
      pdfs[[g]] = generate.pdf(x,y)
      l = min(l,min(x))
      r = max(r,max(x))
    }
    C[[c]] = generate.C(pdfs,l,r,power)
  }
  difference(C)
}

test.size = 100
set.seed(102)
test.answer = rbern(n = 100,prob = 0.5)
test.answer = as.logical(test.answer)
mean(test.answer)
D = vector("list",test.size)
for (t in 1:test.size) {
  d = generate.data(n=1000,anticausal = test.answer[t], vars = c(0.2,0.6))
  D[[t]] = d
}

for (t in 1:test.size) {
  grDevices::pdf(paste("synthetic/",as.character(t),".pdf",sep = ""))
  d = D[[t]]
  plot(d[[2]]$y~d[[2]]$x,cex = 0.3,pch = 16, col= "blue")
  points(d[[1]]$y~d[[1]]$x,cex = 0.3,pch = 16,col = "red")
  dev.off()
}

test.causal = vector(length = test.size)
test.anticausal = vector(length = test.size)
for (t in 1:test.size) {
  print(t)
  d = D[[t]]
  plot(d[[2]]$y~d[[2]]$x,pch = 16, cex = 0.4)
  points(d[[1]]$y~d[[1]]$x,pch = 16, cex = 0.4,col = 'red')
  if (test.answer[t]) print("anticausal")
  else print("causal")
  causal = identification(d)
  print(causal)
  anticausal = identification(reverse.data(d))
  print(anticausal)
  test.causal[t] = causal
  test.anticausal[t] = anticausal
  print("=======")
}

warnings()

hist(test.causal[!test.answer],breaks = 100,xlim = c(min(test.causal),max(test.causal)))
hist(test.causal[test.answer],breaks = 100,xlim = c(min(test.causal),max(test.causal)))

test.tresh = 6
test.classifier = test.causal > test.tresh
table(test.classifier, test.answer)
mean(test.classifier*test.answer + (!test.classifier)*(!test.answer))

test.directions = test.causal > test.anticausal
table(test.directions, test.answer)
mean(test.directions*test.answer + (!test.directions)*(!test.answer))
