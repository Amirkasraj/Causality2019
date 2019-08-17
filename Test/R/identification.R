generate.coditions <- function(Xs, k = 7) {
  # ToDo interval should be small ?! use hdr
  m = mean(Xs[[1]])
  s = sqrt(var(Xs[[1]]))
  seq(from = m-1*s, to = m+1*s, length.out = k)
}

average <-function(functions) {
  function(x) {
    ans = 0
    for (f in functions)
      ans = ans + f(x)
    ans
  }
}

generate.cdf <- function(xs,ys) {
  l=min(xs)
  r=max(xs)
  d = splinefun(xs,ys)
  i = integrate(d,lower = l ,upper = r)$value
  single <- function(x) {
    y = 0
    if (x<l) y = 0
    else if (x>r) y = 1
    else y = integrate(d,lower = l, upper = x)$value
    y/i
  }
  function(x) {
    vapply(FUN = single, X=x, FUN.VALUE = 1)
  }
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

generate.sample <- function(distribution,n=1) {
  u = runif(n,min = 0,max = 1)
  root <- function(a) {
    uniroot(function(x){distribution(x) - a},interval = c(-1000,1000))$root
  }
  vapply(FUN = root, X = u, FUN.VALUE = 1)
}

phi <-function(pdfs,x) {
  single <- function(x) {
    ans = NULL
    for (d in pdfs)
      ans = c(ans, d(x))
    ans
  }
  vapply(FUN = single,X = x,FUN.VALUE = double(length = length(pdfs)))
}

generate.C <- function(pdfs, cdfs, n=1000, ref=average(cdfs)){
  groups = length(pdfs)
  points = generate.sample(ref,n)
  samples = phi(cdfs,points)
  s = 0
  for (g in 1:groups)
    s = s + samples[g,]
  for (g in 1:groups)
    samples[g,] = samples[g,]/s
  #assume groups = 2
  sqrt((1-samples[1,])^2 + samples[2,]^2)/sqrt(2)
}

bootstrap <- function(x,iter = 10*length(x),gridsize=1000) {
  xs = NULL
  for (b in 1:iter) {
    smpl = sample(x,length(x),replace = T)
    l = 0.1
    r = 0.9
    dens = bkde(smpl,range.x = c(l,r),gridsize = gridsize, truncate = T,kernel = 'box')
    ys = dens$y
    xs = dens$x
    ys[ys<0.00000001] = 0
    ys = dens$y
    ys = ys/mean(ys)
    avg = avg + ys
  }
  avg/iter
}

distance <- function(a,b) {
  mean(a*(log(a)-log(b)) + b*(log(b)-log(a)))
}

difference <- function(h) {
  groups = length(h)
  ans = 0
  for (i in 1:groups)
    for (j in 1:groups)
      ans = max(ans,distance(h[[i]],h[[j]]))
  ans
}

library(KernSmooth)
library(hdrcde)

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
  n=300
  # :(
  Cs = data.frame(sample = 1:n)
  Cs = Cs[-1]
  K = length(conditions)
  for (c in 1:K) {
    pdfs = list(groups)
    cdfs = list(groups)
    for (g in 1:groups) {
      pdfs[[g]] = generate.pdf(cdes[[g]]$y, cdes[[g]]$z[c,])
      cdfs[[g]] = generate.cdf(cdes[[g]]$y, cdes[[g]]$z[c,])
    }
    C =  generate.C(pdfs,cdfs,n)
    Cs[[c]] = C
  }
  avg = 0
  # :(
  h = data.frame(sample = 1:1000)
  h = h[-1]
  for (c in 1:K) {
    x = as.vector(Cs[[c]])
    h[[c]] = bootstrap(x)
    avg = avg + h[[c]]
  }
  difference(h)
}

library(Rlab)

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

grDevices::jpeg("grid_of_datasets.jpg",width = 3000,height = 3000)
par(mfrow=c(10,10))
par(mar = c(1,1,1,1))

for (t in 1:test.size) {
  d = D[[t]]
  plot(d[[2]]$y~d[[2]]$x,cex = 0.5,pch = 16, col= "black",ylab = 'Y',xlab = 'X')
  points(d[[1]]$y~d[[1]]$x,cex = 0.4,pch = 16,col = "red")
}

dev.off()

for (t in 1:test.size) {
  grDevices::pdf(paste("synthetic/",as.character(t),".pdf",sep = ""))
  d = D[[t]]
  plot(d[[2]]$y~d[[2]]$x,cex = 0.4,pch = 16, col= "black")
  points(d[[1]]$y~d[[1]]$x,cex = 0.4,pch = 16,col = "red")
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

hist(test.causal[!test.answer],breaks = 100,xlim = c(min(test.causal),max(test.causal)))
hist(test.causal[test.answer],breaks = 100,xlim = c(min(test.causal),max(test.causal)))

test.tresh = 6
test.classifier = test.causal > test.tresh
table(test.classifier, test.answer)
mean(test.classifier*test.answer + (!test.classifier)*(!test.answer))

test.directions = test.causal > test.anticausal
table(test.directions, test.answer)
mean(test.directions*test.answer + (!test.directions)*(!test.answer))
