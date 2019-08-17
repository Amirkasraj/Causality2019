# f : [-0.5,0.5] -> [0,inf]
generate.f <- function(complexity = 5) {
  l = -0.5
  r = 0.5
  x = seq(from = l, to = r,length.out = complexity)
  normal = rnorm(complexity)
  pos = exp(normal)
  f = splinefun(pos~x)
  u = seq(from = l, to=r, length.out = 1000)
  shift = -min(f(u)) + 0.5
  scale = 1/(integrate(f,lower = l,upper = r)$value+shift*(r-l))
  function(x) {
    ans = vapply(FUN = f,X = x,FUN.VALUE = 1)
    (ans + shift)*scale
  }
}

indefinite_integral = function(f) {
  my_integral <- function(x) {
    integrate(f,lower = -0.5,upper = x)$value
  }
  function(x) {
    vapply(FUN = my_integral,X = x,FUN.VALUE = 1)
  }
}

library(sigmoid)

projection <- function(x) {
  #x/(1+abs(x)) * (0.5/1)
  sigmoid(x)-0.5
}

generate.F <- function(...) {
  f = generate.f(...)
  F = indefinite_integral(f)
  #curve(F,xlim = c(-0.5,0.5))
  function(x) {
    F(projection(x))
  }
}

generate.weights <- function(x,l) {
  l = abs(l-x)
  l = 1/(l^2+0.00000001)
  l/sum(l)
}

combination <- function(functions, weights) {
  n = length(functions)
  function(e) {
    ans = 0
    for (i in 1:n)
      ans = ans + functions[[i]](e)*weights[i]
    ans
  }
}

generate.mechanism <- function(stations.x) {
  n = length(stations.x)
  functions = vector()
  for (i in 1:n)
    functions = c(functions,generate.F())
  #for (i in 1:length(stations.x)) {
  #  F  = stations[[i]]
  #  curve(F,add = (i!=1),xlim = c(-3,3))
  #}
  #print(stations.x)
  function(x,e) {
    weights = generate.weights(x,stations.x)
    combination(functions,weights)(e)
  }
}

generate.data <- function(n=10000, groups=2, anticausal = FALSE, vars = 1:groups) {
  f = generate.mechanism(seq(from = -1, to = 1, length.out = 3))
  res = list(groups)
  for (g in 1:groups) {
    X = rnorm(n,sd = 1)
    E = rnorm(n,sd = vars[g],mean = 0)#runif(1,05,2))
    #print(sqrt(var(E)))
    Y = mapply(FUN = f, x = X, e = E)
    res[[g]] = data.frame(x = X,y = Y)
  }
  if (anticausal)
    res = reverse.data(res)
  res
}

reverse.data <- function(D) {
  groups = length(D)
  res = D
  for (g in 1:groups) {
    names(res[[g]])[names(res[[g]]) == "x"] = "z"
    names(res[[g]])[names(res[[g]]) == "y"] = "x"
    names(res[[g]])[names(res[[g]]) == "z"] = "y"
  }
  res
}

library(plotly)

f = generate.mechanism(seq(from = -10, to = 10, length.out = 5))
m = matrix(nrow = 100,ncol = 100)
for (i in 1:100) {
  for (j in 1:100) {
    x = i/100*12 -6
    e = j/100*20 -10
    m[i,j] = f(x,e)
  }
}

plot_ly(z= m)%>% add_surface()
