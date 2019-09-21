# to compute a full Euclidean distance matrix for
dis = function(x)
{
  x = as.matrix(x)
  u = apply(x*x,1,sum) %*% matrix(1.0,1,nrow(x))
  sqrt(abs(u + t(u) - 2 * x %*% t(x)))
}

iorder = function(m)
{
  N = nrow(m) + 1
  iorder = rep(0,N)
  iorder[1] = m[N-1,1]
  iorder[2] = m[N-1,2]
  loc = 2
  for(i in seq(N-2,1))
  {
    for(j in seq(1,loc))
    {
      if(iorder[j] == i)
      {
        iorder[j] = m[i,1]
        if(j==loc)
        {
          loc = loc + 1
          iorder[loc] = m[i,2]
        } else
        {
          loc = loc + 1
          for(k in seq(loc, j+2)) iorder[k] = iorder[k-1]
          iorder[j+1] = m[i,2]
        }
      }
    }
  }
  -iorder
}

hc = function(d, method=c("single","complete","average"))
{
  if(!is.matrix(d)) d = as.matrix(d)
  # Pick a clustering function:
  method_fn = switch(match.arg(method),
                     single   = min,
                     complete = max,
                     average  = mean)
  N = nrow(d)
  diag(d)=Inf
  n = -(1:N)                       # Tracks group membership
  m = matrix(0,nrow=N-1, ncol=2)   # hclust merge output
  h = rep(0,N-1)                   # hclust height output
  for(j in seq(1,N-1))
  {
    # Find smallest distance and corresponding indices
    h[j] = min(d)
    
    i = which(d - h[j] == 0, arr.ind=TRUE)
    i = i[1,,drop=FALSE]
    p = n[i]
    p = p[order(p)]
    m[j,] = p
    grp = c(i, which(n %in% n[i[1,n[i]>0]]))
    n[grp] = j
    r = apply(d[i,],2,method_fn)

    d[min(i),] = d[,min(i)] = r
    d[min(i),min(i)]        = Inf
    d[max(i),] = d[,max(i)] = Inf
  }
  
  structure(list(merge = m, height = h, order = iorder(m),
                 labels = rownames(d), method = method, 
                 call = match.call(), dist.method = "euclidean"), 
            class = "hclust")
}

# Comparing ours and original hclust on single 
actual_hclust_singl = hclust(dist(USArrests),method="single")
defined_hclust_singl = hc(dis(USArrests), method="single")
plot(actual_hclust_singl, main = 'In-built function Single Linkage')
plot(defined_hclust_singl, main = 'Our defined function Single Linkage') # our defined function output

actual_hclust_avg = hclust(dist(USArrests),method="average")
plot(actual_hclust_avg, main = 'In-built function Average Linkage')
defined_hclust_avg = hc(dis(USArrests), method="average")
plot(defined_hclust_avg, main = 'Our defined function Average Linkage')

actual_hclust_compl = hclust(dist(USArrests),method="average")
plot(actual_hclust_compl)
defined_hclust_compl = hc(dis(USArrests), method="complete")
plot(actual_hclust_compl)