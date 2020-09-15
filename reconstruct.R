#gene X sample
SxC = as.matrix(Est.prop.GSE50244$Est.prop.weighted)
SxC = SxC[,-4]
CxS = t(SxC)
dim(CxS)
A = CxS%*%t(CxS) #+ diag(rep(1e-10, 11)) if not invertible
B = solve(A)
GxS = as.matrix(matrixNN)
dim(GxS)
GxC = GxS%*%t(CxS)%*%B
dim(GxC)
write.csv(GxC, 'GeneByCell.csv')
a<-1e4
b<-84
c<-20

A<-matrix(rnorm(a*b,0,1),b,a)
C<-matrix(rnorm(c*b,0,1),b,c)
X<-matrix(rnorm(a*c,0,1),c,a)


A = t(GxS)
C = t(CxS)
X<- t(matrix(rnorm(62449*11,0,1),62449,11))

lambda<-0.01 # adjust lambda according to your "C%*%X - A"

res<- vector()
for (i in 1:1000) {
  X<- X - lambda /sqrt(i) * t(C)%*%( C%*%X - A)
  X <- (X>0)*X
  if(i%/%50==0){
    tmp<-sum( (C%*%X - A)^2 )
    res<-c(res,tmp)
  }
  tmp<-length(res)
  if( abs(res[tmp] - res[tmp-1])<0.1 && tmp>=2    ){ ## convergence
    break
  }
}

plot(res) # check convergence

