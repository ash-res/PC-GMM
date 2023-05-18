PCA_GMM <- function(bx,sx,by,sy,ld,cor.x,nx,ny,r=NULL){
p <- nrow(bx); K <- ncol(bx)

# estimate principal components
Phi <- ((rowSums(abs(bx))/sy)%*%t(rowSums(abs(bx))/sy))*ld
if(missing(r)){r <- which(cumsum(prcomp(Phi,scale=FALSE)$sdev^2/sum((prcomp(Phi,scale=FALSE)$sdev^2)))>0.999)[1]
} else {r <- r}
lambda <- sqrt(p)*prcomp(Phi,scale=FALSE)$rotation[,1:r]
evec <- eigen((t(lambda)%*%lambda))$vectors
eval <- eigen((t(lambda)%*%lambda))$values
lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
dim(lambda) <- c(p,r)

# outcome quantities of interest
ay <- 1/((ny*sy^2)+by^2)
Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
By <- ay*by
Ay.f <- t(lambda)%*%Ay%*%lambda; By.f <- as.vector(t(lambda)%*%By)
SigY <- solve(Ay.f)*(1-as.numeric(t(By.f)%*%solve(Ay.f)%*%By.f)); SigY <- SigY*(1/(ny-r+1))
gamY_est <- as.vector(solve(Ay.f)%*%By.f)

# exposures quantities of interest
ax <- matrix(NA,nrow=p,ncol=K); Ax <- list()
for (k in 1:K){ax[,k] <- 1/((nx[k]*sx[,k]^2)+bx[,k]^2)}
for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
Ax.f <- list()
for (k in 1:K){Ax.f[[k]] <- t(lambda)%*%Ax[[k]]%*%lambda}
Bx.f <- t(lambda)%*%Bx
sqrt.Ax.f <- function(k){
  evec <- eigen(Ax.f[[k]])$vectors; eval <- eigen(Ax.f[[k]])$values
  return((evec%*%diag(sqrt(eval))%*%t(evec)))
}
sqrt.Ax.f <- lapply(1:K,sqrt.Ax.f)

SigX <- function(k,l){
  solve(sqrt.Ax.f[[k]]%*%t(sqrt.Ax.f[[l]]))*(cor.x[k,l]-(as.numeric(t(Bx.f[,k])%*%solve(sqrt.Ax.f[[k]]%*%t(sqrt.Ax.f[[l]]))%*%Bx.f[,l])))*(1/(sqrt(nx[k]-r+1)*sqrt(nx[l]-r+1)))
}

SigX2 <- function(m){
  SigX2a <- list()
  for (m1 in 1:K){SigX2a[[m1]] <- SigX(m1,m)}
  SigX2a <- do.call(rbind, SigX2a)
  return(SigX2a)
}

SigX3 <- list()
for (m1 in 1:K){SigX3[[m1]] <- SigX2(m1)}
SigX3 <- do.call(cbind, SigX3)
# check: SigX3[(2*p+1):(3*p),(1*p+1):(2*p)] == SigX(3,2)
SigX <- SigX3; rm(SigX2, SigX3)
gamX_est <- function(k){as.vector(solve(Ax.f[[k]])%*%Bx.f[,k])}
gamX_est <- sapply(1:K,gamX_est)

# conditional F-test
if(K>2){
condF <- function(j){
SigXX <- SigX[c(((((j-1)*r)+1):(j*r)),((1:(r*K))[-((((j-1)*r)+1):(j*r))])),c(((((j-1)*r)+1):(j*r)),((1:(r*K))[-((((j-1)*r)+1):(j*r))]))]
g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]%*%tet))}
Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
Om.nr <- function(tet){as.matrix(cbind(diag(r),kronecker(t(-tet),diag(r)))%*%SigXX%*%t(cbind(diag(r),kronecker(t(-tet),diag(r)))))}
Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
G <- -gamX_est[,-j]
DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(r-K+1)
return(condF)
}
condF <- sapply(1:K,condF)
}

if(K==2){
  condF <- function(j){
    SigXX <- SigX[c(((((j-1)*r)+1):(j*r)),((1:(r*K))[-((((j-1)*r)+1):(j*r))])),c(((((j-1)*r)+1):(j*r)),((1:(r*K))[-((((j-1)*r)+1):(j*r))]))]
    g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]*tet))}
    Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
    tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
    Om.nr <- function(tet){as.matrix(cbind(diag(r),kronecker(t(-tet),diag(r)))%*%SigXX%*%t(cbind(diag(r),kronecker(t(-tet),diag(r)))))}
    Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
    G <- -gamX_est[,-j]
    DQ.nr <- function(tet){as.numeric(2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet)))}
    condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(r-K+1)
    return(condF)
  }
  condF <- sapply(1:K,condF)
}

# non-robust LIML estimate based on an incorrect weighting matrix (estimate should be consistent)
g <- function(tet){as.vector(gamY_est - (gamX_est%*%tet))}
Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
tet.gg <- nlminb(rep(0,K),objective=Q.gg)$par
Om.nr <- function(tet){as.matrix(SigY+(kronecker(t(tet),diag(r))%*%SigX%*%t(kronecker(t(tet),diag(r)))))}
Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
G <- -gamX_est
DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
liml.nr <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$par
var.liml.nr <- as.matrix(solve(t(G)%*%solve(Om.nr(liml.nr))%*%G))
Q.nr <- Q.nr(liml.nr)

# estimating the overdispersion variance parameter
Om <- function(tet,kappa){as.matrix(SigY+(diag(r)*kappa/ny)+(kronecker(t(tet),diag(r))%*%SigX%*%t(kronecker(t(tet),diag(r)))))}
Q.kap <- function(kappa){as.numeric(t(g(liml.nr))%*%solve(Om(liml.nr,kappa))%*%g(liml.nr))-(r-K)}
kap.sq <- seq(-50,50,0.1)
Q.sq <- vector(,length=length(kap.sq))
for (q in 1:length(kap.sq)){Q.sq[q] <- Q.kap(kap.sq[q])}
if(sum(Q.sq>0)==0){kappa.est <- 0; default <- 1}
if(sum(Q.sq>0)>0){kappa.est <- uniroot(Q.kap,c(kap.sq[max(which(Q.sq>0))],20), extendInt="yes")$root; default <- 0}

# robust LIML estimate based on a plug-in overdispersion parameter estimate
Om.r <- function(tet){as.matrix(SigY+(diag(r)*max(0,kappa.est)/ny)+(kronecker(t(tet),diag(r))%*%SigX%*%t(kronecker(t(tet),diag(r)))))}
Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om.r(tet))%*%g(tet))}
DQ <- function(tet){2*as.matrix(t(G)%*%solve(Om.r(tet))%*%g(tet))}
liml <- nlminb(liml.nr,objective=Q,gradient=DQ)$par
var.liml <- as.matrix(solve(t(G)%*%solve(Om.r(liml))%*%G))

res.list <- list("liml"=liml, "se.liml"=sqrt(diag(var.liml)), "kappa"=kappa.est, "default"=default, "factors"=r, "Q"=Q(liml), "liml.nr"=liml.nr, "se.liml.nr"=sqrt(diag(var.liml.nr)), "Q.nr"=Q.nr, "condF"=condF)
return(res.list)
}