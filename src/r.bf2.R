library(MASS); 
logBF2=function(y,W, Q,D,lam){
   #W contains a column of 1 and has full rank.
   #D is a vector of eigenvalues;
#    y=mph[,1]; 
#    W=cbind(rep(1,ni)); 
#    Q=cbind(ek$vectors[,1]); 
#    D=ek$values[1]; 

    n=length(y)
    y=cbind(y);
    X = cbind(Q, W)
    ytX = t(y) %*% X;
    yty = sum(y*y);
    ytW = t(y) %*% W;
    wtw = t(W) %*% W;
    var0=ytW %*% ginv(wtw) %*% t(ytW)
    #wtw have low dimension.
    lD=lam*D;
    if(length(lD) > 1) {
	mF=diag(lD/(1+lD));
    } else {
	mF=lD/(1+lD);
    }
    QtW = t(Q) %*% W;    #QtW;
    M = ginv(wtw - t(QtW) %*% mF %*% QtW);
    fQw = mF%*%QtW;
    a11=mF+ fQw %*% M %*% t(fQw);
    a12=-fQw %*% M;
    a21=-M %*% t(fQw);
    a22=M;
    aa=rbind(cbind(a11,a12), cbind(a21,a22));
    bhat =  aa %*% t(ytX); 
    yhat = X %*% bhat; 

    res1= 0.5*log10(det(wtw))-0.5*sum(log10(1+lD))+0.5*log10(det(M));
    f1=log10(yty - ytX %*% bhat);
    f2=log10(yty - var0)
    bf = res1 - (n/2)*(f1-f2); 
    h2=var(yhat)/var(y); 


#    xtx = t(X) %*% X; 
#    bxtx = xtx+diag(c(1/lD, 100000))
#    bb=ginv(bxtx);
##    res2= 0.5*log10(det(wtw))-0.5*det(bxtx)-0.5*sum(log(lD)); 
#    res2= 0.5*log10(det(wtw))-0.5*sum(log10(1+lD))+0.5*log10(det(M));
#    f1=log10(yty - ytX %*% bb %*% t(ytX));
#    res2 = res2 - (n/2)*(f1-f2); 

    return(c(bf,h2));
}


library(nloptr); 
myopt=function(ph, W, vec, val) {
    oz = ph; 
    fn=function(x) {
	return(-logBF2(oz, W, vec, val, x)[1]);
    }
    fit=neldermead(0.5, fn, lower=0, upper=10)
#    h=fit$par / (1+fit$par); 
#    bf=-fit$value;
    return(logBF2(oz, W, vec, val, fit$par));
}


