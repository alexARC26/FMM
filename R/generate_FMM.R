# Simulating data from FMM models
#
# Arguments:
#    M: value of the intercept parameter M.
#    A: value of the FMM wave amplitude parameter A.
#       To simulate a FMM model with m components vector of length m.
#    a: value of the FMM wave phase translation parameter alpha.
#       To simulate a FMM model with m components vector of length m.
#    b: value of the FMM wave phase translation parameter beta.
#       To simulate a FMM model with m components vector of length m.
#    w: value of the FMM wave phase translation parameter omega
#       To simulate a FMM model with m components vector of length m.
#    from: initial time point. Default 0.
#    to: final time point. Default 2*pi.
#    length.out: desired length of the simulation. Default 100.
#    timePoints: A numeric vector containing the time points at which each data of one single period will be simulated.
#                By default it is a sequence of equally spaced values from 'from' to 'to' of length 'length.out' arguments.
#    plot: TRUE when the simulated data are plotted.
#    outvalues: TRUE when the simulated data are returned.
#    sigmaNoise: Standard deviation of the gaussian noise to be added.
#                Its default value is 0 equivalent to a simulation set-up without noise.
generate_FMM <- function(M,A,a,b,w,from=0,to=2*pi,length.out=100,timePoints=seq(from,to,length=length.out),
                plot=TRUE,outvalues=TRUE,sigmaNoise=0){

  narg <- max(length(M),length(A),length(a),length(b),length(w))

  if(length(M)>1){
    warning("M parameter should be a vector of length 1.
            The intercept parameter used in the simulation is the sum of the elements of the argument M.")
    M <- sum(M)
  }
  M <- rep(M/narg,length.out=narg)


  A <- rep(A,length.out=narg)
  if(sum(A <= 0) > 0) stop("A parameter must be positive.")

  a <- rep(a,length.out=narg)
  a <- a%%(2*pi) # between 0 and 2*pi

  b <- rep(b,length.out=narg)
  b <- b%%(2*pi) # between 0 and 2*pi

  w <- rep(w,length.out=narg)
  if(sum(w<0)>0 | sum(w>1)>0) stop("w parameter must be between 0 and 1.")

  t <- timePoints

  phi <- list()
  for(i in 1:narg){
    phi[[i]] <- b[i]+2*atan(w[i]*tan((t-a[i])/2))
  }

  ym <- list()
  for(i in 1:narg){
    ym[[i]] <- M[i]+A[i]*cos(phi[[i]])
  }

  y <- rep(0,length(t))
  for(i in 1:narg){
    y <- y + ym[[i]]
  }

  if (sigmaNoise > 0) y <- y + rnorm(length.out,0,sigmaNoise)

  if(plot) {
    type_<-ifelse(sigmaNoise==0,"l","p")
    plot(t,y,type=type_,lwd=2,col=2,xlab="Time",ylab="Response",
                main=paste("Simulated data from FMM model"))

  }


  if(outvalues) return(list(input=list(M = M[1]*narg, A=A,alpha=a,beta=b,omega=w),t=t,y=y))
}
