graddesc_const_step <- function(f,g,x0,tbar,epsilon, Print=FALSE, maxiter=10^4){
  # Gradient method with constant stepsize. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # x0 - initial point
  # t  -  constant stepsize
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  
  
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  
  x=x0
  grad=g(x)
  print(grad)
  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[1,]<-x
  while (sum(grad^2)>epsilon^2 && iter<maxiter){
    iter=iter+1
    x=x-tbar*grad
    grad=g(x);
    trajectory[iter,]<-x
    if(Print){
      fun_val=f(x)  # we don't need this unless we are printing it out
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(fun_val,3), " norm_grad = ", signif(sqrt(sum(grad^2)),3)))
      print("x = ")
      print(x)
      print(grad)
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}



#############################################

graddesc_exact_step_quadratic <- function(A,b,x0,epsilon, maxiter=10^6, Print=FALSE){
  # Gradient method with exactstepsize for quadratic 
  # functions.
  # The method terminates when ||g(x)|| <epsilon
  # INPUT
  # A ....... the positive definite matrix associated with the objective function
  # b ....... a column vector associated with the linear part of the objective function   
  #  x0 - initial point
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  
  x=x0
  iter=1
  grad=2*(A%*%x+b)
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[1,]<-x
  while (sum(grad^2)>epsilon^2 && iter<maxiter){
    iter=iter+1
    tbar=c(sum(grad^2)/(2*t(grad)%*%A%*%grad)) # exact stepsize
    x=x-tbar*grad
    fun_val=t(x)%*%A%*%x+t(b)%*%x
   
    grad=2*(A%*%x+b)
    trajectory[iter,]<-x
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(fun_val,3), " norm_grad = ", signif(sqrt(sum(grad^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}

############################################

graddesc_backtracking <- function(f,g,x0,s, alpha,  beta, epsilon, Print=FALSE, maxiter=10^6){
  # Gradient method with constant stepsize. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # x0 - initial point
  # s ......... initial choice of stepsize
  # alpha ..... tolerance parameter for the stepsize selection
  # beta ...... the constant in which the stepsize is multiplied 
  #              at each backtracking step (0<beta<1)
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  
  stopifnot(beta>0, beta<1, alpha>0, alpha<1)
  
  x=x0
  grad=g(x)
  fun_val = f(x)
  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[1,]<-x
  while (sum(grad^2)>epsilon^2&&iter<maxiter){
    iter=iter+1
    tbar=s
    while (fun_val-f(x-tbar*grad) < alpha*tbar*sum(grad^2)){
      tbar=beta*tbar
    }
    x=x-tbar*grad
    fun_val=f(x)  # we don't need this unless we want to print it out
    grad=g(x);
    trajectory[iter,]<-c(x)
    
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
     
      print(paste("f(x) = ", signif(fun_val,3), " norm_grad = ", signif(sqrt(sum(grad^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}


#########################################################################


scaled_graddesc_const_step <- function(f,g,D, x0,tbar,epsilon, Print=TRUE, maxiter=10^6){
  # Gradient method with constant stepsize. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # x0 - initial point
  # t  -  constant stepsize
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  print("SCALED GRADIENT DESCENT - CONSTANT STEPSIZE")
  x=x0
  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[iter,] <- x0
  while (sum(g(x)^2)>epsilon^2 && iter < maxiter){
    iter=iter+1
    x=x - tbar*D(x)%*%g(x)
    trajectory[iter,] <- x
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(f(x),3), " norm_grad = ", signif(sqrt(sum(g(x)^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}


###################

scaled_graddesc_backtracking <- function(f,g,D, x0,s,alpha, beta, epsilon, Print=FALSE, maxiter=10^6){
  # Gradient method with constant stepsize. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # D - scaling matrix function
  # x0 - initial point
  # t  -  constant stepsize
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  print("SCALED GRADIENT DESCENT - backtracking")
  x=x0
  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[iter,] <- x0
  while (sum(g(x)^2)>epsilon^2&& iter <maxiter){
    iter=iter+1
    tbar=s
    while(f(x)-f(x-tbar*D(x)%*%g(x))<alpha *tbar *t(g(x))%*%D(x)%*%g(x)){
        tbar=beta*tbar
      }
    x=x-tbar*D(x)%*%g(x)
    trajectory=rbind(trajectory, c(x));
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(f(x),3), " norm_grad = ", signif(sqrt(sum(g(x)^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}







##########################################################

pure_Newton <- function(f,g,H, x0, epsilon, Print=TRUE, maxiter=10^4){
  # Pure Newton's method. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # H - a function to compute the Hessian matrix 
  # x0 - initial point
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  # iter = number of iterations
  
  x=x0

  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[1,]<-x
  while (sum(g(x)^2)>epsilon^2&& iter<maxiter){
    iter=iter+1
    x=x-solve(H(x), g(x))
    trajectory[iter,]<-x
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(f(x),3), " norm_grad = ", signif(sqrt(sum(g(x)^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}




##################################


damped_Newton <- function(f,g,H, x0, alpha, beta, epsilon, Print=TRUE, maxiter=10^6){
  
  # Pure Newton's method. The method
  # terminates when ||g(x)|| <epsilon
  # INPUT
  # f - the objective function
  # g - gradient of the objective function
  # H - a function to compute the Hessian matrix 
  # x0 - initial point
  # alpha
  # beta
  # epsilon ... tolerance parameter
  # Print =T/F - whether to print the output each step.
  # maxiter - the maximum number of iterations
  ####
  # OUTPUT
  # x.opt = optimal solution (up to a tolerance) of min f(x)
  # f.opt = optimal function value
  # iter = number of iterations
  
  x=x0
  
  iter=1
  trajectory <- matrix(0, nr=maxiter, nc=length(x0))
  trajectory[iter,] <- x0
  while (sum(g(x)^2)>epsilon^2 &&iter <maxiter){
    iter=iter+1
    d=solve(H(x), g(x))
    tbar=1
    while(f(x)-f(x-tbar*d)< alpha*tbar* t(d)%*%g(x)){
      tbar = beta*tbar
    }
    x=x-tbar*d    
    trajectory[iter,]<-x
    if(Print){
      print(paste('-------------- Iteration ', iter, ' -------------'))
      print(paste("f(x) = ", signif(f(x),3), " norm_grad = ", signif(sqrt(sum(g(x)^2)),3)))
    }
  }
  return(list(x.opt=x, trajectory=trajectory[1:iter,], iter=iter-1))
}


