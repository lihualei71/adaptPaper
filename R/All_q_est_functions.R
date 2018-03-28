## The code is from Rina Barber's website
##   https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R

# this file provides functions for computing q-hat for the SABHA method
# includes:
# Solve_q_step = q-hat of the form (eps,...,eps,1,...,1)
# Solve_q_ordered = q-hat of the form (q1,...,qn) with eps<=q1<=q2<=...<=qn<=1
# Solve_q_block = q-hat that is block-wise constant, and eps<=q<=1
# Solve_q_TV_1dim = q-hat that has bounded 1d total variation norm, and eps<=q<=1
# Solve_q_TV = q-hat that has bounded total variation norm on a specified graph, and eps<=q<=1


##########################################################
## SABHA with ordered q: q = step function
##########################################################

# Solve_q_step: returns a vector q subject to the constraint that (q1, q2, ...,  qn) = (eps, eps, ..., eps, 1, 1, ..., 1), with as many eps's as possible, subject to sum_i 1{P[i]>tau}/q[i](1-tau) <= n
Solve_q_step = function(Pvals, tau, eps){
  n = length(Pvals)
  sum_p_over_threshold = sum(Pvals > tau)
  K = max(which(cumsum(Pvals > tau) <= (n*(1-tau) - sum_p_over_threshold) / (1/eps - 1)), 0)
  
  return (c(rep(eps, K), rep(1, n-K)))
}

##########################################################
## SABHA with ordered q: q satisfies eps<=q1<=q2<=...<=qn=1
##########################################################

Solve_q_ordered = function(Pvals, tau, eps, ADMM_params){
	PAVA_alg = create_PAVA_alg_function()
	M = diag(length(Pvals))
	q = Solve_q_ADMM(Pvals, tau, eps, M, PAVA_alg, ADMM_params)
	q
}

create_PAVA_alg_function = function(){
	function(y){
         # solving: min{1/2 ||x-y||^2_2 : x_1<=..<=x_n}
	       # PAVA algorithm (Barlow et al 1972)
	       n=length(y)
	       thresh = 1e-8
	       groups = 1:n
	       block = 1
	
	       stop = FALSE
	       while(!stop){
		         if(any(groups==block+1)){
			             block_plus = which(groups==block+1)
			             if(mean(y[which(groups==block)])<=mean(y[which(groups==block+1)])+thresh){
				              block = block+1
			             } else{
				              groups[which(groups>block)] = groups[which(groups>block)] - 1
				              stop_inner = FALSE
				              while(!stop_inner)
					               if(any(groups==block-1)){
						               if(mean(y[which(groups==block-1)])>mean(y[which(groups==block)])+thresh){
							                groups[which(groups>=block)] = groups[which(groups>=block)] - 1
							                block = block-1
						               } else{
							                stop_inner=TRUE
						               }
					               } else{
						                stop_inner=TRUE
					               }
				           }
			       } else{
				        stop=TRUE
			       }
	       }
	       x=y
	       for(i in 1:max(groups)){
		         x[which(groups==i)]=mean(y[which(groups==i)])
	       }
	       x
     }
     
}    
     
     
##########################################################
## SABHA with block q: q = constant over blocks
##########################################################
     
Solve_q_block = function(Pvals, tau, eps, blocks, ADMM_params){
	# blocks[i] gives the index of the block to which Pvals[i] belongs
	block_proj = create_block_function(blocks)
	q_init = block_proj((Pvals>tau)/(1-tau))
	if(min(q_init)>=eps & max(q_init)<=1){
		q = q_init
	}else{
		M = diag(length(Pvals))
		q = Solve_q_ADMM(Pvals, tau, eps, M, block_proj, ADMM_params)
	}
	q
}  
     
create_block_function = function(blocks){
	function(y){
         # solving: min{1/2 ||x-y||^2_2 : x is constant over the blocks}
		x = y
		block_inds = sort(unique(blocks))
		for(i in block_inds){
			x[which(blocks==block_inds[i])] = mean(x[which(blocks==block_inds[i])])
		}
		x
	}
}     
     
##########################################################
## SABHA with TV norm constraint on q: ||q||_TV <= TV_bd
##########################################################

Solve_q_TV_1dim = function(Pvals, tau, eps, TV_bd, ADMM_params){
	edges = cbind(1:(length(Pvals)-1),2:length(Pvals))
	Solve_q_TV(Pvals, tau, eps, edges, TV_bd, ADMM_params)
}

Solve_q_TV_2dim = function(Pvals, tau, eps, TV_bd, ADMM_params){
	n1 = dim(Pvals)[1]
	n2 = dim(Pvals)[2]
	edges = NULL
	get_ind = function(i,j){i+(j-1)*n1}
	# horizontal edges
	for(i in 1:n1){for(j in 1:(n2-1)){edges=rbind(edges,c(get_ind(i,j),get_ind(i,j+1)))}}
	# vertical edges
	for(j in 1:n2){for(i in 1:(n1-1)){edges=rbind(edges,c(get_ind(i,j),get_ind(i+1,j)))}}
	Solve_q_TV(c(Pvals), tau, eps, edges, TV_bd, ADMM_params)
}
     
Solve_q_TV = function(Pvals, tau, eps, edges, TV_bd, ADMM_params){
	# edges is a e-by-2 matrix giving the edges of the adjacency graph
	# edges[i,1:2] gives the indices of the nodes on the i-th edge
	# constraint: sum_{i=1,..,e} |q[edges[i,1]] - q[edges[i,2]]| <= TV_bd
	L1_proj = create_L1_function(TV_bd)
	nedge = dim(edges)[1]; n = length(Pvals)
	M = matrix(0,nedge,n); for(i in 1:nedge){M[i,edges[i,1]]=1; M[i,edges[i,2]]=-1}
	q = Solve_q_ADMM(Pvals, tau, eps, M, L1_proj, ADMM_params)
	q
}  
    
create_L1_function = function(bound){
	function(y){
         # solving: min{1/2 ||x-y||^2_2 : ||x||_1 <= bound}
        if(sum(abs(y))<=bound){x=y} else{
			    mu = sort(abs(y), decreasing = TRUE)
    	    xi = max(which(mu - (cumsum(mu)-bound)/(1:length(mu))>0))
        	theta = (sum(mu[1:xi])-bound)/xi
	        tmp = abs(y)-theta
    	    x = rep(0, length(tmp))
        	x[which(tmp>0)] = tmp[which(tmp>0)]
	        x[which(tmp<=0)] = 0
    	    x = x*sign(y)
    	  }
        x
	}
}        
     
     
##########################################################
## SABHA ADMM algorithm
##########################################################

Solve_q_ADMM = function(Pvals, tau, eps, M, projection, ADMM_params){
# min -sum_i (B[i]*log((1-tau) q[i]) + (1-B[i])*log(1-(1-tau) q[i]))
# subject to (1) q \in Qset (characterized by M*q \in Mset)
# and (2) sum_i B[i]/q[i] <= gamma and (3) eps<=q<=1
# introduce auxiliary variables x, y under the constraint Mq = x, q = y
# ADMM optimization:
# minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
# where qt is the previous iteration's q value
  
# ADMM_params are: alpha, beta, eta, max_iters, converge_thr
	alpha_ADMM = ADMM_params[1]
	beta = ADMM_params[2]
	eta = ADMM_params[3]
	max_iters = ADMM_params[4]
	converge_thr = ADMM_params[5]

	n = length(Pvals)
	B = (Pvals > tau) 
  gamma = n*(1-tau) # bound on sum_i (Pvals[i]>tau) / q[i]*(1-tau)
	q = y = rep(1,n)
	v = rep(0,n)
	u = x = rep(0,dim(M)[1])
	
	converge_check = function(q,x,y,u,v,q_old,x_old,y_old,u_old,v_old){
		max(c(sqrt(sum((q-q_old)^2))/sqrt(1+sum(q_old^2)),
          sqrt(sum((x-x_old)^2))/sqrt(1+sum(x_old^2)),
          sqrt(sum((y-y_old)^2))/sqrt(1+sum(y_old^2)),
          sqrt(sum((u-u_old)^2))/sqrt(1+sum(u_old^2)),
          sqrt(sum((v-v_old)^2))/sqrt(1+sum(v_old^2))))
	}
	
	stop = FALSE
	iter = 0
  while(!stop){
    iter = iter+1
    q_old = q; x_old = x; y_old = y; u_old = u; v_old = v
    q = q_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    x = x_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta, projection)
    y = y_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    u = u_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
	  v = v_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
	  if(converge_check(q,x,y,u,v,q_old,x_old,y_old,u_old,v_old)<=converge_thr){stop=TRUE}
	  if(iter>=max_iters){stop=TRUE}
  }
    
  return(q)
    
}


# inverse_sum_prox solves: min{1/2 ||x-y||^2 : x_i>0, sum_i 1/x_i <= bound}
# Used in y-update step of ADMM
inverse_sum_prox = function(y,bound){

	y = pmax(0,y) # the solution will have all positive x_i's now
					# and we can now ignore the constraint x_i>0
	
	if(sum(1/y)<= bound){
		x=y
	}else{ # use Lagrange multipliers
		
		# we should have - lambda * d/dx_j (sum_i 1/x_i) = d/dx_j (1/2 ||x-y||^2)
		# for all j, for some single lambda>0
		# in other words, lambda / x^2 = x-y (this holds elementwise)
		# rearranging, lambda = x^3 - x^2*y
		# let c = log(lambda) so that it's real-valued
		# we need to solve x^3 - x^2*y - exp(c) = 0 (elementwise)
		
		cuberoot = function(c){ # this solves the cubic equation x^3-x^2*y-exp(c)=0
			temp1 = ((y/3)^3 + exp(c)/2 + (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
			temp2 = ((y/3)^3 + exp(c)/2 - (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
			x = sign(temp1)*abs(temp1)^(1/3) + sign(temp2)*abs(temp2)^(1/3) + (y/3)
			x
		}
		
		# now we need to choose c, i.e. choose the lagrange multiplier lambda=exp(c)
		# the right value of c is the one that produces an x satisfying sum_i 1/x_i = bound
		
		c = uniroot(function(c){sum(1/cuberoot(c))-bound},c(-100,100))$root
		x = cuberoot(c)
	}
	x
}

q_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
# minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
# where qt is the previous iteration's q value
# equivalently, -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + (alpha eta + beta)/2 * ||q-w||_2^2
# where w = - (M'(ut + alpha (M qt - xt)) + (vt - beta yt - alpha eta qt))/(alpha eta + beta)
  
	w = - ( t(M)%*%(u + alpha*(M%*%q - x)) + (v - beta*y - alpha*eta*q) )/(alpha*eta + beta)
	
	q[B==1] = (w[which(B==1)]+sqrt(w[which(B==1)]^2+4/(alpha*eta + beta)))/2
	q[B==0] = ((w[which(B==0)]+1/(1-tau))-sqrt((w[which(B==0)]-1/(1-tau))^2+4/(alpha*eta+beta)))/2
	q[q<eps] = eps
	q[q>1] = 1
	q
}

x_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta, projection){
	# Proj_Mset (M q + u/alpha)
	x = projection(M%*%q + u/alpha) 
}

y_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
	# Prof_B (q + v/beta)
	# where B = {sum_i B[i]/y[i]<= gamma}
	y = q + v/beta
	y[which(B==1)] = inverse_sum_prox((q+v/beta)[which(B==1)], gamma)
	y
}

u_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
	u = u + alpha * (M%*%q -x)
	u
}

v_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
	v = v + beta * (q-y)
	v
}







