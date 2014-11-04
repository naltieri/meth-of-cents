function makeTestMatrices(d,n)

	A = zeros(d,d,n);
	B = zeros(d,d,n);
	C = zeros(d,d,n);


	for i = 1:n
		a = rand(d,d);
		b = rand(d,d);
		c = rand(d,d);
		a = a + a' - 1;
		b = b + b' - 1;
		c = c + c' - 1;
		A[:,:,i] = a;
		B[:,:,i] = b;
		C[:,:,i] = c;

	end

	return(A,B,C)

end



function Mx(M,x)
	out = zeros(size(M,1),size(M,2))
	for i = 1:size(M,3)
		out = out + M[:,:,i]*x[i]
	end
	return(out)
end


function analyticCenter(x0,F, alpha, tol)
	# iter = 1;

	xPrev = zeros(size(x0,1),1);
	xPrev[1] = Inf;
	xPrev[2] = Inf;
	xPrev[3] = Inf; 
	x = [x0];



	Fx = zeros(size(F,1),size(F,2));

	for i = 1:size(F,3)
	    Fx = Fx + F[:,:,i]*x[i];
	end 

	if(cond(Fx) == Inf)
		println("Fx is singular...quitting")
		return
	end

	g = zeros(1,size(F,3)-1);
	for i = 2:size(F,3)
	    g[i-1] = -trace(pinv(Fx)*F[:,:,i])
	end
	g = g';

	H = zeros(size(F,3)-1,size(F,3)-1)
	for i = 2:size(F,3)
	    for j = 2:size(F,3)
	        H[i-1,j-1] = trace(pinv(Fx) *F[:,:,i]* pinv(Fx)* F[:,:,j]);
	    end
	end
	objHist = log(det(pinv(Fx)));
	gHist = norm(g,2);

	# Second condition is a hack. Happens when numerical underflow
	# occurs in updating x
	# Replacing 0 with 10^(-10) creates bad results. 

	while (norm(g,2) > tol)

		if (norm(xPrev - x,2) <= 0.0)
			break
		end

		delta = norm(sqrtm(pinv(H))*g,2);

		if delta <= 1/4
			alpha = 1;
		else
			alpha = 1/(1+delta)
		end
		xPrev = copy(x);
		x[2:end] = x[2:end] - alpha * pinv(H) *g;



		Fx = zeros(size(F,1),size(F,2))
		for i = 1:size(F,3)
		    Fx = Fx + F[:,:,i]*x[i];
		end

		chol(Fx)
		if(cond(Fx) == Inf)
			println("Fx is singular...quitting")
			return
		end

		L = chol(Fx);

		invLFL = zeros(size(F,1),size(F,2),size(F,3)-1)

		invL = inv(L);
		for i = 2:size(F,3)
			invLFL[:,:,i-1] = invL * F[:,:,i] *invL';
		end




		# #We will precompute the quantity Fx^(-1) * F_i = L^(-1)F_i L^T and store it.
		# #We will compute L^(-1) F_i by doing back substitution 
		# invLF = zeros(size(F,1),size(F,2),size(F,3)-1)
		# for i = 2:size(F,3)
		# 	invLF[:,j,i-1] = L\F[:,j,i]
		# end




		for i = 2:size(F,3)
		    g[i-1] = - trace(invLFL[:,:,i-1]);
		end
		


		for i = 2:size(F,3)
		    for j = 2:size(F,3)
		        H[i-1,j-1] = trace(invLFL[:,:,i-1] * invLFL[:,:,j-1]);
		    end
	    end
	    gHist = [gHist,norm(g,2)];
		objHist = [objHist, log(det(pinv(Fx)))];


	end

	xOpt = x;

	return(xOpt)
end



function methOfCents( A,B,C, lambda,x, theta)

		tol = 10.0^(-7)
		lambdaPrev = Inf;
		x = [1 ; x ]

		Ax = Mx(A,x)
		Bx = Mx(B,x)
		Cx = Mx(C,x)

	while lambdaPrev - lambda > tol
		lambdaPrev = copy(lambda);
		
		lambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
		println(lambda)

		#Concatenate lambda*B - A with C
	    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
	    F[1:size(A,1),1:size(A,1),:] = lambda*B-A
	    F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C

	    x = analyticCenter(x,F, .01,10.0^(-3))

	    #Recompute Ax,Bx,Cx
	    Ax = Mx(A,x);
		Bx = Mx(B,x);
		Cx = Mx(C,x);



	end
	

	return(x,lambda)
end




	
