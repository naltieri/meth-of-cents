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

	while norm(g,2) > tol
		x[2:end] = x[2:end] - alpha * pinv(H) *g;

		Fx = zeros(size(F,1),size(F,2))
		for i = 1:size(F,3)
		    Fx = Fx + F[:,:,i]*x[i];
		end

		if(cond(Fx) == Inf)
			println("Fx is singular...quitting")
			return
		end


		for i = 2:size(F,3)
		    g[i-1] = - trace(pinv(Fx)*F[:,:,i]);
		end
		


		for i = 2:size(F,3)
		    for j = 2:size(F,3)
		        H[i-1,j-1] = trace(pinv(Fx) *F[:,:,i] * pinv(Fx)* F[:,:,j]);
		    end
	    end
	    gHist = [gHist,norm(g,2)];
		objHist = [objHist, log(det(pinv(Fx)))];


	end

	xOpt = x;

	return(xOpt)
end



function methOfCents( A,B,C, lambda,x, theta)
		x = [1 ; x ]

		Ax = Mx(A,x)
		Bx = Mx(B,x)
		Cx = Mx(C,x)

	for i = 1:25,
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




	
