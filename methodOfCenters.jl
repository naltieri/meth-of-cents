using Plotly

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

function traceTest(n)
	A = rand(n,n);
	B = rand(n,n);
	tic()
	trace(A*B)
	toc()
	tic()
	sum(sum(A.*B'))
	toc()
end

function Mx(M,x)
	out = zeros(size(M,1),size(M,2))
	for i = 1:size(M,3)
		out = out + M[:,:,i]*x[i]
	end
	return(out)
end

function computeDiagAB(A,B)
	#Compute diag(A*B), A,B square
	d = zeros(size(A,1),1)
	for i = 1:size(A,1)
		#JULIAA WHYYYYY
		d[i] = (A[i,:]*B[:,i])[1]
		#d[i][:] = A[i,:]*B[:,i]
	end
	return(d)
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
	try
		chol(Fx)
	catch exc
		println("x is not feasible...quitting analyticCenter ")
		throw(DomainError())
	end

	if(cond(Fx) == Inf)
		println("Fx is singular...quitting")
		throw(DomainError())
	end	
	# pinvFx = Array(Float64,size(Fx,1),size(Fx,2))
	# pinvFxBad = pinv(Fx)
	# for i = 1:size(pinvFx,1)
	# 	for j = 1:size(pinvFx,1)
	# 		pinvFx[i,j] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
	# 		# pinvFx[j,i] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
	# 	end
	# end
	pinvFx = pinv(Fx);

	g = zeros(size(F,3)-1,1);		 
	H = zeros(size(F,3)-1,size(F,3)-1)
	for i = 2:size(F,3)
		pinvFxFi = pinvFx*F[:,:,i];
		g[i-1] = -trace(pinvFxFi)
		pinvFxFiPinvFx  = pinvFxFi*pinvFx;
	    for j = 2:size(F,3)
	    	# 	Fj = F[:,:,j];
	        # H[i-1,j-1] = sum(sum(pinvFxFiPinvFx .* Fj'));
	        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
	        # println("IM HAPPY")
	        # println(issym(F[:,:,i]))
	        # println(issym(Fx))
	        # println(issym(pinvFxFiPinvFx))
	        # println(issym(F[:,:,j]))
	 
	        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(sum(pinvFxFiPinvFx .* Fj')))
	        #H[i-1,j-1] = sum(computeDiagAB(pinvFxFiPinvFx , F[:,:,j]));
	        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(pinvFxFiPinvFx[:] .* F[:,:,j][:]))
	    end
	end
	# pinvHhat = float64((pinv(float32(H))));
	# appErr = norm((pinv(H))*g-pinvHhat*g)/(norm(pinv(H)*g));
	# if appErr > 5
	# 	println("Hessian Too Ill-Conditioned...Quitting YOOOOOOOOO")
	# 	throw(DomainError())
	# end

	# objHist = log(det(pinv(Fx)));
	gHist = norm(g,2);
	hHist = cond(H);
	fHist = cond(Fx);

	# Second condition is a hack. Happens when numerical underflow
	# occurs in updating x
	# Replacing 0 with 10^(-10) creates bad results. 


	 #    sFx = float32(Fx);

	 #    sH = zeros(Float32,size(F,3)-1,size(F,3)-1)
		# for i = 2:size(F,3)
		#     for j = 2:size(F,3)
		#         sH[i-1,j-1] = trace(pinv(sFx) *float32(F[:,:,i])* 	pinv(sFx)* float32(F[:,:,j]));
		#     end
		# end

	alphaG = .1;
	curIter = 1.0;
	maxIter = 15;
	while (norm(g,2) > tol) & (curIter < maxIter)	

		curIter += 1;
		#println(curIter)

		if (norm(xPrev - x,2) <= 0.0)
			break
		end
		# alphaG = alphaG/sqrt(curIter);
		# curIter += 1;
		smallNumber = 0;
		H = H + smallNumber
		pinvH = pinv(H)
		delta = norm(sqrtm(pinvH)*g,2);
		# println(svd(H)[2] )
		# println(g)

		# println(delta)
		if delta <= 1/4
			alpha = 1;
		else
			alpha = 1/(1+delta)
		end

		xPrev = copy(x);
		# println(cond(H))
		#Changing pinv to inv, prevented some solutions from becoming infeasible.
		#Fuck if I know why
		# try 
		# 	L = chol(H,:L)
		# catch 
		# 	println("H is not PSD... :(")
		# 	return(xPrev)
		# end
		# L = chol(H,:L)
		# firstRes = L\g;
		# secRes = L'\firstRes;



		x[2:end] = x[2:end] - alpha * pinvH *g;
		# x[2:end] = x[2:end]-alphaG* g;


		Fx = zeros(size(F,1),size(F,2))
		for i = 1:size(F,3)
		    Fx = Fx + F[:,:,i]*x[i];
		end
		# println("Eigvals of Fx are:")
		# println(eigvals(Fx))

		try
			chol(Fx)
		catch
			println("Analytic Center HAS FAILED (made Fx infeasible)")
			return(xPrev)
		end

		if(cond(Fx) == Inf)
			println("Fx is singular...quitting")
			return(xPrev)
		end
		# # println(cond(Fx))
		# L = chol(Fx);
		# invLFL = zeros(size(F,1),size(F,2),size(F,3)-1)

		# invL = inv(L);
		# for i = 2:size(F,3)
		# 	invLFL[:,:,i-1] = invL * F[:,:,i] *invL';
		# end







		# # #We will precompute the quantity Fx^(-1) * F_i = L^(-1)F_i L^T and store it.
		# # #We will compute L^(-1) F_i by doing back substitution 
		# # invLF = zeros(size(F,1),size(F,2),size(F,3)-1)
		# # for i = 2:size(F,3)
		# # 	invLF[:,j,i-1] = L\F[:,j,i]
		# # end




		# for i = 2:size(F,3)
		#     g[i-1] = - trace(invLFL[:,:,i-1]);
		# end
		


		# for i = 2:size(F,3)
		#     for j = 2:size(F,3)
		#         H[i-1,j-1] = trace(invLFL[:,:,i-1] * invLFL[:,:,j-1]);
		#     end
	 #    end



		# g = zeros(1,size(F,3)-1);
		# for i = 2:size(F,3)
		#     g[i-1] = -trace(pinv(Fx)*F[:,:,i])
		# end
		# g = g';





		# (U,S,V) = svd(Fx);
		# svdTol = maximum(size(Fx))*norm(Fx)*eps();
		# Strunc = S.*(S .> svdTol);
		# FxStab = U*diagm(Strunc)*V';

		# try
		# 	L = chol(FxStab,:L);
		# catch
		# 	println("LOL NOOB")
		# end

		# L = chol(FxStab,:L)


		

		# g2 = zeros(size(F,3)-1,1);
		# for i = 2:size(F,3)
		# 	firstRes = L\F[:,:,i];
		# 	firstRes = L'\firstRes;
		#     g2[i-1] = -trace(firstRes);
		# end

		# g = g2;
		# println("TESTING!")
		# println(issym(Fx))
	
		# pinvFxBad = pinv(Fx)
		# for i = 1:size(pinvFx,1)
		# 	for j = 1:size(pinvFx,1)
		# 		pinvFx[i,j] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
		# 		# pinvFx[j,i] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
		# 	end
		# end
		pinvFx = pinv(Fx)
		# println(issym(pinvFx))
		# pinvFxF = Array(Float64,size(F,1),size(F,2),size(F,3)-1)
		# println(size(pinvFxF))
		# for i = 2:size(F,3)
		# 	println(size(pinvFxF[:,:,i]))
		# 	println(size(pinvFx*F[:,:,i]))
		# 	pinvFxF[:,:,i] = pinvFx*F[:,:,i]
		# end

		# g = zeros(size(F,3)-1,1);
		# for i = 2:size(F,3)
		# 	g[i-1] = -trace(pinvFx*F[:,:,i])
		#     # g[i-1] = -sum(sum(pinvFx.*F[:,:,i]))
		# end
		 

		g = zeros(size(F,3)-1,1);		 
		H = zeros(size(F,3)-1,size(F,3)-1)
		for i = 2:size(F,3)
			pinvFxFi = pinvFx*F[:,:,i];
			g[i-1] = -trace(pinvFxFi)
			pinvFxFiPinvFx  = pinvFxFi*pinvFx;
		    for j = 2:size(F,3)
		    	# 	Fj = F[:,:,j];
		        # H[i-1,j-1] = sum(sum(pinvFxFiPinvFx .* Fj'));
		        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
		        # println("IM HAPPY")
		        # println(issym(F[:,:,i]))
		        # println(issym(Fx))
		        # println(issym(pinvFxFiPinvFx))
		        # println(issym(F[:,:,j]))
		 
		        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(sum(pinvFxFiPinvFx .* Fj')))
		        #H[i-1,j-1] = sum(computeDiagAB(pinvFxFiPinvFx , F[:,:,j]));
		        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(pinvFxFiPinvFx[:] .* F[:,:,j][:]))
		    end
		end

		# pinvHhat = float64((pinv(float32(H))));
		# # println("-------------------------")
		# # println(norm(pinv(H)*g-pinvHhat*g)/(norm(pinv(H)*g)))
		# # println("-------------------------")
		# appErr = norm(pinv(H)*g-pinvHhat*g)/(norm(pinv(H)*g));
		# if appErr > 1
		# 	println("Hessian Too Ill-Conditioned...Analytic Center Failed YOYOYOYOYOYO")
		# 	return(x)
		# end



		# println("Sup with it")

		# println(cond(Fx))


		# H = zeros(size(F,3)-1,size(F,3)-1)
		# for i = 2:size(F,3)
		#     for j = 2:size(F,3)
		#     	firstResI = L\F[:,:,i];
		# 		secResI = L'\firstResI;
		# 		firstResJ = L\F[:,:,j];
		# 		secResJ = L'\firstResJ;
		#         H[i-1,j-1] = trace(secResI*secResJ);
		#     end
		# end

	 #    sFx = float32(Fx);

	 #    sH = zeros(Float32,size(F,3)-1,size(F,3)-1)
		# for i = 2:size(F,3)
		#     for j = 2:size(F,3)
		#         sH[i-1,j-1] = trace(inv(sFx) *float32(F[:,:,i])* 	inv(sFx)* float32(F[:,:,j]));
		#     end
		# end

		# println(typeof(sH))
		# println(norm(float64(sH)-H))
	    gHist = [gHist,norm(g,2)];
	 #    try
		# 	# objHist = [objHist, log(det(pinv(Fx)))];
		# 	chol(Fx)
		# catch exc
		# 	if isa(exc,InterruptException)
		# 		throw(InterruptException())
		# 	end
		# 	println("made infeasible in analytic center")
		# 	return(xPrev)
		# end
		hHist = [hHist,cond(H)]
		fHist = [fHist,cond(Fx)]


	end
	# firstArg = 1:length(fHist)
	# secondArg = fHist;
	# Plotly.plot([firstArg secondArg])

	xOpt = x;

	return(xOpt)
end



function methOfCents( A,B,C, lambda,x, thetaInit, thetaMul)
		tic();
		tol = 10.0^(-20)
		lambdaPrev = Inf;
		xPrev = [1;x]
		x = [1 ; x ]

		theta = thetaInit

		# println(x)

		Ax = Mx(A,x)
		Bx = Mx(B,x)
		Cx = Mx(C,x)

	while lambdaPrev - lambda > tol

		# theta = rand();

		lambdaPrev = copy(lambda);
	
		try	
			lambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
			println(lambda)

			#Concatenate lambda*B - A with C
		    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
		    F[1:size(A,1),1:size(A,1),:] = lambda*B-A
		    F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C
		    xPrev = copy(x);
		    x = analyticCenter(x,F, .01,10.0^(-3))

		    #Recompute Ax,Bx,Cx
		    Ax = Mx(A,x);
			Bx = Mx(B,x);
			Cx = Mx(C,x);

		catch exc
			println(exc)
			if isa(exc,DomainError)		
		    	println("Analytic Center has given up")
		    	theta = theta*thetaMul;
		    	#We could do some other kind of back off, it might increase improvement
		    	if theta > 1
		    		println("THETA GREATER THAN 1!")
		    		return(x,lambda)
		    	end
		    else
		    	throw(exc)
		    end

		end

	end

	if (lambdaPrev - lambda <= tol)
		println("Tolerance Achieved")
	end 
	
	# println(eigvals(Bx))
	# println(eigvals(Cx))
	# println(eigvals(lambda*Bx - Ax))
	toc();
	return(x,lambda)
end




	
