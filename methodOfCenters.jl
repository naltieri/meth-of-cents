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

function genEigSvd(A,B)
	#Computes the maximum generalized eigenvalue of A,B by factorization via SVD
	(D,U) = eig(B);
	# dTol = maximum(size(B))*norm(B)*eps();
	# dInv = 1./D .* [D .> dTol]
	dInv = 1./D;

	invDhalf = diagm(sqrt(dInv));
	newMat = invDhalf*U'*A*U*invDhalf;
	C = makeSymmetric(newMat);
	l = eigmax(C)
	return(l)
end


function genEigChol(A,B)
	#Computes the maximum generalized eigenvalue of A,B by factorization via SVD
	L = chol(B,:L);
	# newMat = inv(L)*A*inv(L');
	tempMat = L\A;
	tempMat = tempMat';
	tempMat2 = L\tempMat;
	newMat = tempMat2';

	C = makeSymmetric(newMat);
	l = eigmax(C)
	return(l)
end

function feasTest(lambda,A,B,x)
	# tol = eps()*100;
	tol = 0;
	Ax = Mx(A,x);
	Bx = Mx(B,x);
	lHat = genEigChol(Ax,Bx)
	if(lHat <= lambda+tol)
		# try
		# 	chol(lambda*Bx - Ax)
		# catch
		# 	println("-----GE (YES) AND CHOL (NO) DISAGREE -----")
		# end
		return(true)
	else
		# try
		# 	chol(lambda*Bx - Ax)
		# 	println("-----GE (NO) AND CHOL (YES) DISAGREE-----")
		# catch
		
		# end
		return(false)
	end
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

function symPinv(A)
	#Compute symmetric pseudo inverse
	(D,V) = eig(A);
	dTol = maximum(size(A))*norm(A)*eps();
	D = abs(D);
	# D = D .*[D .> dTol]
	# dInv = zeros(Float64,length(D))
	# for i = 1:length(D)
	# 	if D[i] != 0
	# 		dInv[i] = 1/D[i];
	# 	end
	# end
	dInv = 1./D .* [D .> dTol]
	pinvA = V*diagm(dInv)*V' 
	#Weirdness for roundoff, should not be more machEps?
	makeSymmetric(pinvA)
	# println("symPinvTest")
	# println(issym(pinvA))
	# pinvA = inv(A)
	# pinvA = inv(A)+1e6*eps()*eye(size(A,1))
	return(pinvA)
end

function makeSymmetric(A)
	#Weirdness for roundoff, should not be more machEps?
	for i = 1:size(A,1)
		for j = i+1:size(A,1)
			avg = .5*(A[i,j] + A[j,i])
			A[i,j] = avg;
			A[j,i] = avg;
			# if(rand(1)[1] > .5)
			# 	val = A[i,j]
			# else
			# 	val = A[j,i]
			# end
			# A[i,j] = val;
			# A[j,i] = val;

		end
	end
	return(A)
end

function projectOntoPSD(A)
	(D,V) = eigs(A)
	# println(D)
	D = max(0,D);
	# println(diagm(D))
	posA = V*diagm(D)*V';
	posA = makeSymmetric(posA)
	return(posA)
end
function analyticCenter(x0,F, alpha, tol,A,B,C,lambda)
	# iter = 1;
	minL = Inf;
	minIter = 0;
	minX = x0;
	if(cond(Mx(F,x0)) == Inf)
		println("........ cond is Inf at start")
	end
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
		# chol(Fx)
		if !feasTest(lambda,A,B,x)
			chol(-eye(1))
		end
		chol(Mx(C,x))

		if(cond(Fx) == Inf)
			throw(DomainError())
		end
	catch exc
		# println("x is not feasible...quitting analyticCenter ")
		# println(x)
		throw(DomainError())
		# println("projecting...")
		# Fx = projectOntoPSD(Fx)
	end

	# if(cond(Fx) == Inf)
	# 	println("Fx is singular...quitting")
	# 	throw(DomainError())
	# end	
	# pinvFx = Array(Float64,size(Fx,1),size(Fx,2))
	# pinvFxBad = pinv(Fx)
	# for i = 1:size(pinvFx,1)
	# 	for j = 1:size(pinvFx,1)
	# 		pinvFx[i,j] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
	# 		# pinvFx[j,i] = .5(pinvFxBad[i,j]+pinvFxBad[j,i])
	# 	end
	# end
	pinvFx = symPinv(Fx);

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

	H = makeSymmetric(H)


	# g = zeros(size(F,3)-1,1);		 
	# H = zeros(size(F,3)-1,size(F,3)-1)
	# for i = 2:size(F,3)
	# 	pinvFxFi = pinvFx*F[:,:,i];
	# 	g[i-1] = -trace(pinvFxFi)
	# 	pinvFxFiPinvFx  = pinvFxFi*pinvFx;
	#     for j = i:size(F,3)
	#     	if(j==i)
	# 	        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
	# 	    else
	# 	    	if(rand(1)[1] > .5)
	# 		        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
	# 		        H[j-1,i-1] = copy(H[i-1,j-1])
	# 		    else
	# 				H[i-1,j-1] = trace(F[:,:,j]*pinvFxFiPinvFx);
	# 		        H[j-1,i-1] = copy(H[i-1,j-1])
	# 		    end
	# 	    end

	#         # println("IM HAPPY")
	#         # println(issym(F[:,:,i]))
	#         # println(issym(Fx))
	#         # println(issym(pinvFxFiPinvFx))
	#         # println(issym(F[:,:,j]))
	 
	#         # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(sum(pinvFxFiPinvFx .* Fj')))
	#         #H[i-1,j-1] = sum(computeDiagAB(pinvFxFiPinvFx , F[:,:,j]));
	#         # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(pinvFxFiPinvFx[:] .* F[:,:,j][:]))
	#     end
	# end



		pinvH = symPinv(H)
		delta = norm(sqrtm(pinvH)*g,2);
		if delta <= 1/4
			alpha = 1;
		else
			alpha = 1/(1+delta)
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
	curIter = 1;
	maxIter = 15;
	firstIter = true;

	while (norm(g,2) > tol) & (curIter < maxIter)


		newL = genEigChol(Mx(A,x),Mx(B,x))
		if newL < minL
			minL = copy(newL);
			minX = copy(x);
			minIter = copy(curIter)
		end
	    # lBA = F[1:5,1:5,:];
	    # println(cond(Mx(lBA,x))) 

		curIter += 1;
		# print("Iterating...")
		# println(curIter)
		# print("Current Gradient is ........")
		# println(norm(g,2))


		if curIter == maxIter-1
			println("=====================")
			print("Condition number of Fx is ")
			println(cond(Fx))
			print("Final Gradient is ........................ ")
			println(norm(g,2))
			println("=====================")
		end


		# print("Current Difference is ........")
		# println(norm(xPrev-x,2))

		# print("Current stepsize is ........")
		# println(alpha)

		if (norm(xPrev - x,2) <= 0.0)
			break
		end
		# alphaG = alphaG/sqrt(curIter);
		# curIter += 1;
		smallNumber = 0;
		H = H + smallNumber
		pinvH = symPinv(H)
		# pinvH = pinv(H);
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
		# try 
		# 	L = chol(H,:L)
		# catch 
		# 	println("H is not PSD... :(")
		# 	return(xPrev)
		# end
		# L = chol(H,:L)
		# firstRes = L\g;
		# secRes = L'\firstRes;
		# if firstIter
		# 	y = copy(x)
		# 	success = false
		# 	while success == false
		# 		y[2:end] = x[2:end] - alpha * pinvH *g;
		# 		Fx = zeros(size(F,1),size(F,2))
		# 		for i = 1:size(F,3)
		# 		    Fx = Fx + F[:,:,i]*y[i];
		# 		end
		# 		try
		# 			chol(Fx)
		# 			if(cond(Fx) == Inf)
		# 				# println("Can't take step, boundary failure")
		# 				throw(DomainError())
		# 			end
		# 			success = true;
		# 		catch
		# 			# println("Can't take step, failure")
		# 			alpha = alpha * .9
		# 		end

		# 	end
		# end
		# firstIter = false;
		x[2:end] = x[2:end] - alpha * pinvH *g;
		# x[2:end] = x[2:end]-alphaG* g;


		Fx = zeros(size(F,1),size(F,2))
		for i = 1:size(F,3)
		    Fx = Fx + F[:,:,i]*x[i];
		end
		# println("Eigvals of Fx are:")
		# println(eigvals(Fx))

		try
			# chol(Fx)
			if !feasTest(lambda,A,B,x)
				chol(-eye(1))
			end
			chol(Mx(C,x))

		catch
			println("=====================")
			print("Analytic Center HAS FAILED (made Fx infeasible) on iteration ")
			println(curIter)
			print("Condition number of Fx is ")
			println(cond(Fx))
			print("Final Gradient is ........................ ")
			println(norm(g,2))
			println("=====================")
			# lBA = F[1:5,1:5,:];
	  		#println(cond(Mx(lBA,x))) 
			# return(xPrev)
			println(curIter)
			println(minIter)
			return(minX)
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
		pinvFx = symPinv(Fx)

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
		 

		# g = zeros(size(F,3)-1,1);		 
		# H = zeros(size(F,3)-1,size(F,3)-1)
		# for i = 2:size(F,3)
		# 	pinvFxFi = pinvFx*F[:,:,i];
		# 	g[i-1] = -trace(pinvFxFi)
		# 	pinvFxFiPinvFx  = pinvFxFi*pinvFx;
		#     for j = 2:size(F,3)
		#     	# 	Fj = F[:,:,j];
		#         # H[i-1,j-1] = sum(sum(pinvFxFiPinvFx .* Fj'));
		#         H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
		#         # println("IM HAPPY")
		#         # println(issym(F[:,:,i]))
		#         # println(issym(Fx))
		#         # println(issym(pinvFxFiPinvFx))
		#         # println(issym(F[:,:,j]))
		 
		#         # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(sum(pinvFxFiPinvFx .* Fj')))
		#         #H[i-1,j-1] = sum(computeDiagAB(pinvFxFiPinvFx , F[:,:,j]));
		#         # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(pinvFxFiPinvFx[:] .* F[:,:,j][:]))
		#     end
		# end

		# H = makeSymmetric(H)

		g = zeros(size(F,3)-1,1);		 
		H = zeros(size(F,3)-1,size(F,3)-1)
		for i = 2:size(F,3)
			pinvFxFi = pinvFx*F[:,:,i];
			g[i-1] = -trace(pinvFxFi)
			pinvFxFiPinvFx  = pinvFxFi*pinvFx;
		    for j = i:size(F,3)
		    	if(j==i)
			        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
			    else
			    	if(rand(1)[1] > .5)
				        H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
				        H[j-1,i-1] = copy(H[i-1,j-1])
				    else
						H[i-1,j-1] = trace(F[:,:,j]*pinvFxFiPinvFx);
				        H[j-1,i-1] = copy(H[i-1,j-1])
				    end
			    end

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

	# return(xOpt)
	println(curIter)
	println(minIter)
	return(minX)
end



function methOfCents( A,B,C, lambda,x, thetaInit, thetaMul,k)
		tic();
		# tol = 10.0^(-20)
		tol = 0;
		lambdaPrev = Inf;
		xPrev = [1;x]
		x = [1 ; x ]
		#Check lambda values
		println("-----------------")
			println(lambda)
			println(maximum(eig(Mx(A,x),Mx(B,x))[1]))
			# println(maximum(eig(float32(Mx(A,x)),float32(Mx(B,x)))[1]))

		println("-----------------")
		#Test feasibility
		try
			chol(Mx(C,x))
		catch
			println("Cx not feasible")
		end
		try
			if !feasTest(lambda,A,B,x)
				chol(-eye(1))
			end
			# chol(lambda*Mx(B,x)-Mx(A,x))
		catch
			println("lambdaB-A not feasible")
		end
		theta = thetaInit

		# println(x)

		Ax = Mx(A,x)
		Bx = Mx(B,x)
		Cx = Mx(C,x)
		# thetaMul = 1;
		# thetaInit = .1;
		thetaList = .1:.1:.9
		newX = NaN;
		newLambda = NaN;

		F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
	    F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C


	while lambdaPrev - lambda > tol

		# if lambda > 1.4
		# 	return
		# end


		# theta = rand();

		lambdaPrev = copy(lambda);
		xPrev = copy(x);



		# success = false;
		# for i = 1:length(thetaList)
		# 	if success == false
		# 		try
		# 			theta = thetaList[i]
		# 			newLambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
		# 			println(newLambda)

		# 			#Concatenate lambda*B - A with C
		# 		    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
		# 		    F[1:size(A,1),1:size(A,1),:] = lambda*B-A
		# 		    F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C
		# 		    newX = analyticCenter(x,F, .01,10.0^(-3))

		# 		    #Recompute Ax,Bx,Cx
		# 		    Ax = Mx(A,x);
		# 			Bx = Mx(B,x);
		# 			Cx = Mx(C,x);
		# 			success = true;
		# 		catch exc
		# 			println(exc)
		# 		end
		# 	end
		# end

		# if success == true
		# 	x = copy(newX)
		# 	lambda = copy(newLambda)
		# 	xPrev = copy(x);
		# end

		# if success == false
		# 	println("analytic is sad")
		# 	return(x,lambda)
		# end

		# success = false;
		# while success == false && theta < 1
		# 	if success == false
		# 		try
		# 			newLambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
		# 			println(newLambda)

		# 			#Concatenate lambda*B - A with C
		# 		    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
		# 		    F[1:size(A,1),1:size(A,1),:] = lambda*B-A
		# 		    F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C
		# 		    newX = analyticCenter(x,F, .01,10.0^(-3))

		# 		    #Recompute Ax,Bx,Cx
		# 		    Ax = Mx(A,x);
		# 			Bx = Mx(B,x);
		# 			Cx = Mx(C,x);
		# 			success = true;
		# 		catch exc
		# 			println(exc)
		# 			theta = theta*1.1
		# 		end
		# 	end
		# end

		# if success == true
		# 	x = copy(newX)
		# 	lambda = copy(newLambda)
		# 	xPrev = copy(x);
		# end

		# if success == false
		# 	println("analytic is sad")
		# 	return(x,lambda)
		# end
		try	


			# lambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
			# success = false
			# while success == false
			# 	try
			# 		println("IM TRYINGGG")
			# 		newLambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
			# 		F[1:size(A,1),1:size(A,1),:] = newLambda*B-A
			# 		try
			# 			Fx = Mx(F,x)
			# 			chol(Fx)
			# 			success = true
			# 		catch
			# 			theta = theta*1.1
			# 			theta = min(theta,1);
			# 			if theta == 1
			# 				success = true;
			# 			end
			# 		end
			# 	end
			# end
			# println("----")
			# println(cond(Ax))
			# println(cond(Bx))
			# println("-----")
						# println(eig(lambda*Bx-Ax)[1])

			# println(maximum(eig(Mx(A,x),Mx(B,x))[1]) - maximum(eig(float32(Mx(A,x)),float32(Mx(B,x)))[1]))
			# lambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
			oldLambda = lambda;
			lambda =  (1-theta)*genEigChol(Ax,Bx)+ theta*lambda;

			print("Improvement is ------------- ")
			println(oldLambda-lambda)
			# if(oldLambda - lambda < 0)
			# 	return(xPrev,lambdaPrev)
			# end

			println(lambda)

			#Concatenate lambda*B - A with C
		    # F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3)
		    # smallNumber = 1e4*eps();
		    # pushFactor = smallNumber*eye(size(Ax,1))
		    F[1:size(A,1),1:size(A,1),:] = copy(lambda*B-A)
		    # F[1:size(A,1),1:size(A,1),end] -= pushFactor

		    # F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C
		    xPrev = copy(x);
			# print("COND OF FX IS ----- ")
			# println(cond(Mx(F,x)))
			print("Cond (lB-A)(x): ")
			println(cond(lambda*Bx-Ax))
			# println(minimum(eig(lambda*Bx-Ax)[1]))
			# println(eig(lambda*Bx-Ax)[1])

  
		    x = analyticCenter(x,F, .01,10.0^(-3),A,B,C,lambda)
		  #   if(cond(lambda*Bx-Ax) > 1e15)
			 #    try
			 #    	x = getMosekX(k,lambda)
			 #    catch
			 #    	println("MOSEK HAS FAILED")
				#     x = analyticCenter(x,F, .01,10.0^(-3))
				# end
		  #   else
		  #  		x = analyticCenter(x,F, .01,10.0^(-3))
		  #  	end
		    #Recompute Ax,Bx,Cx
		    Ax = Mx(A,x);

			Bx = Mx(B,x);
		 #    println(cond(Ax))
			# println(cond(Bx))
			# Cx = Mx(C,x);
			# println(cond(Cx))



		catch exc
			println(exc)
			toc();
			return(x,lambda)
			# if isa(exc,DomainError)		
		 #    	println("Analytic Center has given up")
		 #    	theta = theta*thetaMul;
		 #    	#We could do some other kind of back off, it might increase improvement
		 #    	if theta > 1
		 #    		println("THETA GREATER THAN 1!")
		 #    		return(x,lambda)
		 #    	end
		 #    else
		 #    	throw(exc)
		 #    end

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
