using Plotly
# Plotly.signin("naltieri","829vui2wug")

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
	for i = 1:size(A,1)
		for j = i+1:size(A,2)
			avg = .5*(pinvA[i,j] + pinvA[j,i])
			pinvA[i,j] = avg;
			pinvA[j,i] = avg;
		end
	end
	# println("symPinvTest")
	# println(issym(pinvA))
	return(pinvA)
end

function makeSymmetric(A)
	#Weirdness for roundoff, should not be more machEps?
	for i = 1:size(A,1)
		for j = i+1:size(A,1)
			avg = .5*(A[i,j] + A[j,i])
			A[i,j] = avg;
			A[j,i] = avg;
		end
	end
	return(A)
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

	listOfSV = [];

	Fx = zeros(size(F,1),size(F,2));

	for i = 1:size(F,3)
	    Fx = Fx + F[:,:,i]*x[i];
	end 
	try
		chol(Fx)
	catch exc
		println("x is not feasible...quitting analyticCenter ")
		return(NaN,listOfSV)
	end

	if(cond(Fx) == Inf)
		println("Fx is singular...quitting")
		return(NaN,listOfSV)
	end
	#Compute symmetric pseudo inverse
	# (D,V) = eig(Fx);
	# dTol = maximum(size(Fx))*norm(Fx)*eps();
	# D = abs(D);
	# D = D .*[D .> dTol]
	# dInv = zeros(Float64,length(D))
	# for i = 1:length(D)
	# 	if D[i] != 0
	# 		dInv[i] = 1/D[i];
	# 	end
	# end
	# pinvFx = V*diagm(dInv)*V' 
	# pinvFx = symPinv(Fx);
	pinvFx = pinv(Fx)
	# svdTol = 10.0^(-20)
	# (U,S,V) = svd(Fx);
	# invS = 1./S .* [S .> svdTol]
	# pinvFx = V*diagm(invS)*U';

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

	# H = makeSymmetric(H)
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
	maxIter = 10;
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
		# println(issym(H))
		# #Compute symmetric pseudo inverse
		# (D,V) = eig(H);
		# dTol = maximum(size(H))*norm(H)*eps();
		# D = abs(D);
		# D = D .*[D .> dTol]
		# dInv = zeros(Float64,length(D))
		# for i = 1:length(D)
		# 	if D[i] != 0
		# 		dInv[i] = 1/D[i];
		# 	end
		# end
		# pinvH = V*diagm(dInv)*V' 
		pinvH = pinv(H);
		# pinvH = symPinv(H)
		# svdTol = 10.0^(-5)
		# (U,S,V) = svd(H);
		# invS = 1./S .* [S .> svdTol]
		# pinvH = V*diagm(invS)*U';

		(UH,SH,VH) = svd(H);
		listOfSV = [listOfSV ,SH];
		# print("Difference between symmetric and regular")
		# println(vecnorm(pinvH-pinvH2))
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
			return(xPrev,listOfSV)
		end

		if(cond(Fx) == Inf)
			println("Fx is singular...quitting")
			return(xPrev,listOfSV)
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
		#Compute symmetric pseudo inverse
		# (D,V) = eig(Fx);
		# dTol = maximum(size(Fx))*norm(Fx)*eps();
		# D = abs(D);
		# D = D .*[D .> dTol]
		# dInv = zeros(Float64,length(D))
		# for i = 1:length(D)
		# 	if D[i] != 0
		# 		dInv[i] = 1/D[i];
		# 	end
		# end
		# pinvFx = V*diagm(dInv)*V' 
		pinvFx = pinv(Fx)
		# pinvFx = symPinv(Fx)
		# svdTol = 10.0^(-5)
		# (U,S,V) = svd(Fx);
		# invS = 1./S .* [S .> svdTol]
		# pinvFx = V*diagm(invS)*U';

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
		for i = 2:size(F,3)
			pinvFxFi = pinvFx*F[:,:,i];
			g[i-1] = -trace(pinvFxFi)
			pinvFxFiPinvFx  = pinvFxFi*pinvFx;
		    for j = 2:size(F,3)
		    	# println("sym test 1")
		    	# println(issym(Fx))
		    	# println(issym(F[:,:,i]))
		    	# println(issym(pinvFx))
		    	# iJ = trace(pinvFx*F[:,:,i]*pinvFx*F[:,:,j])
		    	# jI = trace(pinvFx*F[:,:,j]*pinvFx*F[:,:,i])
		    	# println("sym test 2")
		    	# println(iJ == jI)
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

		# H = makeSymmetric(H)



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
	    # gHist = [gHist,norm(g,2)];
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
		# hHist = [hHist,cond(H)]
		# fHist = [fHist,cond(Fx)]


	end
	# println
	# firstArg = 1:length(fHist)
	# secondArg = fHist;
	# Plotly.plot([firstArg secondArg])

	xOpt = x;

	return(xOpt,listOfSV)
end



function methOfCents( A,B,C, lambda,x, thetaInit, thetaMul)
		totalTime = 0;
		lVect = [0 lambda];
		curIter = 1;
		tol = 10.0^(-20)
		# tol = 0;
		lambdaPrev = Inf;
		xPrev = [1;x]
		x = [1 ; x ]

		theta = thetaInit
		bigList = [];

		# println(x)

		Ax = Mx(A,x)
		Bx = Mx(B,x)
		# Cx = Mx(C,x)
		F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3))
		F[(size(A,1)+1):end,(size(A,1)+1):end,:] = C

	while lambdaPrev - lambda > tol

		# theta = rand();
		tic()
		lambdaPrev = copy(lambda);
		# try
			lambda =  (1-theta)*maximum(eig(Ax,Bx)[1])+ theta*lambda;
			println(lambda)

			#Concatenate lambda*B - A with C
		    
		    F[1:size(A,1),1:size(A,1),:] = lambda*B-A

		    xPrev = copy(x);
		    (x,listOfSV) = analyticCenter(x,F, .01,10.0^(-3))
		    # println(isnan(x))
		    if(isnan(x)[1])
		    	println("Analytic Center has given up")
		    	theta = theta*thetaMul;
		    	x = copy(xPrev)
		    	lambda = copy(lambdaPrev)
		    	#We could do some other kind of back off, it might increase improvement
		    	if theta > 1
		    		println("THETA GREATER THAN 1!")
		    		return(x,lambda,lVect)
		    	end
		    else

			    bigList = [bigList ,listOfSV]

			    #Recompute Ax,Bx,Cx
			    Ax = Mx(A,x);
				Bx = Mx(B,x);
				# Cx = Mx(C,x);
			end

		# totalTime += toc();
		# lVect = [lVect ;[totalTime copy(lambda)]]


	end
	toc();

	if (lambdaPrev - lambda <= tol)
		println("Tolerance Achieved")
	end 
	
	# println(eigvals(Bx))
	# println(eigvals(Cx))
	# println(eigvals(lambda*Bx - Ax))
	# toc();
	return(x,lambda,lVect,bigList)
end




	
