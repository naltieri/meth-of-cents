using Mosek
include("feas_point.jl")
include("methodOfCenters.jl")



function makeBasis(n)
	#Creates a basis for nxn symmetric matrices
	basis = zeros(n,n,convert(Int64, n+n*(n-1)/2));
	curK = 1;
	for i = 1:n
		b = zeros(n,n)
		b[i,i] = 1;
		basis[:,:,curK] = b;
		curK += 1;
	end

	for i = 1:n
		for j = i+1:n
			b = zeros(n,n)
			b[i,j] = 1;
			b[j,i] = 1;
			basis[:,:,curK] = b;
			curK += 1;
		end
	end

	return(basis)
end


function getCoeff(P)
	#Takes in a PSD P, and get the coefficients corresponding to the basis vectors in make Basis
	n = size(P,1)
	x = Array(Float64, 1,convert(Int64, n+n*(n-1)/2))
	curK  = 1;
	for i = 1:n
		x[curK] = P[i,i];
		curK += 1;
	end

	for i = 1:n
		for j = i+1:n
			x[curK] = P[i,j] 
			curK += 1;
		end
	end



	return x'
end	

function getMosekX(k,lambdaInit)
	temp = feas_point( k, lambdaInit )
	# println((P,t1,t2))
	if length(temp) == 1
		println("feas_point fail")
		return(NaN)
	end
	(P, t1, t2) = temp;
	x = getCoeff(P);
	return([1;x;t1])
end

function getRho(k,lambdaInit,flag)
	(P, t1, t2) = feas_point( k, lambdaInit )
	# println((P,t1,t2))
	if length(feas_point(k,lambdaInit)) == 1
		println("feas_point fail")
		return
	end

	x = getCoeff(P);
	

	(A,B,C) = formulateGevp(k);


	# G = makeBasis(3);
	# println((Mx(G,x),t1,1.0)) 

	# tX = Mx(A,[1;x;t1])

	# println(tX)



	# trueA = makeTrueA(k,P,t1)

	# 	println(trueA)

	# println(tX - trueA)


	#println((A,B,C))
	if(flag == true)
		(lambdaInit,x) = advBisectionSearchNest(k,0,1.3,1e-5)
		(xOpt,lambdaOpt)= methOfCents( A,B,C, lambdaInit,x[2:end], .1,1.1,k)
	else
		(xOpt,lambdaOpt)= methOfCents( A,B,C, lambdaInit,[x;t1], .1,1.1,k)
		return(xOpt,lambdaOpt)	
	end
end



function createTest(k,lambdaInit)
	(P, t1, t2) = feas_point( k, lambdaInit )
	# println((P,t1,t2))
	if length(feas_point(k,lambdaInit)) == 1
		println("feas_point fail")
		return
	end

	x = getCoeff(P);
	

	(A,B,C) = formulateGevp(k);

	x = [1;x;t1]
	return(x,A,B,C)
end


function formulateGevp(k)

	bet = (sqrt(k)-1)/(sqrt(k)+1);

	xx = [   0  -bet   0  bet
	         1  1+bet  0  -1-bet
	         0    1    0  -1 ];

	x = [ 1 0 0 0
	      0 1 0 0
	      0 0 1 0 ];

	q1 = [ 0  1    0 -1
	       0 -1/k  0  1 ];

	q2 = [ 0  1   -1 -1
	       0 -1/k  0  1 ];

	J = [ 0 1; 1 0];


	# # % parameters of the GEVP.
	# # % note that the optimization variables are {P,t1,t2}
	# # % note also that P is 3x3 here.



	basis = makeBasis(3);

	#Does the GEVP form of this in super tedious painful way
	# A = blkdiag( xx'*P*xx+ t1*q1'*J*q1 + t2*q2'*J*q2, t2 );

	#dimension of first block 
	b1D = size(xx,2)

	#Put terms related to t2 first as the convention of the solver is to 
	#have affine terms first
	lowerOne = zeros(b1D+1,b1D+1);
	lowerOne[b1D+1,b1D+1] = 1;

	temp = zeros(b1D+1,b1D+1);
	temp[1:b1D,1:b1D] = q2'*J*q2;
	A = temp+lowerOne;

	temp = zeros(b1D+1,b1D+1,size(basis,3));
	for i = 1:size(basis,3)
		temp[1:b1D,1:b1D,i] = xx'*basis[:,:,i]*xx;
	end
	A = cat(3,A,temp);
 
	temp = zeros(b1D+1,b1D+1);
	temp[1:b1D,1:b1D] = q1'*J*q1;
	A = cat(3,A,temp);





	#Does the GEVP form of this in super tedious painful way
	# B = blkdiag( x'*P*x, t1 + t2 );



	#dimension of first block
	b1D = size(x,2);



	#Put terms related to t2 first as the convention of the solver is to 
	#have affine terms first

	#HACK ALERT! We put an epsilon in the row of B that has zero for columns
	#and rows to ensure that it can be positive definite.
	smallNumber = eps();
	lowerOneEps = copy(lowerOne);
	lowerOneEps[size(lowerOneEps,1)-1,size(lowerOneEps,2)-1] = smallNumber;
	B = lowerOneEps;



	temp = zeros(b1D+1,b1D+1,size(basis,3));
	for i = 1:size(basis,3)
		temp[1:b1D,1:b1D,i] = x'*basis[:,:,i]*x;
	end
	B = cat(3,B,temp);

	B = cat(3,B,lowerOne);


	#Does the GEVP form of this in super tedious painful way
	# C = blkdiag( P, t1, t2 );

	#dimension of first block
	b1D = 3;

	#Put terms related to t2 first as the convention of the solver is to 
	#have affine terms first
	#Include -eps*I to first affine term to prevent singularities from occuring


	secondLowerOne = zeros(b1D+2,b1D+2);
	secondLowerOne[b1D+2,b1D+2] = 1;

	C = secondLowerOne;

	smallNumber2 = 0; 
	C[1:b1D,1:b1D] = -smallNumber2*eye(b1D)

	temp = zeros(b1D+2,b1D+2,size(basis,3));
	for i = 1:size(basis,3)
		temp[1:b1D,1:b1D,i] = basis[:,:,i];
	end
	C = cat(3,C,temp);

	firstLowerOne = zeros(b1D+2,b1D+2);
	firstLowerOne[b1D+1,b1D+1] = 1;


	C = cat(3,C,firstLowerOne)




	# # % to remove homogeneity, try setting t2 = 1

	return(A,B,C)

end

function makeTrueA(k,P,t1)

	t2 = 1;
	bet = (sqrt(k)-1)/(sqrt(k)+1);

	xx = [   0  -bet   0  bet
	         1  1+bet  0  -1-bet
	         0    1    0  -1 ];

	x = [ 1 0 0 0
	      0 1 0 0
	      0 0 1 0 ];

	q1 = [ 0  1    0 -1
	       0 -1/k  0  1 ];

	q2 = [ 0  1   -1 -1
	       0 -1/k  0  1 ];

	J = [ 0 1; 1 0];

	trueA = full(blkdiag( sparse( xx'*P*xx+ t1*q1'*J*q1 + t2*q2'*J*q2), sparse([1.0]) ))

	return(trueA)
end


function getRhoRepeat(k,lambdaInit)	
	tic()
	#Super hack that results in high precision
	(xOpt,lambdaOpt) = getRho(k,lambdaInit,false)
	minL = lambdaOpt;
	curIter = 1;
	# flag = false;
	while true
		curIter += 1;
		try
			println("=========RESTARTING WITH MOSEK =========")
			(xOpt,lambdaOpt) = getRho(k,lambdaOpt,false)
			minL = min(lambdaOpt,minL);
			# if lambdaOptN > lambdaOpt
			# 	return lambdaOpt
			# else
			# 	lambdaOpt = lambdaOptN;
			# end
		catch exc
			println("Exit Reached")
			println(exc)
			if isa(exc,InterruptException)
				throw(InterruptException())
			end
			if isa(exc,BoundsError)
				print("returning...")
				println(minL)
				toc()
				return(minL)
			end
			
		end
	end
	println(curIter)
	println("end reached")
	toc()
	return(minL)
end

function getRhoSuperRepeat(k,lambdaInit)	
	tic()
	#Super hack that results in high precision
	lambdaPrev = lambdaInit;
	(xOpt,lambdaOpt) = getRho(k,lambdaInit)
	minL = lambdaOpt;
	curIter = 1;
	while true
		curIter += 1;
		try
			println("Hello?")
			xMosek = getMosekX(k,lambdaOpt)
			print("NORM DIFFERENCE IS ... ")
			println(norm(xMosek-xOpt))
							lambdaPrev = copy(lambdaOpt);

			(xOpt,lambdaOpt) = getRho(k,lambdaOpt)

			minL = min(lambdaOpt,minL);
			# if lambdaOptN > lambdaOpt
			# 	return lambdaOpt
			# else
			# 	lambdaOpt = lambdaOptN;
			# end
		catch exc
			println("Exit Reached")
			println(exc)
			if isa(exc,InterruptException)
				throw(InterruptException())
			end
			if isa(exc,BoundsError)
				println("Halving")
				# lambdaHalf = .5*(lambdaOpt+lambdaPrev);
				println(lambdaOpt)
				println(lambdaPrev)
				newOpt = bisectionSearch(k,lambdaOpt,lambdaPrev,1e-10);
				print("NewOpt is ... ")
				println(newOpt)

				nL= getRhoRepeat(k,newOpt);
				biL = bisectionSearch(k,0,1.5,1e-10);
				 biLInit = getRhoRepeat(k,biL)

				print("OldOpt is ... ")
				println(lambdaOpt)
				print("Bisection is ...")

				println(biL)
				print("newOpt is ...")
				println(newOpt)

				print("opt with init at biL is...")

				println(biLInit)

				minL = min(lambdaOpt,newOpt,biLInit);

				print("Returning...  ")
				println(minL)
				toc()
				return(minL)
			end
			throw(exc)
			
		end

	end
	println(curIter)
	println("end reached")
	toc()
	return(minL)
end



function getRhos(kVect,lInit)
	#eg
	# kVect = [1e2, 1e3, 1e4, 1e5]
	# lInit = [1.5,1.5,1.5,1.5]
	#Calls get RhoRepeat for each pair of kVect and lInit.
	maxTrials = 1;
	pVect = zeros(length(kVect))
	for i = 1:length(kVect)
		best = Inf;
		for j = 1:maxTrials
			best = min(best,getRhoRepeat(kVect[i],lInit[i]))
		end

		pVect[i] = best
	end
	return(pVect)
end

function getRhosBasic(kVect,lInit)
	#eg
	# kVect = [1e2, 1e3, 1e4, 1e5]
	# lInit = [1.5,1.5,1.5,1.5]
	#Calls get RhoRepeat for each pair of kVect and lInit.
	pVect = zeros(length(kVect))
	for i = 1:length(kVect)
		pVect[i] = (getRho(kVect[i],lInit[i],true))[2]
	end
	return(pVect)
end
function getRhosCustom(kVect,lInit,A,B,C)
		#eg
	# kVect = [1e2, 1e3, 1e4, 1e5]
	# lInit = [1.5,1.5,1.5,1.5]
	#Calls get RhoRepeat for each pair of kVect and lInit.
	#RETURNS RHO SQUARED <- IMPORTANT!
	pVect = zeros(length(kVect))
	for i = 1:length(kVect)
		pVect[i] = getRhoRepeatCustom(kVect[i],lInit[i],A,B,C)
	end
	return(pVect)
end 


function getRhoRepeatCustom(k,lambdaInit,A,B,C)	
	#Super hack that results in high precision
	lambdaOpt = getRho(k,lambdaInit)
	minL = lambdaOpt;
	curIter = 1;
	while true
		curIter += 1;
		try
			println("== RESTARTING WITH MOSEK ==")
			lambdaOpt = getRhoCustom(k,lambdaOpt,A,B,C)
			minL = min(lambdaOpt,minL);
			# if lambdaOptN > lambdaOpt
			# 	return lambdaOpt
			# else
			# 	lambdaOpt = lambdaOptN;
			# end
		catch exc
			println("Exit Reached")
			println(exc)
			if isa(exc,InterruptException)
				throw(InterruptException())
			end
			if isa(exc,BoundsError)
				print("returning...")
				println(minL)
				return(minL)
			end
			
		end
	end
	println(curIter)
	println("end reached")
	return(minL)
end




function getRhoCustom(k,lambdaInit,A,B,C)

	#You'll need to create your own feas_point function using MOSEK's SDP Solver
	#(P, t1, t2) = feas_point( k, lambdaInit )
	println("Hi Robert, don't forget to create your own initial point finder!")

	#Remember, the convention here is to not include the affine term in the xInit term
	#And the affine terms come first in A,B, and C 
	(xOpt,lambdaOpt)= methOfCents( A,B,C, xInit, lambdaInit, .1)

	return(lambdaOpt)
end


function getRhoTune(k,lambdaInit,tInit,tMul)
	(P, t1, t2) = feas_point( k, lambdaInit )
	# println((P,t1,t2))
	if length(feas_point(k,lambdaInit)) == 1
		println("feas_point fail")
		return
	end

	x = getCoeff(P);
	(A,B,C) = formulateGevp(k);
	# G = makeBasis(3);
	# println((Mx(G,x),t1,1.0)) 

	# tX = Mx(A,[1;x;t1])

	# println(tX)



	# trueA = makeTrueA(k,P,t1)

	# 	println(trueA)

	# println(tX - trueA)


	#println((A,B,C))

	(xOpt,lambdaOpt)= methOfCents( A,B,C, lambdaInit,[x;t1], tInit,tMul)

	return(lambdaOpt)
end

function parameterTune(tInit,tMul)
	paramMat = zeros(length(tInit),length(tMul));
	for i = 1:length(tInit)
		for j = 1:length(tMul)
			paramMat[i,j] = getRhoTune(1e9,1.5,tInit[i],tMul[j]	);
		end
	end
	return(paramMat)
end

function parameterTuneRepeat(tInit,tMul)
	paramMat = zeros(length(tInit),length(tMul));
	for i = 1:length(tInit)
		for j = 1:length(tMul)
			paramMat[i,j] = getRhoRepeatTune(1e9,1.5,tInit[i],tMul[j]	);
		end
	end
	return(paramMat)
end



function getRhoRepeatTune(k,lambdaInit,tInit,tMul)	
	#Super hack that results in high precision
	lambdaOpt = getRhoTune(k,lambdaInit,tInit,tMul)
	minL = lambdaOpt;
	curIter = 1;
	while true
		curIter += 1;
		try
			println("Hello?")
			lambdaOpt = getRhoTune(k,lambdaOpt,tInit,tMul)
			minL = min(lambdaOpt,minL);
			# if lambdaOptN > lambdaOpt
			# 	return lambdaOpt
			# else
			# 	lambdaOpt = lambdaOptN;
			# end
		catch exc
			println("Exit Reached")
			println(exc)
			if isa(exc,InterruptException)
				throw(InterruptException())
			end
			if isa(exc,BoundsError)
				print("returning...")
				println(minL)
				return(minL)
			end
			
		end
	end
	println(curIter)
	println("end reached")
	return(minL)
end

function bisectionSearch(k,lLow,lHigh,tol)
	tic()
	lMed = NaN;
	while lHigh - lLow > tol
		lMed = .5*(lHigh+lLow);
		if (length(feas_point( k, lMed )) == 1)
			lLow = lMed;
		else
			lHigh = lMed;
		end
		lMed
	end
	if(length(feas_point( k, lLow)) > 1)
		#println([feas_point(k, lLow)])
		toc()
		return(lLow)
	end

	if (isnan(lMed))
		println("Uhoh lMed = NaN because tol is too low")
	else
		if(length(feas_point( k, lMed)) > 1)
			#println(feas_point( k, lMed))
			toc()
			return(lMed)
		end
	end
	toc()
	return(lHigh)
end

function advBisectionSearch(A,B,C,k,lLow,lHigh,tol)
	tic()
	lMed = NaN;
	while lHigh - lLow > tol
		lMed = .5*(lHigh+lLow);
		if (feasTest(A,B,C,lMed,getMosekX( k, lMed )) == false)
			lLow = lMed;
		else
			lHigh = lMed;
		end
		lMed
	end
	curX = getMosekX( k, lLow);
	if(feasTest(A,B,C,lLow,curX) == true)
		#println([feas_point(k, lLow)])
		toc()
		return(lLow,curX)
	end

	if (isnan(lMed))
		println("Uhoh lMed = NaN because tol is too low")
	else
		curX = getMosekX( k, lMed);
		if(feasTest(A,B,C,lMed,curX) == true)
			#println(feas_point( k, lMed))
			toc()
			return(lMed,curX)
		end
	end
	toc()
	return(lHigh,getMosekX(k,lHigh))
end

function advBisectionSearchNest(k,lLow,lHigh,tol)
	(A,B,C) = formulateGevp(k);
	return(advBisectionSearch(A,B,C,k,lLow,lHigh,tol))
end

function feasTest(A,B,C,lambda,x)
	res = true;
	if(length(x) == 1 && isnan(x))
		res = false;
		return(res)
	end
	if(maximum(eig(Mx(A,x),Mx(B,x))[1]) > lambda)
		res = false;
		return(res)
	end
	# try
	# 	chol(lambda*Mx(B,x)-Mx(A,x))
	# catch
	# 	res = false;
	# 	return(res)
	# end
	return(res)
end



function bisectionSearchTime(k,lLow,lHigh,tol)
	totalTime = 0;
	lV = [0 lHigh]
	lMed = NaN;
	while lHigh - lLow > tol
		tic()
		lMed = .5*(lHigh+lLow);
		if (length(feas_point( k, lMed )) == 1)
			lLow = lMed;
		else
			lHigh = lMed;
		end
		lMed
		totalTime += toc();
		lV = [lV; [totalTime lHigh]]


	end
	if(length(feas_point( k, lLow)) > 1)
		#println([feas_point(k, lLow)])
		# toc()
		return(lLow,lV)
	end

	if (isnan(lMed))
		println("Uhoh lMed = NaN because tol is too low")
	else
		if(length(feas_point( k, lMed)) > 1)
			#println(feas_point( k, lMed))
			# toc()
			return(lMed,lV)
		end
	end
	
	return(lHigh,lV)
end


function bisectionSearchVect(kVect,lLow,lHigh,tol)
	lBiVect = zeros(length(kVect),1);
	for i = 1:length(kVect)
		lBiVect[i] = bisectionSearch(kVect[i],lLow,lHigh,tol)
	end
	return(lBiVect)
end

