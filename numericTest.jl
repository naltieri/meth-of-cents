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


function getRho(k,lambdaInit)
	(P, t1, t2) = feas_point( k, lambdaInit )
	x = getCoeff(P);
	(A,B,C) = formulateGevp(k);
	(xOpt,lambdaOpt)= methOfCents( A,B,C, lambdaInit,[x;t1], .1)
	return(lambdaOpt)
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
	smallNumber = eps()*100;
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



	secondLowerOne = zeros(b1D+2,b1D+2);
	secondLowerOne[b1D+2,b1D+2] = 1;

	C = secondLowerOne;

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