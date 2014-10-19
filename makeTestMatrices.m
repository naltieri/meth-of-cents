function [ A,B,C ] = makeTestMatrices(d,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1:n,
	a = rand(d);
	b = rand(d);
	c = rand(d);
	a = a+a'-1;
	b = b + b' - 1;
	c = c + c' - 1;
	A(:,:,i) = a;
	B(:,:,i) = b;
	C(:,:,i) = c;
	if i ~= 1,

		A(d,d,i) = -trace(A(:,:,i))+A(d,d,i);
	    B(d,d,i) = -trace(B(:,:,i))+B(d,d,i);
	    C(d,d,i) = -trace(C(:,:,i))+C(d,d,i);
    end

    if i == 1,
    	A(:,:,i) = A(:,:,i)'*A(:,:,i);
    	B(:,:,i) = B(:,:,i)'*B(:,:,i);
    	C(:,:,i) = C(:,:,i)'*C(:,:,i);
end

end


