function [ lambdaOpt,xOpt ] = methOfCents( A,B,C, lambda,x, theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = [1;x];

Ax = Mx(A,x);
Bx = Mx(B,x);
Cx = Mx(C,x);

for i = 1:25,
    
    disp(i)
    
    %Lambda Step
    lambda = (1-theta)*max(eigs(Ax,Bx))+ theta*lambda;
    
     disp(lambda)
%    disp(x)
    %Concatenate lambda*B - A with C
    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3));
    F(1:size(A),1:size(A),:) = lambda*B-A;
    F((size(A)+1):end,(size(A)+1):end,:) = C;
    
    %x step
    x = analyticCenter(x,F, .01,10^(-3));
    
    %disp(x)

    %Recompute Ax,Bx,Cx
    Ax = Mx(A,x);
	Bx = Mx(B,x);
	Cx = Mx(C,x);
end

xOpt = x;
lambdaOpt = lambda;
end


