function [ xSol ] = testValue( A,B,C, lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Concatenate lambda*B-A
D = lambda*B-A;





cvx_begin
	variable x(1)
	variable y(1)
	minimize 1
	D(:,:,1)+D(:,:,2)*x + D(:,:,3)*y  == semidefinite(size(D,1))
	B(:,:,1)+B(:,:,2)*x + B(:,:,3)*y  == semidefinite(size(B,1))
	C(:,:,1)+C(:,:,2)*x + C(:,:,3)*y  == semidefinite(size(C,1))
cvx_end

xSol = [x,y]
end


