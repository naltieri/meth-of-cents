function [ xInit ] = getInit( A,B,C, lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Concatenate lambda*B-A with C
    F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3));
    F(1:size(A),1:size(A),:) = lambda*B-A;
    F((size(A)+1):end,(size(A)+1):end,:) = C;





cvx_begin
	variable x(1)
	variable y(1)
	minimize 1
	F(:,:,1)+F(:,:,2)*x + F(:,:,3)*y  == semidefinite(size(F,1))
cvx_end

xInit = [x,y]
end


