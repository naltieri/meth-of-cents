function [ xInit ] = cvxTest( F)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


cvx_begin
	variable x(1)
	variable y(1)
	minimize 1
	F(:,:,1)+F(:,:,2)*x + F(:,:,3)*y == semidefinite(2)
cvx_end

xInit = [x,y]
end


