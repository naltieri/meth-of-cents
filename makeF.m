function [ F ] = makeF( A,B,C, lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Concatenate lambda*B - A with C
F = zeros(size(A,1)+size(C,1),size(A,1)+size(C,1),size(A,3));
F(1:size(A),1:size(A),:) = lambda*B-A;
F((size(A)+1):end,(size(A)+1):end,:) = C;

end


