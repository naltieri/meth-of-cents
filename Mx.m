function [ out ] = Mx( M,x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



out = zeros(size(M,1),size(M,2));
for i = 1:size(M,3),
   out = out + M(:,:,i)*x(i);
end

end

