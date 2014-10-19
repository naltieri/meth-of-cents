function [Z] = testCenter( F,gridX,gridY)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


	

for i = 1:length(gridX)
	for j = 1:length(gridY)
		x = gridX(i);
		y = gridY(j);
		z = [1;x;y];
		if (isreal(log(det(inv(Mx(F,z)))))),
			Z(i,j) = log(det(inv(Mx(F,z))));
		else
			Z(i,j) = Inf;
		end
	end
end

surf(gridX,gridY,Z);



end


