function [ xOpt ] = analyticCenter(x0,F, alpha, tol )
xPrev = inf;
x = [x0];



Fx = zeros(size(F,1),size(F,2));
for i = 1:size(F,3),
    Fx = Fx + F(:,:,i)*x(i);
end

if(cond(Fx) == Inf),
	disp('Fx is singular...quitting')
	return
end

for i = 2:size(F,3),
    g(i-1) = - trace(pinv(Fx)*F(:,:,i));
end
g = g';

for i = 2:size(F,3),
    for j = 2:size(F,3),
        H(i-1,j-1) = trace(pinv(Fx) *F(:,:,i) * pinv(Fx)* F(:,:,j));
    end
end

objHist = log(det(pinv(Fx)));
gHist = norm(g,2);

%Second condition is a hack. Happens when numerical underflow
%occurs in updating x
%Replacing 0 with 10^(-10) creates infeasible results. 
while norm(g,2) > tol && norm(xPrev-x,2) > 0,
%     if size(gHist,2) > 1
%         disp(gHist(end-1) - gHist(end) )
%     end
%     if size(gHist,2) > 1 && gHist(end-1) - gHist(end) < 10^(-13),
%         disp(norm(g,2))
%     end
    xPrev = x;
	x(2:end) = x(2:end) - alpha * pinv(H) *g;
    %disp(log(det(pinv(Fx))))

	Fx = zeros(size(F,1),size(F,2));
	for i = 1:size(F,3),
	    Fx = Fx + F(:,:,i)*x(i);
	end
	if(cond(Fx) == Inf),
		disp('Fx is singular...quitting')
		return
	end


	for i = 2:size(F,3),
	    g(i-1) = - trace(pinv(Fx)*F(:,:,i));
    end
    
    
	


	for i = 2:size(F,3),
	    for j = 2:size(F,3),
	        H(i-1,j-1) = trace(pinv(Fx) *F(:,:,i) * pinv(Fx)* F(:,:,j));
	    end
    end
    gHist = [gHist,norm(g,2)];
	objHist = [objHist, log(det(pinv(Fx)))];


end

xOpt = x;

