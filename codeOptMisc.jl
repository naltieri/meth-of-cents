
		pinvFx = pinv(Fx)
		# pinvFxF = Array(Float64,size(F,1),size(F,2),size(F,3)-1)
		# println(size(pinvFxF))
		# for i = 2:size(F,3)
		# 	# println(i)
		# 	# println(size(pinvFxF[:,:,i]))
		# 	# println(size(pinvFx*F[:,:,i]))
		# 	pinvFxF[:,:,i-1] = pinvFx*F[:,:,i]
		# end

		g = zeros(size(F,3)-1,1);
		for i = 2:size(F,3)
			# g[i-1] = -trace(pinvFxF[:,:,i])
		    g[i-1] = -sum(sum(pinvFx.*F[:,:,i]))
		end
		 

		H = zeros(size(F,3)-1,size(F,3)-1)
		for i = 2:size(F,3)
			pinvFxFiPinvFx  = pinvFxF[:,:,i]*pinvFx;
		    for j = 2:size(F,3)
		        # H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
		        H[i-1,j-1] = sum(sum(pinvFxFiPinvFx.*F[:,:,j]));
		    end
		end


	pinvFxF = Array(Float64,size(F,1),size(F,2),size(F,3)-1)
		for i = 2:size(F,3)
			pinvFxF[:,:,i-1] = pinvFx*F[:,:,i];
		end

		g = zeros(size(F,3)-1,1);		 
		H = zeros(size(F,3)-1,size(F,3)-1)
		for i = 2:size(F,3)
			pinvFxFi = pinvFx*F[:,:,i];
			# g[i-1] = -trace(pinvFxFi)
			# println(g[i-1])
			g[i-1] = -trace(pinvFxF[:,:,i-1]);
			# println(g[i-1]- (-trace(pinvFxFi)))
			pinvFxFiPinvFx  = pinvFxFi*pinvFx;
			# println(norm(pinvFxF[:,:,i-1]*pinvFx - pinvFxFiPinvFx))
		    for j = 2:size(F,3)
		    	# 	Fj = F[:,:,j];
		        # H[i-1,j-1] = sum(sum(pinvFxFiPinvFx .* Fj'));

		        # H[i-1,j-1] = trace(pinvFxFiPinvFx * F[:,:,j]);
		        H[i-1,j-1] = trace(pinvFxF[:,:,i-1]* pinvFxF[:,:,j-1]);
		        # println(H[i-1,j-1]-trace(pinvFxFiPinvFx * F[:,:,j]))
		        # println("IM HAPPY")
		        # println(issym(F[:,:,i]))
		        # println(issym(Fx))
		        # println(issym(pinvFxFiPinvFx))
		        # println(issym(F[:,:,j]))
		 
		        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(sum(pinvFxFiPinvFx .* Fj')))
		        #H[i-1,j-1] = sum(computeDiagAB(pinvFxFiPinvFx , F[:,:,j]));
		        # println(trace(pinvFxFiPinvFx * F[:,:,j])-sum(pinvFxFiPinvFx[:] .* F[:,:,j][:]))
		    end
		end
