function doTheThing()
	powers = 1:100;
	a = rand(10,10)-.5;
	# a = a'a;

	cList = zeros(length(powers),1)

	for i = 1:length(powers)
		s = a
		for j = 2:i
			s += a^2
		end
		# s = a^i
		cList[i] = cond(s)
	end

	Plotly.plot([1:length(cList) cList])

end
z = zeros(size(rand(3,3)))
b = [rand(3,3) rand(3,3) rand(3,3);rand(3,3) rand(3,3) z;rand(3,3) rand(3,3) rand(3,3)]