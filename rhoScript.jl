include("numericTest.jl")

function plotRho()
	kvals = logspace(.005,8.5,40);
	kvals = kvals[5:end]
	lInit = repmat([1.5],length(kvals),1);

	# lInit[1] = .5;
	# lInit[2] = .7

	# lInit[3:5] = 1;

	lInit[1] = .7;
	lInit[2] = 1

	lInit[3:5] = 1.3;
	tic()
	rhopt = getRhos(kvals,lInit)
	toc()
	# rhopt = bisectionSearchVect(kvals,0,1.5,1e-10);


	# ITERATIONS PLOT
	tracei = [
	    "x" => kvals,
	    "y" => -1./log(sqrt(rhopt)),
	    "mode" => "lines",
	    "type" => "scatter",
	    "line" => [ "color" => "blue" ]
	]
	nesti = [
	    "x" => kvals,
	    "y" => -1./log(1-1./sqrt(kvals)),
	    "name" => "nesterov",
	    "mode" => "lines",
	    "type" => "scatter",
	    "line" => [ "color" => "black" ]
	]
	layouti = [
	    "title" => "Convergence rate",
	    "xaxis" => ["title" => "Condition number L/m","type" => "log",
	                "autorange" => true,"range" => [1e0, 10^(6)]],
	    "yaxis" => ["title" => "Convergence rate","type" => "log",
	                "autorange" => true]
	]
	datai = [tracei,nesti]
	# datai = tracei;
	responsei = Plotly.plot([datai], ["layout" => layouti]);
	println(responsei["url"])

end

