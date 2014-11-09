include("numericTest.jl")

function plotRho()
	kvals = logspace(.005,8,20);
	kvals = kvals[3:end]
	lInit = repmat([1.5],length(kvals),1);

	lInit[1] = .5;
	lInit[2] = .7

	lInit[3:5] = 1;

	rhopt = getRhos(kvals,lInit)


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
	                "autorange" => true,"range" => [1e0, 1e6]],
	    "yaxis" => ["title" => "Convergence rate","type" => "log",
	                "autorange" => true]
	]
	datai = [tracei,nesti]
	responsei = Plotly.plot([datai], ["layout" => layouti, "fileopt" => "overwrite"]);
	println(responsei["url"])

end