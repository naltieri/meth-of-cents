include("numericTest.jl")

function plotRho()
	kvals = logspace(1.5,8.5,60);
	kvals = kvals[5:end]
	lInit = repmat([1.5],length(kvals),1);

	# lInit[1] = .5;
	# lInit[2] = .7

	# lInit[3:5] = 1;

	lInit[1] = 1.3;
	lInit[2] = 1.3

	lInit[3:5] = 1.3;
	tic()
	rhopt = getRhos(kvals,lInit)
	# rhopt = getRhosBasic(kvals,lInit)
	toc()
	# rhopt = bisectionSearchVect(kvals,0,1.5,1e-10)';


	# ITERATIONS PLOT
	tracei = [
	    "x" => kvals,
	    "y" => -1./log(sqrt(rhopt)),
	    "mode" => "lines",
	    "type" => "scatter",
	    "line" => [ "color" => "blue" ]
	]
	# traceiB = [
	#     "x" => kvals,
	#     "y" => -1./log(sqrt(rhoptB)),
	#     "mode" => "lines",
	#     "type" => "scatter",
	#     "line" => [ "color" => "red" ]
	# ]
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
	# matVal = readdlm("gevpVals.csv",',')
	# matlabi = [
	# 	"x" => kvals,
	# 	"y" => -1./log(matVal),
	# 	"mode" => "lines",
	#     "type" => "scatter",
	#     "line" => [ "color" => "red" ]
	# ]
	# datai = [tracei,nesti,matlabi]
	# datai = [tracei,traceiB,nesti]
	datai = [tracei,nesti]

	# datai = tracei;
	# responsei = Plotly.plot([datai]);
	responsei = Plotly.plot([datai], ["layout" => layouti]);
	println(responsei["url"])

end

