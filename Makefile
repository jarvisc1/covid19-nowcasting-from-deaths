
R = Rscript $^ $@

scenarios.csv: scenarios.R scenarios.json
	${R}

parameters.rda: parameters.R parameters.json
	${R}