
R = Rscript $^ $@

scenarios.rds: scenarios.R scenarios.json
	${R}