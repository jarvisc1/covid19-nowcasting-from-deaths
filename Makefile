
R = Rscript $^ $@

RESDIR := results
DATADIR := input

default: $(addprefix ${RESDIR}/,scenarios.ssv parameters.rds)

${RESDIR}:
	mkdir -p $@

${RESDIR}/scenarios.ssv: scenarios.R ${DATADIR}/scenarios.json | ${RESDIR}
	${R}

${RESDIR}/parameters.rda: parameters.R ${DATADIR}/parameters.json | ${RESDIR}
	${R}