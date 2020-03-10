
R = Rscript $^ $@

RESDIR := results
DATADIR := input

default: $(addprefix ${RESDIR}/,scenarios.ssv parameters.rda)

${RESDIR}:
	mkdir -p $@

${RESDIR}/scenarios.ssv: scenarios.R ${DATADIR}/scenarios.json | ${RESDIR}
	${R}

${RESDIR}/parameters.rda: parameters.R ${DATADIR}/parameters.json | ${RESDIR}
	${R}

LEN := $(shell cat ${RESDIR}/scenarios.ssv | wc -l | xargs)
SEQ := $(shell seq 1 ${LEN})

CORES ?= 8
SAMPS ?= 1e4

ARRAYID ?= 1

simtar: ${RESDIR}/bp_${ARRAYID}.rds

${RESDIR}/bp_%.rds: run_sims.R simulator.R | ${RESDIR}
	Rscript $^ `sed -n '$*p' ${RESDIR}/scenarios.ssv` ${CORES} ${SAMPS} $@