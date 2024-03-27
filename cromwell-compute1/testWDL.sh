export WORK_DIR="/storage1/fs1/mgriffit/Active/griffithlab/gc2596/e.schmidt/run_cromwll_compute1/compute1-wdl/rnaseq-wdl-test"


bsub -q oncology -G compute-oncology -g /mgriffit/wdl -M 22000000 -R 'select[mem>22000] rusage[mem=22000]'\
 -oo $WORK_DIR/logs/out.stdout\
 -eo $WORK_DIR/logs/out.stderr\
 -a 'docker(ghcr.io/genome/genome_perl_environment:compute1-58)'\
 /bin/bash $WORK_DIR/runCromwellWDL.sh\
 --cromwell_config $WORK_DIR/cromwell.config.wdl\
 --sample test \
 --wdl $WORK_DIR/analysis-wdls/definitions/rnaseq_star_fusion.wdl \
 --imports $WORK_DIR/analysis-wdls/definitions/workflows.zip\
 --yaml $WORK_DIR/ranseq-test.yaml \
 --results $WORK_DIR/out\
 --cromwell_jar /storage1/fs1/mgriffit/Active/griffithlab/common/cromwell-jars/cromwell-71.jar\
 --cromwell_server_mem 10g --cromwell_submit_mem 10g


cd $WORK_DIR/