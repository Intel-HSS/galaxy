#!/bin/bash
PYTHON_2_7_PATH=/opt/python-2.7/
export http_proxy=http://sparkdmz1:4000
export https_proxy=${http_proxy}
export ftp_proxy=${http_proxy}
export GALAXY_RUN_ALL=1
#Creating env variables for location of exe/jar files, to be used by xml/py files.
GENOMICS_DIR=/opt/Genomics
export BWA_DIR=$GENOMICS_DIR/ohsu/dnapipeline/bwa-0.7.12
export BWA_EXE_PATH=$BWA_DIR/bwa
#export SAMTOOLS_DIR=$GENOMICS_DIR/ohsu/dnapipeline/samtools-0.1.19
#export SAMTOOLS_EXE_PATH=$SAMTOOLS_DIR/samtools
export SAMTOOLS_v020_EXE_PATH=$GENOMICS_DIR/ohsu/dnapipeline/samtools/samtools/samtools
export PICARD_PATH=$GENOMICS_DIR/ohsu/dnapipeline/picard-tools-1.135/dist
export GATK_JAR_PATH=/opt/gsa-unstable/target/GenomeAnalysisTK.jar
export CP_PATH=/cluster_share/cp_pipeline
export LABKEY_USERNAME=ccc@ccc.com
export LABKEY_PASSWORD=Intel123
export OMERO_IMAGE_SELECTION=/cluster_share/omero
export OMICS_QUERY_UI_TOOL=/cluster_share/omics-query-tool
export SEATTLESEQ_GETANNOTATION_JAR_PATH=/opt/Genomics/annotation/getAnnotation/project/
export SEATTLESEQ_WRITEGENOTYPE_JAR_PATH=/opt/Genomics/annotation/writeGenotype/project/
export MUTECT_JAR_PATH=$GENOMICS_DIR/ohsu/dnapipeline/mutect-1.1.7.jar
export BOWTIE_TOOLS_EXE_PATH=$GENOMICS_DIR/ohsu/rnapipeline/bowtie2-2.1.0/
export TOPHAT_EXE_PATH=$GENOMICS_DIR/ohsu/rnapipeline/tophat-2.0.14/
export CUFFLINKS_TOOLS_EXE_PATH=$GENOMICS_DIR/ohsu/rnapipeline/cufflinks-2.2.1/
export CUFFDIFF_SCRIPT_PATH=/opt/Genomics/ohsu/rnapipeline
export SNPEFF_JAR_PATH=$GENOMICS_DIR/annotation/snpEff/
export ONCOTATOR_PATH=$GENOMICS_DIR/annotation/Oncotator/oncotator-1.5.1.0/virt_test_env/bin/
export ONCOTATOR_DB_SOURCE=$GENOMICS_DIR/annotation/Oncotator/oncotator_v1_ds_Jan262014/
export BREAKDANCER_PATH=$GENOMICS_DIR/ohsu/dnapipeline/breakdancer-1.4.5-unstable/
export PINDEL_PATH=$GENOMICS_DIR/ohsu/dnapipeline/pindel/
export SPARK_HOME=/opt/cloudera/parcels/CDH/lib/spark
export ADAM_HOME=/opt/adam/adam-distribution-0.16.0
export SAMBAMBA_PATH=$GENOMICS_DIR/ohsu/dnapipeline/sambamba-0.5.8/
export BAMUTIL_PATH=$GENOMICS_DIR/ohsu/dnapipeline/bamUtil_1.0.13/bamUtil/bin/
export ELPREP_PATH=/home/wshands/elprep/elprep/elprep-v2.34/
export STAR_PATH=/home/aprabh2/STAR-STAR_2.4.2a/bin/Linux_x86_64/
export PATH=$STAR_PATH:$ELPREP_PATH:$BAMUTIL_PATH:$SAMBAMBA_PATH:$PINDEL_PATH:$BREAKDANCER_PATH:$TOPHAT_EXE_PATH:$CUFFLINKS_TOOLS_EXE_PATH:$BOWTIE_TOOLS_EXE_PATH:$GENOMICS_DIR:$BWA_DIR:$SAMTOOLS_v020_EXE_PATH:$PICARD_PATH:$PYTHON_2_7_PATH/bin:$SPARK_HOME/bin:$ADAM_HOME/bin:$PATH
export NSLOTS=16
export PYTHON_EGG_CACHE=.eggs_cache
export BIOBLEND_REPO_DIR=/cluster_share/Galaxy/api/bioblend
export FASTQC_PERL_SCRIPT=/opt/Genomics/ohsu/dnapipeline/FastQC/fastqc
export BCFTOOLS=/home/karthikg/broad/non_variant_db/bcftools/bcftools
export TILEDB_IMPORT_EXE=/home/karthikg/broad/non_variant_db/variantDB/TileDB/example/bin/gt_example_loader
export GCC49_PREFIX_PATH=/opt/gcc-4.9.1
export SCATTER_GATHER_JAR_PATH=$GENOMICS_DIR/ohsu/dnapipeline/ScatterGather-1.0.0.jar
