snakemake -s scBPS/scBPS.smk \
	    --cores 8 \
	    --config zscore_file=data/magma_result.txt  \
	    scfile=data/singlecell/data.h5ad \
	    outdir=result/scRPS_res \
	    anno=data/singlecell/cell_annotation.tsv > log.txt 2>&1
