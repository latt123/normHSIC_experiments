
rule SPARKX:    
    input:
        "results/data/{slide}.h5ad",
    output:
        "results/output/SPARKX/{slide}.csv"
    threads: 10
    conda:
        "SPARKX_env"

    benchmark:
        "results/benchmark/SPARKX/{slide}.csv"
    shell:
        """ 
        Rscript snakemake/methods/run_SPARKX.R -i {input} -o {output}
        """