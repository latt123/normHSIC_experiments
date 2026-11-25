
rule normhsic:    
    input:
        "results/data/{slide}.h5ad",
    output:
        "results/output/normhsic/{slide}.csv"
    conda:
        "knit_env"
    threads: 10

    benchmark:
        "results/benchmark/normhsic/{slide}.csv"
    shell:
        """ 
        python snakemake/methods/run_normhsic.py -i {input} -o {output}
        """