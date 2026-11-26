rule hsic:    
    input:
        "results/data/data_n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}.npy"
    output:
        "results/snakemake/output/hsic/n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}.csv"
    conda:
        "knit_env"
    threads: 10
    benchmark:
        "results/snakemake/benchmark/hsic/n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}.txt"
    shell:
        """ 
        python snakemake/methods/run_hsic.py -i {input} -o {output}
        """