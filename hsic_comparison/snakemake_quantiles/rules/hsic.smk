rule hsic:    
    input:
        "results/data/data_n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}.npy"
    output:
        "results/snakemake_quantiles/output/hsic/n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}_proptrain_{proportion_train}_S_{truncation_X}_T_{truncation_e}.csv"
    conda:
        "knit_env"
    threads: 10
    params:
        proportion_train="{proportion_train}",
        truncation_X="{truncation_X}",
        truncation_e="{truncation_e}"

    benchmark:
        "results/snakemake_quantiles/benchmark/hsic/n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}_proptrain_{proportion_train}_S_{truncation_X}_T_{truncation_e}.txt"
    shell:
        """ 
        python snakemake_quantiles/methods/run_hsic.py -i {input} --proportion_train {params.proportion_train} --truncation_X {params.truncation_X} --truncation_e {params.truncation_e} -o {output}
        """
