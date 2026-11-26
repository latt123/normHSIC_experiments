import pandas as pd
import numpy as np
import scanpy as sc
import argparse
import time


# Functions to run HSIC
from knit.run_test import compute_gram_matrix
from knit.compute_statistic import centre_rows_and_columns
from scipy.stats import gamma


def compute_pval_HSIC_gamma_from_grams(K_X, K_Y):
    # Centre grams
    K_X_centered = centre_rows_and_columns(K_X)
    K_Y_centered = centre_rows_and_columns(K_Y)

    #  Compute HSIC Statistic
    n = K_X.shape[0]
    T = (1/n**2)*np.trace(K_X_centered@K_Y_centered)

    ## Compute expectation T
    trace_Cov_X = (1/n)*np.trace(K_X_centered)
    trace_Cov_Y = (1/n)*np.trace(K_Y_centered)
    expectation_T = (1/n)*trace_Cov_X*trace_Cov_Y

    ## Compute variance T
    norm_Cov_X = (1/n**2)*np.trace(K_X_centered@K_X_centered)
    norm_Cov_Y = (1/n**2)*np.trace(K_Y_centered@K_Y_centered)
    variance_T = ((2*(n-4)*(n-5))/((n-3)*(n-2)*(n-1)*n))*norm_Cov_X*norm_Cov_Y

    ## Compute pvalue from alpha approximation
    alpha = expectation_T**2/variance_T
    beta = n*variance_T/expectation_T
    p = 1 - gamma.cdf(n*T, a = alpha, scale = beta)
    
    return p



def parse_args():
    parser = argparse.ArgumentParser(description="Run hsic with command line")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-o", "--output", type=str, default=None, help="Output filename for results"
    )

    return parser.parse_args()




# Parameters for gram
kernel_X = 'gaussian'
kernel_Y = 'gaussian'
q_X = 0.5
q_Y = 0.5

# For test
alpha = 0.05


def main():
    args = parse_args()
    simulated_data_set = np.load(args.input)
    rep = simulated_data_set.shape[1]
    pvalues  = np.zeros((1, rep))
    truncation_Xs = []
    truncation_es = [] 
    t0 = time.time()
    for r in range(rep):
        X = simulated_data_set[0,r,:,:]
        Y = simulated_data_set[1,r,:,:]
        K_X = compute_gram_matrix(X, X, kernel_X, q_X)  
        K_Y = compute_gram_matrix(Y, Y,  kernel_Y, q_Y)
        
        # Run HSIC
        p_hsic = compute_pval_HSIC_gamma_from_grams(K_X, K_Y)
        pvalues[0,r] = p_hsic

        if (r+1) % 100 == 0:
            print(f' {r + 1}/{rep} repetitions complete')
            print(f'Time for 100 rep = {round((time.time() - t0)/60, 2)} minutes')
            t0 = time.time()


    means = np.mean(pvalues < alpha, axis = 1)
    df = pd.DataFrame({'mean': means, 'rep': rep*np.ones(len(means))})
    df.to_csv(args.output)

# Needs to be at the end of the script
if __name__ == "__main__":
    main()


