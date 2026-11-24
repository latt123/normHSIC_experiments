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
        "--proportion_train",
        type=float,
        default=None,
        help="Proportion of dataset to use in train set",
    )
    
    parser.add_argument(
        "--truncation_X",
        type=int,
        default=None,
        help="Truncation for Cov X",
    )


    parser.add_argument(
        "--truncation_e",
        type=int,
        default=None,
        help="Truncation for Cov E",
    )

    parser.add_argument(
        "--q_X",
        type=float,
        default=None,
        help="Bandwith Quantile for kernel X",
    )

    parser.add_argument(
        "--q_Y",
        type=float,
        default=None,
        help="Bandwith Quantile for kernel Y",
    )

    parser.add_argument(
        "-o", "--output", type=str, default=None, help="Output filename for results"
    )

    return parser.parse_args()

# Parameters for gram
kernel_X = 'gaussian'
kernel_Y = 'gaussian'


# For test
alpha = 0.05

# quantiles to run over
q_X_list = [2**(-1/5), 2**(-1/3), 2**(-1), 2**(-3), 2**(-5)]
q_Y_list = [2**(-1/5), 2**(-1/3), 2**(-1), 2**(-3), 2**(-5)]


def main():
    args = parse_args()
    simulated_data_set = np.load(args.input)
    # q_X = args.q_X
    # q_Y = args.q_Y
    rep = simulated_data_set.shape[1]
    pvalues  = np.zeros((rep, len(q_X_list)*len(q_Y_list)))
    bandwidth = []  

    for q_X in q_X_list:
        for q_Y in q_Y_list:
            bandwidth.append(f'q_X: {q_X}, q_Y: {q_Y}')

    t0 = time.time()
    for r in range(rep):
        X = simulated_data_set[0,r,:,:]
        Y = simulated_data_set[1,r,:,:]
        q =0
        for q_X in q_X_list:
            for q_Y in q_Y_list:

                K_X = compute_gram_matrix(X, X, kernel_X, q_X)  
                K_Y = compute_gram_matrix(Y, Y,  kernel_Y, q_Y)
                
                # Run HSIC
                p_hsic = compute_pval_HSIC_gamma_from_grams(K_X, K_Y)
                pvalues[r, q] = p_hsic
                q += 1
    n = X.shape[0] 
    print(f'{r + 1}/{rep} repetitions (n = {n}) complete in {round((time.time() - t0)/60, 2)} minutes')        
    df = pd.DataFrame(pvalues, columns = bandwidth, index = np.arange(rep))
    df.to_csv(args.output)

# Needs to be at the end of the script
if __name__ == "__main__":
    main()


