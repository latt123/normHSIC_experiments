import pandas as pd
import numpy as np
import scanpy as sc
import argparse
import time

# Functions to run knit
from knit.run_test import compute_gram_matrix
from knit.compute_statistic import compute_test_statistic_from_gram_matrices
from scipy.stats import chi2
def parse_args():
    parser = argparse.ArgumentParser(description="Run knit with command line")
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
    proportion_train = args.proportion_train
    truncation_X = args.truncation_X
    truncation_e = args.truncation_e
    pvalues  = np.zeros(rep)
    t0 = time.time()
    for r in range(rep):
        X = simulated_data_set[0,r,:,:]
        Y = simulated_data_set[1,r,:,:]
        K_X = compute_gram_matrix(X, X, kernel_X, q_X)  
        K_Y = compute_gram_matrix(Y, Y,  kernel_Y, q_Y)
        

        # First need to split gram into train and test sets
        index_train_set = np.random.uniform(0, 1, X.shape[0]) < proportion_train

        K_X_train_train = K_X[index_train_set,:]
        K_X_train_train = K_X_train_train[:, index_train_set]   

        K_X_test_train = K_X[~index_train_set,:]
        K_X_test_train = K_X_test_train[:, index_train_set]

        K_Y_train_train = K_Y[index_train_set,:]
        K_Y_train_train = K_Y_train_train[:, index_train_set]   

        K_Y_test_test = K_Y[~index_train_set,:]
        K_Y_test_test = K_Y_test_test[:, ~index_train_set]

  
        t, A, B = compute_test_statistic_from_gram_matrices(K_X_train_train, K_X_test_train, K_Y_train_train, K_Y_test_test, projection_X='pls', truncation_X = truncation_X, truncation_e = truncation_e)
        p = chi2.sf(t, df = truncation_X*truncation_e)
        pvalues[r] = p

    n = X.shape[0] 
    print(f'{r + 1}/{rep} repetitions (n = {n}) complete in {round((time.time() - t0)/60, 2)} minutes')        
    mean = np.mean(pvalues < alpha)
    df = pd.DataFrame({'truncation_X' : truncation_X, 'truncation_e': truncation_e, 'mean': mean, 'rep': rep, 'proportion_train': proportion_train}, index = [0])
    df.to_csv(args.output)

# Needs to be at the end of the script
if __name__ == "__main__":
    main()



