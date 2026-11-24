import pandas as pd
import numpy as np
import scanpy as sc
import argparse
import time

# Functions to run normHSIC
from knit.run_test import compute_gram_matrix
from knit.compute_statistic import compute_test_statistic_from_gram_matrices
from scipy.stats import chi2
def parse_args():
    parser = argparse.ArgumentParser(description="Run normHSIC with command line")
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



# Parameters for normHSIC
truncation_X_list = [1, 2, 4, 8]
truncation_e_list = [1, 2, 4, 8]
proportion_train = 0.3

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
    pvalues  = np.zeros((len(truncation_X_list)*len(truncation_e_list), rep))
    truncation_Xs = []
    truncation_es = [] 
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

        # Then run for each compina
        i = 0
        for truncation_X in truncation_X_list:
            for truncation_e in truncation_e_list:
                t, A, B = compute_test_statistic_from_gram_matrices(K_X_train_train, K_X_test_train, K_Y_train_train, K_Y_test_test, projection_X='pls', truncation_X = truncation_X, truncation_e = truncation_e)
                p = chi2.sf(t, df = truncation_X*truncation_e)
                pvalues[i,r] = p
                if r==0:
                    truncation_Xs.append(truncation_X)
                    truncation_es.append(truncation_e)
                i+= 1

        if (r+1) % 100 == 0:
            print(f' {r + 1}/{rep} repetitions complete (all truncations)')
            print(f'Time for 100 rep = {round((time.time() - t0)/60, 2)} minutes')
            t0 = time.time()


    means = np.mean(pvalues < alpha, axis = 1)
    df = pd.DataFrame({'truncation_X' : truncation_Xs, 'truncation_e': truncation_es, 'mean': means, 'rep': rep*np.ones(len(means))})
    df.to_csv(args.output)

# Needs to be at the end of the script
if __name__ == "__main__":
    main()



