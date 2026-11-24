import numpy as np
from scipy.stats import t, laplace

# Functions for simulation
def generate_datasets_under_alternative_HSIC(d, m, scale_pi_over_four = 0, dist_1 = 'g', dist_2 = 'g'):

    """
    Using notation constitent with HSIC.
    Parameters
        - m - sample size
        - d - dimension (of X an Y)
        - scale_pi_over_four, rotation angle = scale_pi_over_four * pi/4
        - mixture_mean_1, mixture_mean_2 - means of the symmetric mixtures in source_1 and source_2 (these aren't specified in HSIC paper)
    
        
    Returns
        - Variables X and Y of dimension d where dependence between X and Y is specified by scale_pi_over_four.
    """

    ## STEP 1
    ## Generate univariate random varibales from 'symmetric mixture 2 gauss multimodal' distribution with zero mean and variance 1.
    source_1 = simulate_source(dist_1, m)
    source_2 = simulate_source(dist_2, m)

    ## STEP 2
    ## Mix variables
    theta = (np.pi/4)*scale_pi_over_four
    rotation_matrix = np.array([np.cos(theta), -np.sin(theta), np.sin(theta), np.cos(theta)]).reshape(2,2)
    mix_variables = np.concatenate((source_1.reshape(-1,1),source_2.reshape(-1,1)), axis = 1) @ rotation_matrix
    X = mix_variables[:,0].reshape(-1,1)
    Y = mix_variables[:,1].reshape(-1,1)

    ## STEP 3
    ## Add d-1 dimensional gaussian noise
    if d > 1:
        random_noise = np.random.normal(0, 1, (m, d-1))
        X = np.concatenate((X, random_noise), axis = 1)
        random_noise = np.random.normal(0, 1, (m, d-1))
        Y = np.concatenate((Y, random_noise), axis = 1)

    ## STEP 4
    # Multiply each vector by a random d-dimensional orthogonal matrix
    orthogonal_matrix = np.random.normal(0,1,(d,d))
    X = X@orthogonal_matrix
    Y = Y@orthogonal_matrix
    return X, Y


def simulate_source(distribution, m):
    """
    Simulate the source distribution. Using the labels from table 3 of 'Kernel Methods for Measuring Independence', Gretton 2005.

    Parameters
        - distribution - 'a', 'b', 'c', 'g', 'h', 'i'. 
            - a - t-distribution with 3 degrees of freedom
            - b - laplace distribution
            - c - uniform distribution
            - d - symmetric mixture of multimodal gaussians
            - e - symmetric mixture of 'transitional' gaussians
            - f - symmetric mixture of 'unimodal' gaussians
        - m - sample size

    """

    if distribution == 'a':
        source = t.rvs(df = 3, size = (m))


    elif distribution == 'b':
        source = laplace.rvs(size = (m))

    elif distribution == 'c':
        source = np.random.uniform(size = (m))


    elif distribution == 'g':
        mean = 3
        modes = 2*np.random.binomial(1,1/2,m) - 1
        source = np.random.normal(mean,1,(m))*modes
        
    elif distribution == 'h':
        scale = 3
        mean = 1
        modes = np.random.binomial(1,1/2,m)
        source = np.random.normal(-mean,1,(m))*modes + np.random.normal(mean,scale,(m))*(1- modes)


    elif distribution == 'i':
        scale = 2
        modes = np.random.binomial(1,1/2,m)
        source = np.random.normal(0,1,(m))*modes + np.random.normal(0,scale,(m))*(1- modes)

    else:
        print('Distribution not supported')


    source = (source - np.mean(source))/np.std(source)
    return source


# Simulation parameters
n_list = [128, 512, 2048]
d_list = [1, 2, 4]
dist_list = ['a', 'b', 'c', 'g','i']
scale_pi_over_four_list = [0, 0.2, 0.4, 0.6, 0.8, 1]
rep = 500 # number of simulations to run

def main():
    for n in n_list:
        for d in d_list:
            for dist in dist_list:
                for scale_pi_over_four in scale_pi_over_four_list:
                    generated_data = np.zeros((2,rep, n, d)) # 2 for X and Y, rep for the number of repetitions, n for sample size, d for dimension
                    for r in range(rep):
                        X, Y = generate_datasets_under_alternative_HSIC(d, n, scale_pi_over_four , dist_1 = dist, dist_2 = dist)
                        generated_data[0,r,:,:] = X
                        generated_data[1,r,:,:] = Y
                    np.save(f'results/data/data_n_{n}_d_{d}_dist_{dist}_angle_{scale_pi_over_four}.npy', generated_data)


if __name__ == "__main__":
    main()


