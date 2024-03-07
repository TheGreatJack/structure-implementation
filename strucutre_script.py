import numpy as np
from scipy.stats import dirichlet, multinomial


## Define functions

def get_vals_and_alpha_prior_full_alleles(data_mat,data_general):
    summary_locus_all_samples = np.unique(data_general) 
    summary_locus_pop = np.unique(data_mat, return_counts = True)
    
    # Create new column for unique values to save the population counts
    to_concatenate = np.zeros(summary_locus_all_samples.shape[0])
    summary_locus_all_samples = np.column_stack((summary_locus_all_samples, to_concatenate))

    if summary_locus_pop[0].shape[0] > 0: # If population wasnt sampled, no alleles would be found, this tests if this array is not empty
        # Add to the previous column the appropiate number of allele counts per locus per population
        for allele in range(summary_locus_all_samples.shape[0]):
            if summary_locus_all_samples[allele,0] in summary_locus_pop[0]:
                pop_allele_index = np.where(summary_locus_pop[0] == summary_locus_all_samples[allele,0])

                summary_locus_all_samples[allele,1] = summary_locus_pop[1][pop_allele_index[0][0]]

    # Generate values and prior
    vals = summary_locus_all_samples[:,0]
    prior = summary_locus_all_samples[:,1] + 1/len(vals) ## Definition of prior in original paper?, counts + some alpha that always sums to 1

    return vals,prior

# Create dictionary to find the index (samples) assignated to a population
def get_indexes_of_population(arr):
    unique_values = np.unique(arr)
    indexes_dict = {}
    for value in unique_values:
        #print(np.where(arr == value))
        indexes_dict[value] = list(np.where(arr == value)[0])
    return indexes_dict

# Sample population assignments, as a initial starting point assuming a uniform distr
# Assure that all populations are sampled to not generate problems down the line. Just resample if any pop is missing
def initial_pop_assignment(K,Z):
    all_pops_sampled = False
    while not(all_pops_sampled):
        for j in range(X.shape[0]):
            Z[j] = np.random.multinomial(1, np.ones(K) / K).argmax()
        if not(False in np.isin([x for x in range(K)],Z)):
            all_pops_sampled = True
    return Z

# Sample allele freqs on a per pooulation and per loci basis
def sample_allele_freqs(X,K, indexes_dict):
    # List that will save the dictionary- allele frequencies of the loci per population
    freq_dict_per_pop = []
    for pop in range(K):
        print(pop)
        print(indexes_dict[pop])
        pop_samples = X[indexes_dict[pop]]
        # Define population dictionary
        pop_dict = {}
        # Iterate over each loci within the population
        for locus in range(X.shape[1]):
            values,prior = get_vals_and_alpha_prior_full_alleles(pop_samples[:,locus,:],X)
            
            # Sample allele freqs from dirichlet prior from loci
            sampled_allele_freqs_L = np.random.dirichlet(prior)

            # Save allele freqs in population dictionary
            allele_dict = {}
            for value, freq in zip(values, sampled_allele_freqs_L):
                allele_dict[value] = freq
            pop_dict[locus] = allele_dict

        # Save population dictionary in freq_dict_per_pop
        freq_dict_per_pop.append(pop_dict)
    return freq_dict_per_pop


# Generate population assignment probabilites from estimated allele frequencies. 
# Estimate the "probability" of the full genotype of an individual in a population and ponder 

def generate_pop_probabilities(X, K, freq_dict_per_pop):
    sample_probs = []
    for sample in range(X.shape[0]):
        likelihoods = []
        for pop in range(K):
            pop_dict = freq_dict_per_pop[pop]
            product = 1
            for locus in range(X.shape[1]):
                genotype = X[sample,locus]
                allele_0_freq = pop_dict[locus][genotype[0]]
                allele_1_freq = pop_dict[locus][genotype[1]]
                product *= allele_0_freq*allele_1_freq
            likelihoods.append(product)
        probs = np.array(likelihoods)/sum(likelihoods)
        sample_probs.append(probs)
    sample_probs = np.array(sample_probs)
    return sample_probs


# Assign population from population probabilities estimated before
def assign_pop_from_probs(Z, K, X, sample_probs):
    for sample in range(X.shape[0]):
        Z[j] = np.random.multinomial(1, sample_probs[sample]).argmax()
    return Z


#------------------------------------------------------------------------------------------

# Define the model parameters
K = 3  # Number of populations
L = 5  # Number of loci

# Define the data
X = np.array([[[1, 1],
  [3, 3],
  [4, 4],
  [7, 7],
  [9, 9]],
 [[1, 2],
  [3, 3],
  [4, 6],
  [7, 7],
  [9, 9]],
 [[2, 2],
  [3, 3],
  [4, 5],
  [8, 8],
  [9, 9]],
 [[2, 2],
  [3, 3],
  [4, 5],
  [8, 7],
  [9, 10]]])  # Genotype data

# Initialize the MCMC chain
iterations = 1000

Z_accum = []
P_accum = []

Z = np.zeros(X.shape[0], dtype=int)  # Population assignments
# Sample population assignments, as a initial starting point assuming a uniform distr
Z = initial_pop_assignment(K,Z)
Z_accum.append(Z)

# Run the MCMC chain

for i in range(iterations):
    # Get index from population

    indexes_dict = get_indexes_of_population(Z)

    # Sample allele freqs on a per pooulation and per loci basis

    freq_dict_per_pop = sample_allele_freqs(X,K, indexes_dict)# Sample allele freqs on a per pooulation and per loci basis
    P_accum.append(freq_dict_per_pop)

    # Implement Z assignation from new allele freq matrix

    sample_probs = generate_pop_probabilities(X, K, freq_dict_per_pop)
    Z = assign_pop_from_probs(Z, K, X, sample_probs)

print(Z_accum)
print(P_accum)