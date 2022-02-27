# Sofia García del Barrio
# Foundations of Computational Biology and Bioinformatics
# February 25, 2022
# mutation_caller.py

import argparse
import sys
import pysam
import math

# function to calculate probability p(b|A)
def prob_func(b, A):
    # b=base in pileup
    # A=base in genome

    e = 0.1 # quality score

    # calculate probability p(b|A)
    if (b == A):
        prob = 1-e
    else: # if b!=A
        prob = e/3
    return prob

# function to get list of coverages at each pileup
def get_coverage(samfile):
    coverage_list = [] # empty dictionary to store coverage
    base_dict = {} # empty dictionary to store bases and positions
    
    # iterate through pileup
    for pileupcolumn in samfile.pileup():
        coverage = pileupcolumn.n # coverage at each position
        column = pileupcolumn.pos # pileup column number (position)
        coverage_list.append(coverage) # store coverage in each pileup

        base_list = [] # empty list
        for pileupread in pileupcolumn.pileups: # go through each pileup and get base
            if not pileupread.is_del and not pileupread.is_refskip:
                # get each base in pileup
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_list.append(base)
        
        # add pileup sequence to dictionary with position
        base_dict[column] = base_list 
		
    return coverage_list, base_dict
    
# function to get final positions with coverage. Add zeros until equal length
def equal_coverage(normal_coverage_list,cancer_coverage_list):
	while (len(normal_coverage_list) != len(cancer_coverage_list)): 
		if (len(normal_dict) < len(cancer_dict)): # normal is smaller
			normal_coverage_list.append(0)
		else: # cancer is smaller
			cancer_coverage_list.append(0) 
	return normal_coverage_list, cancer_coverage_list
    
# function to calculate log likelihood from a list of pileups
def log_like_func(normal_coverage_list, normal_bases_dict, cancer_coverage_list, cancer_bases_dict):
    prob_dict_n = {} # initialize dictionary to store log probs and column
    possible_genotypes = ['AA','CC','GG','TT','AC','AG','AT','CG','CT','GT'] # list of possible genotypes
    for col in range(len(normal_coverage_list)):
        # get coverages
        normal_coverage = normal_coverage_list[col]
        cancer_coverage = cancer_coverage_list[col]

	# check coverage for one position if both files
        if (normal_coverage < 20) or (cancer_coverage < 20):
            sys.stdout.write("Insufficient coverage at position %d\n" % col)
        else: # if coverage is enough, calculate probability
            normal_bases = normal_bases_dict[col] # list of normal bases at that pileup
            normal_bases=normal_bases # make it upper for comparison

            # loop through each base and calculate probability for one genotype
            prob_genotype_list = [] # initialize list of probs for each genotype
            for genotype in possible_genotypes: # look at one genotype at a time for entire pileup
                n_prob_DG = 1 # initialize probability 
                for normal_base in normal_bases:
                    normal_base=normal_base.upper()
                    n_prob1 = prob_func(normal_base, genotype[0]) # prob(b|A1)
                    n_prob2 = prob_func(normal_base, genotype[1]) # prob(b|A2)  
                    n_prob_base = 0.5*(n_prob1+n_prob2) # prob(b|G)
                    n_prob_DG = n_prob_DG * n_prob_base # prob(D|G) is product of probs for each base in pileup
                prob_genotype_list.append(n_prob_DG) # store probabilities for each base

            max_prob = max(prob_genotype_list) # max of prob(b|G)
            max_genotype=[]
            # store all the alleles with maximum likelihood
            for num in range(10):
            	if (prob_genotype_list[num]==max_prob):
           	    max_genotype.append(possible_genotypes[num])

            log_prob_normal = math.log(max_prob) # log(likelihood) for normal

            # check probability 
            if (log_prob_normal < -50):
                sys.stdout.write('Position %d has ambiguous genotype\n' % col)  
            else: # if it is non-ambigous, check same position in cancer

                cancer_bases = cancer_bases_dict[col] # get cancer bases at this column
                c_prob_DG = 1 # initialize probability
                log_prob_cancer=[]
                for max_allele in max_genotype: # iterate over all the alleles with maximum likelihood
                    for cancer_base in cancer_bases: # look at each base for the optimal genotype obtained from normal
                        c_prob1 = prob_func(cancer_base, max_allele[0]) # prob(b|A1)
                        c_prob2 = prob_func(cancer_base, max_allele[1]) # prob(b|A2)

                        c_prob_base = 0.5*(c_prob1+c_prob2) # prob(b|G)
                        c_prob_DG = c_prob_DG * c_prob_base # prob(D|G)

                    log_prob_cancer.append(math.log(c_prob_DG)) # log likelihood for cancer
                # check probability for all the alleles
                for prob in log_prob_cancer:
                    if (prob <-75): # if any allele is <-75, there is a possible mutation
                        sys.stdout.write('Position %d has a candidate somatic mutation. (Log-likelihood)=%.3f\n' % (col,prob))
                        break # we only want to know if the position is a candidate. If at least one allele gives a candidate, it is enough

# Our script to compare cancer and normal BAM files

# create the parser with argparse
my_parser = argparse.ArgumentParser(description='Somatic Mutation Caller')

# add arguments
my_parser.add_argument('-n')
my_parser.add_argument('-c')
args = my_parser.parse_args()
# read normal file
try:
    samfile_normal=pysam.AlignmentFile(args.n, "rb")
except FileNotFoundError: # check it is in the same directory
    sys.stdout.write("ERROR: -n file not in the directory\n")
except ValueError: # check it is bam format
    sys.stdout.write("ERROR: -n file is not .bam\n")  
else: # if normal file is correct, check cancer
	try:
	    samfile_cancer=pysam.AlignmentFile(args.c, "rb")
	except FileNotFoundError:
	    sys.stdout.write("ERROR: -c file not in the directory\n")
	except ValueError:
	    sys.stdout.write("ERROR: -c file is not .bam\n") 
	else: # only run program if both files are correct

	    # Calculate coverage for both
	    normal_coverage_list, normal_dict = get_coverage(samfile_normal)
	    cancer_coverage_list, cancer_dict = get_coverage(samfile_cancer)
            
	    # Make their lengths equal
	    normal_coverage_list, cancer_coverage_list = equal_coverage(normal_coverage_list, cancer_coverage_list)
            
	    # Calculate the probability and most likely genotype
	    log_like_func(normal_coverage_list, normal_dict, cancer_coverage_list, cancer_dict)        
            
	    # Close the files
	    samfile_normal.close()
	    samfile_cancer.close()
