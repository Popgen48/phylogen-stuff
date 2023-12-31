import pysam
import pandas as pd

# Get the population of each animal_i
pop_file = './examples/DRZdivAlp1M.uTiere.TxT'
with open(pop_file, 'r') as file:
    pop_lines = file.readlines()
    pop_dict = {}
    n_pop = {}
    for line in pop_lines:
        parts = line.split()
        animal_id = parts[1]
        pop = parts[2]
        if pop not in n_pop:
            n_pop[pop] = 0
        pop_dict[animal_id] = pop
        n_pop[pop] += 1
pop_keys = list(n_pop.keys())
pop_keys = sorted(pop_keys, key=lambda x: int(x))
n_pop = {k: n_pop[k] for k in sorted(n_pop, key=lambda x: int(x))}

sample_alleles = {}
spa_counts = {}

for key in pop_dict.keys():
    sample_alleles[key] = [[], []]
    spa_counts[key] = [0] * len(n_pop.keys())
my_keys = list(sample_alleles.keys())

def compare(dict, spa_criterion):
    # my_keys = list(dict.keys())
    for key, values in dict.items():
        pop1 = pop_dict[key]
        for allele in values:
            count = 0
            idx = []
            pops_considered = []
            for key2, values2 in dict.items():
                pop2 = pop_dict[key2]
                if count > spa_criterion:
                    break
                if pop2 == pop1:
                    continue
                elif allele in values2:
                    if pop2 in pops_considered:
                        idx.append(key2) 
                    else:
                        pops_considered.append(pop2)
                        count += 1
                        idx.append(key2)
            if len(pops_considered) > 0 and len(pops_considered) <= spa_criterion:
                for i in pops_considered:
                    spa_counts[key][pop_keys.index(i)] += 1 
                    
def flush(dict):
    for key, values in dict.items():
        dict[key] = [[], []]

def get_spa(vcf_path, threshold, group_size, spa_criterion):
    '''
    vcf_path(str): path to the vcf file
    threshold(int): The max number of base pairs i.e. max distance between two SNPs to be considered
    group_size(int): Max grouping size of the haplotypes i.e. block size
    spa_criterion(int): The max number of distinct populations (apart from itself) in which an allele can be present to be classified as a semi-private allele
    
    returns a dataframe containing the spa counts (adjusted by the pop size) for each animal along with its population id
    '''

    vcf = pysam.VariantFile(vcf_path)
    previous_record = None
    previous_pos = 0
    group_i = 0
    num_blocks = 0
            
    for i, record in enumerate(vcf):
        gt_map = {0: record.ref, 1: record.alts[0], None: '-'}
        current_pos = record.pos
        if i == 0:
            previous_pos = record.pos   
        if group_i == group_size:
            group_i = 0
            for key, values in sample_alleles.items():
                new_allele1 = ''.join(values[0])
                values[0] = new_allele1
                new_allele2 = ''.join(values[1])
                values[1] = new_allele2
                sample_alleles[key] = values 
            # Add the block to counter
            num_blocks += 1
            # Compare 
            compare(sample_alleles, spa_criterion)
            flush(sample_alleles)
        if current_pos - previous_pos <= threshold and current_pos - previous_pos > 0:
            if previous_record:
                for sample in previous_record.samples:
                    sample_name = sample 
                    if sample_name not in sample_alleles:
                        continue
                    sample_values = previous_record.samples[sample]['GT']
                    sample_alleles[sample_name][0].append(gt_map[sample_values[0]])
                    sample_alleles[sample_name][1].append(gt_map[sample_values[1]])
                group_i += 1
                previous_record = None
                
            for sample in record.samples:
                sample_name = sample 
                if sample_name not in sample_alleles:
                        continue
                sample_values = record.samples[sample]['GT']
                sample_alleles[sample_name][0].append(gt_map[sample_values[0]])
                sample_alleles[sample_name][1].append(gt_map[sample_values[1]])
            group_i += 1
        
        else:
            group_i = 0
            flush(sample_alleles)
            previous_record = record
        previous_pos = current_pos

    print('\n\n====For semi-private alleles====')
    print(f"Total blocks: {num_blocks}")

    spa_df = pd.DataFrame.from_dict(spa_counts, orient='index')
    spa_df.columns = pop_keys
    spa_df.insert(0, 'Animal_ID', spa_df.index)
    spa_df.reset_index(drop=True, inplace=True)
    pop_col = spa_df['Animal_ID'].map(pop_dict)
    spa_df.insert(1, 'Pop_ID', pop_col)

    spa_df_sorted = spa_df.sort_values(by=['Pop_ID'])
    spa_df_sorted.reset_index(drop=True, inplace=True)

    print()
    print('Pop sizes')
    print(n_pop)

    print('\nSPA counts')
    print(spa_df_sorted)
    spa_df_sorted.to_csv('spa_counts.csv', index=False)

    # Get the mean x variance for each population
    max_mean_x_variance = {}
    means_x_variances = {}
    for pop in pop_keys:
        for col in spa_df.columns[2:]:
            spa_df[spa_df['Pop_ID'] == pop][col] = round((spa_df[spa_df['Pop_ID'] == pop][col] / n_pop[col]) * 100, 2) # correct the values by sample size
            means_x_variances[col] = spa_df[spa_df['Pop_ID'] == pop][col].mean() * spa_df[spa_df['Pop_ID'] == pop][col].std()
            # Get the population with max mean x variance
            max_mxv_pop = max(means_x_variances.items(), key=lambda x: x[1])[0]
            max_mean_x_variance[pop] = max_mxv_pop
    
    result_data = spa_df[["Animal_ID", "Pop_ID"]]
    result_data['spa'] = spa_df[max_mean_x_variance[result_data['Pop_ID']]]

    return result_data
