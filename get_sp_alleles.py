# Usage: python3 get_sp_alleles.py <vcf_path> <threshold> <grouping_size> <spa_criterion>

import sys
import pysam
import pandas as pd

vcf_path = sys.argv[1]
threshold = float(sys.argv[2]) # The max number of base pairs i.e. max distance between two SNPs to be considered
group_size = int(sys.argv[3]) # # Max grouping size of the haplotypes
spa_criterion = int(sys.argv[4]) # The max number of distinct groups (apart from itself) in which an allele can be present to be classified as a semi-private allele

def compare(dict):
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
                for i in idx:
                    spa_counts[key][my_keys.index(i)] += 1 
                    
def flush(dict):
    for key, values in dict.items():
        dict[key] = [[], []]

vcf = pysam.VariantFile(vcf_path)

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

sample_alleles = {}
spa_counts = {}

for key in pop_dict.keys():
    sample_alleles[key] = [[], []]
    spa_counts[key] = [0] * len(pop_dict.keys())
my_keys = list(sample_alleles.keys())

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
        compare(sample_alleles)
        flush(sample_alleles)
    
    if current_pos - previous_pos < threshold:
        if previous_record:
            for sample in previous_record.samples:
                sample_name = sample # '_'.join(sample.split('_')[1:])
                if sample_name not in sample_alleles:
                    continue
                # if sample_name not in sample_alleles:
                #     sample_alleles[sample_name] = [[], []]
                sample_values = previous_record.samples[sample]['GT']
                sample_alleles[sample_name][0].append(gt_map[sample_values[0]])
                sample_alleles[sample_name][1].append(gt_map[sample_values[1]])
            group_i += 1
            previous_record = None
            
        for sample in record.samples:
            sample_name = sample # '_'.join(sample.split('_')[1:])
            if sample_name not in sample_alleles:
                    continue
            # if sample_name not in sample_alleles:
            #     sample_alleles[sample_name] = [[], []]
            sample_values = record.samples[sample]['GT']
            sample_alleles[sample_name][0].append(gt_map[sample_values[0]])
            sample_alleles[sample_name][1].append(gt_map[sample_values[1]])
        group_i += 1
    
    else:
        group_i = 0
        flush(sample_alleles)
        previous_record = record
        
    previous_pos = current_pos

print(f"Total blocks: {num_blocks}")

temp_data = []

for key, values in spa_counts.items():
    dict = {}
    for p in list(n_pop.keys()):
        dict[p] = 0
    for val in values:
        if val > 0:
            dict[pop_dict[my_keys[values.index(val)]]] += 1
    
    rec = {
        "Animal_ID": key,
        "Own_Population": pop_dict[key],
    }
    rec = {**rec, **dict}
    temp_data.append(rec)

temp_data = pd.DataFrame(temp_data)

print()
print('Pop sizes')
print(n_pop)

print()
columns_to_check = temp_data.columns[2:]
print(temp_data[(temp_data[columns_to_check] != 0).any(axis=1)])
temp_data.to_csv('temp_data.csv', index=False)

# Get the mean x variance for each population
means_x_variances = {}
for col in temp_data.columns[2:]:
    temp_data[col] = round((temp_data[col] / n_pop[col]) * 100, 2) # correct the values by sample size
    means_x_variances[col] = temp_data[col].mean() * temp_data[col].std()

# Get the population with max mean x variance
max_mxv_pop = max(means_x_variances.items(), key=lambda x: x[1])[0]

result_data = temp_data[["Animal_ID", "Own_Population", max_mxv_pop]].rename(columns={'Own_Population': 'Pop_ID', max_mxv_pop: 'spa_'+max_mxv_pop})

print()
print(result_data[result_data['spa_'+max_mxv_pop] > 0])
result_data.to_csv('result_data.csv', index=False)

# print()
# print(temp_data[(temp_data[columns_to_check] != 0).any(axis=1)])