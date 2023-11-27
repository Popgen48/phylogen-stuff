# Usage: python3 get_sp_alleles.py <vcf_path> <threshold> <grouping_size> <spa_criterion>

import sys
import pysam
import pandas as pd

vcf_path = sys.argv[1]
threshold = float(sys.argv[2]) # The max number of base pairs i.e. max distance between two SNPs to be considered
group_size = int(sys.argv[3]) # # Max grouping size of the haplotypes
spa_criterion = int(sys.argv[4]) # The max number of distinct groups (apart from itself) in which an allele can be present to be classified as a semi-private allele

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
previous_record = None
group_i = 0
flag = False

for i, record in enumerate(vcf):
    gt_map = {0: record.ref, 1: record.alts[0], None: '-'}
    current_pos = record.pos
    if not previous_record:
        previous_pos = record.pos
    else:
        previous_pos = previous_record.pos

    # Continue only if the current position is lesser than previous one by the number of threshold base pairs.
    if current_pos - previous_pos < threshold:
        if group_i == group_size:
            for key, values in sample_alleles.items():
                new_allele1 = [''.join(values[0][-group_i:])]
                values[0] = values[0][:-group_i] + new_allele1

                new_allele2 = [''.join(values[1][-group_i:])]
                values[1] = values[1][:-group_i] + new_allele2

                sample_alleles[key] = values
            group_i = 0
            
        if flag:
            # Also Add the previous genotye to group
            for sample in previous_record.samples:
                sample_name = '_'.join(sample.split('_')[1:])
                if sample_name not in sample_alleles:
                    sample_alleles[sample_name] = [[], []]
                sample_values = previous_record.samples[sample]['GT']
                allele1 = gt_map[sample_values[0]]
                allele2 = gt_map[sample_values[1]]
                sample_alleles[sample_name][0].append(allele1)
                sample_alleles[sample_name][1].append(allele2)
            group_i += 1
            flag = False

        # Add the current genotye to group
        for sample in record.samples:
            sample_name = '_'.join(sample.split('_')[1:])
            if sample_name not in sample_alleles:
                sample_alleles[sample_name] = [[], []]
            sample_values = record.samples[sample]['GT']       
            allele1 = gt_map[sample_values[0]]
            allele2 = gt_map[sample_values[1]]
            sample_alleles[sample_name][0].append(allele1)
            sample_alleles[sample_name][1].append(allele2)
        group_i += 1
    else:
        flag = True
        if group_i == 0:
            continue
        else:
            for key, values in sample_alleles.items():
                new_allele1 = [''.join(values[0][-group_i:])]
                values[0] = values[0][:-group_i] + new_allele1

                new_allele2 = [''.join(values[1][-group_i:])]
                values[1] = values[1][:-group_i] + new_allele2

                sample_alleles[key] = values
        group_i = 0
    previous_record = record

spa_counts = {}
for key, values in sample_alleles.items():
    sample_alleles[key] = values[0] + values[1]
    spa_counts[key] = [0] * len(sample_alleles)


print('APPgoat1')
print(sample_alleles['APPgoat1'])

my_keys = list(sample_alleles.keys())
for key, values in sample_alleles.items():
    pop1 = pop_dict[key]
    for allele in values:
        count = 0
        idx = []
        pops_considered = []
        for key2, values2 in sample_alleles.items():
            pop2 = pop_dict[key2]
            if pop2 == pop1:
                break
            if allele in values2:
                if pop2 in pops_considered:
                    idx.append(key2) 
                else:
                    pops_considered.append(pop2)
                    count += 1
                    idx.append(key2)
            if count > spa_criterion:
                break
        if len(idx) > 0 and len(idx) <= spa_criterion:
            for i in idx:
                spa_counts[key][my_keys.index(i)] += 1

temp_data = []

for key, values in spa_counts.items():
    dict = {}
    for p in list(n_pop.keys()):
        dict[p] = 0
    for val in values:
        if val > 0:
            dict[pop_dict[my_keys[values.index(val)]]] += 1
    
    # Correct the number by sample size
    for p in list(dict.keys()):
        dict[p] = (dict[p] / n_pop[p]) * 100

    # shared_pop, max_count = max(dict.items(), key=lambda x: x[1])

    # if max_count == 0:
        # shared_pop = '-'
    
    rec = {
        "Animal_ID": key,
        "Own_Population": pop_dict[key],
        # "Shared_Population": shared_pop,
        # "SPA_Count": max_count,
    }
    rec = {**rec, **dict}
    temp_data.append(rec)

temp_data = pd.DataFrame(temp_data)

print()
print(temp_data)

# Get the mean x variance for each population
means_x_variances = {}
for col in temp_data.columns[2:]:
    means_x_variances[col] = temp_data[col].mean() * temp_data[col].var()

# Get the population with max mean x variance
max_mxv_pop = max(means_x_variances.items(), key=lambda x: x[1])[0]

result_data = temp_data[["Animal_ID", "Own_Population", max_mxv_pop]].rename(columns={'Own_Population': 'Pop_ID', max_mxv_pop: 'spa'})

print()
print(result_data[result_data['spa'] > 0])
