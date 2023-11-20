import sys
import pysam

vcf_path = sys.argv[1]
threshold = float(sys.argv[2])
group_size = int(sys.argv[3])

vcf = pysam.VariantFile(vcf_path)

sample_alleles = {}

previous_record = None
group_i = 0

for i, record in enumerate(vcf):
    gt_map = {0: record.ref, 1: record.alts[0], None: '-'}
    current_pos = record.pos
    if not previous_record:
        previous_pos = record.pos
    else:
        previous_pos = previous_record.pos
    flag = False
    # Continue only if the current position is greater than previous one by the number of threshold base pairs.
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
                if sample not in sample_alleles:
                    sample_alleles[sample] = [[], []]
                sample_values = previous_record.samples[sample]['GT']
                allele1 = gt_map[sample_values[0]]
                allele2 = gt_map[sample_values[1]]
                sample_alleles[sample][0].append(allele1)
                sample_alleles[sample][1].append(allele2)
            group_i += 1
            flag = False

        # Add the current genotye to group
        for sample in record.samples:
            if sample not in sample_alleles:
                sample_alleles[sample] = [[], []]
            sample_values = record.samples[sample]['GT']       
            allele1 = gt_map[sample_values[0]]
            allele2 = gt_map[sample_values[1]]
            sample_alleles[sample][0].append(allele1)
            sample_alleles[sample][1].append(allele2)
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

print(sample_alleles['AGT05'])

my_keys = list(sample_alleles.keys())
for key, values in sample_alleles.items():
    for allele in values:
        count = 0
        idx = []
        for key2, values2 in sample_alleles.items():
            if allele in values2:
                count += 1
                idx.append(key2) 
            if count > 2:
                break
        if len(idx) == 1:
            spa_counts[key][my_keys.index(idx[0])] += 1

print(spa_counts['AGT05'])

