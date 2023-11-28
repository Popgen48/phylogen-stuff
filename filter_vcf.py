import pysam
import sys

vcf_path = sys.argv[1]
maf_threshold = float(sys.argv[2])

def get_maf(record):
    mafs = [0, 0]
    for val in record.samples:
        if record.samples[val]['GT'][0] != None:
            mafs[record.samples[val]['GT'][0]] += 1
        if record.samples[val]['GT'][1] != None:
            mafs[record.samples[val]['GT'][1]] += 1
    return min(mafs) / sum(mafs)

vcf = pysam.VariantFile(vcf_path)

# Filter the vcf file by MAF
filtered_vcf_path = vcf_path[:-4] + '_filtered.vcf' 
filtered_vcf = pysam.VariantFile(filtered_vcf_path, 'w', header=vcf.header)
for record in vcf:
    maf = get_maf(record)
    if maf >= maf_threshold:
        filtered_vcf.write(record)
        
vcf.close()
filtered_vcf.close()

print(f"Filtered vcf file is saved at '{filtered_vcf_path}'")