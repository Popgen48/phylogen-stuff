import pysam
import sys

vcf_path = sys.argv[1]
vcf = pysam.VariantFile(vcf_path)

map_path = vcf_path[:-4] + '.map'

with open(map_path, 'w') as file:
    for record in vcf:
        file.write(f"{record.contig}\t{record.id}\t0\t{record.pos}\n")

print(f'Map file written to {map_path}')