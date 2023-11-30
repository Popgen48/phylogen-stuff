f1 = 'blocks.txt'
f2 = 'DRZdivA1m025_HapS4_GapS50000.TxT'

with open('blocks1.txt', 'w') as out1, open('blocks2.txt', 'w') as out2:
    with open(f1, 'r') as file1, open(f2, 'r') as file2:
        lines1 = file1.readlines()
        lines2 = file2.readlines()
        for line in lines1:
            for text in line.split('\t'):
                if 'snp' in text:
                    out1.write(text+'\t')
            out1.write('\n')
        for i, line in enumerate(lines2):
            for text in line.split():
                if text.startswith('snp'):
                    out2.write(text+'\t')
            out2.write('\n')