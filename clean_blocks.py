import sys

with open (sys.argv[1], 'r') as i, open (sys.argv[1][:-4] + '_clean' + sys.argv[1][-4:], 'w') as o:
    lines = i.readlines()
    for line in lines:
        line = line.split()
        for i, marker in enumerate(line):
            o.write(marker + '\t')
            if (i+1) % 4 == 0:
                o.write('\n')
