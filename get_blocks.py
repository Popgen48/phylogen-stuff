import sys

map_file = sys.argv[1]
threshold = int(sys.argv[2])
block_size = int(sys.argv[3])
output_file = 'blocks.txt'

with open(map_file, 'r') as file:
    map_lines = file.readlines()

temp = []

previous_pos = 0
previous_marker = None

with open(output_file, 'w') as out_file:
    for line in map_lines:
        parts = line.split()
        current_pos = int(parts[3])
        marker = parts[1]
        
        if current_pos - previous_pos < threshold:
            if previous_marker:
                temp.append(previous_marker)
            temp.append(marker)
        else:
            if len(temp) < block_size:
                temp = []
            else:
                blocks = [temp[x:x+block_size] for x in range(0, len(temp), block_size)]
                blocks = [x for x in blocks if len(x) == block_size]
                out_file.write('-------------------------------------------------------------------\n')
                out_file.write(f'{len(temp)} \t [{temp}]\n')
                out_file.write(f'Blocks: {len(blocks)} \t {blocks}\n')
                out_file.write('-------------------------------------------------------------------\n')
        previous_pos = current_pos
        previous_marker = marker
    