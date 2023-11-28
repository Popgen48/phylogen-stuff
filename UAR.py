import pandas as pd
import math

# Read the file
uar_file_path = './examples/DRZdivAlp1M_DRZdivA1m025_PLINK.uar'  
with open(uar_file_path, 'r') as file:
    lines = file.readlines()

# Get the population of each animal_i
pop_file = './examples/DRZdivAlp1M.uTiere.TxT'
with open(pop_file, 'r') as file:
    pop_lines = file.readlines()
    pop_dict = {}
    for line in pop_lines:
        parts = line.split()
        animal_id = parts[1]
        pop = parts[2]
        pop_dict[animal_id] = pop

# Process the data
data = []
header = lines[0].split()  # Extract column headers

for line in lines[1:]:
    parts = line.split()
    row_label = parts[0]
    values = parts[1:]
    for col_label, value in zip(header, values):
        # if col_label == row_label:
            # continue
        data.append([pop_dict[row_label], row_label, col_label, float(value)])

# Create DataFrame
df = pd.DataFrame(data, columns=['population_animal_id_1_belongs_to', 'animal_id_1', 'animal_id_2', 'UAR score'])
print(df)

sum_within = {}
mean_within = {}
meanIBDip = {}

# Calculate genetic distance for each row
def calculate_genetic_distance(row):
    animal_i = row['animal_id_1']
    pop_i = row['population_animal_id_1_belongs_to']

    interacted_within = df[(df['animal_id_2'] == animal_i) & (df['population_animal_id_1_belongs_to'] == pop_i)]
    interacted_within = interacted_within[animal_i != interacted_within['animal_id_1']]
    n_animals = len(interacted_within)
    
    if (animal_i, pop_i) not in meanIBDip:
        meanIBDip[(animal_i, pop_i)] = interacted_within['UAR score'].sum()
        meanIBDip[(animal_i, pop_i)] = meanIBDip[(animal_i, pop_i)] / n_animals

    if pop_i not in mean_within:
        population_breeds = df[(df['population_animal_id_1_belongs_to'] == pop_i)]
        population_breeds["pop_to_compare"] = population_breeds.apply(lambda row: pop_dict[row['animal_id_2']], axis=1)
        population_breeds = population_breeds[population_breeds['pop_to_compare'] == pop_i]
        population_breeds = population_breeds[population_breeds['animal_id_1'] != population_breeds['animal_id_2']]
        
        sum_within[pop_i] = population_breeds['UAR score'].sum()
        mean_within[pop_i] = sum_within[pop_i] / len(population_breeds)

    mUARi = meanIBDip[(animal_i, pop_i)]
    mUARmp = mean_within[pop_i]
    
    genetic_distance = -math.log(abs(mUARi + mUARmp))
    return genetic_distance

# Calculate maxUAR1(B) and maxUAR2(B) for each row
def calculate_max_uar(row):
    animal_i = row['animal_id_1']
    population_i = row['population_animal_id_1_belongs_to']
    
    # Filter rows where animal_id_1 is different from current animal_i and from a different population
    foreign_breeds = df[(df['animal_id_1'] != animal_i) & (df['population_animal_id_1_belongs_to'] != population_i)]
    
    # For Max_UAR1(B) and Max_UAR2(B) 
    # Get the UAR scores between animal_i i and animals from different populations
    uar_values = foreign_breeds[foreign_breeds['animal_id_2'] == animal_i]['UAR score']
    # Sort and get the highest (maxUAR1(B)) and second highest UAR score (maxUAR2(B))
    sorted_uar_values = uar_values.sort_values(ascending=False)
    max_uar1_b = sorted_uar_values.iloc[0] if len(sorted_uar_values) > 0 else 0.0
    max_uar2_b = sorted_uar_values.iloc[1] if len(sorted_uar_values) > 1 else 0.0

    # For Max_UAR(p) 
    foreign_breeds = foreign_breeds[(foreign_breeds['animal_id_2'] == animal_i)]
    # Group by the breed of animal_id_1 and calculate the mean UAR for each breed
    breed_avg_uar = foreign_breeds.groupby('population_animal_id_1_belongs_to')['UAR score'].mean()
    # Get the maximum average UAR across breeds
    max_uar_p = breed_avg_uar.max()
    max_uar_p = max_uar_p if not pd.isnull(max_uar_p) else 0.0

    return max_uar1_b, max_uar2_b, max_uar_p

# Calculate maxUAR(P)
# def calculate_max_uar_p(row):
#     animal_i = row['animal_id_1']
#     population_i = row['population_animal_id_1_belongs_to']
    
#     # Filter rows where animal_id_1 is different from current animal_i and from a different population
#     foreign_breeds = df[(df['animal_id_1'] != animal_i) & (df['population_animal_id_1_belongs_to'] != population_i)]
#     foreign_breeds = foreign_breeds[(foreign_breeds['animal_id_2'] == animal_i)]
    
#     # Group by the breed of animal_id_1 and calculate the mean UAR for each breed
#     breed_avg_uar = foreign_breeds.groupby('population_animal_id_1_belongs_to')['UAR score'].mean()
    
#     # Get the maximum average UAR across breeds
#     max_uar_p = breed_avg_uar.max()
    
#     return max_uar_p if not pd.isnull(max_uar_p) else 0.0  # If no interactions with foreign breeds, return 0

# Get unique animal_id values
unique_animals = pd.unique(df[['animal_id_1', 'animal_id_2']].values.ravel('K'))

result_data = []

for animal_i in unique_animals:
    temp = df[df['animal_id_1'] == animal_i].head(n=1)
    gd = temp.apply(calculate_genetic_distance, axis=1).values[0]
    max_uar1_b, max_uar2_b, max_uar_p = temp.apply(calculate_max_uar, axis=1).values[0]
    # max_uar_p = temp.apply(calculate_max_uar_p, axis=1).values[0]
    animal_data = {
        'Animal_id': animal_i,
        'Genetic Distance': gd,  
        'maxUAR1(B)': max_uar1_b,  
        'maxUAR2(B)': max_uar2_b,  
        'maxUAR(P)': max_uar_p  
    }
    result_data.append(animal_data)
    
result_df = pd.DataFrame(result_data)
print(result_df)