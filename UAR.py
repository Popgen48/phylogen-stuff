import pandas as pd
import numpy as np
import math

# Sample DataFrame
data = {
    'population_animal_id_1_belongs_to': ['pop1', 'pop1', 'pop2', 'pop2', 'pop2'],
    'animal_id_1': ['A1', 'A1', 'A2', 'A2', 'A3'],
    'animal_id_2': ['A2', 'A3', 'A1', 'A3', 'A1'],
    'UAR score': [0.987, 0.767, 0.987, 0.123, 0.766]
}

df = pd.DataFrame(data)

# Calculate mean UAR for each animal across the population
animal_mean_uar = df.groupby('animal_id_1')['UAR score'].mean().reset_index()
animal_mean_uar.columns = ['animal_id', 'mean_UAR']

# Calculate mean UAR for the entire population
mean_uar_population = df.groupby('population_animal_id_1_belongs_to')['UAR score'].mean().reset_index()
mean_uar_population.columns = ['population', 'mean_UAR']

# Calculate genetic distance for each row
def calculate_genetic_distance(row):
    mUARi = animal_mean_uar.loc[animal_mean_uar['animal_id'] == row['animal_id_1'], 'mean_UAR'].values[0]
    mUARmp = mean_uar_population.loc[mean_uar_population['population'] == row['population_animal_id_1_belongs_to'], 'mean_UAR'].values[0]
    genetic_distance = -math.log(mUARi + mUARmp)
    return genetic_distance

# df['Genetic Distance'] = df.apply(calculate_genetic_distance, axis=1)

# Calculate maxUAR1(B) and maxUAR2(B) for each row
def calculate_max_uar_b(row):
    animal_i = row['animal_id_1']
    population_i = row['population_animal_id_1_belongs_to']
    
    # Filter rows where animal_id_1 is different from current animal and from a different population
    foreign_breeds = df[(df['animal_id_1'] != animal_i) & (df['population_animal_id_1_belongs_to'] != population_i)]
    # print(foreign_breeds)
    
    # Get the UAR scores between animal i and animals from different populations
    uar_values = foreign_breeds[foreign_breeds['animal_id_2'] == animal_i]['UAR score']
    
    # Sort and get the highest (maxUAR1(B)) and second highest UAR score (maxUAR2(B))
    sorted_uar_values = uar_values.sort_values(ascending=False)
    max_uar1_b = sorted_uar_values.iloc[0] if len(sorted_uar_values) > 0 else 0.0
    max_uar2_b = sorted_uar_values.iloc[1] if len(sorted_uar_values) > 1 else 0.0
    
    return max_uar1_b, max_uar2_b
    
# df['maxUAR1(B)'], df['maxUAR2(B)'] = zip(*df.apply(calculate_max_uar_b, axis=1))

# Calculate maxUAR(P)
def calculate_max_uar_p(row):
    animal_i = row['animal_id_1']
    population_i = row['population_animal_id_1_belongs_to']
    
    # Filter rows where animal_id_1 is different from current animal and from a different population
    foreign_breeds = df[(df['animal_id_1'] != animal_i) & (df['population_animal_id_1_belongs_to'] != population_i)]
    foreign_breeds = foreign_breeds[(foreign_breeds['animal_id_2'] == animal_i)]
    
    # Group by the breed of animal_id_1 and calculate the mean UAR for each breed
    breed_avg_uar = foreign_breeds.groupby('population_animal_id_1_belongs_to')['UAR score'].mean()
    
    # Get the maximum average UAR across breeds
    max_uar_p = breed_avg_uar.max()
    
    return max_uar_p if not pd.isnull(max_uar_p) else 0.0  # If no interactions with foreign breeds, return 0

# df['maxUAR(P)'] = df.apply(calculate_max_uar_p, axis=1)

print(df)

# Get unique animal_id values
unique_animals = pd.unique(df[['animal_id_1', 'animal_id_2']].values.ravel('K'))

# result_df = pd.DataFrame(columns=['Animal_id', 'Genetic Distance', 'maxUAR1(B)', 'maxUAR2(B)', 'maxUAR(P)'])
result_data = []

for animal in unique_animals:
    temp = df[df['animal_id_1'] == animal].head(n=1)
    gd = temp.apply(calculate_genetic_distance, axis=1).values[0]
    max_uar1_b, max_uar2_b = temp.apply(calculate_max_uar_b, axis=1).values[0]
    max_uar_p = temp.apply(calculate_max_uar_p, axis=1).values[0]
    animal_data = {
        'Animal_id': animal,
        'Genetic Distance': gd,  
        'maxUAR1(B)': max_uar1_b,  
        'maxUAR2(B)': max_uar2_b,  
        'maxUAR(P)': max_uar_p  
    }
    result_data.append(animal_data)
    
result_df = pd.DataFrame(result_data)
print(result_df)