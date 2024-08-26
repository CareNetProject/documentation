import sys
import os
import csv
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import re
from collections import Counter

def create_csv(list_of_dicts, patientId):
    # Combine dictionaries into a single dictionary
    combined_dict = {}
    for d in list_of_dicts:
        for key, value in d.items():
            if key not in combined_dict:
                combined_dict[key] = []
    for d in list_of_dicts:
        for key in combined_dict:
            if key in d:
                combined_dict[key].append(d[key])
            else:
                combined_dict[key].append(0)
    
    # Specify the CSV file path
    csv_file_path = f"results/{patientId}_centrality_scores.csv"
    
    # Open the CSV file in write mode
    with open(csv_file_path, 'w', newline='') as csvfile:
        # Create a CSV writer object
        csvwriter = csv.writer(csvfile)
    
        # Write header row
        header = ["node", 'man_centrality_icu', 'man_indegree_icu', 'man_outdegree_icu', 'man_centrality_ward', 'man_indegree_ward', 'man_outdegree_ward',
                  'betweenness_centrality_icu', 'betweenness_centrality_ward']  
        csvwriter.writerow(header)
    
        # Write data rows
        for key, values in combined_dict.items():
            csvwriter.writerow([key] + values)
    
    #print(f"Data from dictionaries has been saved to {csv_file_path}")



def count_reverse_pairs(data):
    tuple_set = set()
    for entry in data:
        if len(entry) >= 2:
            pair = (entry[0].strip(), entry[1].strip())
            tuple_set.add(pair)

    pair_counts = {}
    for entry in data:
        if len(entry) >= 2:
            pair = (entry[1].strip(), entry[0].strip())
            if pair in tuple_set:
                pair_counts[pair] = 1

    return pair_counts

def combine_unique_pairs(data):
    unique_pairs = set()

    for entry in data:
        if len(entry) >= 2:
            pair = tuple(sorted([entry[0].strip(), entry[1].strip()]))
            unique_pairs.add(pair)

    return list(unique_pairs)

def combine_unique_values(data):
    unique_values = set()

    for entry in data:
        if len(entry) >= 2:
            unique_values.add(entry[0].strip())
            unique_values.add(entry[1].strip())

    return list(unique_values)

#based on: https://gist.github.com/aldous-rey/e6ee7b0e82e23a686d5440bf3987ee23
#altered to fit project needs
def getCentralization(centrality, c_type, label='none', columns=[], values=[], size=-1):
    c_denominator = float(1)

    if size==-1:
        n_val = float(len(centrality))
    else:
        n_val = size

    if (c_type == "degree" or c_type == 'degree_normalized'):
        c_denominator = (n_val - 1) * (n_val - 2)
        columns.append(label + '_denominator')
        values.append(c_denominator)
        columns.append(label + '_denom_formula')
        values.append('(n-1)*(n-2)')

    if (c_type == "close"):
        c_top = (n_val - 1) * (n_val - 2)
        c_bottom = (2 * n_val) - 3
        c_denominator = float(c_top / c_bottom)

    if (c_type == "between"):
        c_denominator = ((n_val - 1) * (n_val -1) * (n_val - 2))
        columns.append(label + '_denominator')
        values.append(c_denominator)
        columns.append(label + '_denom_formula')
        values.append('(n-1)*(n-1)*(n-2)')

    if (c_type == "eigen"):
        '''
        M = nx.to_scipy_sparse_matrix(G, nodelist=G.nodes(),weight='weight',dtype=float)
        eigenvalue, eigenvector = linalg.eigs(M.T, k=1, which='LR') 
        largest = eigenvector.flatten().real
        norm = sp.sign(largest.sum())*sp.linalg.norm(largest)
        centrality = dict(zip(G,map(float,largest)))
        '''

        c_denominator = sqrt(2) / 2 * (n_val - 2)

    if (c_denominator == 0):
        values.append('null')
        values.append('null')
        values.append('null')
        values.append('null')
        values.append('null')
        return None

    # start calculations

    if len(centrality.values()) == 0: 
        values.append('null')
        values.append('null')
        values.append('null')
        values.append('null')
        values.append('null')
        return 'null'
    c_node_max = max(centrality.values())
    columns.append(label + '_max_centrality')
    values.append(c_node_max)

    c_sorted = sorted(centrality.values(), reverse=True)

    c_numerator = 0

    for value in c_sorted:

        if c_type == "degree_normalized":
            # remove normalisation for each value
            c_numerator += (c_node_max * (n_val - 1) - value * (n_val - 1))
        else:
            c_numerator += (c_node_max - value)

    columns.append(label + '_numerator')
    values.append(c_numerator)
    columns.append(label + '_numberator_formula')
    if c_type == "degree_normalized":
        values.append('sum of ( (max centrality * (n-1)) - (centrality value * (n-1)) )')
    else:
        values.append('sum of (max centrality - centrality value)')

    network_centrality = float(c_numerator / c_denominator)

    return network_centrality



#get all file names from graph directory
directory = 'graphs/'
file_names = [file for file in os.listdir(directory) if file.endswith('.csv')]

#get patientId from file name
patientIds = set()

#add patientid (uncomment below to add one id, and comment out file_name loop
#patientids.add('')

for file_name in file_names:
    patientId = file_name.split('_')[0]
    patientIds.add(patientId)


#prepare output file
final_calc = open('results/auto_calcs.csv', 'w+')
man_calc = open('results/man_calcs.csv', 'w+')
writer = csv.writer(final_calc)
man_writer = csv.writer(man_calc)
column_names = ["deid_patientID", "number_of_nodes_in_icu",	"number_of_nodes_in_ward",	"number_of_edges_in_icu",	"number_of_edges_in_ward",	"network_density_in_icu",	"network_density_in_ward",	"network_reciprocity_in_icu",	"network_reciprocity_in_ward",	"degree_centralization_in_icu",	"degree_centralization_in_ward", "betweenness_centralization_in_icu",	"betweenness_centralization_in_ward"]

#Please note: man_column_names was created this way when this script was used for 1 record at a time
#As long as the first record has the correct columns, this will work as-is.
#'null' and None are currently used interchangably
man_column_names = []
write_man_columns = True
writer.writerow(column_names)

for patientId in patientIds:
    #get each patient's icu and ward files
    icu_file = f"graphs/{patientId}_icu.csv"
    ward_file = f"graphs/{patientId}_ward.csv"

    if not (os.path.exists(icu_file)):
        with open(icu_file, 'w'):
            pass
    if not (os.path.exists(ward_file)):
        with open(ward_file, 'w'):
            pass


    with open(icu_file, 'r', newline='') as file_icu, open(ward_file, 'r', newline='') as file_ward:
        data_icu = list(csv.reader(file_icu, delimiter = ','))
        data_ward = list(csv.reader(file_ward, delimiter = ','))


    # Graph - g
    g_icu = nx.DiGraph()
    g_icu.add_weighted_edges_from(([x, y, int(z)] for (x, y, z) in data_icu), weight='weight')
    g_ward = nx.DiGraph()
    g_ward.add_weighted_edges_from(([x, y, int(z)] for (x, y, z) in data_ward), weight='weight')

    values = []
    man_values = []
    man_column_names.append('deid_patientId')
    values.append(patientId)
    man_values.append(patientId)


    #manual number of nodes
    #first take all the unique providers from the graph (by using a set)
    man_nodes_icu = combine_unique_values(data_icu)
    man_nodes_ward = combine_unique_values(data_ward)
    #then get a count of the unique providers 
    man_column_names.append('man_number_nodes_icu')
    man_values.append(len(man_nodes_icu))
    man_column_names.append('man_number_nodes_ward')
    man_values.append(len(man_nodes_ward))

    if (nx.number_of_nodes(g_icu) != 0):
        number_of_nodes= nx.number_of_nodes(g_icu)
        values.append(number_of_nodes)
    else:
        values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        number_of_nodes= nx.number_of_nodes(g_ward)
        values.append(number_of_nodes)
    else:
        values.append('null')

    #manual number of eddges
    #reasoning: the <patientId>_icu/ward.csv graph files are already constructed
    #as a list of edges (with weights).  The total entries are already the total edges.
    #man_column_names.append('man_number_edges_icu')
    #man_values.append(len(data_icu))
    #man_column_names.append('man_number_edges_ward')
    #man_values.append(len(data_ward))
    

    if (nx.number_of_nodes(g_icu) != 0):
        number_of_edges = nx.number_of_edges(g_icu)
        values.append(number_of_edges)
    else:
        values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        number_of_edges = nx.number_of_edges(g_ward)
        values.append(number_of_edges)
    else:
        values.append('null')


    #manual density
    #just followed the formula:  m/(n * (n-1)) where m is total edges, and n is total nodes
    # from previous calculations
    if (len(man_nodes_icu) > 1):
        man_density_icu = (len(data_icu)) / (len(man_nodes_icu) * (len(man_nodes_icu) - 1))
    else:
        man_density_icu = 'null'
    if (len(man_nodes_ward) > 1):
        man_density_ward = (len(data_ward)) / (len(man_nodes_ward) * (len(man_nodes_ward) - 1))
    else:
        man_density_ward = 'null'
    man_column_names.append('man_density_icu')
    man_values.append(man_density_icu)
    man_column_names.append('man_density_ward')
    man_values.append(man_density_ward)
    
    if (nx.number_of_nodes(g_icu) != 0):
        density = nx.density(g_icu)
        values.append(density)
    else:
        values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        density = nx.density(g_ward)
        values.append(density)
    else:
        values.append('null')

    #manual reciprocated edges
    #function count_reverse_pairs will count all instances of, for example,
    #Node A has an edge to Node B and Node B also has an edge to Node A.
    #Important note: the above scenario will result in 2 edges being added.
    #This behavior does match what networkX will return
    man_recip_edges_icu = count_reverse_pairs(data_icu)
    man_recip_edges_ward = count_reverse_pairs(data_ward)
    man_column_names.append('man_reciprocated_edges_icu')
    man_values.append(len(man_recip_edges_icu))
    man_column_names.append('man_reciprocated_edges_ward')
    man_values.append(len(man_recip_edges_ward))

    #manual reciprocity
    #reciprocating edges over total.
    if (len(data_icu) > 0):
        man_reciprocity_icu = len(man_recip_edges_icu) / len(data_icu)
    else:
        man_reciprocity_icu = 'null'
    if (len(data_ward) > 0):
        man_reciprocity_ward = len(man_recip_edges_ward) / len(data_ward)
    else:
        man_reciprocity_ward = 'null'
    man_column_names.append('man_reciprocity_icu')
    man_values.append(man_reciprocity_icu)
    man_column_names.append('man_reciprocity_ward')
    man_values.append(man_reciprocity_ward)

    
    #arc reciprocity / same as overall reciprocity
    if (nx.number_of_nodes(g_icu) != 0):
        reciprocity = nx.reciprocity(g_icu)
        values.append(reciprocity)
    else:
        values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        reciprocity = nx.reciprocity(g_ward)
        values.append(reciprocity)
    else:
        values.append('null')
    
    #node centrality measures

    #manual calculation for degree centrality

    name_list = []
    out_list = []
    in_list = []
    generator = (x for (x, y, z) in data_icu)
    #take every  'from' name from the list
    for value in generator:
        name_list.append(value)
        out_list.append(value)
    generator = (y for (x, y, z) in data_icu)
    #combine with every 'to' name
    for value in generator:
        name_list.append(value)
        in_list.append(value)
    #this creates a dictionary with each name and count of the number of edges the name has (degree centrality, not normalized)
    man_centrality_icu = Counter(name_list)
    size = len(man_centrality_icu)
    man_indegree_icu = Counter(in_list)
    unique_elements = set(name_list)  # Get unique elements in your data
    for element in unique_elements:
        if element not in man_indegree_icu:
            man_indegree_icu[element] = 0
    man_outdegree_icu = Counter(out_list)
    unique_elements = set(name_list)  # Get unique elements in your data
    for element in unique_elements:
        if element not in man_outdegree_icu:
            man_outdegree_icu[element] = 0

    man_degree_centralization_icu = getCentralization(man_centrality_icu, 'degree', 'man_degree_centralization_icu', man_column_names, man_values, size)
    man_column_names.append('man_degree_centralization_icu')
    man_values.append(man_degree_centralization_icu)
          
    #man_indegree_centralization_icu = getCentralization(man_indegree_icu, 'degree', 'man_indegree_centralization_icu', man_column_names, man_values, size)
    #print('manual in degree centralization icu: ' + str(man_indegree_centralization_icu))
    #man_column_names.append('man_indegree_centralization_icu')
    #man_values.append(man_indegree_centralization_icu)

    #man_outdegree_centralization_icu = getCentralization(man_outdegree_icu, 'degree', 'man_outdegree_centralization_icu', man_column_names, man_values, size)
    #print('manual out degree centralization icu: ' + str(man_outdegree_centralization_icu))
    #man_column_names.append('man_outdegree_centralization_icu')
    #man_values.append(man_outdegree_centralization_icu)


    name_list = []
    out_list = []
    in_list = []
    generator = (x for (x, y, z) in data_ward)
    for value in generator:
        name_list.append(value)
        out_list.append(value)
    generator = (y for (x, y, z) in data_ward)
    for value in generator:
        name_list.append(value)
        in_list.append(value)
    man_centrality_ward = Counter(name_list)
    size = len(man_centrality_ward)
    man_indegree_ward = Counter(in_list)
    unique_elements = set(name_list)  # Get unique elements in your data
    for element in unique_elements:
        if element not in man_indegree_ward:
            man_indegree_ward[element] = 0
    man_outdegree_ward = Counter(out_list)
    unique_elements = set(name_list)  # Get unique elements in your data
    for element in unique_elements:
        if element not in man_outdegree_ward:
            man_outdegree_ward[element] = 0

    man_degree_centralization_ward = getCentralization(man_centrality_ward, 'degree', 'man_degree_centralization_ward', man_column_names, man_values, size)
    man_column_names.append('man_degree_centralization_ward')
    man_values.append(man_degree_centralization_ward)
          
    #man_indegree_centralization_ward = getCentralization(man_indegree_ward, 'degree', 'man_indegree_centralization_ward', man_column_names, man_values, size)
    #print('manual in degree centralization ward: ' + str(man_indegree_centralization_ward))
    #man_column_names.append('man_indegree_centralization_ward')
    #man_values.append(man_indegree_centralization_ward)

    #man_outdegree_centralization_ward = getCentralization(man_outdegree_ward, 'degree', 'man_outdegree_centralization_ward', man_column_names, man_values, size)
    #print('manual out degree centralization ward: ' + str(man_outdegree_centralization_ward))
    #man_column_names.append('man_outdegree_centralization_ward')
    #man_values.append(man_outdegree_centralization_ward)

    
    
    if (nx.number_of_nodes(g_icu) != 0):
        #networkX is returning the degree_centrality as a normalized value
        degree_centrality = nx.degree_centrality(g_icu)
        indegree_centrality = nx.in_degree_centrality(g_icu)
        outdegree_centrality = nx.out_degree_centrality(g_icu)
        #so we changed the getCentralization to account for this
        values.append(getCentralization(degree_centrality, 'degree_normalized'))
        #values.append(getCentralization(indegree_centrality, 'degree_normalized'))
        #values.append(getCentralization(outdegree_centrality, 'degree_normalized'))
    else:
        values.append('null')
        #values.append('null')
        #values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        degree_centrality = nx.degree_centrality(g_ward)
        indegree_centrality = nx.in_degree_centrality(g_ward)
        outdegree_centrality = nx.out_degree_centrality(g_ward)
        values.append(getCentralization(degree_centrality, 'degree_normalized'))
        #values.append(getCentralization(indegree_centrality, 'degree_normalized'))
        #values.append(getCentralization(outdegree_centrality, 'degree_normalized'))
    else:
        values.append('null')
        #values.append('null')
        #values.append('null')


    if (nx.number_of_nodes(g_icu) != 0):
        betweenness_centrality_icu = nx.betweenness_centrality(g_icu, normalized=False)

        betweenness_centralization = getCentralization(betweenness_centrality_icu, 'between', 'man_betweenness_centralization_icu', man_column_names, man_values) 
        values.append(betweenness_centralization)

        man_column_names.append('man_betweenness_centralization_icu')
        man_values.append(betweenness_centralization)
    else:
        values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
    if (nx.number_of_nodes(g_ward) != 0):
        betweenness_centrality_ward = nx.betweenness_centrality(g_ward, normalized=False)
        betweenness_centralization = getCentralization(betweenness_centrality_ward, 'between', 'man_betweenness_centralization_ward', man_column_names, man_values) 
        values.append(betweenness_centralization)

        man_column_names.append('man_betweenness_centralization_ward')
        man_values.append(betweenness_centralization)
    else:
        values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')
        man_values.append('null')

    create_csv([man_centrality_icu, man_indegree_icu, man_outdegree_icu, man_centrality_ward, man_indegree_ward, man_outdegree_ward,
                betweenness_centrality_icu, betweenness_centrality_ward], patientId)
    
    #in_degree_centrality = nx.in_degree_centrality(g)
    #prints every node
    #file.write('In degree centrality: ' + str(in_degree_centrality) + '\n')
    
    #out_degree_centrality = nx.out_degree_centrality(g)
    #prints every node
    #file.write('Out degree centrality: ' + str(out_degree_centrality) + '\n')
    
    #component analysis
    
    #if (nx.number_of_nodes(g_icu) != 0):
    #    is_strongly_connected = nx.is_strongly_connected(g_icu)
    #    values.append(is_strongly_connected)
    #else:
    #    values.append('null')
    #if (nx.number_of_nodes(g_ward) != 0):
    #    is_strongly_connected = nx.is_strongly_connected(g_ward)
    #    values.append(is_strongly_connected)
    #else:
    #    values.append('null')
    #
    #if (nx.number_of_nodes(g_icu) != 0):
    #    number_strongly_connected_components = nx.number_strongly_connected_components(g_icu)
    #    values.append(number_strongly_connected_components)
    #else:
    #    values.append('null')
    #if (nx.number_of_nodes(g_ward) != 0):
    #    number_strongly_connected_components = nx.number_strongly_connected_components(g_ward)
    #    values.append(number_strongly_connected_components)
    #else:
    #    values.append('null')
    
    #strongly_connected_components = nx.strongly_connected_components(g)
    #file.write('Strongly connected components: ' + str(strongly_connected_components) + '\n')
    
    #if (nx.number_of_nodes(g_icu) != 0):
    #    is_weakly_connected = nx.is_weakly_connected(g_icu)
    #    values.append(is_weakly_connected)
    #else:
    #    values.append('null')
    #if (nx.number_of_nodes(g_ward) != 0):
    #    is_weakly_connected = nx.is_weakly_connected(g_ward)
    #    values.append(is_weakly_connected)
    #else:
    #    values.append('null')
    #
    #if (nx.number_of_nodes(g_icu) != 0):
    #    number_weakly_connected_components = nx.number_weakly_connected_components(g_icu)
    #    values.append(number_weakly_connected_components)
    #else:
    #    values.append('null')
    #if (nx.number_of_nodes(g_ward) != 0):
    #    number_weakly_connected_components = nx.number_weakly_connected_components(g_ward)
    #    values.append(number_weakly_connected_components)
    #else:
    #    values.append('null')
    
    #weakly_connected_components = nx.weakly_connected_components(g)
    #file.write('Weakly connected components: ' + str(weakly_connected_components) + '\n')

    # Get the edge weights
    #weights = nx.get_edge_attributes(g, 'weight').values()
    
    # Count the occurrences of each weight
    #weight_counts = Counter(weights)

    # Sort the weight_counts dictionary by weight values
    #sorted_counts = sorted(weight_counts.items(), key=lambda x: x[0])
    
    # Print the count of edges with each weight
    #for weight, count in sorted_counts:
    #    file.write(f"Weight {weight}: {count} edge(s)\n")

    writer.writerow(values)
    if (write_man_columns):
        man_writer.writerow(man_column_names)
        write_man_columns = False
    man_writer.writerow(man_values)


final_calc.close()
print("Results saved in 'results' folder")

# Now to split the centrality scores into separate files
"""
import pandas as pd

# Load the original CSV file
original_df = pd.read_csv('centrality_scores.csv')
columns = ["man_centrality_icu","man_indegree_icu","man_outdegree_icu","man_centrality_ward","man_indegree_ward","man_outdegree_ward","betweenness_centrality_icu","betweenness_centrality_ward"]

# Iterate over each of the 8 columns (assuming they are named col2, col3, ..., col9)
for col in columns:
    # Create a new DataFrame with the first column and the i-th column
    new_df = original_df[['node', col]]
    
    if 'icu' in col:
        new_df = new_df[original_df['man_centrality_icu'] != 0]
    elif 'ward' in col:
        new_df = new_df[original_df['man_centrality_ward'] != 0]

    
    # Sort the DataFrame by the value in the i-th column
    new_df = new_df.sort_values(by=[col], ascending=False)
    
    # Define the filename for the new CSV file
    filename = f'{col}.csv'
    
    # Save the new DataFrame as a CSV file
    new_df.to_csv(filename, index=False)

"""
