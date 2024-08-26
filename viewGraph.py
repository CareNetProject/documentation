import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import csv

def normalize(data):
    # Min-max scaling
    min_val = np.min(data)
    max_val = np.max(data)
    normalized_data = (data - min_val) / (max_val - min_val)
    #adding .1 to all values so the min is not zero
    normalized_data = normalized_data + .1
    return normalized_data


if len(sys.argv) < 2:
    print("Usage: python3 viewGraph.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, 'r') as file:
        data = list(csv.reader(file, delimiter=','))
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    sys.exit(1)

file.close()

# Graph - g
g = nx.DiGraph()

# ignores first row (header)
g.add_weighted_edges_from(([x, y, int(z)] for (x, y, z) in data[1:]), weight='weight')# data [1:] means we start at the second row.

#arc reciprocity / same as overall reciprocity
reciprocity = nx.reciprocity(g)
arrayReciprocity = []
for node in g.nodes():
    nReciprocity = nx.reciprocity(g, node)
    arrayReciprocity.append(nReciprocity)



###############
# Everything that follows is for drawing the graph

# Adjusts size of figure
plt.figure(figsize=(15, 9))

# positions for the drawing- pos
# can uncomment whichever layout is desired

pos = nx.fruchterman_reingold_layout(g,  scale = 1, k=1)
#pos = nx.circular_layout(g)
#pos = nx.spectral_layout(g)
#pos = nx.shell_layout(g)
#pos = nx.spring_layout(g, scale = 1, k=1)
#pos = nx.random_layout(g)

nodeSize = arrayReciprocity
#nodeSize = arrayDegCentrality
#nodeSize = arrayBetCentrality
nodeSize= normalize(nodeSize)
# Node size is based on degree (with weight)
node_sizes = [(size * 500) for (size) in nodeSize]#the size will be the actual weight.

# edge colors/transparency assigned by weight
#edge_colors = [int(z) for (x, y, z) in data[1:]]#edge_colors is a list of weights
edge_colors = [1 for (x, y, z) in data[1:]]#edge_colors is a list of weights
maxColor = max(edge_colors)
#edge_alphas = [int(z) / maxColor for (x, y, z) in data[1:]]
edge_alphas = [1 / maxColor for (x, y, z) in data[1:]]
#cmap = plt.cm.plasma

# determines node colors
node_colors = []
for x, y in zip(g.in_degree(), g.out_degree()):
    node_colors.append('red')
    #inD = x[1]
    #outD = y[1]
    #if inD > 0 and outD > 0:
    #    node_colors.append("green")
    #elif inD > 0:
    #    node_colors.append("red")
    #else:
    #    node_colors.append("blue")

## a separate check for reciprocal nodes
#nodes = list(g)
#for x in range(len(nodes)):
#    for compare_node in g.nodes(): #for every node compare it to every node, but
#        if (compare_node != nodes[x]): #don't compare node a with node a.
#            if g.has_edge(nodes[x], compare_node) and g.has_edge(compare_node, nodes[x]): #if node a when compared to node b have reciprocal edges then color them purple.
#                node_colors[x] = "purple"

# draws each node in graph g with given positions pos, node sizes node_sizes, and node_colors
nodes = nx.draw_networkx_nodes(g, pos, node_size=node_sizes, node_color=node_colors)
# draws each edge in graph g with positions pos, node sizes node_sizes, etc
edges = nx.draw_networkx_edges(
    g,
    pos,
    node_size=node_sizes,
    arrowstyle="->",
    arrowsize=10,
    #edge_color=edge_colors,
    edge_color='black',
    #edge_cmap=cmap,
    width=.5,
)

# set alpha value for each edge
m = g.number_of_edges()
for i in range(m):
    edges[i].set_alpha(edge_alphas[i])#setting edges_alpha[i] to the edges_alpha list we created previously

# setup the colorbar and draw everything
#pc = mpl.collections.PatchCollection(edges, cmap=cmap)
pc = mpl.collections.PatchCollection(edges)
pc.set_array(edge_colors)

ax = plt.gca()
ax.set_axis_off()
#plt.colorbar(pc, ax=ax)
plt.show()
