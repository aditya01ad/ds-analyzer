import networkx as nx

# Create the Petersen graph
G = nx.petersen_graph()

# Save it to a file named 'petersen.gml'
# The file will be created in the same directory where you run this script
nx.write_gml(G, "petersen.gml")

print("Successfully created petersen.gml")