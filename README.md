# ds-analyzer

A Python tool to analyze whether a graph is **Determined by its Spectrum (DS)**.

A graph is *determined by its spectrum* if no other non-isomorphic graph shares the same set of eigenvalues. `ds-analyzer` computes structural and spectral invariants of a graph and searches for cospectral mates using the **Godsil–McKay switching** technique to assess the DS property.

## Features

- **Structural invariants** – number of vertices and edges, degree sequence, connectivity, diameter, radius, girth, vertex/edge connectivity, and automorphism group size.
- **Spectral invariants** – adjacency spectrum, Laplacian spectrum, and signless Laplacian spectrum, together with derived quantities (number of components, spanning trees, triangles, etc.).
- **Godsil–McKay switching search** – exhaustive search for a non-isomorphic cospectral mate by enumerating valid switching partitions.
- **CLI interface** – analyze any graph file from the command line with a single command.

## Requirements

- Python ≥ 3.8
- [NetworkX](https://networkx.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [Click](https://click.palletsprojects.com/)

## Installation

```bash
# Clone the repository
git clone https://github.com/aditya01ad/ds-analyzer.git
cd ds-analyzer

# Install in editable mode (recommended for development)
pip install -e .

# Or install directly
pip install .
```

## Usage

### Command-line interface

```bash
ds-analyzer GRAPH_FILE [OPTIONS]
```

**Arguments**

| Argument     | Description                                        |
|--------------|----------------------------------------------------|
| `GRAPH_FILE` | Path to a graph file (GML, GraphML, or edge list). |

**Options**

| Option                   | Default | Description                                              |
|--------------------------|---------|----------------------------------------------------------|
| `--max-set-size INTEGER` | `8`     | Maximum size of the Godsil–McKay switching set to check. |
| `--time-limit INTEGER`   | `300`   | Time limit in seconds for the switching search.          |

**Example**

```bash
# Generate the example Petersen graph file
python create_test_graph.py

# Analyze it
ds-analyzer petersen.gml --max-set-size 6 --time-limit 60
```

Sample output:

```
Loading graph from petersen.gml...
Initializing DS Analyzer...
Step 1: Computing structural invariants...
Step 2: Computing spectral invariants...
Step 3: Searching for Godsil-McKay switching partitions...
Analysis complete.

--- DS Analysis Report ---

[Structural Properties]
num_vertices: 10
num_edges: 15
...

[Spectral Properties]
adjacency_spectrum: [-2.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 3.0]
...

[Godsil-McKay Switching Analysis]
  - Starting GM switching search ...
  - Search completed. No non-isomorphic cospectral mate found via GM switching.

Conclusion: The graph is LIKELY a DS graph.
No cospectral mate was found via Godsil-McKay switching.
```

### Python API

```python
import networkx as nx
from ds_analyzer.analyzer import DSAnalyzer

# Load or create a graph
G = nx.petersen_graph()

# Initialize and run the full analysis
analyzer = DSAnalyzer(G)
analyzer.run_full_analysis(max_set_size=6, time_limit=60)

# Access results
print(analyzer.results["graph_properties"])
print(analyzer.results["spectral_properties"])
print(analyzer.results["gm_switching_analysis"]["found_cospectral_mate"])
```

You can also run individual steps:

```python
analyzer.compute_structural_invariants()
analyzer.compute_spectral_invariants()
analyzer.find_gm_switching_partitions(max_set_size=8, time_limit=300)
```

## Supported Graph File Formats

| Extension           | Format    | NetworkX reader    |
|---------------------|-----------|--------------------|
| `.gml`              | GML       | `nx.read_gml`      |
| `.graphml` / `.xml` | GraphML   | `nx.read_graphml`  |
| `.edgelist` / `.txt`| Edge list | `nx.read_edgelist` |

## How It Works

1. **Structural invariants** – basic graph properties are computed using NetworkX algorithms.
2. **Spectral invariants** – eigenvalues of the adjacency, Laplacian, and signless Laplacian matrices are computed. For graphs with fewer than 500 nodes, dense linear algebra (NumPy) is used; for larger graphs, sparse iterative solvers (SciPy) are used.
3. **Godsil–McKay switching** – the algorithm searches for a subset *X* of vertices such that:
   - The induced subgraph *G[X]* is regular.
   - Every vertex outside *X* has exactly 0, |X|/2, or |X| neighbours in *X*.

   If such a partition is found and the switched graph *G'* is non-isomorphic to *G*, then *G* is **not** DS. Otherwise the graph is reported as *likely* DS (no guarantee without an exhaustive search over all possible cospectral mates).

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.
