import networkx as nx
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import eigs, eigsh
from itertools import combinations
from collections import Counter
import time

class DSAnalyzer:
    """
    A comprehensive toolkit for analyzing the structural and spectral properties
    of a graph G to assess its likelihood of being Determined by its Spectrum (DS).
    """

    def __init__(self, graph: nx.Graph):
        """
        Initializes the analyzer with a NetworkX graph object.

        Args:
            graph (nx.Graph): The input graph to be analyzed. Must be a simple,
                              undirected graph.
        """
        # This is the correct line
        if not isinstance(graph, nx.Graph) or graph.is_directed() or graph.is_multigraph():
            raise TypeError("Input must be a simple, undirected NetworkX Graph object.")
        
        self.graph = graph
        self.n = graph.number_of_nodes()
        self.m = graph.number_of_edges()
        self.adj_matrix = nx.to_scipy_sparse_array(self.graph, format='csr')
        
        self.results = {
            "graph_properties": {},
            "spectral_properties": {},
            "gm_switching_analysis": {
                "found_cospectral_mate": False,
                "mate_graph": None,
                "switching_partition": None,
                "analysis_log": []
            }
        }
        self._is_connected = None

    def _compute_eigenvalues(self, matrix, k, symmetric=False):
        """Helper to compute eigenvalues, handling small or dense cases."""
        if self.n == 0 or matrix.shape[0] == 0:
            return np.array([])
        # For small matrices, dense computation is faster
        if self.n < 500:
            dense_matrix = matrix.toarray()
            if symmetric:
                return np.linalg.eigvalsh(dense_matrix)
            else:
                return np.linalg.eigvals(dense_matrix)
        
        # For large sparse matrices, use iterative solvers
        try:
            if symmetric:
                # eigsh is for symmetric matrices
                vals = eigsh(matrix, k=k, which='LM', return_eigenvectors=False)
            else:
                # eigs is for general matrices
                vals = eigs(matrix, k=k, which='LM', return_eigenvectors=False)
            return np.real(vals)
        except Exception as e:
            self.results['gm_switching_analysis']['analysis_log'].append(f"Eigenvalue computation failed: {e}")
            return np.array()

    def compute_structural_invariants(self):
        """
        Computes fundamental structural properties of the graph.
        These properties are not necessarily determined by the spectrum.
        """
        props = self.results["graph_properties"]
        props["num_vertices"] = self.n
        props["num_edges"] = self.m
        
        degrees = [d for _, d in self.graph.degree()]
        props["degree_sequence"] = sorted(degrees, reverse=True)
        
        self._is_connected = nx.is_connected(self.graph)
        props["is_connected"] = self._is_connected
        
        if self._is_connected:
            props["num_connected_components"] = 1
            props["diameter"] = nx.diameter(self.graph)
            props["radius"] = nx.radius(self.graph)
            try:
                props["girth"] = nx.girth(self.graph)
                #props["girth"] = self._compute_girth(self.graph)
            except nx.NetworkXNoCycle:
                props["girth"] = float('inf')
            props["vertex_connectivity"] = nx.node_connectivity(self.graph)
            props["edge_connectivity"] = nx.edge_connectivity(self.graph)
        else:
            props["num_connected_components"] = nx.number_connected_components(self.graph)
            # Metrics below are undefined for disconnected graphs
            props["diameter"] = float('inf')
            props["radius"] = float('inf')
            props["girth"] = min([nx.girth(c) for c in (self.graph.subgraph(comp).copy() for comp in nx.connected_components(self.graph)) if c.number_of_nodes() > 2], default=float('inf'))
            props["vertex_connectivity"] = 0
            props["edge_connectivity"] = 0

        # Automorphism group size - computationally intensive
        # Note: This relies on an external library (e.g., bliss) if installed,
        # otherwise falls back to a slower VF2-based algorithm.
        try:
            # Using a placeholder for direct computation as it can be very slow.
            # In a real scenario, this would call a dedicated isomorphism backend.
            # For demonstration, we assume a function `_get_aut_group_size` exists.
            # This is a known hard problem (GI-complete).
            is_isomorphic = nx.isomorphism.GraphMatcher(self.graph, self.graph)
            props["automorphism_group_size"] = sum(1 for _ in is_isomorphic.isomorphisms_iter())
        except Exception:
             props["automorphism_group_size"] = "Computation failed or too complex"

    def compute_spectral_invariants(self):
        """
        Computes spectral properties from Adjacency, Laplacian, and other matrices.
        """
        spec_props = self.results["spectral_properties"]
        
        # Adjacency Spectrum
        adj_eigenvalues = self._compute_eigenvalues(self.adj_matrix, k=self.n - 1, symmetric=True)
        spec_props["adjacency_spectrum"] = np.round(adj_eigenvalues, 8).tolist()
        
        # Properties from Adjacency Spectrum
        spec_props["num_vertices_from_spec"] = len(adj_eigenvalues)
        spec_props["num_edges_from_spec"] = 0.5 * np.sum(adj_eigenvalues**2)
        spec_props["num_triangles_from_spec"] = (1/6) * np.sum(adj_eigenvalues**3)
        
        # Check for regularity
        #is_regular = all(d == self.results["graph_properties"]["degree_sequence"] for d in self.results["graph_properties"]["degree_sequence"])
        #self.results["graph_properties"]["is_regular"] = is_regular
        deg_seq = self.results["graph_properties"]["degree_sequence"]
        is_regular = all(d == deg_seq[0] for d in deg_seq)

        
        # Check for bipartiteness
        sorted_spec = sorted(np.round(adj_eigenvalues, 8))
        is_bipartite = np.allclose(sorted_spec, -np.array(sorted_spec)[::-1])
        self.results["graph_properties"]["is_bipartite"] = is_bipartite

        # Laplacian Spectrum
        L = nx.laplacian_matrix(self.graph)
        lap_eigenvalues = self._compute_eigenvalues(L, k=self.n - 1, symmetric=True)
        spec_props["laplacian_spectrum"] = np.round(lap_eigenvalues, 8).tolist()
        
        # Properties from Laplacian Spectrum
        spec_props["num_components_from_spec"] = np.sum(np.isclose(lap_eigenvalues, 0))
        if self._is_connected and self.n > 1:
            spec_props["num_spanning_trees"] = np.prod(lap_eigenvalues[lap_eigenvalues > 1e-8]) / self.n
        else:
            spec_props["num_spanning_trees"] = 0

        # This is the corrected block
        # Signless Laplacian Spectrum
        degrees = [d for _, d in self.graph.degree()]
        D = scipy.sparse.diags(degrees, format='csr')
        A = self.adj_matrix
        Q = D + A
        q_eigenvalues = self._compute_eigenvalues(Q, k=self.n - 1, symmetric=True)
        spec_props["signless_laplacian_spectrum"] = np.round(q_eigenvalues, 8).tolist()

    def find_gm_switching_partitions(self, max_set_size=8, time_limit=300):
        """
        Searches for a Godsil-McKay switching partition to find a cospectral mate.

        Args:
            max_set_size (int): Maximum size of the switching set X to check.
            time_limit (int): Time limit in seconds for the search.
        """
        log = self.results['gm_switching_analysis']['analysis_log']
        start_time = time.time()
        log.append(f"Starting GM switching search with max_set_size={max_set_size} and time_limit={time_limit}s.")

        nodes = list(self.graph.nodes())
        
        for k in range(4, min(max_set_size, self.n // 2) + 1, 2):
            if time.time() - start_time > time_limit:
                log.append("Time limit exceeded. Terminating search.")
                return

            log.append(f"Checking for switching sets of size k={k}.")
            
            for X_nodes in combinations(nodes, k):
                if time.time() - start_time > time_limit:
                    log.append("Time limit exceeded. Terminating search.")
                    return

                X = set(X_nodes)
                
                # Condition 1: Induced subgraph G[X] must be regular.
                subgraph_X = self.graph.subgraph(X)
                degrees_in_X = [d for _, d in subgraph_X.degree()]
                #if not all(d == degrees_in_X for d in degrees_in_X):
                if not all(d == degrees_in_X[0] for d in degrees_in_X):
                    continue
                
                log.append(f"Found candidate set X={X} (induces a regular subgraph).")

                # Condition 2: Vertices in Y must have 0, k/2, or k neighbors in X.
                Y = set(nodes) - X
                is_valid_partition = True
                switching_vertices_in_Y = []
                
                for y in Y:
                    neighbors_in_X = len(set(self.graph.neighbors(y)) & X)
                    if neighbors_in_X not in {0, k / 2, k}:
                        is_valid_partition = False
                        break
                    if neighbors_in_X == k / 2:
                        switching_vertices_in_Y.append(y)
                
                if is_valid_partition and switching_vertices_in_Y:
                    log.append(f"Found valid GM partition with X={X}.")
                    
                    # Perform the switch
                    g_prime = self.graph.copy()
                    for y in switching_vertices_in_Y:
                        neighbors_in_X = set(self.graph.neighbors(y)) & X
                        non_neighbors_in_X = X - neighbors_in_X
                        
                        edges_to_remove = [(y, x) for x in neighbors_in_X]
                        edges_to_add = [(y, x) for x in non_neighbors_in_X]
                        
                        g_prime.remove_edges_from(edges_to_remove)
                        g_prime.add_edges_from(edges_to_add)

                    # Check for non-isomorphism
                    if not nx.is_isomorphic(self.graph, g_prime):
                        log.append(f"SUCCESS: Found non-isomorphic cospectral mate G' with switching set X={X}.")
                        self.results['gm_switching_analysis']['found_cospectral_mate'] = True
                        self.results['gm_switching_analysis']['mate_graph'] = g_prime
                        self.results['gm_switching_analysis']['switching_partition'] = (X, Y)
                        return # Terminate on first find

        log.append("Search completed. No non-isomorphic cospectral mate found via GM switching.")

    def run_full_analysis(self, **kwargs):
        """
        Runs the complete analysis pipeline.
        """
        print("Step 1: Computing structural invariants...")
        self.compute_structural_invariants()
        print("Step 2: Computing spectral invariants...")
        self.compute_spectral_invariants()
        print("Step 3: Searching for Godsil-McKay switching partitions...")
        self.find_gm_switching_partitions(**kwargs)
        print("Analysis complete.")

    def generate_report(self):
        """
        Prints a summary report of the analysis.
        """
        #... (Implementation to format and print self.results)...
        # This method would format the data from the self.results dictionary
        # into a human-readable report similar to the structure of this document.
        print("\n--- DS Analysis Report ---")
        # Print structural properties
        # Print spectral properties
        # Print GM switching results and final conclusion
        pass