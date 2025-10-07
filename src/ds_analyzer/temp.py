

""" import click
import networkx as nx
from .analyzer import DSAnalyzer

@click.command()
@click.argument('graph_file', type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option('--max-set-size', default=8, help='Maximum size of the GM switching set to check.')
@click.option('--time-limit', default=300, help='Time limit in seconds for the GM search.')
def main(graph_file, max_set_size, time_limit):
    """
   # Analyzes a graph to assess its likelihood of being Determined by its Spectrum (DS).

   # GRAPH_FILE: Path to the graph file (e.g., in GML, GraphML, or edge list format).
    """
    try:
        # NetworkX can infer the file type from the extension for many formats
        click.echo(f"Loading graph from {graph_file}...")
        G = nx.read_gml(graph_file) # Or another appropriate nx.read_* function
    except Exception as e:
        click.echo(f"Error: Could not read graph file. {e}", err=True)
        return

    click.echo("Initializing DS Analyzer...")
    analyzer = DSAnalyzer(G)

    # Run the full analysis
    analyzer.run_full_analysis(max_set_size=max_set_size, time_limit=time_limit)

    # Print a summary report
    click.echo("\n--- DS Analysis Report ---")
    
    # Print structural properties
    click.echo("\n")
    for key, value in analyzer.results["graph_properties"].items():
        click.echo(f"{key}: {value}")

    # Print GM switching results and final conclusion
    click.echo("\n")
    gm_results = analyzer.results['gm_switching_analysis']
    if gm_results['found_cospectral_mate']:
        click.secho("\nConclusion: The graph is NOT a DS graph.", fg='red', bold=True)
        click.echo("A non-isomorphic cospectral mate was found.")
    else:
        click.secho("\nConclusion: The graph is LIKELY a DS graph.", fg='green', bold=True)
        click.echo("No cospectral mate was found via Godsil-McKay switching.")

if __name__ == '__main__':
    main() 
""" 
