import os
import click
import networkx as nx
from .analyzer import DSAnalyzer  # or 'from analyzer import DSAnalyzer' if run as standalone script

@click.command()
@click.argument('graph_file', type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option('--max-set-size', default=8, help='Maximum size of the GM switching set to check.')
@click.option('--time-limit', default=300, help='Time limit in seconds for the GM search.')
def main(graph_file, max_set_size, time_limit):
    """
    Analyzes a graph to assess its likelihood of being Determined by its Spectrum (DS).

    GRAPH_FILE: Path to the graph file (e.g., in GML, GraphML, or edge list format).
    """
    try:
        click.secho(f"Loading graph from {graph_file}...", fg='cyan')
        ext = os.path.splitext(graph_file)[1].lower()
        if ext == '.gml':
            G = nx.read_gml(graph_file)
        elif ext in ['.graphml', '.xml']:
            G = nx.read_graphml(graph_file)
        elif ext in ['.edgelist', '.txt']:
            G = nx.read_edgelist(graph_file)
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    except Exception as e:
        click.secho(f"Error: Could not read graph file. {e}", fg='red', err=True)
        return

    click.secho("Initializing DS Analyzer...", fg='yellow')
    analyzer = DSAnalyzer(G)

    # Run the full analysis
    analyzer.run_full_analysis(max_set_size=max_set_size, time_limit=time_limit)

    # Print report
    click.echo("\n--- DS Analysis Report ---")

    click.echo("\n[Structural Properties]")
    for key, value in analyzer.results["graph_properties"].items():
        click.echo(f"{key}: {value}")

    click.echo("\n[Spectral Properties]")
    for key, value in analyzer.results["spectral_properties"].items():
        click.echo(f"{key}: {value}")

    click.echo("\n[Godsil–McKay Switching Analysis]")
    gm_results = analyzer.results['gm_switching_analysis']
    for entry in gm_results['analysis_log']:
        click.echo(f"  - {entry}")

    if gm_results['found_cospectral_mate']:
        click.secho("\nConclusion: The graph is NOT a DS graph.", fg='red', bold=True)
        click.echo("A non-isomorphic cospectral mate was found.")
    else:
        click.secho("\nConclusion: The graph is LIKELY a DS graph.", fg='green', bold=True)
        click.echo("No cospectral mate was found via Godsil–McKay switching.")

if __name__ == '__main__':
    main()
