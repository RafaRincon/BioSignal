"""
BioSignal Discovery Engine
CLI Principal
==============
Interfaz de línea de comandos para ejecutar el pipeline completo
o agentes individuales.

Uso:
    python -m biosignal.cli run --disease 'pancreatic cancer'
    python -m biosignal.cli run --disease 'alzheimer disease' --max-datasets 10
    python -m biosignal.cli status
    python -m biosignal.cli list-diseases
"""

import click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint

console = Console()


@click.group()
@click.version_option(version="1.0.0", prog_name="biosignal")
def cli():
    """
    🧬 BioSignal Discovery Engine v1.0.0

    Plataforma de descubrimiento biológico basada en transcriptómica
    multi-dataset e inteligencia artificial.
    """
    pass


@cli.command()
@click.option("--disease", "-d", required=True,
              help="Nombre de la enfermedad a analizar (ej. 'pancreatic cancer')")
@click.option("--max-datasets", default=15, show_default=True,
              help="Número máximo de datasets a analizar")
@click.option("--min-samples", default=10, show_default=True,
              help="Mínimo de muestras por dataset")
@click.option("--parallel", default=4, show_default=True,
              help="Workers paralelos para descarga")
@click.option("--output-dir", default="data/", show_default=True,
              help="Directorio base de outputs")
@click.option("--report-format", default="pdf,json,csv", show_default=True,
              help="Formatos de reporte (pdf,json,csv,markdown)")
@click.option("--start-agent", default=1, show_default=True,
              help="Número de agente desde donde iniciar (1-8)")
@click.option("--stop-agent", default=8, show_default=True,
              help="Número de agente donde detener (1-8)")
@click.option("--config", default="config/settings.yaml", show_default=True,
              help="Path al archivo de configuración")
def run(disease, max_datasets, min_samples, parallel, output_dir,
        report_format, start_agent, stop_agent, config):
    """
    Ejecuta el pipeline completo de análisis para una enfermedad.

    Ejemplo:

        python -m biosignal.cli run --disease 'pancreatic cancer' --max-datasets 15

        python -m biosignal.cli run \\
            --disease 'alzheimer disease' \\
            --max-datasets 20 \\
            --parallel 6 \\
            --report-format pdf,json,csv
    """
    console.print(Panel.fit(
        f"[bold green]🧬 BioSignal Discovery Engine[/bold green]\n"
        f"Analizando: [bold cyan]{disease}[/bold cyan]",
        border_style="green"
    ))

    from biosignal.pipeline import BioSignalPipeline

    pipeline = BioSignalPipeline(config_path=config)
    results = pipeline.run(
        disease=disease,
        max_datasets=max_datasets,
        min_samples=min_samples,
        parallel=parallel,
        output_dir=output_dir,
        report_format=report_format,
        start_from_agent=start_agent,
        stop_after_agent=stop_agent,
    )

    # Mostrar resumen final
    table = Table(title="Resumen del Pipeline", show_header=True, header_style="bold magenta")
    table.add_column("Agente", style="cyan")
    table.add_column("Estado", style="green")
    table.add_column("Tiempo (s)")

    for agent_key, agent_result in results.get("agents", {}).items():
        status = agent_result.get("status", "unknown")
        elapsed = str(agent_result.get("elapsed_s", "—"))
        icon = "✓" if status == "success" else ("⚠" if status.startswith("skipped") else "○")
        status_style = "green" if status == "success" else "yellow"
        table.add_row(agent_key, f"[{status_style}]{icon} {status}[/{status_style}]", elapsed)

    console.print(table)
    console.print(f"\n[bold]Tiempo total:[/bold] {results.get('total_elapsed_s', 0)/60:.1f} minutos")
    console.print(f"[bold]Outputs en:[/bold] {output_dir}")


@cli.command()
def list_diseases():
    """Lista las enfermedades con aliases predefinidos."""
    import json
    from pathlib import Path

    aliases_path = Path("config/disease_aliases.json")
    if not aliases_path.exists():
        console.print("[red]No encontrado: config/disease_aliases.json[/red]")
        return

    with open(aliases_path) as f:
        aliases = json.load(f)

    for category, diseases in aliases.items():
        if category.startswith("_"):
            continue
        table = Table(title=f"📂 {category.upper()}", show_header=True, header_style="bold blue")
        table.add_column("Nombre coloquial", style="cyan")
        table.add_column("Término GEO")

        if isinstance(diseases, dict):
            for name, term in diseases.items():
                table.add_row(name, term)
        console.print(table)
        console.print()


@cli.command()
@click.option("--output-dir", default="data/", show_default=True)
def status(output_dir):
    """Muestra el estado de los datos procesados en el directorio de outputs."""
    from pathlib import Path

    base = Path(output_dir)
    subdirs = ["discovery", "raw", "processed", "dea", "meta", "pathways", "insights", "reports"]

    table = Table(title="Estado del Pipeline", show_header=True)
    table.add_column("Directorio")
    table.add_column("Archivos")
    table.add_column("Tamaño")

    for subdir in subdirs:
        path = base / subdir
        if path.exists():
            files = list(path.rglob("*"))
            file_count = sum(1 for f in files if f.is_file())
            total_size = sum(f.stat().st_size for f in files if f.is_file())
            size_str = f"{total_size / 1024 / 1024:.1f} MB" if total_size > 0 else "—"
            table.add_row(f"data/{subdir}/", str(file_count), size_str)
        else:
            table.add_row(f"data/{subdir}/", "[dim]no existe[/dim]", "—")

    console.print(table)


@cli.command()
def doctor():
    """Verifica que todas las dependencias estén instaladas correctamente."""
    console.print("\n[bold]🔍 Verificando dependencias...[/bold]\n")

    checks = [
        ("Python ≥ 3.11", "python", lambda: __import__("sys").version_info >= (3, 11)),
        ("pandas", "pandas", lambda: bool(__import__("pandas"))),
        ("numpy", "numpy", lambda: bool(__import__("numpy"))),
        ("scipy", "scipy", lambda: bool(__import__("scipy"))),
        ("biopython", "Bio", lambda: bool(__import__("Bio"))),
        ("GEOparse", "GEOparse", lambda: bool(__import__("GEOparse"))),
        ("pydeseq2", "pydeseq2", lambda: bool(__import__("pydeseq2"))),
        ("gseapy", "gseapy", lambda: bool(__import__("gseapy"))),
        ("langchain", "langchain", lambda: bool(__import__("langchain"))),
        ("pydantic v2", "pydantic", lambda: int(__import__("pydantic").VERSION.split(".")[0]) >= 2),
        ("rpy2", "rpy2", lambda: bool(__import__("rpy2"))),
        ("scikit-learn", "sklearn", lambda: bool(__import__("sklearn"))),
        ("loguru", "loguru", lambda: bool(__import__("loguru"))),
        ("rich", "rich", lambda: bool(__import__("rich"))),
    ]

    all_ok = True
    for name, module, check_fn in checks:
        try:
            ok = check_fn()
            icon = "[green]✓[/green]" if ok else "[red]✗[/red]"
            console.print(f"  {icon} {name}")
            if not ok:
                all_ok = False
        except ImportError:
            console.print(f"  [red]✗[/red] {name} [dim](no instalado)[/dim]")
            all_ok = False

    console.print()
    if all_ok:
        console.print("[bold green]✓ Todas las dependencias OK![/bold green]")
    else:
        console.print("[bold red]✗ Faltan dependencias. Ejecuta:[/bold red]")
        console.print("  pip install -r requirements.txt")
        console.print("  Rscript install_r_packages.R")


if __name__ == "__main__":
    cli()