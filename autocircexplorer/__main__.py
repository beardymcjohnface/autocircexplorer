"""
Entrypoint for autoCircExplorer

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("autocircexplorer.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("autocircexplorer.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="autocircexplorer.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=40, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--apptainer-dir",
            default="apptainer",
            show_default=False,
            callback=default_to_output,
            help="Custom apptainer directory [default: (outputDir)/apptainer/]",
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
                "--use-singularity",
                "--singularity-args \"-B `pwd`\""
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="autocircexplorer.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Run CIRCexplorer2 but better
    \b
    For more options, run:
    autoce2 command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
autoce2 run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           autoce2 run --input [file]
Specify threads:    autoce2 run ... --threads [threads]
Disable conda:      autoce2 run ... --no-use-conda 
Change defaults:    autoce2 run ... --snake-default="-k --nolock"
Add Snakemake args: autoce2 run ... --dry-run --keep-going --touch
Specify targets:    autoce2 run ... all print_targets
Available targets:
    all             Run everything (default)
    star            Align with STAR
    salmon          Count with salmon
    rmats           Run rmats-turbo
    ce2             Run circexplorer2
    ciri2           Run CIRI2
    multiqc         Make QC of rawreads and STAR alignments
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input1", help="Input file/directory sample group 1", type=str, required=True)
@click.option("--input2", help="Input file/directory sample group 2", type=str, required=False)
@click.option("--ref", help="STAR reference directory", type=str, required=False)
@click.option("--gtf", help="Reference GTF file", type=str, required=False)
@click.option("--fa", help="Reference FASTA file", type=str, required=False)
@click.option("--readlen", help="Read length", type=int, required=False, default=150)
@common_options
def run(**kwargs):
    """Run autoCircExplorer"""
    # Config to add or update in configfile
    merge_config = {
        "autocircexplorer": {
            "args": kwargs
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@common_options
def config(configfile, system_config, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=system_config)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
