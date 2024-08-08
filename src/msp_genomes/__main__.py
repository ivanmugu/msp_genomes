from msp_genomes.utils.get_cli import parse_command_line_input
from msp_genomes.utils.consolidate_compilations import consolidate
from msp_genomes.utils.config import OUTPUT_NAMES_COMPILATIONS
from msp_genomes.find_resistances.run_resfinder import find_resistances
from msp_genomes.find_plasmids.run_plasmidfinder import find_plasmids
from msp_genomes.find_st.run_mlst import find_st
from msp_genomes.find_virulence.run_virulencefinder import find_virulence
from msp_genomes.find_mge.run_mge_finder import find_mge


def run_find_resistances():
    cli = parse_command_line_input("find_resistances")
    find_resistances(cli)


def run_find_plasmids():
    cli = parse_command_line_input("find_plasmids")
    find_plasmids(cli)


def run_find_st():
    cli = parse_command_line_input("find_st")
    find_st(cli)


def run_find_virulence():
    cli = parse_command_line_input("find_virulence")
    find_virulence(cli)


def run_find_mge():
    cli = parse_command_line_input("find_mge")
    find_mge(cli)


def find_all():
    cli = parse_command_line_input("find_all")
    # Find resistances
    find_resistances(cli)
    # Find plasmids
    find_plasmids(cli)
    # Find sequence types
    find_st(cli)
    # Find virulence genes
    find_virulence(cli)
    # Consolidate compilations.
    consolidate(
        cli["compilation_output"], OUTPUT_NAMES_COMPILATIONS, cli["output_folder"]
    )
