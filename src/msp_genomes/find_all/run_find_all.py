"""Run resfinder, plasmidfinder, mlst, and virulencefinder."""

from msp_genomes.find_plasmids import run_plasmidfinder
from msp_genomes.find_resistances import run_resfinder
from msp_genomes.find_st import run_mlst
from msp_genomes.find_virulence import run_virulencefinder
from msp_genomes.utils.get_cli import parse_command_line_input


def main():
    pass


if __name__ == "__main__":
    main()
