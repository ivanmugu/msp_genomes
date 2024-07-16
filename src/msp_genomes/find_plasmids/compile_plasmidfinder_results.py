"""Compile plasmidfinder results iterating over folders with plasmidfinder outputs.

This script uses `results_tab.tsv` files stored in a directory tree structured as shown
below.

```
assemblies/
├── SW0001/
│   └── SW0001_plasmidfinder/
│       └── results_tab.tsv
└── SW0002/
    └── SW0002_plasmidfinder/
        └── results_tab.tsv
`
"""

from pathlib import Path
import pandas as pd
from pandas import DataFrame

from msp_genomes.utils.get_cli import parse_command_line_input
from msp_genomes.utils.miscellaneous import (
    get_assemblies_info,
    compile_info_from_assemblies_into_dataframe,
)

_plasmidfinder_results_file = "results_tab.tsv"
_script_name = "plasmidfinder"


def extract_plasmidfinder_results_by_molecule_size(
    plasmidfinder_results_file_path: Path,
) -> dict:
    """Organize plasmid types from plasmidfinder results into a dictionary.

    This function uses the `results_tab.tsv` file.

    Return example:
    {
        3034: 'Col(pHAD28)',
        176653: 'IncFIB(pB171),IncFII(Yp)',
        121178: 'IncFII(Cf)'
    }
    """
    # Convert plasmid finder results file into DataFrame
    plasmidfinder_df = pd.read_csv(plasmidfinder_results_file_path, sep="\t")
    # Initialize a dictionary to store the results.
    results = {}
    # Iterate over DataFrame rows to get plasmmid types
    for _, row in plasmidfinder_df.iterrows():
        plasmid = row["Plasmid"]
        contig = row["Contig"]
        # Get molecule length.
        molecule = int(contig.split(" ")[1].split("=")[1])
        # If results doesn't have any molecule record, start it
        if not results.get(molecule):
            results[molecule] = plasmid
        # If results already has a molecule record, append plasmid type.
        else:
            results[molecule] += f",{plasmid}"
    return results


def compile_plasmidfinder_results_into_dataframe(strains_info: dict) -> DataFrame:
    """Compile plasmidfinder results from folders into a DataFrame."""
    results = {}  # To compile information
    counter = 0  # to use it as keys in results. It will help to make the DataFrame.
    # Iterate over plasmidfinder results files to compile the information.
    for strain, info in strains_info.items():
        # Get plasmids types from molecules.
        plasmids = extract_plasmidfinder_results_by_molecule_size(
            info["output_folder"] / _plasmidfinder_results_file
        )
        # Iterate over the plasmids dictionary.
        # key is molecule size and value is plasmid type.
        for key, value in plasmids.items():
            results[counter] = {
                "Run info": info["run_info"],
                "Isolate ID": strain,
                "Molecule size": key,
                "Plasmids": value,
            }
            counter += 1
    # Covert results into a DataFrame
    results = pd.DataFrame.from_dict(results, orient="index")

    return results


if __name__ == "__main__":
    # Get command line input
    cli = parse_command_line_input("find_plasmids")

    # Collect information to compile plasmidfinder results.
    strains_info = get_assemblies_info(
        cli["input_folder"],
        cli["output_folder"],
        cli["compilation_output"],
        _script_name,
    )

    # Compile plasmidfinder results
    compiled_results = compile_plasmidfinder_results_into_dataframe(strains_info)

    # print(compiled_results)

    # Compile info from the assembly files into a DataFrame
    df_assemblies = compile_info_from_assemblies_into_dataframe(cli["input_folder"])
    df_assemblies.sort_values(
        by=["Isolate ID", "Molecule size"], ascending=[True, False]
    )
    # print(df_assemblies)

    merged_df = pd.merge(
        df_assemblies, compiled_results, on=["Isolate ID", "Molecule size"], how="left"
    )
    print(merged_df)
