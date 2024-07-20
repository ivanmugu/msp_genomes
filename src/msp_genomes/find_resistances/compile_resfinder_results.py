"""Compile resfinder results iterating over folders with resfinder outputs.

This script uses `Resfinder_results_tab.txt` files stored in a directory tree structured
as shown below.

```
assemblies/
├── SW0001/
│   └── SW0001_resfinder/
│       └── ResFinder_results_tab.txt
└── SW0002/
    └── SW0002_resfinder/
        └── ResFinder_results_tab.txt
```
"""

import pandas as pd
from pandas import DataFrame
import sys
from pathlib import Path

from msp_genomes.utils.get_cli import parse_command_line_input
from msp_genomes.utils.miscellaneous import (
    get_assemblies_info,
    compile_info_from_assemblies_into_dataframe,
    rm_folder,
)
from msp_genomes.utils.config import PHENOTYPES_TABLE

_resfinder_results_file = "ResFinder_results_tab.txt"


def split_gene_accession(gene_accession: str) -> tuple[str]:
    """Split Gene_accession no. column from phenotypes.txt"""
    underscores = gene_accession.count("_")
    items = gene_accession.split("_")
    if underscores == 3:
        gene = items[0]
        acc_no = "_".join(items[-2:])
        return (gene, acc_no)
    elif underscores == 2:
        gene = items[0]
        acc_no = items[-1]
        return (gene, acc_no)
    else:
        sys.exit(f"Error with the Gene_accession no.: {gene_accession}")


def get_unique_antibiotic_classes(df: pd.DataFrame) -> list[str]:
    """Get unique antibiotic classes from class column."""
    antibiotic_classes = set()
    # Get Class column from df.
    for antibiotic_class in df["Class"]:
        # Remove spaces and convert string to list.
        antibiotic_class = antibiotic_class.replace(", ", ",").split(",")
        # Update the set.
        antibiotic_classes.update(antibiotic_class)
    # Convert set to list.
    antibiotic_classes = list(antibiotic_classes)
    # Return the sorted list.
    return sorted(antibiotic_classes)


def get_unique_and_composed_antibiotic_classes(df: pd.DataFrame) -> list[str]:
    """Get unique and composed antibiotic classes from class column."""
    antibiotic_classes = set()
    for antibiotic_class in df["Class"]:
        antibiotic_classes.update([antibiotic_class])
    antibiotic_classes = list(antibiotic_classes)
    # # I will add Beta-lactam in a different position when making the DataFrame
    # antibiotic_classes.remove("Beta-lactam")
    return sorted(antibiotic_classes)


def make_phenotypes_df(phenotypes_table: Path = PHENOTYPES_TABLE) -> DataFrame:
    """Make a DataFrame using the phenotypes.txt from resfinder_db.

    This DataFrame will help to get the antibiotic classes of ARGs.
    """
    phenotypes_df = pd.read_csv(phenotypes_table, delimiter="\t")
    # Add Accession no. and Gene to phenotype_df
    genes = []
    accessions = []
    for _, row in phenotypes_df.iterrows():
        genes.append(split_gene_accession(row["Gene_accession no."])[0])
        accessions.append(split_gene_accession(row["Gene_accession no."])[1])
    phenotypes_df["Accession no."] = accessions
    phenotypes_df["Gene"] = genes
    return phenotypes_df


def extract_resfinder_results_by_molecule_size(
    resfinder_results_file_path: Path,
    phenotypes_df: DataFrame,
) -> dict:
    """Organize genes from resfinder results by classes.

    This function uses the `ResFinder_results_tab.txt` file. Returns a dictionary with
    molecule sizes as keys and a dictionary with classes as value.

    Return example:
    {
        4784016: {
            'Aminoglycoside': "aadA1,aac(6')-Ib3",
            'Beta-lactam': 'blaMOX-6,blaOXA-427,blaOXA-10',
            'Folate pathway antagonist': 'dfrA14'
        },
        42792: {
            'Beta-lactam': 'blaKPC-2'
        }
    }
    """
    # Convert `ResFinder_results_tab.txt` into DataFrame
    resfinder_df = pd.read_csv(resfinder_results_file_path, sep="\t")
    # Initiate a dictionary to store the results.
    results = {}
    # Iterate over DataFrame rows to get ARGs
    for _, row in resfinder_df.iterrows():
        gene = row["Resistance gene"]
        contig = row["Contig"]
        # Get molecule length.
        molecule = int(contig.split(" ")[1].split("=")[1])
        # Get ARG class using phenotypes_df
        gene_class = phenotypes_df.loc[phenotypes_df["Gene"] == gene, "Class"].iloc[0]
        # If results doesn't have any molecule record, start it.
        if not results.get(molecule):
            results[molecule] = {gene_class: gene}
        # If results has a molecule record but not a class, add class.
        elif not results[molecule].get(gene_class):
            results[molecule][gene_class] = gene
        # If results has a molecule record, a class, but not a gene, add gene.
        elif results[molecule].get(gene_class) and (
            gene not in results[molecule][gene_class]
        ):
            results[molecule][gene_class] += f",{gene}"

    return results


def compile_resfinder_results_into_dataframe(
    strains_info: dict,
    extended_output: bool = True,
    include_all_phenotypes_headers: bool = True,
) -> DataFrame:
    """Compile resfinder results from folders into a DataFrame."""
    # Make a dataframe using phenotypes.txt from resfinder_db to find ARGs classes
    phenotypes_df = make_phenotypes_df()

    results = {}  # To compile information
    counter = 0  # To use it as keys in results. It will help to make the DataFrame.
    # Iterate over resfinder results files to compile the information.
    for info in strains_info.values():
        # Get ARGs from molecules.
        resistances = extract_resfinder_results_by_molecule_size(
            info["output_folder"] / _resfinder_results_file, phenotypes_df
        )
        # Iterater over the resistances dictionary.
        # key is molecule size and value is a dictionary with classes as keys and ARGs as
        # values.
        for key, value in resistances.items():
            results[counter] = {
                "Isolate ID": info["strain"],
                "Run info": info["run_info"],
                "Molecule size": key,
            }
            results[counter].update(value)
            counter += 1
        # If not extended extended output
        if not extended_output:
            rm_folder(info["output_folder"])

    # Convert results into a DataFrame.
    results = pd.DataFrame.from_dict(results, orient="index")

    if include_all_phenotypes_headers:
        headers = ["Isolate ID"]
        headers.extend(get_unique_and_composed_antibiotic_classes(phenotypes_df))
        results_all_headers = pd.DataFrame(columns=headers)
        results = pd.concat([results_all_headers, results], ignore_index=True)

    return results


if __name__ == "__main__":
    # Get command line input
    cli = parse_command_line_input("find_resistances")

    # Collect information to compile resfinder results.
    strains_info = get_assemblies_info(
        cli["input_folder"], cli["output_folder"], cli["compilation_output"]
    )

    # Compile resfinder results.
    compiled_results = compile_resfinder_results_into_dataframe(strains_info)

    # Compile info from the assembly files into a DataFrame
    df_assemblies = compile_info_from_assemblies_into_dataframe(cli["input_folder"])

    merged_df = pd.merge(
        df_assemblies, compiled_results, on=["Isolate ID", "Molecule size"], how="left"
    )
    merged_df.to_excel(Path("../../test.xlsx"), index=False)
    print(merged_df)
