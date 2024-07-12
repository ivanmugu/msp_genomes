"""Consolidate compilations databases."""

import pandas as pd
from pandas import DataFrame
from pathlib import Path
from msp_genomes.utils.config import OUTPUT_NAMES_COMPILATIONS


def consolidate(
    path_to_compilations: Path,
    output_names_compilations: dict,
    output_folder: Path,
    output_name: str = "consolidated.xlsx",
    export_consolidation: bool = True,
) -> DataFrame:
    counter = 0
    for value in output_names_compilations.values():
        df = pd.read_excel(path_to_compilations / value)
        if counter == 0:
            merged_df = df
        else:
            merged_df = pd.merge(
                merged_df,
                df,
                on=["Isolate ID", "Molecule size", "Species", "Score", "Topology"],
                how="left",
            )
        counter += 1

    if export_consolidation:
        merged_df.to_excel(output_folder / output_name, index=False)

    return merged_df


if __name__ == "__main__":
    consolidated = consolidate(
        Path("/Users/msp/Documents/Coding/python_projects/MSP_genomes_project/results"),
        OUTPUT_NAMES_COMPILATIONS,
    )
    print(consolidated)
    consolidated.to_excel(
        Path(
            "/Users/msp/Documents/Coding/python_projects/MSP_genomes_project/results/consolidated.xlsx"
        ),
        index=False,
    )
