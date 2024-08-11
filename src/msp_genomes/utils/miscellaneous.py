"""Miscellaneious functions."""

from Bio import SeqIO
from pathlib import Path
import pandas as pd
from pandas import DataFrame
import importlib.resources as resources
import shutil
import subprocess

from msp_genomes.utils.config import MSP_COLLECTIONS, TMP_DIR

_msp_collection = pd.read_excel(MSP_COLLECTIONS)


def get_molecules_info_from_assembly(fasta_file: Path) -> dict:
    """Extract molecules information from fasta headers created by Unicycler."""
    molecules_info = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Check if molecule is circular.
        header = record.description
        if "circular=true" in header:
            circular = True
        else:
            circular = False
        # make header into list.
        header = header.split(" ")
        molecule = int(header[0])
        length = int(header[1].split("=")[1])
        depth = header[2].split("=")[1]
        # Store info in dictionary.
        molecules_info[molecule] = {
            "length": length,
            "depth": depth,
            "circular": circular,
        }

    return molecules_info


def compile_info_from_assemblies_into_dataframe(assemblies: Path) -> DataFrame:
    """Compile information from assembly files into a DataFrame."""
    # Dictionary to hold assemblies information.
    assemblies_info = {}
    counter = 0
    for assembly in assemblies.iterdir():
        assembly_path = assembly / "assembly.fasta"
        if assembly.is_dir() and assembly_path.exists():
            # Get strain ID from folder's name
            folder_name = assembly.name.split("_", 1)
            strain = folder_name[0]
            # Check if folder name has run info
            if len(folder_name) > 1:
                run_info = folder_name[1]
            else:
                run_info = "Not found"
            # Get strain's info
            strain_info = get_strain_info(strain, _msp_collection)
            # Get molecules info from assembly.fasta
            molecules = get_molecules_info_from_assembly(assembly_path)
            for value in molecules.values():
                assemblies_info[counter] = {
                    "Isolate ID": strain,
                    "Run info": run_info,
                    "Species": strain_info[0],
                    "Score": strain_info[1],
                    "Molecule size": value.get("length"),
                    "Topology": "circular" if value.get("circular") else "linear",
                }
                counter += 1
    assemblies_info = pd.DataFrame.from_dict(assemblies_info, orient="index")
    return assemblies_info


def make_tmp_directory(destiny_path: Path = TMP_DIR) -> None:
    """Make tmp directory if it doesn't exists"""
    if not destiny_path.exists():
        subprocess.run(["mkdir", destiny_path])


def clear_folder(folder_path: Path) -> None:
    """Removes all files and directories within the specified directory"""
    for item in folder_path.iterdir():
        try:
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        except Exception as e:
            print(f"Failed to delete {item}. Reason: {e}")


def rm_folder(folder_path: Path) -> None:
    """Remove folder"""
    if folder_path.is_dir():
        try:
            shutil.rmtree(folder_path)
        except Exception as e:
            print(f"Failed to delete {folder_path}. Reason: {e}")


def get_pkg_path(package: str = "msp_genomes") -> Path:
    """Get path to package directory in a src-layout"""
    return resources.files(f"{package}")


def get_strain_info(strain_id: str, msp_collection: DataFrame) -> tuple:
    """Find the species and MALDI-TOF score from a strain ID.

    This scrip uses the filtered_collections.xlsx database created from the MSP
    collections database.
    """
    species = msp_collection.loc[
        msp_collection["Isolate ID"] == strain_id, "Very best match"
    ]
    score = msp_collection.loc[
        msp_collection["Isolate ID"] == strain_id, "Very best score"
    ]
    if not species.empty:
        species = species.values[0]
        score = score.values[0]
    else:
        species = "Not found"
        score = "Not found"

    return (species, score)


def get_assemblies_info(
    assemblies: Path,
    output_folder: Path,
    compilation_output: Path,
    cge_script_name: str,
) -> dict:
    """Compile information from the assembly files and output folder.

    Makes a dictionary using isolate IDs keys and another dictionary as values. The keys
    of a dictionary associated with an isolate ID are 'assembly_path', 'species',
    'score', and 'output_folder'.

    Return example:

    {
        'SW3902_n1565-n1566_R129AB-NB03': {
            'strain': 'SW3902',
            'run_info': 'n1565-n1566_R129AB-NB03',
            'assembly_path': PosixPath('assemblies/SW3902_n1565/assembly.fasta'),
            'species': 'Escherichia coli',
            'score': 2.351,
            'output_folder': PosixPath(
                'assemblies/SW3902_n1565/SW3902_n1565-n1566_R129AB-NB03_resfinder'
            ),
        },
        'SW1327R_n2219_R129AB-NB11': {
            'strain': 'SW1327R',
            'run_info': 'n2219_R129AB-NB11',
            'assembly_path': PosixPath('assemblies/SW1327R_n1610/assembly.fasta')
            'species': 'Not found',
            'score': 'Not found',
            'output_folder': PosixPath(
                'assemblies/SW1327R_n1610/SW1327R_n2219_R129AB-NB11_resfinder'
            ),
        }
    }
    """
    # Dictionary to hold assemblies information.
    assemblies_info = {}
    for assembly in assemblies.iterdir():
        assembly_path = assembly / "assembly.fasta"
        if assembly.is_dir() and assembly_path.exists():
            folder_name = assembly.name.split("_", 1)
            strain = folder_name[0]
            # Check if folder name has run info
            if len(folder_name) > 1:
                run_info = folder_name[1]
            else:
                run_info = "Not found"
            strain_info = get_strain_info(strain, _msp_collection)
            assemblies_info[assembly.name] = {
                "strain": strain,
                "run_info": run_info,
                "assembly_path": assembly_path,
                "species": strain_info[0],
                "score": strain_info[1],
                "output_folder": make_output_path(
                    assemblies, assembly, output_folder, assembly.name, cge_script_name
                ),
            }
        elif assembly.is_dir() and not assembly_path.exists():
            message = f"`{assembly}` folder doesn't have `assembly.fasta`\n"
            log_writer(message, compilation_output)

    return assemblies_info


def log_writer(message: str, compilation_output: Path) -> None:
    """Write a message to the log file."""
    log_file = compilation_output / "log_file.txt"
    with open(log_file, "a") as log:
        log.write(message)


def make_output_path(
    assemblies: Path, assembly: Path, output_folder: Path, strain: str, script_name: str
) -> Path:
    if assemblies == output_folder:
        output_path = assembly / f"{strain}_{script_name}"
    else:
        output_path = output_folder / f"{strain}_{script_name}"
    return output_path


def make_output_folders(assemblies_info: dict) -> None:
    for strain in assemblies_info:
        subprocess.run(["mkdir", assemblies_info[strain]["output_folder"]])


if __name__ == "__main__":
    # print(get_pkg_path())
    # clear_folder(get_pkg_path() / "tmp")
    info = get_assemblies_info(
        Path(
            "/Users/msp/Documents/Coding/python_projects/MSP_genomes_project/assemblies"
        ),
        Path(
            "/Users/msp/Documents/Coding/python_projects/MSP_genomes_project/assemblies"
        ),
        Path(
            "/Users/msp/Documents/Coding/python_projects/MSP_genomes_project/assemblies"
        ),
        "test",
    )
    print(info)
    # for isolale, info in info.items():
    #     print(f"{isolale} -> {info}\n")
