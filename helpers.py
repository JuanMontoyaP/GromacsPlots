"""
Helper functions
"""
from pathlib import Path
import numpy as np

def read_xvg_files(file_path: str):
    """
    Read data from xvg files.
    """
    data = []

    with open(file_path, 'r', encoding='utf8') as file:
        for line in file:

            if line.startswith('#') or line.startswith('@'):
                continue

            values = [float(val) for val in line.strip().split()]
            data.append(values)

    data_array = np.array(data, dtype=float)
    return data_array

def create_folder(path:str, folder_name:str="images"):
    """
    Create a folder in the given path
    """
    # Specify the path of the new folder
    folder_path = Path(f'{path}/{folder_name}')

    # Use the mkdir method to create the folder
    folder_path.mkdir(parents=True, exist_ok=True)

def find_folders_with_same_pattern(path, pattern):
    """
    Folder with the same pattern
    """
    search_directory = Path(path)
    search_pattern = sorted(search_directory.glob(pattern))

    return [folder_path for folder_path in search_pattern if folder_path.is_dir()]
