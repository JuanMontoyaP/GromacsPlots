"""
This file generates the data for plotting
"""
import logging
from pathlib import Path
import tarfile
import re
import docker
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms

from helpers import find_files_with_same_pattern, create_folder

logging.basicConfig(
    format='[%(asctime)s] - %(levelname)s - %(filename)s - %(funcName)s:%(lineno)d - %(message)s',
    level=logging.INFO
)


class GromacsData:
    """
    This class manages all the functions to generate Gromacs data.
    """
    CONTAINER_WORK_DIR = "/container/data"
    XTC_FILE = "final.xtc"
    XTC_NO_PBC_FILE = "fixed.xtc"

    def __init__(self, path: Path):
        self.path: Path = Path(path)
        self.volume: dict = {
            self.path: {
                'bind': self.CONTAINER_WORK_DIR,
                'mode': 'rw'
            }
        }

        create_folder(self.path)
        self.results_folder: Path = Path(f"{self.path}/results")

    def _run_gromacs_container(
        self,
        command,
        working_dir=CONTAINER_WORK_DIR,
        **kwargs
    ):
        """
        Runs a Gromacs container with the specified command.

        Parameters:
            command (str): The command to be executed inside the container.
            working_dir (str, optional): The working directory inside the container. Defaults
                to CONTAINER_WORK_DIR.
            **kwargs: Additional keyword arguments to be passed to the container run method.

        Returns:
            docker.models.containers.Container: The container object representing the
                running Gromacs container.
        """
        gromacs_image = "jpmontoya19/gromacs:latest"
        client = docker.from_env()
        return client.containers.run(
            gromacs_image,
            command,
            working_dir=working_dir,
            **kwargs
        )

    def _find_last_tpr_file(self) -> Path:
        """
        Finds the last TPR file in the given path.

        Returns:
            A tuple containing the name and path of the last TPR file.
        """
        last_tpr = find_files_with_same_pattern(self.path, "md_0_*.tpr")[-1]
        return (last_tpr.name, last_tpr)

    def generate_minimization_data(self):
        """
        Function to generate the minimization energy data
        """
        command = [
            "sh",
            "-c",
            f"echo 10 0 | gmx energy -f em.edr -o {self.results_folder.name}/potential.xvg"
        ]

        logging.info(command[-1])
        energy_minimization = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")
        return energy_minimization

    def generate_nvt_data(self):
        """
        Function to generate the temperature data
        """
        command = [
            "sh",
            "-c",
            f"""
                echo 16 0 | \
                gmx energy -f nvt.edr -o {self.results_folder.name}/temperature.xvg
            """
        ]

        logging.info(command[-1])
        temperature = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")
        return temperature

    def generate_npt_data(self):
        """
        Function to generate the pressure & density data
        """
        command_pressure = [
            "sh",
            "-c",
            f"""
                echo 18 0 | \
                gmx energy -f npt.edr -o {self.results_folder.name}/pressure.xvg
            """
        ]

        logging.info(command_pressure[-1])
        pressure = self._run_gromacs_container(
            command_pressure,
            volumes=self.volume
        ).decode("utf-8")

        command_density = [
            "sh",
            "-c",
            f"""
                echo 24 0 | \
                gmx energy -f npt.edr -o {self.results_folder.name}/density.xvg
            """
        ]

        logging.info(command_density[-1])
        density = self._run_gromacs_container(
            command_density,
            volumes=self.volume
        ).decode("utf-8")

        return f"{pressure}, \n, {density}"

    def generate_final_xtc_file(self):
        """
        Generates the final xvg file joining all the parts of the trajectory
        """
        steps = find_files_with_same_pattern(
            self.path, "my_job.output_*.tar.gz")

        xtc_files = []
        for step in steps:
            with tarfile.open(step, "r:gz") as tar:
                for member in tar.getmembers():
                    if member.name.endswith(".xtc"):
                        tar.extract(member, path=f"{self.path}")
                        xtc_files.append(member.name)

        command = [
            "sh",
            "-c",
            f"gmx trjcat -f {' '.join(xtc_files)} -o {self.results_folder.name}/{self.XTC_FILE}"
        ]

        logging.info(command[-1])
        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def fix_periodicity(self):
        """
        Function to fix the periodicity of the trajectory.
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 1 0 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_FILE} \
                    -o {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -pbc mol -center
            """
        ]

        logging.info(command[-1])
        periodicity = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return periodicity

    def generate_video(self):
        """
        Generate a video from the trajectory
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 1 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -o {self.results_folder.name}/reduced.xtc \
                    -dt 250
            """
        ]

        logging.info(command[-1])
        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def generate_radius_gyration(self):
        """
        Generate radius of gyration        
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 1 | \
                gmx gyrate \
                -s {filename} \
                -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                -o {self.results_folder.name}/gyrate.xvg
            """
        ]

        logging.info(command[-1])
        gyr = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode('utf-8')

        return gyr

    def get_initial_configuration(self):
        """
        Retrieves the initial configuration of the system.

        Returns:
            None
        """
        filename, _ = self._find_last_tpr_file()

        command = [
            "sh",
            "-c",
            f"""
                echo 1 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -o {self.results_folder.name}/initial_conf.pdb \
                    -dump 0
            """
        ]

        logging.info(command[-1])

        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def get_md_frames(self) -> str:
        """
        Generate a doc with the number of frames in the trajectory
        """
        command = [
            "sh",
            "-c",
            f"""
                echo 1 | \
                gmx check \
                -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                > {self.results_folder.name}/trajectory_info.txt 2>&1
            """
        ]

        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

        with open(
                f"{self.results_folder}/trajectory_info.txt",
                "r",
                encoding='utf-8'
            ) as file:
            lines = file.readlines()

        last_frame = None
        for line in lines:
            match = re.search(r'Last frame\s+(\d+(\.\d+)?)', line)

            if match:
                last_frame = match.group(1)
                break

        logging.info("Last frame: %s", last_frame)
        return last_frame
    
    def generate_rmsd(self):
        """"
        Generate the RMSD data
        """
        print(self.path)

        u = mda.Universe(
            f"{self.path}/npt.gro",
            f"{self.path}/results/fixed.xtc",
            topology_format="GRO",
            trajectory_format="XTC"
        )

        total_residues = u.select_atoms('protein').residues.n_residues

        logging.info(f"Total residues: {total_residues}")

        selections = {
            'Backbone': f"backbone and not (resid {total_residues-200} to {total_residues}) and not (resid 1 to 211)",
            'Head': "resid 1 to 211",
            'CRD1': "resid 225 to 341",
            'CRD2': "resid 369 to 487",
            'CRD3': "resid 511 to 626",
            'CRD4': "resid 655 to 778",
            'CRD5': "resid 807 to 923",
            'CRD6': "resid 951 to 1079",
            'CRD7': "resid 1101 to 1212",
            'CRD8': "resid 1240 to 1355",
            'Tail': f"resid {total_residues-200} to {total_residues}",
        }

        df = pd.DataFrame()

        for domain, selection in selections.items():
            ref = u.select_atoms(selection)
            mobile = u.select_atoms(selection)

            R = rms.RMSD(u, ref, select=selection, ref_frame=0)
            R.run()

            print(f"Average RMSD for {domain}: {np.mean(R.rmsd[:, 2]):.2f} Å")

            df['Time (ps)'] = R.rmsd[:, 1]
            df[domain] = R.rmsd[:, 2]

        df.to_csv(f"{self.results_folder}/rmsd.csv", index=True)

        with open(f"{self.results_folder}/rmsd_summary.txt", "w") as file:
            file.write(df.describe().to_string())

        return df
