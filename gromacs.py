"""
This file generates the data for plotting
"""
import logging
import docker

from helpers import find_folders_with_same_pattern

logging.basicConfig(
    format='[%(asctime)s] - %(levelname)s - %(filename)s - %(funcName)s:%(lineno)d - %(message)s',
    level=logging.INFO
)

class GromacsData:
    """
    This class manages all the functions to generate Gromacs data.
    """
    CONTAINER_WORK_DIR = "/container/data"

    def __init__(self, path: str):
        self.path = path
        self.volume = {
            self.path : {
                'bind': self.CONTAINER_WORK_DIR,
                'mode': 'rw'
            }
        }

    def _get_last_tpr_file(self):
        last_folder = find_folders_with_same_pattern(self.path, "my_output_*")[-1]
        tpr_file = sorted(last_folder.glob('*.tpr'))[0].name
        return [last_folder.name, tpr_file]

    def _run_gromacs_container(
            self,
            command,
            working_dir=CONTAINER_WORK_DIR,
            **kwargs
        ):
        """
        Run the Gromacs container with the specified command and arguments.
        """
        gromacs_image = "gromacs/gromacs:2022.2"
        client = docker.from_env()
        return client.containers.run(
            gromacs_image,
            command,
            working_dir=working_dir,
            **kwargs
        )

    def generate_minimization_data(self):
        """
        Function to generate the minimization energy data
        """
        command = [
            "sh",
            "-c",
            "echo 10 0 | gmx energy -f em.edr -o potential.xvg"
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
            "echo 16 0 | gmx energy -f nvt.edr -o temperature.xvg"
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
            "echo 18 0 | gmx energy -f npt.edr -o pressure.xvg"
        ]

        logging.info(command_pressure[-1])
        pressure = self._run_gromacs_container(
            command_pressure,
            volumes=self.volume
        ).decode("utf-8")

        command_density = [
            "sh",
            "-c",
            "echo 24 0 | gmx energy -f npt.edr -o density.xvg"
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
        steps = find_folders_with_same_pattern(self.path, "my_output_*")

        xtc_files = ""
        for step in steps:
            filename = sorted(step.glob('*.xtc'))
            if not filename:
                raise ValueError("Could not find the xtc file")

            xtc_files += f"{step.name}/{filename[0].name} "

        command = ["sh", "-c", f"gmx trjcat -f {xtc_files} -o combined.xtc"]
        logging.info(command[-1])
        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def fix_periodicity(self):
        """
        Function to fix the periodicity of the trajectory.
        """
        folder, tpr = self._get_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"echo 1 0 | gmx trjconv -s {folder}/{tpr} -f combined.xtc -o combinedNoPBC.xtc -pbc mol -center" #pylint: disable=line-too-long
        ]

        logging.info(command[-1])
        periodicity = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return periodicity

    def generate_rmsd(self):
        """
        Generate rmsd data
        """
        folder, tpr = self._get_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"echo 4 4 | gmx rms -s {folder}/{tpr} -f combinedNoPBC.xtc -o rmsd.xvg -tu ns"
        ]

        logging.info(command[-1])
        rmsd = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return rmsd

    def generate_radius_gyration(self):
        """
        Generate radius of gyration        
        """
        folder, tpr = self._get_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"echo 1 | gmx gyrate -s {folder}/{tpr} -f combinedNoPBC.xtc -o gyrate.xvg"
        ]

        logging.info(command[-1])
        gyr = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode('utf-8')

        return gyr
