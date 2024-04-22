"""
This file contains the calls to graph the gromacs data
"""
from matplotlib import pyplot as plt

from helpers import read_xvg_files, create_folder

class GromacsPlot:
    """
    Class for plotting all the data
    """

    def __init__(self, path, protein):
        self.path = path
        self.protein = protein
        self.images_path = f"{self.path}/images/"

        create_folder(self.path)

    def plot_data(self, x, y, titles):
        """
        Plot the specified data
        """
        fig, ax = plt.subplots()

        ax.plot(x, y, color='black', linewidth=0.7)

        # Add labels and title
        plt.xlabel(titles[0])
        plt.ylabel(titles[1])
        plt.title(titles[2])

        if len(titles) == 4:
            plt.legend(titles[-1])

        plt.xlim(left=x[0], right=x[-1])

        plt.minorticks_on()

        plt.tick_params(axis='both', which='both',
                        direction='in', right=True, top=True)

        return fig

    def plot_energy_minimization(self):
        """
        Plot the energy minimization function
        """
        energy_data = read_xvg_files(f"{self.path}/potential.xvg")

        _ = self.plot_data(
            energy_data[:, 0]/1000,
            energy_data[:, 1],
            [
                'Time $(ns)$',
                'Potential Energy $(kJ/mol)$',
                self.protein
            ]
        )

        plt.suptitle('Energy Minimization', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "minimization.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_temperature(self):
        """
        Plot temperature in the NVT equilibration
        """
        temperature = read_xvg_files(f"{self.path}/temperature.xvg")

        _ = self.plot_data(
            temperature[:, 0],
            temperature[:, 1],
            [
                'Time $(ps)$',
                'Potential Energy $(K)$',
                f"{self.protein}, NVT Equilibration"
            ]
        )

        plt.suptitle('Temperature', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "temperature.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_pressure(self):
        """
        Plot pressure in NPT equilibration
        """
        pressure = read_xvg_files(f"{self.path}/pressure.xvg")

        _ = self.plot_data(
            pressure[:, 0],
            pressure[:, 1],
            [
                'Time $(ps)$',
                'Pressure $(bar)$',
                f"{self.protein}, NPT Equilibration"
            ]
        )

        plt.suptitle('Pressure', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "pressure.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_density(self):
        """
        Plot density in NPT equilibration
        """
        density = read_xvg_files(f"{self.path}/density.xvg")

        _ = self.plot_data(
            density[:, 0],
            density[:, 1],
            [
                'Time $(ps)$',
                'Density $(kg/m^{3})$',
                f"{self.protein}, NPT Equilibration"
            ]
        )

        plt.suptitle('Density', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "density.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_rmsd(self):
        """
        Plot RMSD
        """
        rmsd = read_xvg_files(f"{self.path}/rmsd.xvg")

        _ = self.plot_data(
            rmsd[:, 0],
            rmsd[:, 1],
            [
                'Time $(ns)$',
                'RMSD $(nm)$',
                f"{self.protein}, Backbone"
            ]
        )

        plt.suptitle('RMSD', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "rmsd.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_radius_gyration(self):
        """
        PLot radius of gyration
        """
        gyr = read_xvg_files(f"{self.path}/gyrate.xvg")

        _ = self.plot_data(
            gyr[:, 0]/1000,
            gyr[:, 1],
            [
                'Time $(ns)$',
                '$R_{g}$ $(nm)$',
                f"{self.protein}, Unrestrained MD"
            ]
        )

        plt.suptitle('Radius of gyration', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "gyration.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )
