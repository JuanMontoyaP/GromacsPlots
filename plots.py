"""
This file contains the calls to graph the gromacs data
"""
import pandas as pd
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

        create_folder(self.path, "images")

    def plot_data(self, x, y, titles):
        """
        Plot the specified data
        """
        fig, ax = plt.subplots()

        ax.plot(x, y, linewidth=0.7)

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

        return fig, ax

    def plot_energy_minimization(self):
        """
        Plot the energy minimization function
        """
        energy_data = read_xvg_files(f"{self.path}/results/potential.xvg")

        _, ax = self.plot_data(
            energy_data[:, 0]/1000,
            energy_data[:, 1],
            [
                'Time $(ns)$',
                'Potential Energy $(kJ/mol)$',
                self.protein
            ]
        )

        ax.lines[0].set_color('black')

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
        temperature = read_xvg_files(f"{self.path}/results/temperature.xvg")

        _, ax = self.plot_data(
            temperature[:, 0],
            temperature[:, 1],
            [
                'Time $(ps)$',
                'Potential Energy $(K)$',
                f"{self.protein}, NVT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

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
        pressure = read_xvg_files(f"{self.path}/results/pressure.xvg")

        _, ax = self.plot_data(
            pressure[:, 0],
            pressure[:, 1],
            [
                'Time $(ps)$',
                'Pressure $(bar)$',
                f"{self.protein}, NPT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

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
        density = read_xvg_files(f"{self.path}/results/density.xvg")

        _, ax = self.plot_data(
            density[:, 0],
            density[:, 1],
            [
                'Time $(ps)$',
                'Density $(kg/m^{3})$',
                f"{self.protein}, NPT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

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
        rmsd = pd.read_csv(f"{self.path}/results/rmsd.csv")

        fig, ax = self.plot_data(
            rmsd['Time (ps)'].values/1000,
            list(zip(
                rmsd['Backbone'].values/10,
                rmsd['CRD1'].values/10,
                rmsd['CRD2'].values/10,
                rmsd['CRD3'].values/10,
                rmsd['CRD4'].values/10,
                rmsd['CRD5'].values/10,
                rmsd['CRD6'].values/10,
                rmsd['CRD7'].values/10,
                rmsd['CRD8'].values/10
            )),
            [
                'Time $(ns)$',
                'RMSD $(nm)$',
                f"{self.protein}, Backbone",
            ]
        )

        plt.legend(
            [
                'Backbone',
                'CRD1',
                'CRD2',
                'CRD3',
                'CRD4',
                'CRD5',
                'CRD6',
                'CRD7',
                'CRD8'
            ],
            loc='center left',
            bbox_to_anchor=(1, 0.5),
            ncol=1,
            fancybox=True,
            fontsize='small'
        )
        
        plt.ylim(0, (rmsd['Backbone'].max()+2)/10)

        plt.suptitle('RMSD', fontsize=20, y=1)
        
        plt.tight_layout()
    

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
        gyr = read_xvg_files(f"{self.path}/results/gyrate.xvg")

        _, ax = self.plot_data(
            gyr[:, 0]/1000,
            gyr[:, 1],
            [
                'Time $(ns)$',
                '$R_{g}$ $(nm)$',
                f"{self.protein}, Unrestrained MD"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            min(gyr[:, 1]) - 0.2,
            max(gyr[:, 1]) + 0.2
        )

        plt.suptitle('Radius of gyration', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "gyration.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )
