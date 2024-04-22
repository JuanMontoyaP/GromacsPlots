"""
This script graphs the results of MD simulations
"""
#!/usr/bin/env python3 -u
import sys
import logging
import click

from gromacs import GromacsData
from plots import GromacsPlot

logging.basicConfig(
    format='[%(asctime)s] - %(levelname)s - %(filename)s - %(funcName)s:%(lineno)d - %(message)s',
    level=logging.DEBUG
)

class GromacsPlots:
    """
    This class generates the plots
    """
    def __init__(self, path: str, protein: str):
        self.path = path
        self.protein = protein

    def generate_gromacs_data(self):
        """
        This function generates the gromacs data needed for plotting
        """
        gromacs = GromacsData(self.path)

        logging.info('Generating minimization data')
        energy_minimization = gromacs.generate_minimization_data()
        logging.info(energy_minimization)

        logging.info('Generating NVT data')
        temperature_data = gromacs.generate_nvt_data()
        logging.info(temperature_data)

        logging.info('Generating NPT data')
        npt_data = gromacs.generate_npt_data()
        logging.info(npt_data)

        try:
            logging.info('Generating final xtc file')
            gromacs.generate_final_xtc_file()
        except ValueError as e:
            logging.error('Error generating final xtc file: %s', e)
            sys.exit(1)

        logging.info('Fixing periodicity')
        periodicity = gromacs.fix_periodicity()
        logging.info(periodicity)

        logging.info('Writing RMSD')
        rmsd = gromacs.generate_rmsd()
        logging.info(rmsd)

        logging.info('Writing Radius of gyration')
        gyr = gromacs.generate_radius_gyration()
        logging.info(gyr)

    def generate_gromacs_plots(self):
        """
        Generate plots.
        """
        plots = GromacsPlot(self.path, self.protein)

        logging.info("Plotting energy minimization")
        plots.plot_energy_minimization()

        logging.info("Plotting temperature")
        plots.plot_temperature()

        logging.info("Plotting pressure")
        plots.plot_pressure()

        logging.info("Plotting density")
        plots.plot_density()

        logging.info("Plotting RMSD")
        plots.plot_rmsd()

        logging.info("Plotting radius of gyration")
        plots.plot_radius_gyration()

    def Run(self, run_gromacs: bool):
        """
        This class run the necessary commands to generate the plots
        """

        if run_gromacs:
            self.generate_gromacs_data()

        logging.info('Generating gromacs plots')
        self.generate_gromacs_plots()

@click.command()
@click.option(
    '--path',
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True
    ),
    help="Path where the data files should be.",
    required=True
)
@click.option(
    '--protein',
    type=str,
    required=True,
    help="Name of the protein to be plotted"
)
@click.option(
    '--run_gromacs',
    is_flag=True,
    help='False if you dont need to generate gromacs data'
)
def cli(path: str, protein: str, run_gromacs: bool = False):
    """
    This is the cli for managing the command line interface
    """
    gromacs_plots = GromacsPlots(path, protein)
    gromacs_plots.Run(run_gromacs)

if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
