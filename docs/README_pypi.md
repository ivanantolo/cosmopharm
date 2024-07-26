# COSMOPharm

Welcome to the COSMOPharm package, accompanying [our paper in *Molecular Pharmaceutics*](https://doi.org/10.1021/acs.molpharmaceut.4c00342). This project and its associated publication offer insights and a practical toolkit for researching drug-polymer and drug-solvent systems, aiming to provide the scientific community with the means to reproduce our findings and further the development of COSMO-SAC-based models.

<p align="center">
  <!-- <img src="https://github.com/usnistgov/COSMOSAC/raw/master/JCTC2020.PNG" alt="TOC Figure" width="500"> -->
  <img src="https://github.com/ivanantolo/cosmopharm/raw/main/TOC.png" alt="TOC Figure">
</p>

## About 

COSMOPharm is a Python package designed to streamline the predictive modeling of drug-polymer compatibility, crucial for the development of pharmaceutical amorphous solid dispersions. Apart from that, it can also be used for the miscibility/solubility of drugs with/in common solvents. Leveraging the COSMO-SAC (Conductor-like Screening Model Segment Activity Coefficient) model, COSMOPharm offers a robust platform for scientists and researchers to predict solubility, miscibility, and phase behavior in drug formulation processes.

## Features

- **Compatibility Prediction**: Utilize open-source COSMO-SAC model for prediction of drug-polymer compatibility.
- **Solubility Calculation**: Calculate drug-polymer solubilities to guide the selection of suitable polymers for drug formulations.
- **Miscibility and Phase Behavior Analysis**: Analyze the miscibility of drug-polymer pairs and understand their phase behavior under various conditions.
- **User-friendly Interface**: Easy-to-use functions and comprehensive documentation to facilitate research and development in pharmaceutical sciences.

## Installation

Install COSMOPharm with pip:

`pip install cosmopharm`

Ensure you have installed the `cCOSMO` library as per instructions on the [COSMOSAC GitHub page](https://github.com/usnistgov/COSMOSAC).

## Quick Start

Get started with COSMOPharm using the minimal example below, which demonstrates how to calculate the solubility of a drug in a polymer. This example succinctly showcases the use of COSMOPharm for solubility calculations:


```python
import matplotlib.pyplot as plt
import cCOSMO
from cosmopharm import SLE, COSMOSAC
from cosmopharm.utils import create_components, read_params

# Define components
names = ['SIM','PLGA50']
params_file = "data/sle/table_params.xlsx"

# Load parameters and create components
parameters = read_params(params_file)
mixture = create_components(names, parameters)

# Initialize COSMO-SAC model - replace paths with your local paths to COSMO profiles
db = cCOSMO.DelawareProfileDatabase(
    "./profiles/_import_methods/UD/complist.txt",
    "./profiles/_import_methods/UD/sigma3/")

for name in names:
    iden = db.normalize_identifier(name)
    db.add_profile(iden)
COSMO = cCOSMO.COSMO3(names, db)

# Setup the COSMO-SAC model with components
actmodel = COSMOSAC(COSMO, mixture=mixture)

# Calculate solubility (SLE)
sle = SLE(actmodel=actmodel)
solubility = sle.solubility(mix='real')

# Output the solubility
print(solubility[['T', 'w', 'x']].to_string(index=False))

# Plot results
plt.plot(*solubility[['w','T']].values.T,'.-', label='Solubility (w)')

# Settings
plt.xlim(0,1)
plt.ylim(300,500)
# Adding title and labels
plt.title('Solubility vs. Temperature')
plt.ylabel("T / K")
xlabel = {'w':'Weight', 'x':'Mole'}
plt.xlabel(f"Weight fraction {mixture[0].name}")
plt.legend()
# Save the figure to a PNG or PDF file
plt.savefig('solubility_plot.png')  # Saves the plot as a PNG file
# plt.savefig('solubility_plot.pdf')  # Saves the plot as a PDF file
plt.show()
```

For a more comprehensive demonstration, including advanced functionalities and plotting results, please see the [example_usage.py](https://github.com/ivanantolo/cosmopharm/blob/main/example_usage.py) script in this repository. This detailed example walks through the process of setting up COSMOPharm, initializing models, and visualizing the results of solubility and miscibility calculations.

## Contributing

Contributions are welcome! Please refer to our [GitHub repository](https://github.com/ivanantolo/cosmopharm) for more information.

## Citation

We appreciate citations to our work as they help acknowledge and spread our research contributions. If you use COSMOPharm in your research, please cite the associated paper as follows:

```bibtex
@article{Antolovic2024COSMOPharm,
  title={COSMOPharm: Drug--Polymer Compatibility of Pharmaceutical Amorphous Solid Dispersions from COSMO-SAC},
  author={Antolovic, Ivan and Vrabec, Jadran and Klajmon, Martin},
  journal={Molecular Pharmaceutics},
  year={2024},
  volume={1}, # Will be adjusted accordingly
  issue={1}, # Will be adjusted accordingly
  month={3}, # Will be adjusted accordingly
  pages={1--10},  # Will be adjusted accordingly
  doi={10.1021/acs.molpharmaceut.3c12345} # Will be adjusted accordingly
}
```

## License

COSMOPharm is released under the MIT License. For more details, see the [LICENSE](https://github.com/ivanantolo/cosmopharm/LICENSE) file.
