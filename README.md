# COSMOPharm

<p align="center">
  <img src=https://github.com/ivanantolo/cosmopharm/raw/main/TOC.png alt="TOC Figure" width="500"/>
</p>


Welcome to the COSMOPharm repository, accompanying [our paper in *Molecular Pharmaceutics*](https://doi.org/10.1021/acs.molpharmaceut.4c00342). This project and its associated publication offer insights and a practical toolkit for researching drug-polymer and drug-solvent systems, aiming to provide the scientific community with the means to reproduce our findings and further the development of COSMO-SAC-based models.

## About 

COSMOPharm is a Python package designed to streamline the predictive modeling of drug-polymer compatibility, crucial for the development of pharmaceutical amorphous solid dispersions. Apart from that, it can also be used for the miscibility/solubility of drugs with/in common solvents. Leveraging the COSMO-SAC (Conductor-like Screening Model Segment Activity Coefficient) model, COSMOPharm offers a robust platform for scientists and researchers to predict solubility, miscibility, and phase behavior in drug formulation processes.

## Features

- **Compatibility Prediction**: Utilize open-source COSMO-SAC model for prediction of drug-polymer compatibility.
- **Solubility Calculation**: Calculate drug-polymer solubilities to guide the selection of suitable polymers for drug formulations.
- **Miscibility and Phase Behavior Analysis**: Analyze the miscibility of drug-polymer pairs and understand their phase behavior under various conditions.
- **User-friendly Interface**: Easy-to-use functions and comprehensive documentation to facilitate research and development in pharmaceutical sciences.

## Installation

### Quick Installation
For most users, the quickest and easiest way to install COSMOPharm is via pip, which will manage all dependencies for you. Ensure you have already installed the `cCOSMO` library by following the instructions on the [COSMOSAC GitHub page](https://github.com/usnistgov/COSMOSAC).

Once `cCOSMO` is installed, you can install COSMOPharm directly from [PyPI](https://pypi.org/project/cosmopharm/):

```
pip install cosmopharm
```

### Advanced Installation Options

For users who need more control over the installation process (e.g., for development purposes or when integrating with other projects), COSMOPharm can also be installed by cloning the repository and installing manually. 

#### Step 1: Clone the Repository

First, clone the COSMOPharm repository. The project includes a submodule named "pharmaceuticals" that stores essential data files (.sigma, .cosmo, and .xyz files) for the pharmaceutical components. Use the --recurse-submodules option to ensure that the "pharmaceuticals" submodule is correctly initialized and updated along with the main project.
```
git clone --recurse-submodules https://github.com/ivanantolo/cosmopharm
```

#### Step 2: Navigate to the Repository Directory

```
cd cosmopharm
```

#### Option 1: Using pip to Install from Local Source

This method installs COSMOPharm and manages all dependencies efficiently:

```
pip install .
```


#### Option 2: Using setup.py for Installation

Alternatively, you can run the setup script directly:

```
python setup.py install
```

While this method is straightforward, using `pip` is generally preferred for its dependency management capabilities.

Please note: Before proceeding with either advanced installation option, ensure the `cCOSMO` library is installed as described at the beginning of this section.

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

## Contributing / Getting Help

Contributions to COSMOPharm are welcome! We accept contributions via pull requests to the [GitHub repository](https://github.com/ivanantolo/cosmopharm). 

For bugs, feature requests, or other queries, please [open an issue](https://github.com/ivanantolo/cosmopharm/issues) on GitHub.


## Citation

We appreciate citations to our work as they help acknowledge and spread our research contributions. If you use COSMOPharm in your research, please cite the associated paper. Citation details are provided in the [`CITATION.cff`](https://github.com/ivanantolo/cosmopharm/CITATION.cff) file, and GitHub generates APA or BibTeX formats accessible under the "Cite this repository" dropdown on our repository page.

For convenience, here's the citation in BibTeX format:

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

### Gaussian Citation
**Permission from Gaussian, Inc.** was obtained to make the [.cosmo files](https://github.com/ivanantolo/cosmopharm/tree/main/profiles/polymers/cosmo) for the oligomers included in [our paper in *Mol. Pharm.*](https://doi.org/10.1021/acs.molpharmaceut.4c00342) available for academic and research (noncommercial) use. This is to enhance research transparency and facilitate the validation process.

When using COSMOPharm with the sigma-profiles provided therein, please also cite Gaussian as follows:

```bibtex
@misc{g16,
author={M. J. Frisch and G. W. Trucks and H. B. Schlegel and G. E. Scuseria and M. A. Robb and J. R. Cheeseman and G. Scalmani and V. Barone and G. A. Petersson and H. Nakatsuji and X. Li and M. Caricato and A. V. Marenich and J. Bloino and B. G. Janesko and R. Gomperts and B. Mennucci and H. P. Hratchian and J. V. Ortiz and A. F. Izmaylov and J. L. Sonnenberg and D. Williams-Young and F. Ding and F. Lipparini and F. Egidi and J. Goings and B. Peng and A. Petrone and T. Henderson and D. Ranasinghe and V. G. Zakrzewski and J. Gao and N. Rega and G. Zheng and W. Liang and M. Hada and M. Ehara and K. Toyota and R. Fukuda and J. Hasegawa and M. Ishida and T. Nakajima and Y. Honda and O. Kitao and H. Nakai and T. Vreven and K. Throssell and Montgomery, {Jr.}, J. A. and J. E. Peralta and F. Ogliaro and M. J. Bearpark and J. J. Heyd and E. N. Brothers and K. N. Kudin and V. N. Staroverov and T. A. Keith and R. Kobayashi and J. Normand and K. Raghavachari and A. P. Rendell and J. C. Burant and S. S. Iyengar and J. Tomasi and M. Cossi and J. M. Millam and M. Klene and C. Adamo and R. Cammi and J. W. Ochterski and R. L. Martin and K. Morokuma and O. Farkas and J. B. Foresman and D. J. Fox},
title={Gaussian˜16 {R}evision {C}.01},
year={2016},
note={Gaussian Inc. Wallingford CT}
} 
```

## License

COSMOPharm is released under the MIT License. See the [LICENSE](https://github.com/ivanantolo/cosmopharm/LICENSE) file for more details.

