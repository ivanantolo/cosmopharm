# DirectImport for COSMO-SAC

## Introduction
The `DirectImport` feature is a custom enhancement to the [COSMO-SAC package](https://github.com/usnistgov/COSMOSAC), designed to provide more flexibility and organization in handling sigma profiles for different components. Unlike the standard import methods (`VirginiaTechProfileDatabase`, `DelawareProfileDatabase`) in COSMO-SAC, which require all sigma profiles to be stored in one folder and listed in a TXT file, `DirectImport` allows users to store sigma profiles in separate directories. This feature enables specifying the path and name of the sigma file for each component individually, offering a significant improvement in directory structure and accessibility.

## Features
- **Enhanced Organization:** With `DirectImport`, users can separate sigma profiles into different folders, such as one for polymers and another for APIs (Active Pharmaceutical Ingredients), and provide their paths separately. This organization is particularly beneficial for maintaining a clear and structured directory.
- **Flexible Testing:** `DirectImport` facilitates the testing of different sigma profiles for the same component by allowing users to select which sigma profile to use explicitly. This flexibility is invaluable for research and development purposes, where multiple iterations of sigma profiles may need to be evaluated.

## Installation
To incorporate the `DirectImport` feature into your COSMO-SAC package, follow these steps:
1. Clone the repository to your local machine: 'git clone --recursive --shallow-submodules https://github.com/usnistgov/COSMOSAC'
2. Navigate to the Repository Directory: 'cd COSMOSAC'
3. Replace the original `include` and `src` folders in the [COSMO-SAC package](https://github.com/usnistgov/COSMOSAC) with the ones provided in this `directimport` directory.
4. Proceed with the standard installation process of the COSMO-SAC package. The addition of `DirectImport` does not alter the remaining installation steps.

## Usage
To utilize the `DirectImport` feature, specify the path and name of the sigma file for each component in your simulations. This approach allows for a more organized and efficient management of sigma profiles.

### Example
Here, an example code snippet will be provided to illustrate the difference between the old and new methods of importing sigma profiles. (Note: insert example code here comparing the standard and `DirectImport` methods.)

## Contribution
This enhancement is currently not a part of the official COSMO-SAC package. It has been developed as an add-on to offer a more convenient and organized way to handle sigma profiles. A pull request (PR) will be made to the COSMO-SAC repository for potential inclusion. Until then, and until the publication associated with this work is released, the source code and instructions will be available here for anyone interested in adopting this feature.

### Note
Please note that `DirectImport` is a custom feature created specifically for the COSMO-SAC package by [Your Name/Your Institution]. It is not yet an official part of the COSMO-SAC project. Feedback and contributions are welcome to improve and integrate this feature into the main COSMO-SAC package.
