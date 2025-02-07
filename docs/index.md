# Cosmic-Ray Neutron Python (CRNPy) Library

![CRNPy logo](img/logo/crnpy-logo.png)

[![GitHub Workflow Status (building)](https://img.shields.io/github/actions/workflow/status/soilwater/crnpy/python-package.yml)](https://github.com/soilwater/crnpy/actions/workflows/python-package.yml)
[![GitHub Workflow Status (publish)](https://img.shields.io/github/actions/workflow/status/soilwater/crnpy/python-publish.yml?label=publish)](https://github.com/soilwater/crnpy/actions/workflows/python-publish.yml)
[![PyPI - Status](https://img.shields.io/pypi/v/crnpy)](https://pypi.org/project/crnpy/)
[![GitHub commits since latest release (by SemVer including pre-releases)](https://img.shields.io/github/commits-since/soilwater/crnpy/latest/main)](https://github.com/soilwater/crnpy)
[![JOSS submission status](https://joss.theoj.org/papers/e65c1bb5fee58c39289efc4547d1fd10/status.svg)](https://joss.theoj.org/papers/e65c1bb5fee58c39289efc4547d1fd10)

Welcome to the homepage of the CRNPy (Cosmic-Ray Neutron Python) library, an open-source Python library designed for the processing and conversion of raw neutron counts from cosmic-ray neutron probes (CRNP) into soil moisture data.

This library has been developed with the intent of providing a comprehensive yet easy-to-use workflow for processing raw data from a variety of CRNP, encompassing multiple manufacturers and models.

## Statement of Need
CRNPs are a valuable tool for non-invasive soil moisture estimation at the hectometer scale (e.g., typical agricultural fields), filling the gap between point-level sensors and large-scale (i.e., several kilometers) remote sensors onboard orbiting satellites. However, cleaning, processing, and analyzing CRNP data involves multiple corrections and filtering steps spread across multiple peer-reviewed manuscripts. CRNPy simplifies these steps by providing a complete, user-friendly, and well-documented library with minimal dependencies that includes examples to convert raw CRNP data into soil moisture. The library is designed to be accessible to both researchers and instrument manufacturers. Unlike other similar libraries, CRNPy does not require any specific naming convention for the input data or large external data sources, or reanalysis data.

## Key Features
- Versatile and instrument agnostic: CRNPy can handle data from various CRNP manufacturers and models. It has been successfully tested on both roving and stationary CRNP.

- Modular: The library is designed to be modular, allowing users to easily customize the processing workflow to their needs.



## Installation
To install the CRNPy library, you can use Python's package manager. Open a terminal and type:

```pip install crnpy```

from the Jupyter notebook, type:

```!pip install crnpy```

Ideally dependencies should be installed automatically. If not, you can install them manually by typing:

```pip install -r requirements.txt```

The CRNPy library is compatible with Python 3.7 and above.
See [requirements.txt](https://github.com/soilwater/crnpy/blob/main/requirements.txt) for a list of dependencies.

## Authors

The CRNPy library was developed at the Kansas State University Soil Water Processes Lab by:

- Joaquin Peraza

- Andres Patrignani

The Soil Water Processes Lab at Kansas State University combines a range of experimental and computational approaches to tackle pressing issues in soil and water research. The development of the CRNPy library is a step forward to creating reproducible data processing workflows across the scientific community using cosmic-ray neutrons probes for soil moisture sensing. 


## Community Guidelines

### Contributing
To contribute to the software, please first fork the repository and create your own branch from `main`. Ensure your code adheres to our established code structure and includes appropriate test/examples coverage. CRNPy source code is located in the `/src/crnpy/` folder, and tests implemented using `pytest` are stored in the `/src/tests/` folder. Submit a pull request with a clear and detailed description of your changes to include them in the main repository.

### Reporting Issues
If you encounter any issues or problems with the software, please report them on our [issues page](https://github.com/soilwater/crnpy/issues). Include a detailed description of the issue, steps to reproduce the problem, any error messages you received, and details about your operating system and software version.

### Seeking Support
If you need support, please first refer to the documentation. If you still require assistance, post a question on the [issues page](https://github.com/soilwater/crnpy/issues) with the `question` tag. For private inquiries, you can reach us via email at jperaza@ksu.edu or andrespatrignani@ksu.edu.

