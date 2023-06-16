# Cosmic-Ray Neutron Python (CRNPy) Library

<img src="/docs/img/logo/crnpy-logo.png" alt="CRNPY logo" width="250"/>

## Overview

Welcome to the homepage of the CRNPy (Cosmic-Ray Neutron Python) library, an open-source Python library designed for the processing and conversion of raw neutron counts from Cosmic-Ray Neutron Probes (CRNP) into accurate field-level soil moisture data.

This library has been developed with the intent of providing a comprehensive yet easy-to-use workflow for processing raw data from a variety of CRNP, encompassing multiple manufacturers and models.

## Key Features
- Versatility: CRNPy can handle data from various CRNP manufacturers and models. It has been successfully tested on both roving and stationary CRNP.
- Modularity: The library is designed to be modular, allowing users to easily customize the processing workflow to their needs.
- Accuracy: The library shows excellent agreement with ground-based measurements, making it a reliable tool for scientific research and practical applications.
- Correction Routines: The correction functions of the library have been effectively used to represent field and watershed scale soil moisture conditions.

![CRNPy Processing Workflow](docs/img/workflow.png)
Overview of the proposed CRNPy processing workflow. Final user can choose to use the entire workflow, part of it, or build functions on top of it depending on their needs, dashed lines indicate optional steps.


## Installation

To install the CRNPy library, you can use Python's package manager, pip. Simply open your command line or terminal and type:

```pip install crnpy```

## Authors
The CRNPy library was conceived and developed by Joaquin Peraza and Andres Patrignani in the Soil Water Processes Lab at Kansas State University. The team's passion for making soil moisture data more accessible and easier to process culminated in this powerful tool.

The Soil Water Processes Lab at Kansas State University is a leading research group focused on understanding and modeling soil water processes. The lab combines a range of experimental and computational approaches to tackle some of the most pressing issues in soil and water research. The development of the CRNPy library is a testament to the lab's commitment to pushing the boundaries of soil science through the innovative use of technology. like to acknowledge the contributions of the wider scientific community in testing and providing feedback on the library, which has been instrumental in its ongoing development.
