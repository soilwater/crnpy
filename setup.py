# contents of setup.py
import setuptools

setuptools.setup(
    name="crnpy",
    version="0.3.2",
    packages=['crnpy'],
    package_dir = {"": "src"},
    description="A Python package for the estimation and processing of soil moisture data from cosmic-ray neutron counts.",
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/soilwater/crnpy",
    author="Joaquin Peraza, Andres Patrignani",
    author_email="jperaza@ksu.edu, andrespatrignani@ksu.edu",
    license="MIT"
)