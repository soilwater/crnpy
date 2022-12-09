# contents of setup.py
import setuptools

setuptools.setup(
    name="crnpy",
    version="0.2-18",
    packages=setuptools.find_packages(),
    package_data={'': ['*.csv']},
    include_package_data=True,
)