# contents of setup.py
import setuptools

setuptools.setup(
    name="crnpy",
    version="0.2-28",
    packages=['crnpy'],
    package_dir = {"": "src"},
    package_data={'crnpy': ['*.csv']},
    include_package_data=True,
)