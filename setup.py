# contents of setup.py
import setuptools

setuptools.setup(
    name="crnpy",
    version="0.2-29",
    packages=['crnpy'],
    package_dir = {"": "src"},
    include_package_data=True,
)