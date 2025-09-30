from setuptools import setup, find_packages

setup(
    packages=find_packages(where="src"),
    name='gmx-top4py',
    version='0.0.0',
    include_package_data=True,
)