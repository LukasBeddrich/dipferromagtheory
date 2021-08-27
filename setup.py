from setuptools import setup, find_packages

setup(
    name="dipferromagtheory",
    version="0.1.0",
    packages=find_packages(),
    author="Lukas Beddrich",
    author_email="lukas.beddrich@frm2.tum.de",
    description="package to calculate a ferromagnets magnetic response (observed in neutron scattering) including dipolar interactions.",
    package_data={"dipferromagtheory" : ["res/*.json"]},
    include_package_data=True
)