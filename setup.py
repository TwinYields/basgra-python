#!/usr/bin/env python3
import setuptools
import site
from numpy.distutils.core import setup, Extension

site.ENABLE_USER_SITE = True

setup(
    name = "basgra",
    version = "0.0.1",
    packages = setuptools.find_packages(),
    author = "Matti Pastell",
    ext_modules=[
        Extension(name="bglib", sources=["src/parameters_plant.f90", "src/parameters_site.f90",
                                        "src/resources.f90", "src/environment.f90",
                                        "src/soil.f90", "src/plant.f90", "src/set_params.f90",
                                        "src/bglib.f90"]),
     ],
    include_package_data=True
)