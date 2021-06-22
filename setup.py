# -*- coding: utf-8 -*-
import sys

from setuptools import find_packages
from setuptools import setup

if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must be using Python 3")

    setup(
        name="ablc",
        python_requires=">=3.9.0",
        version=open("version.txt", "r").read(),
        description="Analysis of Light Chain-Antigen Binding Interface",
        url="https://github.com/Dariusrussellkish/ablc-binding-interface",
        author="Darius Russell Kish",
        packages=find_packages(),
        install_requires=["biopython", "numpy", "pandas"],
        setup_requires=["pytest-runner"],
        tests_require=["pytest"],
        include_package_data=True,
        zip_safe=False,
        scripts=[],
    )