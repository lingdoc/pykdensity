import os
from setuptools import setup, find_packages

# read the readme file safely if it exists on disk
readme_path = "README.md"
long_description = ""
if os.path.exists(readme_path):
    with open(readme_path, encoding="utf-8") as f:
        long_description = f.read()

setup(
    name="pykdensity",
    version="0.1.0",
    author="Hiram Ring",
    description="Network adjacency engine for measuring spatial and structural data connectivity",
    long_description=open("README.md", encoding="utf-8").read() if open("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/lingdoc/pykdensity",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.22.0",
        "pandas>=1.4.0",
        "scipy>=1.8.0",
        "biopython>=1.79",
        "openpyxl>=3.0.0",  # necessary for loading excel spreadsheets in the culture track
    ],
)
