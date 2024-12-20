from setuptools import find_packages, setup
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = [
    "attrs",
    "requests",
    "fqfa>=1.2.1",
    "mavehgvs",
    "numpy",
    "pandas",
    "idutils",
    "metapub"
]

# fqfa requires backported dataclasses in Python 3.6
if sys.version_info.major == 3 and sys.version_info.minor == 6:
    requirements.append("dataclasses")

setup(
    name="mavetools",
    version="0.2.0",
    author="MaveDB Developers",
    author_email="alan.rubin@wehi.edu.au",
    description=(
        "Useful functions for manipulating Multiplex Assay of Variant Effect datasets."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VariantEffect/mavetools",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=requirements,
    test_suite="tests",
)
