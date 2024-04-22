#!/usr/bin/env python

import codecs
from distutils.core import setup
import os.path
from pathlib import Path

# read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# Read the requirements from requirements.txt
with open("requirements.txt", "rt") as handle:
    requirements = [
        line.rstrip("\n")
        for line in handle
    ]

setup(
    name="mudata_explorer",
    version=get_version("mudata_explorer/__init__.py"),
    description="A tool for exploring MuData objects",
    long_description_content_type='text/markdown',
    long_description=long_description,
    author="Samuel Minot",
    author_email="sminot@fredhutch.org",
    url="https://github.com/CirroBioApps/mudata-explorer",
    packages=["mudata_explorer"],
    license="MIT",
    entry_points={},
    install_requires=requirements,
    include_package_data=True
)
