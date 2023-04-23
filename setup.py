#!/usr/bin/env python
import setuptools
import glob
import os

setuptools.setup(
    name="spacipy",
    version="1.0.1",
    author="Kai Han",
    author_email="hankai@genomics.cn",
    url="https://github.com/lskfs/spacipy",
    #packages=setuptools.find_packages(),
    packages=setuptools.find_packages(where="spacipy"),
    package_dir={"": "spacipy"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("spacipy/*.py")
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    license="GPL-3.0+",
    description="Collection of in-house python modules for some basic STOmic data processing and visualization. It write mainly for format converting and visualization.",
    platforms='any'
)
