#!/usr/bin/env python
import setuptools
import glob
import os

install_requires = [
    'numpy>=1.22.4',
    'scipy>=1.9.0',
    'pandas>=1.4.3',
    'matplotlib>=3.6.2',
    'opencv-python>=4.6.0.66',
    'scikit-image>=0.19.2',
    'alphashape>=1.3.1',
    'anndata>=0.8.0',
]
extras_require = {}

setuptools.setup(
    name="spacipy",
    version="0.0.1",
    author="Kai Han",
    author_email="hankai@genomics.cn",
    url="https://github.com/lskfs/spacipy",

    packages=setuptools.find_packages(where="spacipy"),
    package_dir={"": "spacipy"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("spacipy/*.py")
    ],
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    license="GPL-3.0+",
    description="Collection of in-house python modules for some basic STOmic data processing and visualization.",
    install_requires=install_requires,
    extras_require=extras_require,
)
