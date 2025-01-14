#!/usr/bin/env python
"""Shim to allow Github to detect the package, build is done with hatch."""

# !/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="popv",
    version="0.1.0",
    packages=find_packages(include=['popv', 'popv.*']),  # 只包含popv包
    python_requires='>=3.9',
    install_requires=[
        'scvi-tools>=1.0.0',
        'scimilarity',
        'scanpy',
        'numpy',
        'torch',
        'scikit-learn',
        'pandas',
        # Other dependencies
    ],
)
