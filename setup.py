#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme_file:
    readme = readme_file.read()


with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = []

test_requirements = []

setup(
    author="PengRan",
    author_email="21112030023@m.fudan.edu.cn",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="A Python library that presents a standardized dataset-based algorithm designed to reduce variation in large-scale data-independent acquisition (DIA) mass spectrometry data.",
    entry_points={
        "console_scripts": [
            "staver=staver.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="staver",
    name="staver",
    packages=find_packages(include=["staver", "staver.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/Ran485/STAVER",
    version="0.1.2",
    zip_safe=False,
)
