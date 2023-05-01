#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="PengRan",
    author_email='2502388440@qq.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A standardized dataset-based approach designed to reduce variation in large-scale Data-Independent Acquisition (DIA) mass spectrometry data. By utilizing a reference dataset to standardize mass spectrometry signals, STAVER effectively reduces noise and enhances the accuracy of protein quantification, particularly in the context of multi-library searches.",
    entry_points={
        'console_scripts': [
            'staver=staver.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='staver',
    name='staver',
    packages=find_packages(include=['staver', 'staver.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Ran485/staver',
    version='0.1.0',
    zip_safe=False,
)
