"""
A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from alexov import __version__, _logo

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-alexov',  # Required
    version=__version__,  # Required
    description='Scipion plugin for Alexov group in http://compbio.clemson.edu/lab/',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-chem/scipion-chem-alexov',  # Optional
    author='Natalia del Rey',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion chemoinformatics',  # Optional
    classifiers=[  # Optional
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'
    ],
    packages=find_packages(),
    install_requires=[requirements],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/scipion-chem/scipion-chem-alexov/issues',
        'Source': 'https://github.com/scipion-chem/scipion-chem-alexov/',
        'Pull Requests': 'https://github.com/scipion-chem/scipion-chem-alexov/pulls'
    },
    entry_points={'pyworkflow.plugin': 'alexov = alexov'},
    package_data={  # Optional
       'alexov': [_logo, 'protocols.conf'],
    }
)
