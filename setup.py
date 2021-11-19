# based on https://realpython.com/pypi-publish-python-package
import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

# This call to setup() does all the work
setup(
    name='ncbi_blast',
    version='0.1.0',
    description='Wrapper for NCBI Blast',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/MrTomRod/ncbi_blast',
    author='Thomas Roder',
    author_email='roder.thomas@gmail.com',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    packages=['ncbi_blast'],
    install_requires=['biopython'],
)
