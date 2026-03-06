
from setuptools import setup, find_packages

with open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='scipion-em-biometall',
    version='1.0.0',
    description='Scipion plugin for BioMetAll metal-binding site prediction',
    long_description=long_description,
    url='',
    author='',
    author_email='',
    keywords='scipion biometall metal-binding structural-biology scipion-3.0',
    packages=find_packages(),
    install_requires=['scipion-em', 'biometall', 'freesasa', 'pandas'],
    package_data={
        'scipionbiometall': ['protocols.conf'],
    },
    entry_points={
        'pyworkflow.plugin': ['scipionbiometall=scipionbiometall']
    },
)
