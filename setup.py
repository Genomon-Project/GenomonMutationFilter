#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='genomon_mutation_filter',
    version='0.2.9',
    description='Python tools to verify the validity of somatic mutations.',
    author='Ken-ichi Chiba',
    author_email='kchiba@hgc.jp',
    url='https://github.com/Genomon-Project/GenomonMutationFilter',
    # package_dir = {'': 'lib'},
    packages = find_packages(exclude = ['tests']),
    test_suite = 'unit_tests.suite',
    license='GPL-3'
)
