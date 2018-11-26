
#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_mutation_filter',
    version='0.2.3',
    description='Python tools to verify the validity of somatic mutations.',
    author='Ken-ichi Chiba',
    author_email='kchiba@hgc.jp',
    url='https://github.com/ken0-1n/GenomonMutationFilter',
    package_dir = {'': 'lib'},
    packages=['mutfilter'],
    scripts=['mutfilter'],
    license='GPL-3'
)
