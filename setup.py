from setuptools import setup, find_packages
from os import path
here = path.abspath(path.dirname(__file__))

def get_version():
    with open(path.join(here, "genomon_mutation_filter/version.py")) as hin:
        for line in hin:
            if line.startswith("__version__"):
                version = line.partition('=')[2]
                return version.strip().strip('\'"')
    raise ValueError('Could not find version.')

setup(
      name='fisher',
      version=get_version(),
      description="Python programs for filtering mutation results.",
      long_description="""""",

      classifiers=[
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 5 - Stable',
          # Indicate who your project is intended for
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      
      keywords='Bio-informatics',
      author='Ken-ichi Chiba',
      author_email='kchiba@hgc.jp',
      url='https://github.com/Genomon-Project/GenomonMutationFilter.git',
      license='GPL-3',
      
      packages = find_packages(exclude = ['tests']),
      install_requires=[
      ],
      entry_points = {'console_scripts': ['mutfilter = genomon_mutation_filter:main']},
      test_suite = 'unit_tests.suite'
)


