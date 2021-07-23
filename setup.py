from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from CIS import __version__, _program

setup(name='CIS',
      version=__version__,
      packages=find_packages(),
      scripts=['CIS/bin/Snakefile'
                ],
      package_data={"CIS":["bin/*"]},
      install_requires=[
            'pandas>=1.0.1',
        ],
      description='COVID Illumina Pipeline from SAVID: Snakemake edition',
      url='https://github.com/charlesfoster/covid-illumina-snakemake',
      author='Dr Charles Foster',
      author_email='charles.foster@unsw.edu.au',
      entry_points="""
      [console_scripts]
      {program} = CIS.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
