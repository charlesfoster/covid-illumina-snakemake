from setuptools import setup, find_packages
import glob
import os
import pkg_resources
import pathlib

from CIS import __version__, _program

# glob_fix command thanks to https://stackoverflow.com/a/64789489/15793575
def glob_fix(package_name, glob):
    # this assumes setup.py lives in the folder that contains the package
    package_path = pathlib.Path(f'./{package_name}').resolve()
    return [str(path.relative_to(package_path))
            for path in package_path.glob(glob)]

setup(name='CIS',
      version=__version__,
      packages=find_packages(),
      scripts=['CIS/bin/Snakefile'
                ],
#      package_data={"CIS":["bin/*"]},
      package_data={"CIS":[*glob_fix("CIS", "bin/**/*")]},
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
