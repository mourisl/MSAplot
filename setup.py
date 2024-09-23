#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
import sys

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))

SRC_DIR = 'msaplot'

import msaplot
version = msaplot.__version__
AUTHOR = 'Li Song'
EMAIL = 'mourisl@hotmail.coom'
URL = 'https://github.com/mourisl/msaplot'
DESCRIPTION = 'https://github.com/mourisl/msaplot'
PACKAGES = [SRC_DIR,
            ]
setup(name='msaplot',
        version=version,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        packages=PACKAGES,
        license='MIT',
        install_requires=["seaborn", "matplotlib", "pandas"],
        package_data={
        },
        entry_points={
        },
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        keywords= ['Data Visualization'], 
        classifiers=[ "Programming Language :: Python :: 3"],
)
