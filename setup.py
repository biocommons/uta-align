import importlib
import logging
import subprocess

from setuptools import find_packages, setup, Extension


# N.B. log messages typically hidden; use `pip install -v` to see
_logger = logging.getLogger()

try:
    from Cython.Build import cythonize, build_ext
    has_cython = True
    _logger.warning("# cython available; will cythonize .pyx sources")
except ImportError:
    has_cython = False
    _logger.warning("# cython unavailable; using included .c sources")



# ###########################################################################
# extnsion setup
# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#configuring-the-c-build
# http://docs.cython.org/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules

import pysam
import glob

glob_ext = "pyx" if has_cython else "c"
src_glob = "uta_align/align/*." + glob_ext

if has_cython:
    sources = [src_glob]
else:
    sources = glob.glob(src_glob)

extensions = [
    Extension("uta_align", sources,
              include_dirs = pysam.get_include(),
              libraries = [pysam.libchtslib.__file__])
    ]

if has_cython:
    compiler_directives = {
        'c_string_encoding': 'ascii',
        'embedsignature': True,
        'language_level': 3
        }
    extensions = cythonize(extensions, compiler_directives=compiler_directives)



install_requires = [
    "pysam>=0.15.0",
]

setup_requires = [
    "pytest-runner",
    "setuptools > 41"
    "setuptools_scm",
    "wheel",
]

tests_require = [
    "coverage",
    "pytest",
]




setup(
    author='Kevin Jacobs',
    author_email='reecehart+biocommons@gmail.com',
    description='C-based alignment for Python',
    maintainer='Reece Hart',
    maintainer_email='reecehart+biocommons@gmail.com',
    url='https://github.com/biocommons/uta-align/',
    name='uta-align',

    #cmdclass={'build_ext': build_ext},
    ext_modules=extensions,
    install_requires=install_requires,
    license='Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    packages=find_packages(),
    setup_requires=setup_requires,
    tests_require=tests_require,
    use_scm_version=True,
    zip_safe=False,
)


## <LICENSE>
## Copyright 2014 uta-align Contributors (https://bitbucket.org/biocommons/uta-align)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
