import glob
import importlib
import logging
import subprocess

from setuptools import find_packages, setup, Extension


# N.B. log messages typically hidden; use `pip install -v` to see
_logger = logging.getLogger()


try:
    import pysam
except ModuleNotFoundError:
    _logger.critical("!! pysam must be installed prior to running setup.py (pip install pysam)")
    raise

try:
    from Cython.Build import cythonize, build_ext
    has_cython = True
    _logger.info("# cython available; will cythonize .pyx sources")
except ImportError:
    has_cython = False
    _logger.info("# cython unavailable; using included .c sources")



# ###########################################################################
# extnsion setup
# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#configuring-the-c-build
# http://docs.cython.org/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules

glob_ext = "pyx" if has_cython else "c"
src_glob = "uta_align/align/*." + glob_ext

if has_cython:
    sources = [src_glob]
else:
    sources = glob.glob(src_glob)

extensions = [
    Extension("*", sources, include_dirs = pysam.get_include())
    ]

if has_cython:
    compiler_directives = {
        'c_string_encoding': 'ascii',
        'embedsignature': True,
        'language_level': 3
        }
    extensions = cythonize(extensions, compiler_directives=compiler_directives)


setup(
    ext_modules=extensions,
    use_scm_version=True,
    )


## <LICENSE>
## Copyright 2014-2020 uta-align Contributors (https://github.com/biocommons/uta-align)
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
