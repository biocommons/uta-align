import importlib
import subprocess

from setuptools import setup, find_packages

# build_requires is not a setuptools keyword. These packages are
# required just to process this setup.py, and therefore are installed
# by setup.py prior to invoking setup().
build_requires = [
    'Cython>=0.20',             # for extension hooks below
    'pysam>=0.8',               # for include file inference below
]

install_requires = build_requires + [
    'pysam>=0.10.0',
]

setup_requires = [
    'setuptools_scm',
    'wheel',
]

tests_require = [
    'coverage',
    'nose',
]


# ###########################################################################
# install prerequisites that aren't well-handled by setuptools

def install_build_requires(pkg_targets):
    """Iterate through build_requires list and pip install if package is not present
    accounting for version"""

    def pip_install(pkg_name, pkg_vers=None):
        pkg_name_version = '%s==%s' % (pkg_name, pkg_vers) if pkg_vers else pkg_name
        print '[WARNING] %s not found, attempting to install using a raw "pip install" call!' % pkg_name_version
        subprocess.Popen('pip install %s' % pkg_name_version, shell=True).communicate()

    def get_pkg_info(pkg):
        """Get package name and version given a build_requires element"""
        pkg_name, pkg_vers = None, None
        if '==' in pkg:
            pkg_name, pkg_vers = pkg.split('==')
        else:
            pkg_name = pkg.replace('>', '').replace('<', '').split('=')[0]
        return pkg_name, pkg_vers

    for pkg in pkg_targets:
        pkg_name, pkg_vers = get_pkg_info(pkg)
        try:
            pkg_name_version = '%s==%s' % (pkg_name, pkg_vers) if pkg_vers else pkg_name
            if pkg_vers:
                version = getattr(importlib.import_module(pkg_name), '__version__')
                if version != pkg_vers:
                    pip_install(pkg_name, pkg_vers)
            else:
                importlib.import_module(pkg_name)
        except ImportError:
            pip_install(pkg_name, pkg_vers)

install_build_requires(build_requires)


# ###########################################################################
# cython dependencies/setup
#
from cy_distribute import CyExtension, cy_build_ext

def pysam_incl(cy_ext):
    import pysam
    cy_ext.extend_includes(pysam.get_include())
    cy_ext.extend_macros(pysam.get_defines())
    cy_ext.extend_extra_objects([pysam.libchtslib.__file__])


ext_modules = [
    CyExtension('uta_align.align.cigar_utils',      ['uta_align/align/cigar_utils.pyx'],      init_func=pysam_incl),
    CyExtension('uta_align.align.algorithms',       ['uta_align/align/algorithms.pyx'],       init_func=pysam_incl),
]



setup(
    author='Kevin Jacobs',
    author_email='reecehart+biocommons@gmail.com',
    cmdclass={'build_ext': cy_build_ext},
    description='C-based alignment for Python',
    ext_modules=ext_modules,
    install_requires=install_requires,
    license='Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    maintainer='Reece Hart',
    maintainer_email='reecehart+biocommons@gmail.com',
    name='uta-align',
    packages=find_packages(),
    setup_requires=setup_requires,
    test_suite='nose.collector',
    tests_require=tests_require,
    url='https://github.com/biocommons/uta-align/',
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
