``uta-align`` provides C-based Needleman-Wunsch and Smith-Waterman
alignment algorithms with a Python interface.

|build_status| |pypi_badge|


.. important:: Python 3.6+ only



Installation
@@@@@@@@@@@@

From source::

  $ python setup.py install



For Development
###############

* Install developer environment::

  $ make devready

* Activate the environment::

  $ source venv/3.7/bin/activate

(or other python version as appropriate for your system)


* During development, you'll need to recompile any code changes to the
  .pyx files, like this::

  $ pip install -e .


* Test it (in current Python environment)::

  $ make test
  ======================================== test session starts ========================================
  platform linux -- Python 3.7.5, pytest-5.3.5, py-1.8.1, pluggy-0.13.1
  rootdir: /home/reece/projects/biocommons/uta-align, inifile: setup.cfg, testpaths: tests
  collected 5 items                                                                                   
   
  tests/test_align_algorithms.py .                                                              [ 20%]
  tests/test_align_gap_open.py ..                                                               [ 60%]
  tests/test_cigar_utils.py .                                                                   [ 80%]
  tests/test_nw.py .                                                                            [100%]
   
  ========================================= 5 passed in 0.12s =========================================


* Tox test it (in multiple environments::
  
  $ make tox
  
  GLOB sdist-make: /home/reece/projects/biocommons/uta-align/setup.py
  py37 inst-nodeps: /home/reece/projects/biocommons/uta-align/.tox/.tmp/package/1/uta-align-1.1.13a2.dev7+g98b7cc9.d20200303.zip
  ...
  py37: commands succeeded
  congratulations :)


.. |pypi_badge| image:: https://badge.fury.io/py/uta-align.png
  :target: https://pypi.python.org/pypi?name=uta-align
  :align: middle

