``uta-align`` provides C-based Needleman-Wunsch and Smith-Waterman
alignment algorithms with a Python interface.

|build_status| |pypi_badge|


.. important:: Python 3.6+ only



Installation
@@@@@@@@@@@@


Virtual Environment and Prerequisites
#####################################

::

   $ python3 -m venv venv
   $ source venv/bin/activate
   $ pip install cython ipython setuptools_scm pytest wheel
   $ pip install pysam

Installing pysam can take a few minutes. Be patient.


For Development
###############

::

   $ pip install -e .

Because cython (.pyx) files are compiled, you'll need to rerun this
command after code changes.

Test it::

  $ pytest



.. |pypi_badge| image:: https://badge.fury.io/py/uta-align.png
  :target: https://pypi.python.org/pypi?name=uta-align
  :align: middle
