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

::

  $ make devready

Test it::

  $ pytest

Because cython (.pyx) files are compiled, you'll need to reinstall the
development package after code changes, like this::
  
  $ pip install -e .




.. |pypi_badge| image:: https://badge.fury.io/py/uta-align.png
  :target: https://pypi.python.org/pypi?name=uta-align
  :align: middle
