.. highlight:: shell

============
Installation
============


Stable release
--------------

To install STAVER, run this command in your terminal:

.. code-block:: console

    $ pip install staver

This is the preferred method to install STAVER, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. note::
    To avoid potential dependency conflicts, installing within a
    `conda environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__
    is recommended.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for STAVER can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/Ran485/staver

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/Ran485/staver/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/Ran485/staver
.. _tarball: https://github.com/Ran485/staver/tarball/master


Environmental dependencies
--------------------------
The `STAVER` package use some functions for the data processing, leveraging the open-source Python packages.

You may install additional environmental dependencies:

.. code-block:: console

    $ pip install -r requirements_dev.txt
    $ pip install -r requirements.txt

.. note::
    To avoid potential dependency conflicts, installing within a
    `conda environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__
    is recommended.

Now you are all set. Proceed to tutorials for how to use the `STAVER` package.