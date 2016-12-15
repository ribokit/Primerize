Installation
------------------

**Primerize** can be downloaded for non-commercial use at the `Primerize Server <https://primerize.stanford.edu/license/>`_. Commercial users, please `contact us <https://primerize.stanford.edu/about/#contact>`_.

* Download the zip or tar file of the repository and unpack; or 

.. code-block:: bash

    git clone https://github.com/DasLab/Primerize.git

* To install ``Primerize``, simply:

.. code-block:: bash

    cd path/to/Primerize/
    python setup.py install

For system-wide installation, you must have permissions and use with ``sudo``.

``Primerize`` requires the following *Python* packages as dependencies, all of which can be installed through `pip <https://pip.pypa.io/>`_:

.. code-block:: js

    matplotlib >= 1.5.0
    numpy >= 1.10.1
    xlwt >= 1.0.0

----------

Loop Optimization with ``numba`` *(Optional)*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To speed up **Primerize** code, we take advantage of ``@jit`` (`see here <http://numba.pydata.org/numba-doc/0.23.1/user/jit.html>`_) decorator of `numba <http://numba.pydata.org/>`_ on loop optimization. **This is totally optional.** Enabling such feature may speed up the run for up to *10x*.

``numba`` requires `llvm <http://llvm.org/>`_, which can be installed through `apt-get <https://help.ubuntu.com/lts/serverguide/apt-get.html>`_ on *Linux* or `brew <http://brew.sh/>`_ on Mac *OSX*. It also requires ``llvmlite``, which can be installed through ``pip``. The compatibility between ``numba``, ``llvmlite``, and ``llvm`` needs to pay special attention to. The below specified ``numba`` and ``llvmlite`` versions have been tested to work with ``llvm 3.6.2`` on *Linux* machines.

.. code-block:: js

    llvmlite == 0.8.0
    numba == 0.23.1

Or the following works with ``llvm 3.7.1``.

.. code-block:: js

    llvmlite == 0.12.1
    numba == 0.27.0

Newer version combinations may work, but we haven't test since.

----------

Test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test if **Primerize** is functioning properly in local installation, run the *unit test* scripts:

.. code-block:: bash

    cd path/to/Primerize/tests/
    python -m unittest discover

All test cases should pass.

