.. Model Balancing documentation master file, created by
   sphinx-quickstart on Mon Jun  7 19:46:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Model balancing
===============

`Model balancing <https://www.metabolic-economics.de/model-balancing/index.html>`_
is a computational method to determine plausible kinetic constants and metabolic
states in kinetic metabolic models. It integrates flux, metabolite, protein,
and kinetic constant data, using prior distributions for all these variables, and
computes the joint posterior mode.

Model balancing can be run in *Matlab* or *Python*. Data tables can be provided in
`SBtab <https://www.sbtab.net>`_ format, models can be provided in 
`SBML <http://sbml.org>`_ or `SBtab <https://www.sbtab.net>`_ format.

This documentation is for the *Python* version of Model Balancing only.
For balancing your model, generating a JSON input file, or running model balancing in Maltab -
see instructions `here <https://github.com/liebermeister/model-balancing>`_.

Installation
============
 
Clone the repository::

    git clone https://github.com/liebermeister/model-balancing.git

Install using the package in a new Virtual Environment using::

    cd python
    python -m venv venv
    source venv/bin/activate
    pip install -e .

Obtain a license for `MOSEK <https://www.mosek.com/>`_, for example you might qualify for a
`free academic license <https://www.mosek.com/products/academic-licenses/>`_.

Now you can then try the example script::

    python examples/comparison_with_matlab.py

which runs model balancing on a list of JSON examples and for a fixed set of values for alpha.

The code was tested with Python 3.9 on Ubuntu Linux 21.04.

License
=======
This package is released under the GNU General Public License.

Contact
=======
Please contact `Wolfram Liebermeister <mailto:wolfram.liebermeister@gmail.com>`_
and `Elad Noor <mailto:elad.noor@weizmann.ac.il>`_ with any questions or comments.

References
==========
Liebermeister W. and Noor E. (2021),
*Model balancing: in search of consistent metabolic states and in-vivo kinetic constants*
`bioRxiv doi:10.1101/2019.12.23.887166v2 <https://www.biorxiv.org/content/10.1101/2019.12.23.887166v2>`_

`www.metabolic-economics.de/model-balancing/ <https://www.metabolic-economics.de/model-balancing/index.html>`_


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   autoapi/index
   autoapi/python/src/model_balancing/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
