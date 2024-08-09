=======================
Scipion Alexov plugin
=======================

This is a plugin for **scipion** wrapping some algorithms developed by `Professor Emil Alexov 
Group <http://compbio.clemson.edu/lab/>`_.

These tools will make it possible to carry out different functions for computational modeling 
of biological macromolecules and their assemblages as well as predicting biophysical quantities 
associated with them (e.g: prediction of the binding free energy change (ΔΔG) for protein-protein 
interactions).

Therefore, this plugin allows to use programs from the SAAMBE-3D software suite
within the Scipion framework.


==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.


1. **Install the plugin in Scipion**

External software `SAAMBE-3D <https://github.com/delphi001/SAAMBE-3D>`_ is installed automatically by 
Scipion.

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or by running command:

.. code-block::

    scipion3 installp -p scipion-chem-alexov


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-alexov.git

    2. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-alexov --devel

