
Scipion plugin for BioMetAll
============================

.. image:: https://img.shields.io/badge/status-beta-yellow
   :target: https://github.com/
.. image:: https://img.shields.io/badge/scipion-3.8-red
   :target: https://scipion.i2pc.es/

A `Scipion <https://scipion.i2pc.es/>`_ plugin to predict metal-binding
sites in proteins using `BioMetAll <https://github.com/insilichem/biometall>`_.

BioMetAll identifies candidate metal coordination sites from backbone
preorganization geometry, without requiring a metal ion in the input
structure. It generates virtual probes at positions that are
geometrically compatible with metal coordination and clusters them to
rank candidate binding sites.


Installation
------------

You will need a working Scipion 3 installation. Then install this plugin
in development mode::

    scipion3 installp -p /path/to/scipion-em-biometall --devel

BioMetAll and its dependencies (freesasa, pandas) are installed
automatically via pip.


Protocols
---------

The plugin provides one protocol, available in Scipion under
**scipionbiometall > predict metal binding sites**.


predict metal binding sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Predicts metal-binding sites in a protein structure by placing virtual
probes on geometrically favourable positions and clustering them.

**Input**

- *Input protein structure*: an AtomStruct object (PDB/mmCIF) imported
  into the Scipion project.

**Basic options**

- *Coordinating residues*: residue types considered as potential metal
  coordinators. Comma-separated three-letter codes.
  Default: ``HIS,CYS,ASP,GLU``

- *Minimum coordinators*: minimum number of residues that must
  simultaneously coordinate a probe for the site to be considered valid.
  Default: ``3``

- *Cutoff fraction*: keep only the top fraction of sites by probe count.
  ``0.0`` = show all results, ``0.5`` = top 50%.
  Default: ``0.0``

**Backbone**

- *Include backbone oxygens*: if enabled, backbone carbonyl oxygens are
  also considered as potential coordination donors in addition to
  sidechains.

- *Backbone residues*: residues whose backbone oxygens are included.
  ``ALL`` or a comma-separated list. Only active when the above option
  is enabled.

**Search zone**

- *Restrict search to a sphere*: if enabled, restricts the probe grid
  to a sphere instead of covering the whole protein.

- *Sphere centre*: XYZ coordinates of the sphere centre.
  Format: ``x,y,z`` (no spaces). Example: ``84.98,42.82,16.04``

- *Sphere radius*: radius in Ångstroms of the search sphere.
  Default: ``10.0``

- *Grid spacing*: distance in Ångstroms between probe grid points.
  Smaller values give higher detail but longer computation time.
  Default: ``1.0``

**Motif**

- *Search for a specific motif*: if enabled, searches for a specific
  coordination pattern instead of any valid site.

- *Motif*: coordination motif to search for.
  Format: three-letter codes comma-separated.
  Examples: ``HIS,HIS,ASP``, ``HIS/CYS,HIS,ASP``

**Output**

- *outputProbes*: an AtomStruct object containing the predicted
  metal-binding probe positions in PDB format. Can be opened directly
  in ChimeraX or PyMOL and superimposed on the input structure.


Usage example
-------------

1. Import a protein structure using **pwem > import atomic structure**.
2. Open **scipionbiometall > predict metal binding sites**.
3. Select the imported structure as input.
4. Set coordinating residues and minimum coordinators.
5. Optionally restrict the search zone or specify a motif.
6. Click **Execute**.
7. When finished, click **Analyze results** to open the probes in
   ChimeraX.

The output ``probes.pdb`` file is also saved at::

    ScipionUserData/projects/<project>/Runs/<run_id>_ProtBioMetAll/extra/probes.pdb

The full BioMetAll text output (cluster table with residues, coordinates
and probe counts) is available in the **Output Log** tab of the protocol.


Reference
---------

If you use this plugin in your work, please cite:

    Sanchez-Aparicio, J.-E., Tiessler-Sala, L., Velasco-Carneros, L.,
    Roldán-Martín, L., Sciortino, G., & Maréchal, J.-D. (2021).
    BioMetAll: Identifying Metal-Binding Sites in Proteins from Backbone
    Preorganization. *Journal of Chemical Information and Modeling*,
    61(1), 311–323.
    https://doi.org/10.1021/acs.jcim.0c00827