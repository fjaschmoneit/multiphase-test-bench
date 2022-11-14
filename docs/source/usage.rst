Usage
=====


Flow models are developed and executed through Python script.
The workflow below indicates the neccessary steps for the simulation set-up and execution.
A simulation instance can be defined within a single script, see for example VALIDATION test,
or distributed across several scripts, as for example in APPLICATION.
This flexible simulation set up allows you to



Workflow
--------

- :ref:`importMTB`
- :ref:`defineGeomMesh`
- :ref:`selectModels`
- :ref:`definePassiveFields`
- :ref:`initializeSimInstance`
- :ref:`defBCs`
- :ref:`setTranspModelParams`
- :ref:`execSim`




.. _importMTB:

Import the Multiphase Test-Bench
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This project is not available as a pip package.
After the sources are downloaded you can access the MPTB with:

>>> import sys
>>> sys.path.append("PATH_TO_SOURCE/MultiphaseTestBench")
>>> import Manager an mptb




.. _defineGeomMesh:

Define the geometry and mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation instance needs both, a geometry object, and a mesh object.


For example, to generate a rectangular 2D geometry with default boundary names,
and to create a mesh mesh from the geometry, simply write:

>>> geom = mptb.createGeometry( 'rectangle', [0.4, 0.3] )
>>> mesh = mptb.createMesh( geom, res=resolution )

You first create or import a geometry using the ``mptb.createGeometry`` function:

    .. autofunction:: Manager.createGeometry




.. _selectModels:

Select your flow models
^^^^^^^^^^^^^^^^^^^^^^^


.. _definePassiveFields:

Define passive fields
^^^^^^^^^^^^^^^^^^^^^


.. _initializeSimInstance:

Initialize the simulation instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _defBCs:

Define boundary conditions and initial field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _setTranspModelParams:

Adjust transport models
^^^^^^^^^^^^^^^^^^^^^^^

.. _execSim:

Execute the simulations
^^^^^^^^^^^^^^^^^^^^^^^






..
    Creating recipes
    ----------------

    To retrieve a list of random ingredients,
    you can use the ``lumache.get_random_ingredients()`` function:

    .. autofunction:: lumache.get_random_ingredients

    The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
    or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
    will raise an exception.

    .. autoexception:: lumache.InvalidKindError

    For example:

    >>> import lumache
    >>> lumache.get_random_ingredients()
    ['shells', 'gorgonzola', 'parsley']

