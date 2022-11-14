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

The simulation is conducted from a global Manager.
For a high-level application of the test-bench, it is sufficient to only import the Manager.
After the sources are downloaded you can access the MPTB with:

>>> import sys
>>> sys.path.append("PATH_TO_SOURCE/MultiphaseTestBench")
>>> import Manager an mptb


.. _defineGeomMesh:

Define the geometry and mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation instance needs both, a geometry object, and a mesh object.
The geometry holds information on the dimension, the boundary names and possible regions.
The mesh is contains information about the discretization and fundamental mesh operations.

For example, to generate a rectangular 2D geometry with default boundary names,
and to create a mesh mesh from the geometry, simply write:

>>> geom = mptb.createGeometry( 'rectangle', [0.4, 0.3] )
>>> mesh = mptb.createMesh( geom, res=resolution )

See :py:func:`Manager.createGeometry` and :py:func:'Manager.createMech' for details.


.. _selectModels:

Select your flow models
^^^^^^^^^^^^^^^^^^^^^^^

Every flow variable is governed by its flow model.
Here you define a name for every variable and link it to a flow model type.
Neither the fields, nor the flow models are initialized here.
They will be brought together with the mesh in the intialize step (see below).

A simple heat transfer problem with only one variable could be set up like this:

>>> myFlowModels = {
>>>    'T' : mptb.TransportModels.scalarTransport
>>> }

For a 2D incompressible flow problem, we could set up the flow model like this:

>>> myFlowModels = {
>>>    'u' : Odin.TransportModels.staggeredTransport_u,
>>>    'v' : Odin.TransportModels.staggeredTransport_v,
>>>    'p' : Odin.PressureModels.Pressure
>>> }

Note that the different velocity components have a different transport model.


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

