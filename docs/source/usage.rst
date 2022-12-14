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
>>>    'u' : mptb.TransportModels.staggeredTransport_u,
>>>    'v' : mptb.TransportModels.staggeredTransport_v,
>>>    'p' : mptb.PressureModels.Pressure
>>> }

Note that the different velocity components have a different transport model.


.. _definePassiveFields:

Define passive fields
^^^^^^^^^^^^^^^^^^^^^

In some cases we would want to prescribe a constant passive field, without solving it.
Depending on the chosen transport model, we have to provide neccessary fields here,
if they don't appear in the flow model dictionary.

For example, in a heat transfer problem the constant advection,
we prescribe a passive velocity field, which affects the heat conduction:

>>> passiveFields = {
>>>    'u' : 'faces_u',
>>>    'v' : 'faces_v'
>>> }

These fields are not solved by the CFD algorithm.
Hence, they don't have boundary conditions.
However, they can still be accessed or modified.

.. _initializeSimInstance:

Initialize the simulation instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step couples the mesh, the fields and the flow models together.
It creates a global object registry, which will be the entry point for later
field executions. If passive fields are required the initializaiton will also
create fields for these and couple them with the respective flow models.
It requires the flow model dictionary, the mesh, geometry and the (possibly emtpy)
passive fields dictionary to be passed:

>>> mptb.initialize(flowmodels=myFlowModels, mesh=mesh, geometry=geom, passiveFields=passiveFields )

After having initialized the fields, you can get a handle to the fields, which gives you access
to its data and its flow model:

>>> T = mptb.getField('T')

See :py:func:`Manager.initialize`.


.. _defBCs:

Define boundary conditions and initial field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the fields are initialized and tied to a flow model, they must be subjected to
boundary conditions. Use the afore-defined handle to the fields and specify a boundary
condition for every boundary in the domain. You can use the function
:py:func:`Manager.listAvailableBoundaryModels` to get a list of supported boundary conditions
for the field's flow model.

For example, a 'fixed value' boundary condition for the scalar field T can be defined as:

>>> mptb.defineBoundaryCondition(field=T, boundaryName='left', type='fixedValue', value=100 )

See :py:func:`Manager.defineBoundaryCondition` for details.

.. _setTranspModelParams:

Adjust transport models
^^^^^^^^^^^^^^^^^^^^^^^

Set the diffusion coefficent of the respective field. For momentum fields, this relates to
the viscosity. This value governs the diffusive transport of the respective flow variable.

In the example of a heat-conduction problem, this defines the thermal conductivity:

>>> T.govModel.setDiffusionCoefficient(value=1e-6)


.. _execSim:

Execute the simulations
^^^^^^^^^^^^^^^^^^^^^^^

The simulation is set up and is now ready for being executed.
The :py:func:`Manager.solve` function takes a field and solves it according to the field's
flow model. The function returns the raw data of one iteration step:

>>> T.data = mptb.solve(T)

From here, you can build iteration schemes combining various solve commands.

