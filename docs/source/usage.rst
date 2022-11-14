Usage
=====


Flow models are developed and executed through Python script.
The workflow below indicates the neccessary steps for the simulation set-up and execution.
A simulation instance can be defined within a single script, see for example VALIDATION test,
or distributed across several scripts, as for example in APPLICATION.
This flexible simulation set up allows you to



Workflow
--------

:ref:`importMTB`

:ref:`defineGeomMesh`

:ref:`selectModels`

:ref:`definePassiveFields`

:ref:`initializeSimInstance`

:ref:`defBCs`





.. _importMTB:
Import the Multiphase Test-Bench
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



.. _defineGeomMesh:
Define the geometry and mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

