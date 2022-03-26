Input File
==========

QPAD requires an input file for initializing the simulation. When QPAD starts running, it will look for the input file named "qpinput.json" located in the current working directory and read parameters from this file. If the file does not exist in the same folder where the executable file is located, the program will stop with an error message.

The Format of the Input File
----------------------------

The input file is written in JSON format. It should be noted that JSON does not allow trailing commas and comment. However, in QPAD we use `JSON-FORTRAN <https://github.com/jacobwilliams/json-fortran>`__ which ignore all the comments before processing the file. In "qpinput.json", you can add any comment after an exclamation mark "!", and any text following the "!" will be ignored until the end of the line is reached.

Sessions of Input File
----------------------

The input file "qpinput.json" consists of ``simulation``, ``beam``, ``species``, ``neutrals`` and ``field`` sessions.

.. toctree::
   :maxdepth: 1

   session_simulation
   session_beam
   session_species
   session_neutral
   session_field