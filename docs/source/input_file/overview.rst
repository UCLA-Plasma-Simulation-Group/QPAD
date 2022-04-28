Input File
==========

QPAD requires an input file for initializing the simulation. When QPAD starts running, it will look for the input file named "qpinput.json" located in the current working directory and read parameters from this file. If the file does not exist in the same folder where the executable file is located, the program will stop with an error message.

Format of Input File
----------------------------

The input file is written in JSON format. It should be noted that JSON does not allow trailing commas and comment. However, in QPAD we use `JSON-FORTRAN <https://github.com/jacobwilliams/json-fortran>`__ which ignore all the comments before processing the file. In "qpinput.json", you can add any comment after an exclamation mark "!", and any text following the "!" will be ignored until the end of the line is reached.

Parameter Unit
--------------

Unless otherwise specified, all input parameters and output data are in the normalized units based on the reference plasma density :math:`n_p`. The normalized units of physical quantities of different dimensions are listed below

================ ====================== 
**dimension**    **normalized unit**
Length           :math:`c/\omega_p`
Time             :math:`1/\omega_p`
Mass             :math:`m_e`
Charge           :math:`e`
Momentum         :math:`m_e c`
Velocity         :math:`c`
Density          :math:`n_p`
Scalar potential :math:`m_ec^2/e`
Vector potential :math:`m_ec/e`
Electric field   :math:`\omega_pm_ec/e`
Magnetic field   :math:`\omega_pm_e/e`
================ ======================

where :math:`\omega_p\equiv\sqrt{n_pe^2/\epsilon_0m_e}` is the reference plasma frequency, :math:`c` the speed of light, :math:`e` the elemental charge and :math:`m_e` the rest mass of electron.

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
