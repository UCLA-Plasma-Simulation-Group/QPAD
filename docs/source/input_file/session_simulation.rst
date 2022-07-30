Simulation Session
==================

* ``algorithm``, string
    The simulation mode. Available modes include ``"standard"`` and ``"popas"``. The ``"popas"`` mode provides the interfaces and emittance/energy spread diagnostics for the beam parameter optimization code developed by ANL.

* ``nodes``, integer array (2)
    The first component define the MPI processors used in each stage (in r-direction); the second component defines the number of pipeline stages used in the simulation (in :math:`\xi`-direction).

* ``grid``, integer array (2)
    The first and second components define the number of cells used along r- and :math:`\xi`-direction respectively.

* ``box``, session
    Define the size of simulation box in r- and :math:`\xi`-direction. This session includes two parameters

    * ``r``, real array(2)
        Starting and end coordinates in r-direction. The lower limit should be always 0.
    * ``z``, real array(2)
        Starting and end coordinates in :math:`\xi`-direction.

* ``field_boundary``, string
    Boundary condition of fields for the simulation. The current available options only include ``"open"``.

* ``max_mode``, integer
    Maximum harmonic number of Fourier azimuthal mode.

* ``interpolation``, string
    Order of interpolation of the simulation. Only ``"linear"`` is available currently.

* ``n0``, real
    Reference plasma density (in the unit of :math:`\text{cm}^{-3}`). **This parameter is required when the ionization is turned on**.

* ``time``, real
    Total duration of the simulation. The initial moment defaults to 0.

* ``dt``, real
    Time step of pushing relativistic particle beams.

* ``nbeams``, integer
    Total number of particle beams.

* ``nspecies``, integer
    Total number of plasma species.

* ``nneutrals``, integer
    Total number of neutral gas species.

* ``nlasers``, integer
    Total number of laser pulses.

.. warning::

    The laser module is still not available in this development version of QPAD, but ``nlasers`` parameter still needs to be provided just as a placeholder.

* ``dump_restart``, logical
    If true, the code will dump the restart files during the simulation.

* ``ndump_restart``, integer
    Frequency of dumping restart files. QPAD will dump the restart files every ``ndump_restart`` steps.

* ``read_restart``, logical
    If true, QPAD will read the restart files when it starts.

* ``restart_timestep``, integer
    The number of time step from which the code will restart when ``read_restart`` is true.

* ``iter_reltol``, real
    The tolerance of the relative error for the transverse magnetic field in the predictor-corrector iterations. A typical choice might be 1.0e-2 to 1.0e-3.

* ``iter_abstol``, real
    The tolerance of the absolute error for the transverse magnetic field in the predictor-corrector iterations.

.. note::
    The predictor-corrector iteration terminates when the relative error reaches ``iter_reltol`` OR the absolute error reaches ``iter_abstol``. The ``iter_abstol`` parameter is needed since the predictor-corrector iteration is sometimes found to converge slowly where the plasma response is perturbative if only relying on the relative error setting by ``iter_reltol``. To avoid a waste of computation time at these relatively unimportant regions, ``iter_abstol`` is often set to be on the same order of magnitude of (or 1 order of magnitude smaller than) the perturbative magnitude. Taking the nonlinear PWFA simulation for example, assuming the maximum transverse B field is around 1.0 [simulation unit] and you want to iterate quickly where the transverse B field is smaller than 0.001, setting ``iter_abstol`` to be 0.001 will make the code only iterate a very few times (typically 1 to 2) to terminate. However, for those situations that these perturbative regions are also important and need to be modeled very accurately, a smaller ``iter_abstol`` is often needed.

* ``iter_max``, integer
    The maximum number of predictor-corrector iterations the code will use when solving the plasma response.

* ``smooth_type``, string
    Smoother used to smooth the source terms. The available options include ``"binomial"`` and ``"compensated"``. Set ``"none"`` to switch off the smoother. 

* ``smooth_order``, integer
    Order of smooth. It only takes effect for "binomial" and "compensated" types.

.. warning::

    The smoother is not currently available now. Please turn it off by setting ``"smooth_type": "none"`` until this issue is fixed in the future updates.

* ``verbose``, integer
    Level of verbosity. This parameter should be always set as 0 unless for debugging purpose. A larger ``verbose`` will output more debugging information into the log files. The maximum value is 5.

* ``random_seed``, integer
    Seed of the pseudo-random number generator. Set 0 to use the operating-system-specified seeds, which varies from run to run. Specify any integer greater than 0 to use the user-specified seed which is fixed for each run.