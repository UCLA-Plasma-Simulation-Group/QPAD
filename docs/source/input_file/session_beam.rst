Beam Session
============

In the input file, each component in the ``beam`` session is a session that defines the parameters related to each particle beam. Currently, there are three types of beam sources configured through ``profile_type`` parameter in the ``beam`` session. They are ``"standard"``, ``"random"`` and ``"file"`` which correspond to

* "Lattice-like" initialization of macro-particles with varying charge
* Initialization using probability distribution function with fixed charge
* Read macro-particles from files

Lattice-like Initialization
---------------------------

Setting ``profile_type`` to be ``"standard"`` enables the "Lattice-like" initialization of beam particles. This method initially inject a fixed number of macro-particles into the simulation box with even spacing in r-, :math:`\phi`- and :math:`\xi`-directions. The charge carried by each macro-particle varies according the local beam density.

* ``profile_type``, string
    Source type of the beam particles. Available options include ``"standard"``, ``"random"`` and ``"file"``. Here it should be set as ``"standard"``.

* ``geometry``, string, optional, default ``"cartesian"``
    Geometry for the particle initialization, including ``"cartesian"`` and ``"cylindrical"``. The initialization assumes the beam density distribution is variable-separable in three directions, i.e., 1-, 2- and 3-direction. Note that the choice of the value will affect the configuration of some of the following parameters. Setting ``geometry`` to be ``"cartesian"`` means the 1-, 2- and 3-directions correspond to x-, y- and :math:`\xi`-directions, or otherwise r-, :math:`\phi`- and :math:`\xi`-directions if ``"cylindrical"`` is selected.

* ``profile``, string array(3)
    Profile types for the beam density distribution. The three components correspond to the profile types in 1-, 2- and 3-direction. The available options for the transverse profile type include ``"uniform"``, ``"gaussian"``, ``"parabolic"``, ``"rational"`` and ``"piecewise-linear"``. The ``"uniform"`` option does not need extra parameters while the other types do.

    The ``"gaussian"`` type defines a Gaussian beam profile, and the following characteristic parameters need to be provided
    
    * ``gauss_center``, real array (3)
        The centers in 1-, 2- and 3-direction.
    * ``gauss_sigma``, real array (3)
        The rms sizes in 1-, 2- and 3-direction.

    The ``"parabolic"`` type defines a parabolic beam profile like :math:`1-(x-x_0)^2/a^2` when :math:`|x|<a` (a is the radius). The following characteristic parameters need to be provided.

    * ``parabolic_center``, real array(3)
        The center :math:`x_0` in 1-, 2- and 3-direction.
    * ``parabolic_radius``, real array(3)
        The radius a in 1-, 2- and 3-direction.

    The ``"rational"`` type defines a rational beam profile like :math:`\sum_k a_k(x-x_0)^k / \sum_k b_k(x-x_0)^k`. The following characteristic parameters need to be provided

    * ``poly_numerator1``, ``poly_numerator2``, ``poly_numerator3``, real array (\*)
        The arrays of polynomial coeffcients :math:`(a_0, a_1, a_2, \cdots)` in 1-, 2- and 3-direction respectively.
    * ``poly_denominator1``, ``poly_denominator2``, ``poly_denominator3``, real array (\*)
        The arrays of polynomial coeffcients :math:`(b_0, b_1, b_2, \cdots)` in 1-, 2- and 3-direction respectively.
    * ``poly_x0``: real array (3)
        The center :math:`x_0` in 1-, 2- and 3-direction.

    The ``"piecewise-linear"`` defines a piecewise linear function according to which the beam density will be initialized. The following parameters are needed
    
    * ``piecewise_x1``, ``piecewise_x2``, ``piecewise_x3``: real array(\*)
        The arrays of position of the piecewise linear function in 1-, 2- and 3-direction respectively. They must be a monotonically increasing array.
    * ``piecewise_fx1``, ``piecewise_fx2``, ``piecewise_fx3``: real array(\*)
        The density defined on each points defined by ``piecewise_x1``, ``piecewise_x2`` and ``piecewise_x3``.

* ``evolution``, logical, optional, default ``"true"``
    If it is true, the code will update the momentum of beam particles every time step.

* ``push_type``, string, optional, default ``"reduced"``
    Type of particle pusher. The available options are ``"reduced"`` and ``"boris"``. The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is slightly less efficient than the reduced pusher.

* ``ppc``, integer array (3)
    The numbers of particles per "virtual" cell in the cylindrical coordinate, i.e. :math:`(\Delta r, \Delta\phi, \Delta z)`. Note it is not affected by the choice of ``geometry``.

* ``num_theta``, integer
    Numbers of "virtual" cells distributed azimuthally, i.e. :math:`\Delta\phi=2\pi`/``num_theta``.

* ``npmax``, integer, optional
    Number of particles allowed for this MPI partition. If not given, the program will automatically calculate an initial guess for this parameter. If ``npmax`` is not large enough during the initialization, the program will automatically resize the particle buffers.
    
.. note::

    Manually setting ``npmax`` to be sufficiently large value can avoid frequent buffer reallocation which may slow down the simulation. However, for some memory-intense tasks, letting the code automatically determine the initial ``npmax`` may more effciently exploit the memory.

* ``den_min``, real, optional, default ``1.0d-10``
    It specifies the minimum density for injecting particles. Particles are only injected when the specified density is above this threshold.

* ``range1``, ``range2``, ``range3``, real array(2)
    The three arrays specifies the lower and upper boundaries in 1-, 2- and 3-direction within which the particles are injected. The particles beyond this region will not be initialized.

* ``has_spin``, logical, optional, default ``"false"``
    Switch of spin dynamics. When this parameter is true, extra coordinates of spin :math:`(s_x, s_y, s_z)` will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles, i.e., ``profile_type`` is ``"file"``. For other beam source types, this parameter must be set as ``"false"``.

* ``q``, real
    Charge for each beam particle. For example, it is ``-1.0`` for an electron and ``1.0`` for a proton or positron.

* ``m``, real
    Rest mass for each beam particle. For example, it is ``1.0`` for an electron and ``1836.15267389`` for a proton.

* ``gamma``, real
    The Lorentz factor for the average energy of the particle beam.

* ``density``, real
    Global multiplication factor for the density profile. Regardless of which profile type you choose the final density value will be product of ``density`` and the value set in the specific beam profile.

* ``quiet_start``, logical, optional, default ``"false"``
    Switch of initializing the beam particles using the "quiet start" method. If it is turned on, a set of image particles will be added to suppress the statistic noise. Note that with this feature on, the total particle number will be doubled.

* ``uth``, real array (3), optional, default ``[0, 0, 0]``
    The thermal proper velocity in x-, y- and z-direction (note it is not affected by ``geometry``.). The thermal distribution is subject to Gaussian distribution. This is usually used to initialize rms beam divergence.

* ``alpha``, real array (2), optional, default ``0``
    The Twiss parameter :math:`\alpha` in x- and y-directions. This is used to initialize a tilt phase-space ellipse. Note that this parameter is **only available** for ``"geometry": "cartesian"`` and Gaussian profile in x- and/or y-directions. The Twiss parameter :math:`\beta` and the emittance will be automatically calculated from ``gauss_sigma``, ``gamma`` and ``uth``, so only :math:`\alpha` needs to be given explicitly.

* ``perp_offset_x``, ``perp_offset_y``, real array(\*), optional
    These two parameters are used to set the transverse position offset in x- and y-directions as a function of :math:`\xi`. Taking ``perp_offset_x`` for example, its form looks like :math:`[\xi_0, P_0, P_1, \cdots]` where :math:`\xi_0` is the reference position and :math:`P_i` are the coefficients of a polynomial. The transverse offset is given by :math:`\Delta x=\sum_{k=0} P_k(\xi-\xi_0)^k`. The configuration in y-direction is similar.

* ``diag``, session array (\*), optional
    Every type of diagnostics must be provided as a session. The parameters of each session include

    * ``name``, string array (\*)
        Available options include ``"charge_cyl_m"`` for dumping beam charge density, and ``"raw"`` for dumping beam particle raw data.
    * ``ndump``, integer
        The code will dump the data every ``ndump`` time steps. The data dump is turned off for ``"ndump": 0``.
    * ``psample``, integer
        Only needed by ``"raw"``diagnostic. The code will dump one particle raw data from every ``psample`` particles.

Initialization using probability distribution function
------------------------------------------------------

Setting ``profile_type`` to be ``"random"`` enables this type of initialization of beam particles. This method initially inject macro-particles into the simulation box using the probability distribution functions of various density profiles. The charge carried by each macro-particle is the same.

* ``profile_type``, string
    The source type of the beam particles. Available options include ``"standard"``, ``"random"`` and ``"file"``. Here it should be set as ``"random"``.

* ``total_num``, integer
    The total number of particles of the entire beam.

* ``total_charge``, real
    The total charge of the beam in the unit of :math:`en_pc^3\omega_p^{-3}`.

.. note::

    For simple beam density profiles, it is easy to connect the total charge with the peak density. Taking the tri-Gaussian profile for example, the total charge :math:`Q=en_b(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z`.
    
.. However, for arbitrary beam profiles it is usually impossible to exactly know the total charge from the peak density, or vice versa. In some special situations where one only knows the peak density but needs to use "random" initialization, a useful trick to know the total charge is:

.. - First, set arbitrary total charge with sufficient number of macro-particles and run the simulation one time step to obtain the initial beam density distribution.
.. - Second, read the peak density from the datasets. Since the number of particles is large enough to suppress the statistic fluctuation in the peak density, the reading should be very accurate. 
.. - Third, scale the total charge according to the desired and reading values of the peak density (the total charge is proportional to the peak density).

Other available parameters for ``"random"`` type beam profile include ``geometry``, ``profile``, ``evolution``, ``push_type``, ``npmax``, ``range1``, ``range2``, ``range3``, ``has_spin``, ``q``, ``m``, ``gamma``, ``quiet_start``, ``uth``, ``alpha``, ``perp_offset_x``, ``perp_offset_y``, ``diag``, and their definitions and configuration are identical to those of the ``"standard"`` profile type.

Importing particles from a HDF5 file
------------------------------------

Setting ``profile_type`` to be ``"file"`` will import macro-particles from a HDF5 file. This file should contains seven datasets named ``"x1"``, ``"x2"``, ``"x3"``, ``"p1"``, ``"p2"``, ``"p3"`` and ``"q"`` which corresponds to the beam positions and momenta in x-, y- and z-direction (not :math:`\xi`-direction), and the charge per particle.

* ``profile_type``, string
    Source type of the beam particles. Here it should be set as "file".

* ``filename``, string
    Name of the HDF5 file.

* ``anom_mag_moment``, real
    Anomalous magnet moment of the particle. Used for spin dynamics.

* ``beam_center``, real array(3)
    Cartesian coordinates :math:`(x, y, \xi)` of the beam center.

* ``file_center``, real array(3)
    Cartesian coordinates :math:`(x, y, z)` of beam center in the HDF5 file.

* ``length_conv_fac``, real, optional, default ``1.0``
    The scaling factor of the quantities with a length dimension. This is often used when the beam defined in the HDF5 file and the QPAD simulation have different reference density. With this parameter, the beam size will be scaled by ``length_conv_fac`` times.

* ``charge_conv_fac``, real, optional, default ``1.0``
    The scaling factor of the charge per particle. This is often used when the beam defined in the HDF5 file and the QPAD simulation have different reference density, or when the beam defined in the HDF5 file is extracted from other simulation (e.g. `OSIRIS <http://epp.tecnico.ulisboa.pt/osiris/>`__) with different cell volume. With this parameter, the charge per particle will be multiplied by ``charge_conv_fac``.

Other available parameters for ``"file"`` type beam profile include  ``evolution``, ``push_type``, ``npmax``, ``has_spin``, ``q``, ``m``, ``diag``, and their definitions and configuration are identical to those of the ``"standard"`` profile type.

Examples
--------

The following example shows the initialization of a beam with Gaussian transverse profile and a sawtooth longitudinal profile using the cylindrical geometry.

.. code-block:: json

  "beam" :
  [
      {
      "profile_type" : "standard",
      "geometry" : "cylindrical",
      "profile" : ["gaussian", "uniform", "piecewise-linear"],
      "evolution" : true,
      "push_type" : "reduced",
      "has_spin" : false,
      "ppc" : [2, 2, 2],
      "num_theta" : 16,
      "npmax" : 20000000,
      "q" : -1.0,
      "m" : 1.0,
      "gamma" : 20000,
      "density" : 4.0,
      "quiet_start" : true,
      "gauss_center" : [0.0, "none", "none"],
      "gauss_sigma" : [0.25, "none", "none"],
      "piecewise_x3" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
      "piecewise_fx3" : [0.0, 1.0, 0.1, 1.0, 0.2, 1.0, 0.3, 1.0, 0.4, 1.0, 0.0],
      "range1" : [0, 1.25],
      "range2" : [0, 6.283185307179586],
      "range3" : [-2.5, 2.5],
      "uth" : [0.0, 0.0, 0.0],
      "den_min" : 1e-10,
      "diag" :
      [
          {
          "name" : ["charge_cyl_m"],
          "ndump" : 1
          },
          {
          "name" : ["raw"],
          "ndump" : 1,
          "psample" : 10
          }
      ]    
      }
  ],

This can also be realized by using the Cartesian geometry.

.. code-block:: json

  "beam" :
  [
      {
      "profile_type" : "standard",
      "geometry" : "cartesian",
      "profile" : ["gaussian", "gaussian", "piecewise-linear"],
      "evolution" : true,
      "push_type" : "reduced",
      "has_spin" : false,
      "ppc" : [2, 2, 2],
      "num_theta" : 16,
      "npmax" : 20000000,
      "q" : -1.0,
      "m" : 1.0,
      "gamma" : 20000,
      "density" : 4.0,
      "quiet_start" : true,
      "gauss_center" : [0.0, 0.0, "none"],
      "gauss_sigma" : [0.25, 0.25, "none"],
      "piecewise_x3" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
      "piecewise_fx3" : [0.0, 1.0, 0.1, 1.0, 0.2, 1.0, 0.3, 1.0, 0.4, 1.0, 0.0],
      "range1" : [-1.25, 1.25],
      "range2" : [-1.25, 1.25],
      "range3" : [-2.5, 2.5],
      "uth" : [0.0, 0.0, 0.0],
      "den_min" : 1e-10,
      "diag" :
      [
          {
          "name" : ["charge_cyl_m"],
          "ndump" : 1
          },
          {
          "name" : ["raw"],
          "ndump" : 1,
          "psample" : 10
          }
      ]    
      }
  ]
