Neutral Session
===============

The ``neutrals`` session is an array and each component is a session that defines the parameters related to each neutral gas species. 

* ``profile``, string array (2)
    Profile types for the plasma density distribution. The first and second components are the profile types of the transverse (in r-direction) and longitudinal (in :math:`\xi`-direction) directions. The available options for the transverse profile type include ``"uniform"``, ``"parabolic-channel"`` and ``"hollow-channel"``. The "uniform" option does not need extra parameters while the other types do.

    The ``"parabolic-channel"`` defines a parabolic plasma channel, and the following characteristic parameters need to be provided

    * ``channel_n0``, real
        On-axis channel density.
    * ``channel_depth``, real
        Channel density difference between r=0 and r=``channel_r0``.
    * ``channel_r0``, real
        Characteristic channel radius.
    * ``channel_width``, real
        Truncated radius of the plasma channel, beyond which the distribution becomes uniform.

    The ``"hollow-channel"`` defines plasma channel that has a uniform distribution between the inner and outer radii, and the following characteristic parameters need to be provided

    * ``channel_rmin``, real
        Inner channel radius.
    * ``channel_rmax``, real
        Outer channel radius.
    * ``channel_depth``, real
        Channel density difference between the inner and outer radii.

    The available options for the longitudinal profile type include ``"uniform"`` and ``"piecewise-linear"``. The ``"uniform"`` option does not need extra parameters while the other types do.

    The ``"piecewise-linear"`` defines a piecewise linear function according to which the plasma density will be updated for each 3D time step. The following parameters are needed

    * ``piecewise_s``, real array(\*)
        Time points of the piecewise linear function. They must be a monotonically increasing array.
    * ``piecewise_fs``, real array(\*) 
        Density defined on each time point. The length should be the same with ``piecewise_s``.

* ``ppc``, integer array (2)
    Numbers of macro-particles per cell in polar coordinate system, i.e. :math:`(\Delta r, \Delta \phi)`.

* ``num_theta``, integer
    Numbers of "virtual" cells distributed azimuthally, i.e. :math:`\Delta\phi=2\pi`/``num_theta``.

* ``npmax``, integer, optional
    Maximum particle number allowed for each partition. If not given, the code will try to determine this value automatically. The default value is twice of the minimum buffer size required to initialize the particles, i.e., ``ppc(1)`` * ``ppc(2)`` * ``num_theta`` * (number of cells in this partition) * 2.0. If the number of particles on a certain node exceeds the provided value, the code will attempt to reallocate the buffer to accommodate the extra particles. However, frequent buffer reallocation will impact the performance, which one should avoid. Warnings of buffer reallocation will be printed out so one can see how frequently this process is happening.

* ``q``, real
    Charge for each plasma particle. For example, it is ``-1.0`` for an electron.

* ``m``, real
    Rest mass for each plasma particle. For example. it is ``1.0`` for an electron and ``1836.15267389`` for a proton.

* ``density``, real
    Global density of the plasma. This factor will be multiplied to the density values defined in ``profile``.

* ``element``, integer
    Atom number of neutral gas species. The available include: H(1), He(2), Li(3), C(6), N(7), O(8), Ar(18), K(19), Rb(37), Xe(54) and Cs(55).

* ``ion_max``, integer
    Maximum ionization status allowed in the simulation.

* ``push_type``, string
    Species particle pusher type. Currently, the valid options are ``"robust"`` and ``"clamp"``. The ``"robust"`` push type is generally the first choice for most cases. The ``"clamp"`` push type will clamp the value of :math:`\gamma/(\gamma-p_z)` (specified by ``fac_clamp`` parameter) of the plasma particles. This can mitigate the code crashing or numerical instability that emerges at the back of the wake for highly nonlinear blowout regime. Note that the clamping is an artificial treatment and will lead to unphysical result. This method is currently provided only for experimental purposes.

* ``fac_clamp``, real, optional, default ``10.0``
    Clamped value of :math:`\gamma/(\gamma-p_z) used by ``"clamp"`` push type.

* ``den_min``, real, optional, default ``1.0d-10``
    It specifies the minimum density for injecting particles. Particles are only injected when the specified density is above this threshold.

* ``uth``, real(3), optional, default ``[0.0, 0.0, 0.0]``
    Initial thermal velocity (proper velocity) in (x,y,z) directions.
    
.. warning::

    This is an experimental functionality. Setting non-zero values to ``uth`` will lead to unphysical plasma oscillation, and the cause is unclear currently. If not necessary, one should always set this parameter zero.

* ``diag``, session array(\*), optional
    For the species, every type of diagnostics must be provided as a session. The parameters of each session include

    * ``name``, string array(\*)
        Available options include ``"charge_cyl_m"`` for dumping species charge density, ``"ion_cyl_m"`` for dumping ion density, and ``"raw"`` for dumping species particle raw data.
    * ``ndump``, integer
        The code will dump the data every ``ndump`` time steps. If ``ndump`` is zero, the dumping is turned off.
    * ``psample``, integer
        Only needed by ``"raw"`` diagnostic. The code will dump one particle raw data from every ``psample`` particles.

Example
-------

This example shows the settings for a uniform Lithum gas.

.. code-block:: json

  "neutrals" :
  [
      {
      "profile" : ["uniform", "uniform"],
      "ppc" : [8, 8],
      "num_theta" : 32,
      "q" : -1.0,
      "m" : 1.0,
      "density" : 1.0,
      "element" : 3,
      "ion_max" : 3,
      "push_type" : "robust",
      "den_min" : 1.0e-10,
      "diag" :
      [
          {
          "name" : ["charge_cyl_m", "ion_cyl_m"],
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
