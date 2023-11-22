Laser Session
===============

.. warning::

    Note that the laser module is still an experimental feature. The functionality might be incomplete and unstable.

    The current laser module DOES NOT include the laser-beam interaction and photoionization features.

The "lasers" is an array and each component is a session defines the parameters of laser pulses which are described by ponderomotive guiding center (PGC) model.

* ``profile``, string array (2)
    Profile types for laser field envelope. The first and second components are the profile types of the transverse (in r-direction) and longitudinal (in :math:`\xi`-direction) directions. The available options for the transverse profile type include ``"gaussian"``, ``"laguerre"``, and ``"astrl_analytic"``.

    The ``"gaussian"`` profile defines a tri-Gaussian laser intensity distribution, and the following characteristic parameters need to be provided

    * ``w0``, real
        Radius (:math:`1/e^2`) of laser pulse.
    * ``focal_distance``, real
        Distance of focal plane from the laser pulse. If negative, the focal plane is behind the laser pulse.

    The ``"laguerre"`` defines a Laguerre-Gaussian laser pulse, and the following characteristic parameters need to be provided

    * ``w0``, real
        Radius of laser pulse.
    * ``focal_distance``, real
        Distance of focal plane from the laser pulse. If negative, the focal plane is behind the laser pulse.
    * ``radial_index``, real
        Radial index of Laguerre-Gaussian mode.
    * ``phi_index``, real
        Azimuthal index of Laguerre-Gaussian mode.

    The ``"astrl_analytic"`` defines an Arbitrarily Structured Laser (ASTRL) pulse, and the following characteristic parameters need to be provided

    * ``w0_math_func``, string
        Radius of laser pulse as a function of :math:`\xi`.
    * ``s0_math_func``, string
        Distance of focal plane from the laser pulse as a function of :math:`\xi`. If negative, the focal plane is behind the laser pulse.
    * ``a0_math_func``, string
        Strength parameter of the laser pulse as a function of :math:`\xi`. Ignores ``a0`` input parameter (not required). 

    The available options for the longitudinal profile type include ``"sin2"``, ``"polynomial"``, and ``"astrl_analytic"``. The ``"sin2"`` defines a profile like :math:`sin^2(\pi\xi/2)` and the ``"polynomial"`` like :math:`10\xi^3-15\xi^4+6\xi^5` for :math:`0<\xi<1`. Both profile types have the same parameters:

    * ``t_rise``, real
        Length of rising edge.
    * ``t_flat``, real
        Length of intensity plateau.
    * ``t_fall``, real
        Length of falling edge.

    The ``"astrl_analytic"`` longitudinal profile defines a constant profile equal to :math: `1`. The ``"astrl_analytic"`` profile should be selected for the longitudinal direction if already selected for the transverse direction.
    .. The ``"piecewise-linear"`` defines a piecewise linear function according to which the plasma density will be updated for each 3D time step. The following parameters are needed

    .. * ``piecewise_s``, real array(\*)
    ..     Time points of the piecewise linear function. They must be a monotonically increasing array.
    .. * ``piecewise_fs``, real array(\*) 
    ..     Density defined on each time point. The length should be the same with ``piecewise_s``.

* ``iteration``, integer
    Interation number of laser PGC solver. 3 or more is recommended to ensure convergence.

* ``k0``, real
    Wavenumber corresponding to the central frequency of the laser pulse.

* ``a0``, real
    Strength parameter (maximum normalized vector potential) of the laser pulse. Parameter skipped if ``profile`` contains ``"astrl_analytic"`` (see ``a0_math_func``).

* ``lon_center``, real
    Center of laser pulse in :math:`\xi`-direction. 

*  ``chirp_coefs``, real array (*), optional
    Frequency chirp coefficients :math:`C` of laser pulse. The frequency chirp distribution is described by :math:`k(\xi)=k_0+C(1)(\xi-\xi_0)+C(2)(\xi-\xi_0)^2+...` where :math:`k_0` and :math:`\xi_0` are the central wavenumber and longitudinal center defined by ``k0`` and ``lon_center`` respectively. The default is [0.0].

* ``diag``, session array(\*), optional
    For lasers, every type of diagnostics must be provided as a session. The parameters of each session include

    * ``name``, string array(\*)
        Currently the available option only includes ``"a_cyl_m"`` for dumping local normalized vector potential.
    * ``ndump``, integer
        The code will dump the data every ``ndump`` time steps. If ``ndump`` is zero, the dumping is turned off.

Example
-------

This example shows the settings for a Gaussian Laser pulse with a `sin^2` longitudinal profile.

.. code-block:: json

    "laser" :
    [
        {
            "profile" : ["gaussian", "sin2"],
            "iteration" : 3,
            "k0" : 20.0,
            "a0" : 2.0,
            "w0" : 2.828427,
            "focal_distance" : 0.0,
            "lon_center" : 0.0,
            "t_rise" : 2.0,
            "t_flat" : 0.0,
            "t_fall" : 2.0,
            "diag" :
            [
                {
                    "name" : ["a_cyl_m"],
                    "ndump" : 1
                }
            ]
        }
    ],


This example shows the settings for an analytic ASTRL laser pulse using mathematical functions for the focal distance :math: `s0(\xi)`, the vector potential :math: `a_0(\xi)`, and the spot size :math: `w_0(\xi)` along the :math:`\xi`-direction.

.. code-block:: json

    "laser": [
        {
            "profile": [
                "astrl_analytic",
                "astrl_analytic"
            ],
            "iteration": 3,
            "k0": 10.0,
            "lon_center": 0.0,
            "s0_math_func": "200 - 20 * xi",
            "a0_math_func": "if( xi < 0, 0, if( xi < (1), sin(1.5707 * xi / 1)^2, if( xi  < (9), 1, if( xi < (10), sin(1.5707 * (xi-10) / 1)^2, 0))))",
            "w0_math_func": "2.0",
            "diag": [
                {
                    "name": [
                        "a_cyl_m"
                    ],
                    "ndump": 5
                }
            ]
        }
    ]

