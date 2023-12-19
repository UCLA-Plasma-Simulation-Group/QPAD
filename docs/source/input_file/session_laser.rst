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

    The ``"astrl_discrete"`` profile defines a discretization of the ASTRL integral defined in Pierce et al., Physical Review Research 5, 013085 (2023). This is accomplished by defining a grid of pulselets with centroids linearly spaced in :math:`\xi` The properties of each pulselet are determined by the specified math functions (``a0_math_func``, ``s0_math_func``, and ``w0_math_func``) evaluated at the central :math:`\xi` of the pulselet. The pulselets may overlap, in which case they are added coherently. The following  parameters need to be provided:

    * ``a0_math_func``, string
        Strength parameter of the pulselets as a function of :math:`\xi`. Ignores ``a0`` input parameter (not required).
    * ``s0_math_func``, string
        Distance of focal plane from the pulselets as a function of :math:`\xi`. If negative, the focal plane is behind the laser pulselet.
    * ``w0_math_func``, string
        Spot size of pulselets as a function of :math:`\xi`.
    * ``pulselet_math_func``, string
        Longitudinal profile of each pulselet in the discrete sum as a function of :math:`\xi`. This is the same as the function :math:`B(\xi)` used in Pierce et al., 2023. In the current implementation, each pulselet uses the same longitudinal profile, but this is not necessary in general.
    * ``pulselet_range``, real
        At a given :math:`xi` value, only pulselets whose centroid lies within ``pulselet_range`` of :math:`\xi` will be added to compute the field. This reduces the cost of the laser field initialization. ``pulselet_range`` should be several times the duration of the pulselets determined by ``pulselet_math_func``. 
    * ``pulselet_delay``, real
        Delay between the central :math:`\xi` of each of the pulselets. 
    * ``pulselet_offset``, real
        This specifies an offset to the grid of :math:`\xi` values over which the pulselet centroids in :math:`\xi` are defined. For example, if ``pulselet_offset = 0`` and ``pulselet_delay = 2``, then the centroids occur at :math:`\xi = \hdots, -4, -2, 0, 2, 4, \hdots`. If ``pulselet_offset = 1`` and ``pulselet_delay = 2``, then the centroids occur at :math:`\xi = \hdots -5, -3, -1, 1, 3, 5, \hdots`. 
    * ``if_norm_a0``, logical
        Whether to normalize the total amplitude of the total laser field to a particular value. If ``if_norm_a0 = True``, then the laser field is normalized such that the total normalized vector potential equals the deck parameter ``a0`` at the specified coordinates (``r_norm``, ``xi_norm,``, ``z_norm``) for the specific mode ``mode_norm``. Normalization is generally useful because the integration leads to the total laser amplitude being larger than the amplitude of each individual pulselet. 
    * ``r_norm``, real
        The transverse coordinate used for optional normalization of an ``astrl_discrete`` profile.
    * ``xi_norm``, real
        The ``xi`` coordinate used for optional normalization. This corresponds to the same :math:`xi` variable used in the math functions. 
    * ``z_norm``, real
        The ``z`` coordinate used for optional normalization. Setting ``z_norm = 0`` corresponds to performing the normalization at simulation initialization. Setting ``z_norm`` to a value greater than 0 corresponds to normalizing the field according to the approximate amplitude the pulse would have at the specified ``z_norm`` if the initial pulse propagated in vacuum. 
    * ``mode_norm``, integer
        Mode used for the normalization. 

    The ``"astrl_analytic"`` profile is the result of evaluation of the ASTRL integral in the above paper on ASTRL pulses when :math:`B(\xi)` is a Dirac delta function. This results in an ASTRL pulse whose transverse properties continuously vary as a function of :math:`\xi`. This general form has not yet been published but will be included in a forthcoming publication. The ``"astrl_analytic"`` profile does not use any of the pulselet profiles or normalization parameters used by the ``"astrl_discrete"`` class. The following must be specified:

    * ``a0_math_func``, string
        Strength parameter of the laser pulse as a function of :math:`\xi`. Ignores ``a0`` input parameter (not required).
    * ``s0_math_func``, string
        Distance of focal plane from the laser pulse as a function of :math:`\xi`. If negative, the focal plane is behind the laser pulse.
    * ``w0_math_func``, string
        Spot size of laser pulse as a function of :math:`\xi`.

    The available options for the longitudinal profile type include ``"sin2"``, ``"polynomial"``, and ``"astrl_analytic"``. The ``"sin2"`` defines a profile like :math:`sin^2(\pi\xi/2)` and the ``"polynomial"`` like :math:`10\xi^3-15\xi^4+6\xi^5` for :math:`0<\xi<1`. Both profile types have the same parameters:

    * ``t_rise``, real
        Length of rising edge.
    * ``t_flat``, real
        Length of intensity plateau.
    * ``t_fall``, real
        Length of falling edge.

    The ``"const"`` longitudinal profile defines a constant profile equal to :math:`1`. The ``"const"`` profile should be selected for the longitudinal direction if the ``"astrl_analytic"`` or ``"astrl_discrete"`` profile is chosen for the transverse direction.

    The ``"piecewise-linear"`` defines a piecewise linear function to describe the laser's longitudinal profile.

    * ``piecewise_t``, real array(\*)
        Time points of the piecewise linear function. They must be a monotonically increasing array.
    * ``piecewise_ft``, real array(\*) 
        Density defined on each time point. The length should be the same as ``piecewise_t``.

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


This example shows the settings for an analytic ASTRL laser pulse using mathematical functions for the focal distance ``s0_math_func``, the vector potential ``a0_math_func``, and the spot size ``w0_math_func`` along the :math:`\xi`-direction.

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

