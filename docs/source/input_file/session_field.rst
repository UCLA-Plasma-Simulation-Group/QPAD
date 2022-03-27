Field Session
=============

In the ``field`` session, the ``diag`` sessions are set to diagnose field quantities.

* ``diag``, session array(\*), optional 
    Every type of field diagnostics must be provided as a session. The parameters of each session include

    * ``name``, string array(\*)
        Available options include
        
        * ``"er_cyl_m"``, ``"ephi_cyl_m"``, ``"ez_cyl_m"``. Diagnostics for :math:`E_r, E_\phi, E_z` respectively.
        * ``"br_cyl_m"``, ``"bphi_cyl_m"``, ``"bz_cyl_m"``. Diagnostics for :math:`B_r, B_\phi, B_z` respectively.
        * ``"spec_er_cyl_m"``, ``"spec_ephi_cyl_m"``. Diagnostics for plasma-induced :math:`E_r, E_\phi` respectively.
        * ``"spec_br_cyl_m"``, ``"spec_bphi_cyl_m"``. Diagnostics for plasma-induced :math:`B_r, B_\phi` respectively.
        * ``"jr_cyl_m"``, ``"jphi_cyl_m"``, ``"jz_cyl_m"``. Diagnostics for :math:`J_r, J_\phi, J_z` respectively.
        * ``"charge_cyl_m"``. Diagnostics for total charge.
        * ``"ar_cyl_m"``, ``"aphi_cyl_m"``, ``"az_cyl_m"``. Diagnostics for vector potential (induced by plasma) :math:`A_r, A_\phi, A_z` respectively.
        * ``"psi_cyl_m"``. Diagnostics for pseudo-potential.

    * ``ndump``, integer
        The code will dump the data every ``ndump`` time steps. If ``ndump`` is zero, the dumping is turned off.

.. note::

    Note that :math:`E_z` induced by the plasma is identical to the complete field, because the beams don not contribute to :math:`E_z`, therefore, there is no specific diagnostics for plasma-induced :math:`E_z`.

.. note::

    Diagnostics of vector potential need to call extra Poisson-like solvers and thus turning on these diagnostics will slightly affect the performance.

Example
-------

.. code-block:: json

  "field":
  [
      "diag":
      [
          {
          "name" : ["er_cyl_m", "ez_cyl_m", "bphi_cyl_m"],
          "ndump" : 1
          }
      ]
  ],