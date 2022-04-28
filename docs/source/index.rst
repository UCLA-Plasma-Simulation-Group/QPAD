.. QPAD documentation master file, created by
   sphinx-quickstart on Fri Mar 25 21:17:08 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QPAD documentation!
==============================

Introduction
------------

QuickPIC with Azimuthal Decomposition (QPAD) is a quasi-3D parallel (MPI-based) quasi-static PIC code written in modern Fortran, which is developed based on the framework of full 3D PIC code `QuickPIC`_. QPAD is a very efficient PIC code for modeling short-pulse laser or relativistic charged particle beam–plasma interactions. QPAD is the first PIC code that combines both quasi-static approximation and azimuthal harmonic decomposition in the community of plasma based accelerators. This can speed up the simulation by orders of magnitude compared with the 3D fully explicit electromagnetic PIC code or full 3D quasi-static code. QPAD utilizes synchronous parallelism in the transverse direction and asynchronous parallelism (pipeline) in the longitudinal direction, which maximizes the parallel scalability. The current version implements many useful features for simultion of plasma based accelerators, including

.. _QuickPIC: https://github.com/UCLA-Plasma-Simulation-Group/QuickPIC-OpenSource

- Various particle pushers
- Field-induced ionization (ADK model)
- Spin dynamics for beam particles
- Multiple initialization for beams and plasma species
- Import beam particles from file
- Parallel cyclic reduction Poisson solver

More simulation features and hardware supports are being developed, including

- GPU support
- Dynamic load balance
- Adaptive time step
- Ponderomotive guiding center model for laser
- Radiation reaction
- Others

Contents
--------

.. toctree::
   :maxdepth: 1

   overview/overview
   installation/installation
   input_file/overview
   post_process/post_process
   developer/contribution

References & Attribution
------------------------

QPAD was originally developed by Fei Li (UCLA) and Weiming An (BNU) and is currently being actively developed and maintained by `PICKSC group <https://picksc.physics.ucla.edu/>`__ at UCLA. If you install QPAD, we ask that you *please contact* **Fei Li** (lifei11@ucla.edu) or **Weiming An** (anweiming@bnu.edu.cn). The development of QPAD relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies. We are pleased to see QPAD can help your researches. For scientific publications, please consider to cite the original paper of QPAD:

- F\. Li et al., Comput. Phys. Commun. (2021), 261, 107784. (`Link <https://doi.org/10.1016/J.CPC.2020.107784>`__)

Other possible relavant papers includes:

- QuickPIC original paper 1: W. An et al., J. Comput. Phys. (2013), 250, 165–177. (`Link <https://doi.org/10.1016/j.jcp.2013.05.020>`__)
- QuickPIC original paper 2: C. Huang et al., J. Comput. Phys. (2006), 217(2), 658–679. (`Link <https://doi.org/10.1016/j.jcp.2006.01.039>`__)
- Pipelining algorithm: B. Feng et al., J. Comput. Phys. (2009), 228(15), 5340–5348. (`Link <https://doi.org/10.1016/J.JCP.2009.04.019>`__)

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
