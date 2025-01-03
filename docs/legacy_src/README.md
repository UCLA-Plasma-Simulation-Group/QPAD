## Documentation

* [**Installation**](./Install-QPAD.md)
* **Input File**
  * [Overview](./Input-File-Overview.md)
  * [Session: Simulation](./Session-Simulation.md)
  * [Session: Beam](./Session-Beam.md)
  * [Session: Species](./Session-Species.md)
  * [Session: Neutrals](./Session-Neutrals.md)
  * [Session: Field](./Session-Field.md)
  * [Post-processing](./Post-processing.md)
* [**Developer guidelines**](./Developer-Guide.md)

## Change log

**Update on 06/18/2021**
- Added new parameters "alpha", "perp_offset_x" and "perp_offset_y" in the "beam" session. See the manual for details.

**Update on 06/17/2021**
- Rewrote the "beam" session. See the manual for details.

**Update on 04/24/2021**
- Added a new parameter "npmax" into "species" and "neutral" sections. See the manual for details.

**Update on 04/22/2021**
- Added the switch "neutralized" in "species" sessions. This parameter is used to control whether to generate the neutralized background.

**Update on 02/22/2021**
- Added a new push type "clamp" in "species" and "neutrals" sessions.
- Change the parameters "uth" and "den_min" in the "species" and "neutrals" sessions to be optional.

**Update on 12/17/2020**
- Changed usage of input file sessions "species" and "neutrals". The previous users should see the wiki pages for the description.

**Update on 12/12/2020**
- Added input file parameter descriptions for "neutrals" session.
- Added a new parameter "nneutrals" into the "simulation" session.
- Modified the definition of "ppc" in the "species" session.
- Added a new parameter "num_theta" into the "species" session.

**Update on 11/16/2020**
- Add a new entry "random_seed" into the "simulation" session of the input file.

**Update on 10/21/2020**
- Add a new entry "algorithm" into the "simulation" session of the input file.
- Modified the instruction of creating users' configuration file.

**Update on 08/23/2020**
- Modified the guidance of installing JSON-FORTRAN library. Now the latest version 8.2.0 is recommended.

**Update on 04/26/2020**
- Add diagnostics for the plasma-induced vector potential.
- Add external beam particle importing from OSIRIS.
- Add spin dynamics (currently only available for external particle import).

**Update on 03/02/2020**
- Add diagnostics for the plasma-induced fields.

**Update on 02/28/2020**
- Fixed some errors in the beam parameter definitions.

