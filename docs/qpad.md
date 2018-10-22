src_dir: ../
output_dir: ./website
project: QuickPIC wit Azimuzal Decomposition
title: QPAD
summary: ![QuickPIC](|media|/quickpic_logo.png)
         {: style="text-align: center" }
author: Fei Li and Weiming An
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: by-nc
extra_filetypes: sh #
include: /usr/local/openmpi/include

QuickPIC is a 3D parallel (MPI & OpenMP Hybrid) Quasi-Static PIC code, which is developed based on the framework UPIC. This is the UCLA Plasma Simulation Group's official open-source repository for QuickPIC.