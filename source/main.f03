program quickpic_quasi3d

use simulation_class
use sys

implicit none

type( simulation ) :: sim

call sim%new()
call sim%run()
call sim%del()

end program quickpic_quasi3d
