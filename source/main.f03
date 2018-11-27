program quickpic

use simulation_class
use system

implicit none

type( simulation ) :: sim

call sim%new()
call sim%run()
call sim%end()

end program quickpic
