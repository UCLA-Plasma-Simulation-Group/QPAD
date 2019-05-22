program test_timing

use parallel_pipe_class
use sys
use mpi

implicit none

type( parallel_pipe ), pointer :: pp => null()

allocate(pp)
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call init_tprof( enable_timing=.true. )
call write_dbg( 'main', 'test_timing', 0, 'starts' )

call start_tprof( 'field solve' )
call sleep(3) ! pretend to be solving the fields
call stop_tprof( 'field solve' )

call start_tprof( 'particle push' )
call sleep(2) ! pretend to be pushing particles
call stop_tprof( 'particle push' )

call start_tprof( 'arithmetics' )
call sleep(1) ! pretend to be performing arithmetics
call stop_tprof( 'arithmetics' )

call write_tprof()

call write_dbg( 'main', 'test_timing', 0, 'ends' )

call end_errors()
call pp%del()

end program test_timing