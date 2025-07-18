! hdf5io module for QuickPIC Open Source 1.0
! update: 04/18/2016

      module hdf5io_class

!! Quick QuickPIC to QPAD conversion note (08)
!! 
!! replaces parallel_pipe_class and perrors_class with parallel module, and 
!! the pp variable is no longer necessary (as parallel is no longer a data type)
!!
!! apply the following substitutions:
!!
!!    pp%getlgrp() --> comm_loc()
!!    pp%getlnvp() --> num_procs_loc()
!!    pp%getkstrt() --> id_proc()
!!    pp%getmint() --> p_dtype_int
!!    pp%getlworld() --> comm_world()
!!    pp%gernvp() --> num_procs()
!!

      use parallel_module
      use sysutil_module
      use HDF5
      use mpi
      use param
         
      implicit none
      
      private
      
      public :: hdf5file, pwfield, pwfield_pipe, wfield_pipe, pwpart_pipe, pwpart,&
      &wpart,rpart, prfield
      
      public :: detect_precision
      
      type hdf5file
         
         private
         
         character(len=100) :: filename = 'file.h5', timeunit = 'a.u.',&
         &dataname = 'Data', units = 'a.u.', label = 'Data'
         character(len=100) :: ty = 'grid',openpmd = '1.1.0', &
         & basepath='/data/%T/',filenamebase = 'data', records='Ex',&
         & dataorder = 'F', software='QPAD',base = 'data',&
         & meshespath = 'fields/', particlespath = 'particles/', &
         & iterationencoding = 'fileBased', iterationformat = 'filename_%T.h5',&
         & filepath ='./', spectype = 'electron'

         integer :: n = 1, rank = 2
         real :: t = 1.0, dt = 1.0
         character(len=100), dimension(3) :: axisname  = (/'x1','x2','x3'/),&
         &axislabel = (/'x1','x2','x3'/), axisunits = (/'a.u.','a.u.','a.u.'/)

         character(len=100), dimension(3) :: partlabel = (/'x','y','z'/)

         real, dimension(3) :: axismax = (/1.0,1.0,1.0/), axismin = (/0.0,0.0,0.0/)
         real :: n0 = 1e26 ! normalizing density in m^{-3}
         real :: wp = 5.6414602e13 ! plasma freq in s^{-1}
         real :: kp = 1.881788e5 ! plasma skin depth in m^{-1}
         ! ------------------------------------------------------------------
         ! ------------------------------------------------------------------
         ! ------------------   OpenPMD attributes   ------------------------
         ! ------------------------------------------------------------------
         ! -----------------------------------------------------------------
         integer :: openpmdextension = 1
         integer :: iter = 0
         character(len=8) :: chiter = '00000000'
         real :: unitsi = 1.0, timeunitsi = 1.0, gridunitsi=1.0
         real :: timeoffset = 0.0
         
         contains
         
         generic :: new => init_hdf5file
         procedure, private :: init_hdf5file
         
      end type hdf5file

      interface add_h5_atribute
         module procedure add_h5_atribute_str
         module procedure add_h5_atribute_str_v1
         module procedure add_h5_atribute_single
         module procedure add_h5_atribute_v1_single
         module procedure add_h5_atribute_int
         module procedure add_h5_atribute_v1_int
      end interface

      interface pwfield
        module procedure pwfield_3d
        module procedure pwfield_2d
      end interface

      interface prfield
        module procedure prfield_2d
      end interface

      interface pwfield_pipe
        module procedure pwfield_3d_pipe
        module procedure pwfield_2d_pipe
      end interface
      
      interface wfield_pipe
        module procedure wfield_2d_pipe
      end interface      

!
! new QPAD pwpart_3d uses x,p,q and maybe spin, in QuickPIC the particle array is stored in a single array called part
      interface pwpart_pipe
        module procedure pwpart_3d_pipe
      end interface
!

      interface pwpart
        module procedure pwpart_2d
        module procedure pwpart_2d_r
      end interface
      
      contains

      subroutine init_hdf5file(this,filename,timeunit,ty,n,t,dt,axisname,&
      & axislabel,axisunits,axismax,axismin,dataname,units,label,rank, &
      & openpmd,openpmdextension,iter,base,records,basepath,meshespath, &
      & particlespath,filenamebase,timeunitsi,iterationencoding, &
      & iterationformat,software, n0)

         implicit none
         
         class(hdf5file), intent(inout) :: this
         character(len=*), intent(in), optional :: filename, timeunit
         character(len=*), intent(in), optional :: ty, dataname, units, label
         integer, intent(in), optional :: n
         real, intent(in), optional :: n0
         integer, optional :: rank
         real, intent(in), optional :: t, dt
         character(len=*), dimension(3), intent(in), optional :: axisname, &
         &axislabel, axisunits
         real, dimension(3), intent(in), optional :: axismax, axismin

         character(len=*), intent(in), optional :: openpmd
         character(len=*), intent(in), optional :: openpmdextension
         integer, intent(in), optional :: iter
         character(len=*), intent(in), optional :: base
         character(len=*), intent(in), optional :: records
         character(len=*), intent(in), optional :: basepath
         character(len=*), intent(in), optional :: meshespath
         character(len=*), intent(in), optional :: particlespath
         character(len=*), intent(in), optional :: filenamebase
         real, intent(in), optional :: timeunitsi
         character(len=*), intent(in), optional :: iterationencoding
         character(len=*), intent(in), optional :: iterationformat
         character(len=*), intent(in), optional :: software




         character(len=8) :: chiter
         integer i

          
         ! DEBUG
         ! write(*,*) 'in init_hdf5file ',trim(this%filename),'  ',trim(this%timeunit)
         ! write(*,*) 'this%records = ', trim(this%records)
         ! write(*,*) 'this%filenamebase = ', trim(this%filenamebase)
         ! write(*,*) 'this%base = ', trim(this%base)
         ! DEBUG

         if(.not. present(rank)) this%rank=2
         
         if (present(filename)) then
            this%filename = filename
         end if

         if (present(dataname)) then
            this%dataname = dataname
         end if

         if (present(units)) then
            this%units = units
         end if

         if (present(label)) then
            this%label = label
         end if

         if (present(timeunit)) then
            this%timeunit = timeunit
         end if

         if (present(ty)) then
            this%ty = ty
         end if
         
         if (present(n)) then
            this%n = n
         end if

         if (present(rank)) then
            this%rank = rank
         end if

         if (present(t)) then
            this%t = t
         end if

         if (present(dt)) then
            this%dt = dt
         end if

         if (present(axisname)) then
             this%axisname = axisname
         end if

         if (present(axislabel)) then
            this%axislabel = axislabel
         end if


         if (present(axisunits)) then
            this%axisunits = axisunits
         end if

         if (present(axismax)) then
            this%axismax = axismax
         end if

         if (present(axismin)) then
            this%axismin = axismin
         end if

         if (present(openpmd)) then    
             this%openpmd = trim(openpmd)
         end if

         if (present(base)) then
             this%base = trim(base)
         end if
         if (present(iter)) then    
             this%iter = iter
         end if

         !write(chiter,'(I8.8)') this%iter
         write(chiter,'(I8.1)') this%iter
         chiter = trim(chiter)
         this%chiter = chiter
         
         


         if (present(records)) then
            this%records = trim(records)
         end if

         if (present(software)) then
             this%software = trim(software)
         end if

         

         if(present(n0)) then
            this%n0 = n0 * 1e6
            this%wp = sqrt(q_e**2 * this%n0/(ep0 * m_e))
            this%kp = this%wp/c
            this%timeunitsi = 1.0/this%wp
            this%gridunitsi = 1.0/this%kp
         end if

         ! if (present(timeunitsi)) then
         !    this%timeunitsi = timeunitsi
         ! else
         !    this%timeunitsi = 1.0
         ! end if

         ! DEBUG
         ! write(*,*) 'out iniit_hdf5file'
         ! write(*,*) 'filename = ', trim(this%filename)
         ! write(*,*) 'ty = ', trim(this%ty)
         ! write(*,*) 'openpmd = ', trim(this%openpmd)
         ! write(*,*) 'meshespath = ', trim(this%meshespath)
         ! write(*,*) 'software = ', trim(this%software)
         ! DEBUG

      end subroutine init_hdf5file
!
      
      subroutine add_h5_atribute_str( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         character( len = * ), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, typeID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer(size_t) :: maxlen
          
         integer :: ierr
          
         dims(1) = 1
         call h5screate_simple_f(0, dims, dataspaceID, ierr ) 
         ! call h5tcopy_f(H5T_FORTRAN_S1, typeID, ierr)
         call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
          
         maxlen = len(attribute)
         call h5tset_size_f(typeID, maxlen, ierr)
          
         call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, typeID, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5tclose_f( typeID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_str        
!
      subroutine add_h5_atribute_str_v1( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         character( len = * ), dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, typeID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer(size_t) :: maxlen
         integer :: i, ierr
          
         dims(1) = size(attribute)
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
          
         maxlen = 0
         do i = 1, size(attribute)-1
            if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
         enddo
          
         call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
         ! call h5tcopy_f(H5T_FORTRAN_S1, typeID, ierr)
         call h5tset_size_f(typeID, maxlen, ierr)
         
         call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, typeID, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5tclose_f( typeID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_str_v1
!
      subroutine add_h5_atribute_single( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         real, intent(in) :: attribute
         
         integer(hid_t) :: dataspaceID, attrID
         integer(hid_t) :: d_float 
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         d_float = detect_precision()
         dims(1) = 1
         call h5screate_simple_f(0, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, d_float, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, d_float, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_single
!
      subroutine add_h5_atribute_v1_single( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         real, dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hid_t) :: d_float
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         d_float = detect_precision()
         dims(1) = size(attribute)
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, d_float, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, d_float, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_v1_single
!
      subroutine add_h5_atribute_int( objID, name, attribute )

         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         integer, intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         dims(1) = 1
         call h5screate_simple_f(0, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID,&
         &attrID, ierr )
         call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_int
!
      subroutine add_h5_atribute_v1_int( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         integer, dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         dims(1) = size(attribute)
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID,&
         &attrID, ierr )
         call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_v1_int
!
! In case real is defaulted to 4 bytes and real and double become 2 
! distinct types these routines are needed
!
!       subroutine add_h5_atribute_double(objID, name, attribute)

!          implicit none
          
!          integer(hid_t), intent(in) :: objID
!          character(len=*), intent(in) :: name
!          double precision, intent(in) :: attribute
          
!          integer(hid_t) :: dataspaceID, attrID
!          integer(hsize_t), dimension(1) :: dims
!          integer :: ierr
          
!          dims(1) = 1
!          call h5screate_simple_f(0, dims, dataspaceID, ierr) 
!          call h5acreate_f(objID, name, H5T_NATIVE_DOUBLE, dataspaceID,&
!          &attrID, ierr)
!          call h5awrite_f(attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)
!          call h5aclose_f(attrID, ierr)
!          call h5sclose_f(dataspaceID, ierr)
        
!       end subroutine add_h5_atribute_double
! !
!       subroutine add_h5_atribute_v1_double(objID, name, attribute)
        
!          implicit none
          
!          integer(hid_t), intent(in) :: objID
!          character(len=*), intent(in) :: name
!          double precision, dimension(:), intent(in) :: attribute
          
!          integer(hid_t) :: dataspaceID, attrID
!          integer(hsize_t), dimension(1) :: dims
!          integer :: ierr
          
!          dims(1) = size(attribute)
!          call h5screate_simple_f(1, dims, dataspaceID, ierr) 
!          call h5acreate_f(objID, name, H5T_NATIVE_DOUBLE, dataspaceID,&
!          &attrID, ierr)
!          call h5awrite_f(attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)
!          call h5aclose_f(attrID, ierr)
!          call h5sclose_f(dataspaceID, ierr)
        
!       end subroutine add_h5_atribute_v1_double
! !
      subroutine wrattr_file(this,file_id,xferID, type)
      
         implicit none
         
         class(hdf5file), intent(in) :: this
         integer(hid_t), intent(in) :: file_id
! local data
          integer(hid_t) :: rootID, aid, dspace_id, dset_id, treal, xferID
          integer(hsize_t), dimension(1) :: dims
          integer, parameter :: zero = ichar('0')          
          integer :: i, ierr
          logical, intent(in), optional :: type
          logical :: mesh

          treal = detect_precision()
          mesh = .true.
           if (present(type)) then
               mesh = type
           endif

          call h5gopen_f(file_id, '/', rootID, ierr)

          call add_h5_atribute(rootID, 'NAME', trim(this%dataname) )
          call add_h5_atribute(rootID, 'TYPE', this%ty) 
          call add_h5_atribute(rootID, 'TIME', this%t)
          call add_h5_atribute(rootID, 'ITER', this%n) 
          call add_h5_atribute(rootID, 'DT', this%dt)
          call add_h5_atribute(rootID, 'TIME UNITS', trim(this%timeunit) )
! OpenPMD attributes

          
          call add_h5_atribute(rootID, 'iterationEncoding',trim(this%iterationencoding) )
          call add_h5_atribute(rootID, 'basePath', trim(this%basepath))  

          if(mesh) then
            call add_h5_atribute(rootID, 'meshesPath', trim(this%meshespath))
            call add_h5_atribute(rootID, 'openPMDextension',this%openpmdextension)
            call add_h5_atribute(rootID, 'openPMD', trim(this%openpmd) )
          else
            ! particles dumped using OpenPMD Beam Physics Ext
            call add_h5_atribute(rootID, 'openPMD', '2.0.0' )
            call add_h5_atribute(rootID, 'particlesPath',trim(this%particlespath))
            call add_h5_atribute(rootID, 'speciesType',trim(this%spectype))
            call add_h5_atribute(rootID, 'openPMDextension','BeamPhysics;SpeciesType')
         endif
          ! call add_h5_atribute(rootID, 'particlesPath',trim(this%particlespath)
          ! write(*,*) 'DATANAME: ',this%dataname
          ! write(*,*) 'meshesPATH: ',this%meshespath
          ! write(*,*) 'openPMD: ', this%openpmd
          call add_h5_atribute(rootID, 'software',trim(this%software))
          call add_h5_atribute(rootID,'iterationFormat','data%T.h5')
! OpenPMD attributes     
!         if (this%ty == 'grid') then
!            call h5gcreate_f(rootID, 'AXIS', aid, ierr) 
!
!            dims(1) = 2 
!            call h5screate_simple_f(1, dims, dspace_id, ierr ) 
!            do i = 1, this%rank
!               call h5dcreate_f(aid, 'AXIS'//char(zero+i), treal, dspace_id,&
!               &dset_id, ierr )
!
!               call add_h5_atribute(dset_id, 'TYPE', 'linear') 
!               call add_h5_atribute(dset_id, 'UNITS', trim(this%axisunits(i))) 
!               call add_h5_atribute(dset_id, 'NAME', trim(this%axisname(i))) 
!               call add_h5_atribute(dset_id, 'LONG_NAME', trim(this%axislabel(i))) 
!
!               call h5dwrite_f(dset_id, treal, (/this%axismin(i),this%axismax(i)/),&
!               &dims, ierr, xfer_prp=xferID)
!
!               call h5dclose_f(dset_id, ierr)
!            enddo
!   
!            call h5sclose_f(dspace_id, ierr) 
!            call h5gclose_f(aid, ierr) 
!         endif
             call h5gclose_f(rootID, ierr)

      end subroutine wrattr_file
!
!
      subroutine wrattr_dataset(this,dset_id,unit,name, type)
      
         implicit none
         
         class(hdf5file), intent(in) :: this
         integer(hid_t), intent(in) :: dset_id
         character(len=*), intent(in), optional :: unit,name

         logical, intent(in), optional :: type
         logical :: mesh
         ! real(8) :: gridtoSI = 1.0
         

         ! if (present(unit)) then
            ! call add_h5_atribute(dset_id, 'UNITS', unit)
         ! else
            ! call add_h5_atribute(dset_id, 'UNITS', this%units)
         ! endif 
          
         ! if (present(name)) then
            ! call add_h5_atribute(dset_id, 'LONG_NAME', name) 
         ! else
            ! call add_h5_atribute(dset_id, 'LONG_NAME', this%label) 
         ! endif 
! OPENPMD

         mesh = .true.
         if(present(type)) then
            mesh = type
         endif
         if(mesh) then
            call add_h5_atribute(dset_id,'dataOrder','F')
            call add_h5_atribute(dset_id,'geometry','cartesian')
            call add_h5_atribute(dset_id,'geometryParameters','cartesian')
            call add_h5_atribute(dset_id,'gridUnitSI',this%gridunitsi)
            call add_h5_atribute(dset_id,'axisLabels',this%axisname(this%rank:1:-1))
            call add_h5_atribute(dset_id,'unitSI',1.0)
         endif

         call add_h5_atribute(dset_id,'timeOffset',0.0)
         
         call add_h5_atribute(dset_id,'time',this%t)
         ! call add_h5_atribute(dset_id,'axisLabels',this%axisname(1:this%rank))
         
         
         if(present(name)) then
            if(name(1:1) == 'x') then
               call add_h5_atribute(dset_id,'unitDimension',[ 1., 0., 0., 0., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
               call add_h5_atribute(dset_id,'unitSI',this%gridunitsi)
               call add_h5_atribute(dset_id,'unitSymbol','m')
            else if (name(1:1) == 'p') then
               call add_h5_atribute(dset_id,'unitDimension',[ 1., 1., -1., 0., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
               call add_h5_atribute(dset_id,'unitSI',5.344285992678308E-28)
               call add_h5_atribute(dset_id,'unitSymbol','eV/c')
            else if(name(1:1) == 'q') then
               call add_h5_atribute(dset_id,'unitDimension',[ 0., 0., 1., 1., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
               call add_h5_atribute(dset_id,'unitSI',1.0)
               call add_h5_atribute(dset_id,'unitSymbol','C')
            else if(name(1:1) == 't') then
               call add_h5_atribute(dset_id,'unitDimension',[ 0., 0., 1., 0., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
               call add_h5_atribute(dset_id,'unitSymbol','s')
               call add_h5_atribute(dset_id,'unitSI',1.0)
            else if(name(:) == 'weight') then
               call add_h5_atribute(dset_id,'unitDimension',[ 0., 0., 0., 0., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
            else if(name(:) == 'mass') then
               call add_h5_atribute(dset_id,'unitDimension',[ 0., 1., 0., 0., 0., 0., 0.])
               call add_h5_atribute(dset_id,'weightingPower',0)
               call add_h5_atribute(dset_id,'macroWeighted',0)
               call add_h5_atribute(dset_id,'unitSymbol','kg')
            else
               call add_h5_atribute(dset_id,'unitDimension',[ 1., 0., 1., 1., 0., 0., 0.])
            endif
         else
            call add_h5_atribute(dset_id,'unitDimension',[ 1., 0., 1., 1., 0., 0., 0.])
         endif
         
! OPENPMD
  
      end subroutine wrattr_dataset
!
!
      subroutine pwfield_3d(file,fd,gs,ls,noff,ierr)

         implicit none
        
         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:,:), intent(in) :: fd
         integer, dimension(3), intent(in) :: gs, ls
         integer, dimension(2), intent(in) :: noff
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id

         integer(hsize_t), dimension(3) :: start
         integer(hsize_t), dimension(3) :: gsize, lsize
         integer(hsize_t), dimension(2) :: lnoff
         integer :: info
         character(len=:), allocatable :: filename
         character(len=8) :: st
! OPENPMD hierachy
         integer(hid_t) :: iterID, meshID, dataID
! OPENPMD local  variables
         real, dimension(3) :: local_gridspacing
         real, dimension(3) :: local_GlobalOffset
         real, dimension(3) :: local_position
         integer i
         character(len=20) :: iter_str

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
         ! write(*,*) 'in pwfield_3D',filename

         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff



         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
      
         call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,&
         &access_prp=flplID) 

         call wrattr_file(file,file_id,xferID)

         call h5screate_simple_f(3, gsize, dspace_id, ierr)
         call h5screate_simple_f(3, lsize, memspaceID, ierr )
         ! In OSIRIS convention, the group ID is simply "/"
         ! OSIRIS
         call h5gopen_f(file_id, '/', rootID, ierr)
         ! OSIRIS
         ! In OpenPMD, the group ID is more complicated because
         ! the file is designed to hold the entire simulation, so
         ! you have to specify whether it is "data" or any word
         ! specified by the keyword
         ! "basePath" that indicates this is the beginning of the data
         !
         ! <ROOT> 
         !
         ! underneath "basePath" contains the timesteps
         ! at each timestep the OpenPMD file can contain
         ! particle or mesh data,
         ! (in OpenPMD) there is a keyword for the beginning 
         ! of the file, and also a keyword for mesh data
         ! (we have always used the keyword "mesh")
         ! and also a keyword for particle data
         ! (we have always used the keyword "particle")
         !
         ! OPENPMD
         call h5gcreate_f(rootID, 'data' ,dataID , ierr)
         if(ierr .ne. 0) then
             write(*,*) 'error open data'
         end if
         ! write(iter_str,'(I0.8)') file%n
         write(iter_str,'(I8)') file%n
         call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
         if(ierr .ne. 0) then
             write(*,*) 'error open iter'
         end if
         ! here we need to write dt, time, and also
         ! time to SI connversion 
         call add_h5_atribute(iterID,'dt',file%dt)
         call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
         call h5gcreate_f(iterID,trim(file%meshespath),meshID,ierr)
         ! Finally we create the mesh data
         !
         call h5dcreate_f(meshID, file%dataname, treal, dspace_id, dset_id,&
         &ierr, dcplID)

         start(1) = 0
         start(2) = lnoff(1)
         start(3) = lnoff(2)
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2),1:lsize(3)),&
         &lsize, ierr, memspaceID, dspace_id, xfer_prp=xferID)

         call wrattr_dataset(file,dset_id)

         do i=1,3 
             local_gridspacing(i) = (file%axismax(i)-file%axismin(i))/gsize(i)
             local_GlobalOffset(i) = real(i)
             local_position(i) = 0.0
         end do
         
         call add_h5_atribute(dset_id,'gridSpacing',local_gridspacing)
         call add_h5_atribute(dset_id,'gridGlobalOffset',local_GlobalOffset)
         call add_h5_atribute(dset_id,'position',local_position)

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         call h5gclose_f(meshID, ierr)
         call h5gclose_f(iterID, ierr)
         call h5gclose_f(dataID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_3d
!
      subroutine pwfield_2d(file,fd,gs,ls,noff,ierr)

         implicit none
        
         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: fd
         integer, dimension(2), intent(in) :: gs, ls
         integer, intent(in) :: noff
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id
         integer(hsize_t), dimension(2) :: start
         integer(hsize_t), dimension(2) :: gsize, lsize
         integer(hsize_t) :: lnoff
         integer :: info
         character(len=:), allocatable :: filename
         character(len=8) :: st
! OPENPMD hierachy
         integer(hid_t) :: iterID, meshID, dataID
! other OPENPMD variables
         real, dimension(2) :: local_gridspacing
         real, dimension(2) :: local_GlobalOffset
         real, dimension(2) :: local_position
         integer i
         character(len=20) :: iter_str

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
         ! write(*,*) ' in pwrite_2d',filename
                  
         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff
         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
      
         call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,&
         &access_prp=flplID) 
         call wrattr_file(file,file_id,xferID)

         call h5screate_simple_f(2, gsize, dspace_id, ierr)
         call h5screate_simple_f(2, lsize, memspaceID, ierr )
         call h5gopen_f(file_id, '/', rootID, ierr)
         ! write(*,*) '(pwfield_2d) open data'
         call h5gcreate_f(rootID,'data',dataID,ierr)
         ! write(iter_str,'(I0.8)') file%n
         write(iter_str,'(I8)') file%n
         ! write(*,*) '(pwfield_2d) open iter  --- ',iter_str
         call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
         ! here we need to write dt, time, and also
         ! time to SI connversion 
         call add_h5_atribute(iterID,'dt',file%dt)
         call add_h5_atribute(iterID,'time',file%t)
         call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
         write(*,*) 'open mesh'
         call h5gcreate_f(iterID,trim(file%meshesPath),meshID,ierr)
         ! Finally we create the mesh data
         !
         call h5dcreate_f(meshID, file%dataname, treal, dspace_id, dset_id,&
         &ierr, dcplID)

         start(1) = 0
         start(2) = lnoff
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)


         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2)),lsize, ierr,&
         &memspaceID, dspace_id, xfer_prp=xferID)

         call wrattr_dataset(file,dset_id)

         do i=1,2 
             ! local_gridspacing(i) = (file%axismax(3-i)-file%axismin(3-i))/gsize(3-i)
             ! local_GlobalOffset(i) = file%axismin(3-i)
             local_gridspacing(i) = (file%axismax(i)-file%axismin(i))/gsize(i)
             local_GlobalOffset(i) = file%axismin(i)
             local_gridspacing(i) = 0.0
         end do
         call add_h5_atribute(dset_id,'gridSpacing',local_gridspacing)
         call add_h5_atribute(dset_id,'gridGlobalOffset',local_GlobalOffset)
         call add_h5_atribute(dset_id,'position',local_position)

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5gclose_f(meshID, ierr)
         call h5gclose_f(iterID, ierr)
         call h5gclose_f(dataID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_2d
!
      subroutine pwfield_3d_pipe(file,fd,gs,ls,noff,rtag,&
      &stag,id,ierr)

         implicit none
        
         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:,:), intent(in) :: fd
         integer, dimension(3), intent(in) :: gs, ls
         integer, dimension(2), intent(in) :: noff
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id, ierr
         
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aid
         integer(hid_t) :: tstring
         integer(hsize_t), dimension(3) :: gsize, lsize
         integer(hsize_t), dimension(2) :: lnoff
         character(len=80) :: string
         integer(hsize_t) :: lstr
         integer(hsize_t), dimension(3) :: start
         integer :: ori, des, nvyp, stageid, mid, message, info
         integer, dimension(10) :: istat
         integer(hsize_t), dimension(1) :: dims 
         character(len=:), allocatable :: filename
         character(len=8) :: st
! OPENPMD hierachy
         integer(hid_t) :: dataID, iterID, meshID
! other OPENPMD variables
         real, dimension(3) :: local_gridspacing
         real, dimension(3) :: local_GlobalOffset
         real, dimension(3) :: local_position
         integer i
         character(len=20) :: iter_str
                  
         ! write(*,*) 'in pwfield_3D_pipe'
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
         
         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff
         nvyp = num_procs_loc()
         ori = id_proc() - nvyp 
         des = id_proc() + nvyp 
         dims = 1
                  
         if (ori >= 0) then
            call MPI_IRECV(message,1,p_dtype_int,ori,rtag,comm_world(),&
            &mid,ierr)
            call MPI_WAIT(mid,istat,ierr)
         endif
         if (ori < 0 ) then
             ! write(*,*) 'in pwfield_3D_pipe',filename
         endif
             
         
         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
         
         if (ori >= 0) then
            call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr,&
            &access_prp=flplID)
            call h5gopen_f(file_id,'/',rootID,ierr)
            ! write(iter_str,'(I0.8)') file%n
            write(iter_str,'(I8)') file%n
            call h5gopen_f(rootID,'data',dataID,ierr)
            call h5gopen_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
            call h5gopen_f(iterID,trim(file%meshespath),meshID,ierr)
            call h5dopen_f(meshID, trim(file%dataname), dset_id, ierr, H5P_DEFAULT_F)
            if (ierr .ne. 0) then
                write(*,*)'pwfield_3d_pipe: error opening dataset'
            endif 
            ! call h5aopen_by_name_f(file_id, "/", "NAME", aid, ierr)
            ! lstr = len(string)
            ! call h5tcopy_f(H5T_NATIVE_CHARACTER, tstring, ierr)
            ! call h5tset_size_f(tstring, lstr, ierr)
            ! call h5aread_f(aid, tstring, string, dims, ierr)  
            ! call h5dopen_f(file_id, string, dset_id, ierr, H5P_DEFAULT_F)
            ! call h5aclose_f(aid, ierr)
            ! call h5tclose_f(tstring, ierr)
         else
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,&
            &access_prp=flplID) 
            call wrattr_file(file,file_id,xferID)
            call h5gopen_f(file_id, '/', rootID, ierr)
            call h5gcreate_f(rootID,'data',dataID,ierr)
            ! write(iter_str,'(I0.8)') file%n
            write(iter_str,'(I8)') file%n
            call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
            ! here we need to write dt, time, and also
            ! time to SI connversion 
            call add_h5_atribute(iterID,'dt',(file%dt))
            call add_h5_atribute(iterID,'time',(file%t))
            call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
            call h5gcreate_f(iterID,file%meshesPath,meshID,ierr)
         endif
         
         call h5screate_simple_f(3, gsize, dspace_id, ierr)
         call h5screate_simple_f(3, lsize, memspaceID, ierr )
         
         if (ori < 0) then
            ! call h5dcreate_f(rootID, file%dataname, treal, dspace_id, dset_id,&
            ! &ierr, dcplID)

            ! Finally we create the mesh data
            !
            call h5dcreate_f(meshID, file%dataname, treal, dspace_id, dset_id,&
            &ierr, dcplID)
            do i=1,3 
                local_gridspacing(i) = (file%axismax(i)-file%axismin(i))/gsize(i)
                local_GlobalOffset(i) = 0.0
                local_position(i) = 0.0
            end do
         
         endif
         
         start(1) = 0
         start(2) = lnoff(1)
         start(3) = lnoff(2)
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2),1:lsize(3)),&
         &lsize, ierr, memspaceID, dspace_id, xfer_prp=xferID)
         if (ori < 0) then
             call wrattr_dataset(file,dset_id)
             call add_h5_atribute(dset_id,'gridSpacing',local_gridspacing)
             call add_h5_atribute(dset_id,'gridGlobalOffset',local_GlobalOffset)
             call add_h5_atribute(dset_id,'position',local_position)
         endif

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5gclose_f(meshID, ierr)
         call h5gclose_f(iterID, ierr)
         call h5gclose_f(dataID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
         if (ierr .ne. 0) write(*,*)'pwfield_3d_pipe: error in fclose()'

         if (des < num_procs()) then
            call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
            &id,ierr)
         else
            id = MPI_REQUEST_NULL         
         endif
                        
      end subroutine pwfield_3d_pipe
!
      subroutine pwfield_2d_pipe(file,fd,gs,ls,noff,rtag,&
      &stag,id,ierr)

         implicit none
        
         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: fd
         integer, dimension(2), intent(in) :: gs, ls
         integer, dimension(2), intent(in) :: noff
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id, ierr
         
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aid
         integer(hid_t) :: tstring
         integer(hsize_t), dimension(2) :: gsize, lsize
         integer(hsize_t), dimension(2) :: lnoff
         character(len=80) :: string
         integer(hsize_t) :: lstr
         integer(hsize_t), dimension(2) :: start
         integer :: ori, des, nvyp, stageid, mid, message, info
         integer, dimension(10) :: istat
         integer(hsize_t), dimension(1) :: dims 
         character(len=:), allocatable :: filename
         character(len=8) :: st
! OPENPMD hierachy
         integer(hid_t) :: dataID, iterID, meshID
! other OPENPMD variables
         real, dimension(2) :: local_gridspacing
         real, dimension(2) :: local_GlobalOffset
         real, dimension(2) :: local_position
         integer i
         character(len=20) :: iter_str
         logical gexist
                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
         
         ierr = 0
         ! gsize = gs
         ! lsize = ls
         ! lnoff = noff
         gsize(1) = gs(2)
         gsize(2) = gs(1)
         lsize(1) = ls(2)
         lsize(2) = ls(1)
         lnoff(1) = noff(2)
         lnoff(2) = noff(1)
         nvyp = num_procs_loc()
         ori = id_proc() - nvyp 
         des = id_proc() + nvyp 
         dims = 1
         ! print *,'ori',ori
         
         if (ori >= 0) then
            call MPI_IRECV(message,1,p_dtype_int,ori,rtag,comm_world(),&
            &mid,ierr)
            call MPI_WAIT(mid,istat,ierr)
         endif
         if (ori < 0 ) then
             ! write(*,*) ' in pwrite_2d_pipe',filename
         end if
         
         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)   

         if (ori >= 0) then ! slave nodes, here we open the file
            call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr,&
            &access_prp=flplID)
            call h5gopen_f(file_id,'/',rootID,ierr)
            ! write(iter_str,'(I0.8)') file%n
            write(iter_str,'(I8)') file%n
            call h5gopen_f(rootID,'data',dataID,ierr)
            call h5gopen_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
            call h5gopen_f(iterID,trim(file%meshespath),meshID,ierr)
            ! call h5aopen_by_name_f(file_id, "/", "NAME", aid, ierr)
            ! lstr = len(string)
            ! call h5tcopy_f(H5T_NATIVE_CHARACTER, tstring, ierr)
            ! call h5tset_size_f(tstring, lstr, ierr)
            ! call h5aread_f(aid, tstring, string, dims, ierr)  
            call h5dopen_f(meshID, trim(file%dataname), dset_id, ierr, H5P_DEFAULT_F)
            if (ierr .ne. 0) then
                write(*,*)'pwfield_2d_pipe: error opening dataset'
            endif 
            ! call h5aclose_f(aid, ierr)
            ! call h5tclose_f(tstring, ierr)
         else ! orig < 0 --> master node.  Here we write the attributes
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,&
            &access_prp=flplID) 
            call wrattr_file(file,file_id,xferID)
         endif 
         
         call h5screate_simple_f(2, gsize, dspace_id, ierr)
         call h5screate_simple_f(2, lsize, memspaceID, ierr )
         if(ori < 0) then
             call h5gopen_f(file_id, '/', rootID, ierr)
             call h5gcreate_f(rootID,'data', dataID, ierr)
             ! write(iter_str,'(I0.8)') file%n
             write(iter_str,'(I8)') file%n
             call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
             if (ierr .ne. 0) then
                 write(*,*) 'pwfield_2d_pipe: error in iterID ', ori
             end if
             ! here we need to write dt, time, and also
             ! time to SI connversion 
             call add_h5_atribute(iterID,'dt',file%dt)
             call add_h5_atribute(iterID,'time',file%t)
             call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
             call h5lexists_f(iterID, trim(file%meshespath),gexist, ierr)
            if (gexist) then
                call h5gopen_f(iterID,trim(file%meshespath),meshID,ierr)
            else
                call h5gcreate_f(iterID,trim(file%meshespath),meshID,ierr)
            end if
         endif
         
         if (ierr .ne. 0) then
             write(*,*) 'pwfield_2d_pipe: error in meshID ', ori
         end if
         if (ori < 0) then
            call h5dcreate_f(meshID, trim(file%dataname), treal, dspace_id, dset_id,&
            &ierr, dcplID)
            do i=1,2 
                local_gridspacing(i) = (file%axismax(i)-file%axismin(i))/gsize(3-i)
                local_GlobalOffset(i) = file%axismin(i)
                local_position(i) = 0.0
            end do
         
         endif
         start = lnoff
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, transpose(fd(1:lsize(2),1:lsize(1))),&
         &lsize, ierr, memspaceID, dspace_id, xfer_prp=xferID)
         
         if (ori < 0 ) then
            call wrattr_dataset(file,dset_id)
            call add_h5_atribute(dset_id,'gridSpacing',local_gridspacing)
            call add_h5_atribute(dset_id,'gridGlobalOffset',local_GlobalOffset)
            call add_h5_atribute(dset_id,'position',local_position)
         end if

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         
         call h5gclose_f(meshID, ierr)
         call h5gclose_f(iterID, ierr)
         call h5gclose_f(dataID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
         if (ierr .ne. 0) write(*,*)'pwfield_2d_pipe: error in fclose()'
         if (des < num_procs()) then
            call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
            &id,ierr)
         else
            id = MPI_REQUEST_NULL         
         endif
                        
      end subroutine pwfield_2d_pipe
!
      subroutine wfield_2d_pipe(file,fd,gs,ls,noff,rtag,&
      &stag,id,ierr)

         implicit none
        
         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: fd
         integer, dimension(2), intent(in) :: gs, ls
         integer, dimension(2), intent(in) :: noff
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id, ierr
         
! local data
         integer(hid_t) :: treal, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aid
         integer(hid_t) :: tstring
         integer(hsize_t), dimension(2) :: gsize, lsize
         integer(hsize_t), dimension(2) :: lnoff
         character(len=80) :: string
         integer(hsize_t) :: lstr
         integer(hsize_t), dimension(2) :: start
         integer :: ori, des, nvyp, stageid, mid, message, info
         integer, dimension(10) :: istat
         integer(hsize_t), dimension(1) :: dims 
         character(len=:), allocatable :: filename
         character(len=8) :: st
! OPENPMD hierachy
         integer(hid_t) :: dataID, iterID, meshID
! other OPENPMD variables
         real, dimension(2) :: local_gridspacing
         real, dimension(2) :: local_position
         real, dimension(2) :: local_GlobalOffset
         integer i
         character(len=20) :: iter_str
                  
                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
         write(*,*) ' in write_2d_pipe',filename
         
         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff
         nvyp = num_procs_loc()
         ori = id_proc() - nvyp 
         des = id_proc() + nvyp 
         dims = 1
         
         if (ori >= 0) then
            call MPI_IRECV(message,1,p_dtype_int,ori,rtag,comm_world(),&
            &mid,ierr)
            call MPI_WAIT(mid,istat,ierr)
         endif
         
         call h5open_f(ierr)
         treal = detect_precision()
         info = MPI_INFO_NULL
         
         if (ori >= 0) then
            call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr)
            call h5gopen_f(file_id,'/',rootID,ierr)
            ! write(iter_str,'(I0.8)') file%n
            write(iter_str,'(I8)') file%n
            call h5gopen_f(rootID,'data',dataID,ierr)
            call h5gopen_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
            call h5gopen_f(iterID,trim(file%meshespath),meshID,ierr)
            ! call h5aopen_by_name_f(file_id, "/", "NAME", aid, ierr)
            ! lstr = len(string)
            ! call h5tcopy_f(H5T_NATIVE_CHARACTER, tstring, ierr)
            ! call h5tset_size_f(tstring, lstr, ierr)
            ! call h5aread_f(aid, tstring, string, dims, ierr)  
            call h5dopen_f(meshID, trim(file%dataname), dset_id, ierr, H5P_DEFAULT_F)
            if (ierr .ne. 0) then
                write(*,*)'wfield_2d_pipe: error opening dataset'
            endif
            ! call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr)
            ! call h5aopen_by_name_f(file_id, "/", "NAME", aid, ierr)
            ! lstr = len(string)
            ! call h5tcopy_f(H5T_NATIVE_CHARACTER, tstring, ierr)
            ! call h5tset_size_f(tstring, lstr, ierr)
            ! call h5aread_f(aid, tstring, string, dims, ierr)  
            ! call h5dopen_f(file_id, string, dset_id, ierr, H5P_DEFAULT_F)
            ! call h5aclose_f(aid, ierr)
            ! call h5tclose_f(tstring, ierr)
         else
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)
            call wrattr_file(file,file_id,H5P_DEFAULT_F)
         endif
         
         call h5screate_simple_f(2, gsize, dspace_id, ierr)
         call h5screate_simple_f(2, lsize, memspaceID, ierr)
         if ( ori < 0 ) then
             call h5gopen_f(file_id, '/', rootID, ierr)
             call h5gcreate_f(rootID,'data',dataID,ierr)
             ! write(iter_str,'(I0.8)') file%n
             write(iter_str,'(I8)') file%n
             call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
         ! here we need to write dt, time, and also
         ! time to SI connversion 
             call add_h5_atribute(iterID,'dt',file%dt)
             call add_h5_atribute(iterID,'time',file%t)
             call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
         call h5gcreate_f(iterID,file%meshesPath,meshID,ierr)
         endif
         
         if (ori < 0) then
            call h5dcreate_f(meshID, file%dataname, treal, dspace_id, dset_id,&
            &ierr)
            call wrattr_dataset(file,dset_id)
            do i=1,2 
                local_gridspacing(i) = (file%axismax(i)-file%axismin(i))/gsize(i)
                local_GlobalOffset(i) = file%axismin(3-i)
                local_position(i) = 0.0
            end do
            call add_h5_atribute(dset_id,'gridSpacing',local_gridspacing)
            call add_h5_atribute(dset_id,'gridGlobalOffset',local_GlobalOffset)
            call add_h5_atribute(dset_id,'position',local_position)
         endif
         
         start(1) = 0
         start(2) = lnoff(2)
         
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2)),&
         &lsize, ierr, memspaceID, dspace_id)
         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5gclose_f(meshID, ierr)
         call h5gclose_f(iterID, ierr)
         call h5gclose_f(dataID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
         if (ierr .ne. 0) write(*,*)'wfield_2d_pipe: error in fclose()'

         if (des < num_procs()) then
            call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
            &id,ierr)
         else
            id = MPI_REQUEST_NULL         
         endif
                        
      end subroutine wfield_2d_pipe

      subroutine prfield_2d(file, fd, gs, ls, noff, ierr)

   implicit none
  
   class(hdf5file), intent(in) :: file
   real, intent(inout), dimension(:,:) :: fd
   integer, intent(in), dimension(2) :: gs, ls
   integer, intent(in), dimension(2) :: noff
   integer, intent(inout) :: ierr
  
   integer(hid_t) :: treal,fapl_id, xfer_id, mspace_id
   integer(hid_t) :: file_id, root_id, dset_id, fspace_id, dtype_id
   integer(hsize_t), dimension(2) :: gsize, lsize, offset
   character(len=:), allocatable :: filename
   character(len=8) :: st
  
   filename = trim(file%filename) // trim(file%dataname) // '_' // num2str(file%n, width=8) // '.h5'
   ierr = 0
   gsize = gs
   lsize = ls
   treal = detect_precision()

   call h5open_f(ierr)

   ! set property for file accessing and data transfer
   call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
   call h5pset_fapl_mpio_f(fapl_id, comm_world(), MPI_INFO_NULL, ierr)
   call h5pcreate_f(H5P_DATASET_XFER_F, xfer_id, ierr)
   call h5pset_dxpl_mpio_f(xfer_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  
   ! open hdf5 file
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr, fapl_id)
  
   ! create dataspaces for dataset in the memory
   call h5screate_simple_f(2, lsize, mspace_id, ierr)

   ! open root group and dataset
   call h5gopen_f(file_id, '/', root_id, ierr)
   call h5dopen_f(root_id, file%dataname, dset_id, ierr)
   call h5dget_space_f(dset_id, fspace_id, ierr)
   call h5dget_type_f(dset_id, dtype_id, ierr)
  
   ! read dataset from the file
   offset = noff
   call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, lsize, ierr)
   call h5dread_f(dset_id, dtype_id, fd(1:lsize(1), 1:lsize(2)), lsize, ierr, &
      mem_space_id=mspace_id, file_space_id=fspace_id, xfer_prp=xfer_id)

   ! close all the objects
   call h5sclose_f(mspace_id, ierr)
   call h5sclose_f(fspace_id, ierr)
   call h5tclose_f(dtype_id, ierr)
   call h5pclose_f(xfer_id, ierr)
   call h5pclose_f(fapl_id, ierr)
   call h5gclose_f(root_id, ierr)
   call h5dclose_f(dset_id, ierr)
   call h5fclose_f(file_id, ierr)
   call h5close_f(ierr)
  
end subroutine prfield_2d
!
      subroutine pwpart_2d_orig(file,part,npp,dspl,delta,ierr)
     
         implicit none

         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: part
         real, dimension(2), intent(in) :: delta
         integer, intent(in) :: npp,dspl
         integer, intent(inout) :: ierr
! local data
         integer :: tnpp, tp, color, pgrp, pid, pnvp, i, j
         integer(hsize_t), dimension(1) :: ldim
         integer, dimension(:), pointer :: np
         integer, dimension(:,:), pointer:: dims
         real, dimension(:), pointer :: buff
         integer(hsize_t), dimension(1) :: start,maxdim
         integer(hid_t) :: treal
         integer(hid_t) :: flplID, xferID, memspaceID, aid
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
         integer :: info
         integer, dimension(10) :: istat
         character(len=:), allocatable :: filename
         character(len=8) :: st

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
      
         ierr = 0
         ldim(1) = 1
         call h5open_f(ierr)
         treal = detect_precision()
         
         tnpp = int(npp/dspl)
         tp = 0
         call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,comm_loc(),ierr)

         if (tp == 0) then
            call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
            call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
            info = MPI_INFO_NULL
            call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
            call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
            &access_prp=flplID) 
            call wrattr_file(file,file_id,xferID)
            call h5gopen_f(file_id, '/', rootID, ierr)
            call h5screate_simple_f(1, ldim, aspace_id, ierr)
            call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
            &aid, ierr )
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
            call h5aclose_f(aid, ierr)
            call h5sclose_f(aspace_id, ierr)
            call h5pclose_f(xferID, ierr)
            call h5pclose_f(flplID, ierr)
            call h5gclose_f(rootID, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
            return
         else 
            if (tnpp > 0) then 
               color = 1
            else
               color = MPI_UNDEFINED
            endif
            call MPI_COMM_SPLIT(comm_loc(), color, 0, pgrp, ierr )

            if (tnpp > 0) then
               call MPI_COMM_RANK(pgrp, pid, ierr)
               call MPI_COMM_SIZE(pgrp, pnvp, ierr)
               allocate(np(pnvp), dims(2,pnvp), stat = ierr)
               call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
               &pgrp, ierr)
               dims(1, 1) = 1
               dims(2, 1) = np(1) 
               do i = 2, pnvp
                  dims(1,i) = dims(2,i-1) + 1
                  dims(2,i) = dims(1,i) + np(i) - 1
               enddo
               allocate(buff(tnpp), stat = ierr)

               call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
               call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
               info = MPI_INFO_NULL
               call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
               call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
               &access_prp=flplID) 
               call wrattr_file(file,file_id,xferID)
               call h5gopen_f(file_id, '/', rootID, ierr)
               call h5screate_simple_f(1, ldim, aspace_id, ierr)
               call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
               &aid, ierr )
               call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
               call h5aclose_f(aid, ierr)
               call h5sclose_f(aspace_id, ierr)

               do i = 1, 2
                  buff(1:tnpp) = part(i,1:(1+(tnpp-1)*dspl):dspl)*delta(i)
                  ldim(1) = tp
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), treal,&
                  &dspace_id, dset_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr)
                  start = dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                  &ldim, ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
                  call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                  &name='x_'//char(iachar('0')+i))
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               do i = 1, 3
                  buff(1:tnpp) = part((i+2),1:(1+(tnpp-1)*dspl):dspl)*delta(1) 
                  ldim(1) = tp
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), treal,&
                  &dspace_id, dset_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr)
                  start = dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                  &ldim, ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
                  call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                  &name='p_'//char(iachar('0')+i))
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               buff(1:tnpp) = part(8,1:(1+(tnpp-1)*dspl):dspl) 
               ldim(1) = tp
               call h5screate_simple_f(1, ldim, dspace_id, ierr)
               call h5dcreate_f(rootID, 'q', treal,&
               &dspace_id, dset_id, ierr)
               ldim(1) = tnpp
               call h5screate_simple_f(1, ldim, memspaceID, ierr)
               start = dims(1,pid+1) - 1
               call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
               &ldim, ierr)
               call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
               &dspace_id, xfer_prp=xferID)
               call wrattr_dataset(file,dset_id,unit='a.u.',&
               &name='q')
               call h5sclose_f(memspaceID, ierr)
               call h5sclose_f(dspace_id, ierr)
               call h5dclose_f(dset_id, ierr)


               call h5pclose_f(xferID, ierr)
               call h5pclose_f(flplID, ierr)
               call h5gclose_f(rootID, ierr)
               call h5fclose_f(file_id, ierr)
               deallocate(np,dims,buff)
            endif            
            if (pgrp /= MPI_COMM_NULL) then
               call MPI_COMM_FREE(pgrp, ierr)
            endif
            call h5close_f(ierr)
         endif
         
      end subroutine pwpart_2d_orig
!

subroutine pwpart_2d(file,x,p,q,npp,dspl,delta,ierr)

 implicit none

 class(hdf5file), intent(in) :: file
 real, dimension(:,:), intent(in) :: x, p
 real, dimension(:), intent(in) :: q
 real, dimension(2), intent(in) :: delta
 integer, intent(in) :: npp,dspl
 integer, intent(inout) :: ierr
! local data
 integer :: tnpp, tp, color, pgrp, pid, pnvp, i
 integer(hsize_t), dimension(1) :: ldim
 integer, dimension(:), pointer :: np => null()
 integer, dimension(:,:), pointer:: dims => null()
 real, dimension(:), pointer :: buff => null()
 integer(hsize_t), dimension(1) :: start
 integer(hid_t) :: treal
 integer(hid_t) :: flplID, xferID, memspaceID, aid
 integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
 integer :: info
 character(len=:), allocatable :: filename
 character(len=8) :: st


 allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
 write (st,'(I8.8)') file%n
 filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
 ! write(*,*) ' in pwpart_2d'
 ierr = 0
 ldim(1) = 1
 call h5open_f(ierr)
 treal = detect_precision()

 tnpp = int(npp/dspl)
 tp = 0
 call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,comm_loc(),ierr)

 if (tp == 0) then
    call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
    info = MPI_INFO_NULL
    call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
    call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
    &access_prp=flplID)
    call wrattr_file(file,file_id,xferID)
    call h5gopen_f(file_id, '/', rootID, ierr)
    call h5screate_simple_f(1, ldim, aspace_id, ierr)
    call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
    &aid, ierr )
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
    call h5aclose_f(aid, ierr)
    call h5sclose_f(aspace_id, ierr)
    call h5pclose_f(xferID, ierr)
    call h5pclose_f(flplID, ierr)
    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
    return
 else
    if (tnpp > 0) then
       color = 1
    else
       color = MPI_UNDEFINED
    endif
    call MPI_COMM_SPLIT(comm_loc(), color, 0, pgrp, ierr )

    if (tnpp > 0) then
       call MPI_COMM_RANK(pgrp, pid, ierr)
       call MPI_COMM_SIZE(pgrp, pnvp, ierr)
       allocate(np(pnvp), dims(2,pnvp), stat = ierr)
       call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
       &pgrp, ierr)
       dims(1, 1) = 1
       dims(2, 1) = np(1)
       do i = 2, pnvp
          dims(1,i) = dims(2,i-1) + 1
          dims(2,i) = dims(1,i) + np(i) - 1
       enddo
       allocate(buff(tnpp), stat = ierr)

       call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
       info = MPI_INFO_NULL
       call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
       call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
       &access_prp=flplID)
       call wrattr_file(file,file_id,xferID)
       call h5gopen_f(file_id, '/', rootID, ierr)
       call h5screate_simple_f(1, ldim, aspace_id, ierr)
       call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
       &aid, ierr )
       call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
       call h5aclose_f(aid, ierr)
       call h5sclose_f(aspace_id, ierr)

       do i = 1, 2
          buff(1:tnpp) = x(i,1:(1+(tnpp-1)*dspl):dspl)*delta(i)
          ldim(1) = tp
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), treal,&
          &dspace_id, dset_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr)
          start = dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
          &name='x_'//char(iachar('0')+i))
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       do i = 1, 3
          buff(1:tnpp) = p(i,1:(1+(tnpp-1)*dspl):dspl)
          ldim(1) = tp
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), treal,&
          &dspace_id, dset_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr)
          start = dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
          &name='p_'//char(iachar('0')+i))
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       buff(1:tnpp) = q(1:(1+(tnpp-1)*dspl):dspl)
       ldim(1) = tp
       call h5screate_simple_f(1, ldim, dspace_id, ierr)
       call h5dcreate_f(rootID, 'q', treal,&
       &dspace_id, dset_id, ierr)
       ldim(1) = tnpp
       call h5screate_simple_f(1, ldim, memspaceID, ierr)
       start = dims(1,pid+1) - 1
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
       &ldim, ierr)
       call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
       &dspace_id, xfer_prp=xferID)
       call wrattr_dataset(file,dset_id,unit='a.u.',&
       &name='q')
       call h5sclose_f(memspaceID, ierr)
       call h5sclose_f(dspace_id, ierr)
       call h5dclose_f(dset_id, ierr)


       call h5pclose_f(xferID, ierr)
       call h5pclose_f(flplID, ierr)
       call h5gclose_f(rootID, ierr)
       call h5fclose_f(file_id, ierr)
       deallocate(np,dims,buff)
    endif
    if (pgrp /= MPI_COMM_NULL) then
       call MPI_COMM_FREE(pgrp, ierr)
    endif
    call h5close_f(ierr)
 endif

end subroutine pwpart_2d

subroutine pwpart_2d_r(file,x,p,q,npp,dspl,ierr)

 implicit none

 class(hdf5file), intent(in) :: file
 real, dimension(:,:), intent(in) :: x, p
 real, dimension(:), intent(in) :: q
 integer(kind=LG), intent(in) :: npp
 integer, intent(in) :: dspl
 integer, intent(inout) :: ierr
! local data
 integer :: tnpp, tp, color, pgrp, pid, pnvp, i
 integer(hsize_t), dimension(1) :: ldim
 integer, dimension(:), pointer :: np => null()
 integer, dimension(:,:), pointer:: dims => null()
 real, dimension(:), pointer :: buff => null()
 integer(hsize_t), dimension(1) :: start
 integer(hid_t) :: treal
 integer(hid_t) :: flplID, xferID, memspaceID, aid
 integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
 integer :: info
 character(len=:), allocatable :: filename
 character(len=8) :: st


 allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
 write (st,'(I8.8)') file%n
 filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'
! write(*,*) ' in pwpart_2d_r'
 ierr = 0
 ldim(1) = 1
 call h5open_f(ierr)
 treal = detect_precision()

 tnpp = int(npp/dspl)
 tp = 0
 call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,comm_loc(),ierr)

 if (tp == 0) then
    call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
    info = MPI_INFO_NULL
    call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
    call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
    &access_prp=flplID)
    call wrattr_file(file,file_id,xferID)
    call h5gopen_f(file_id, '/', rootID, ierr)
    call h5screate_simple_f(1, ldim, aspace_id, ierr)
    call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
    &aid, ierr )
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
    call h5aclose_f(aid, ierr)
    call h5sclose_f(aspace_id, ierr)
    call h5pclose_f(xferID, ierr)
    call h5pclose_f(flplID, ierr)
    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
    return
 else
    if (tnpp > 0) then
       color = 1
    else
       color = MPI_UNDEFINED
    endif
    call MPI_COMM_SPLIT(comm_loc(), color, 0, pgrp, ierr )

    if (tnpp > 0) then
       call MPI_COMM_RANK(pgrp, pid, ierr)
       call MPI_COMM_SIZE(pgrp, pnvp, ierr)
       allocate(np(pnvp), dims(2,pnvp), stat = ierr)
       call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
       &pgrp, ierr)
       dims(1, 1) = 1
       dims(2, 1) = np(1)
       do i = 2, pnvp
          dims(1,i) = dims(2,i-1) + 1
          dims(2,i) = dims(1,i) + np(i) - 1
       enddo
       allocate(buff(tnpp), stat = ierr)

       call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
       info = MPI_INFO_NULL
       call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
       call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
       &access_prp=flplID)
       call wrattr_file(file,file_id,xferID)
       call h5gopen_f(file_id, '/', rootID, ierr)
       call h5screate_simple_f(1, ldim, aspace_id, ierr)
       call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
       &aid, ierr )
       call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
       call h5aclose_f(aid, ierr)
       call h5sclose_f(aspace_id, ierr)

       do i = 1, 2
          buff(1:tnpp) = x(i,1:((tnpp-1)*dspl+1):dspl)

          ! if (i == 1) then
          !    buff(1:tnpp) = part(1,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else if (i == 2) then
          !    buff(1:tnpp) = part(1,1:((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! end if
          ldim(1) = tp
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), treal,&
          &dspace_id, dset_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr)
          start = dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
          &name='x_'//char(iachar('0')+i), type = .false.)
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       do i = 1, 3
          buff(1:tnpp) = p(i,1:((tnpp-1)*dspl+1):dspl)

          ! if (i == 1) then
          !    buff(1:tnpp) = part(3,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl))-part(4,1:&
          !    &((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else if (i == 2) then
          !    buff(1:tnpp) = part(4,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl))+part(3,1:&
          !    &((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else
          !    buff(1:tnpp) = part(5,1:(1+(tnpp-1)*dspl):dspl)
          ! end if
          ldim(1) = tp
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), treal,&
          &dspace_id, dset_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr)
          start = dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
          &name='p_'//char(iachar('0')+i),type = .false.)
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       buff(1:tnpp) = q(1:(1+(tnpp-1)*dspl):dspl)
       ldim(1) = tp
       call h5screate_simple_f(1, ldim, dspace_id, ierr)
       call h5dcreate_f(rootID, 'q', treal,&
       &dspace_id, dset_id, ierr)
       ldim(1) = tnpp
       call h5screate_simple_f(1, ldim, memspaceID, ierr)
       start = dims(1,pid+1) - 1
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
       &ldim, ierr)
       call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
       &dspace_id, xfer_prp=xferID)
       call wrattr_dataset(file,dset_id,unit='a.u.',&
       &name='q', type = .false.)
       call h5sclose_f(memspaceID, ierr)
       call h5sclose_f(dspace_id, ierr)
       call h5dclose_f(dset_id, ierr)


       call h5pclose_f(xferID, ierr)
       call h5pclose_f(flplID, ierr)
       call h5gclose_f(rootID, ierr)
       call h5fclose_f(file_id, ierr)
       deallocate(np,dims,buff)
    endif
    if (pgrp /= MPI_COMM_NULL) then
       call MPI_COMM_FREE(pgrp, ierr)
    endif
    call h5close_f(ierr)
 endif

end subroutine pwpart_2d_r
!
      subroutine pwpart_3d_pipe_orig(file,part,npp,dspl,delta,origin,rtag,stag,&
      &id,ierr)
     
         implicit none

         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: part
         real, dimension(3), intent(in) :: delta,origin
         integer, intent(in) :: npp,dspl
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id, ierr
! local data
         integer :: tnpp, tp, tpo, color, pgrp, pid, pnvp, i, j
         integer(hsize_t), dimension(1) :: ldim
         integer, dimension(:), pointer :: np
         integer, dimension(:,:), pointer:: dims
         real, dimension(:), pointer :: buff
         integer(hsize_t), dimension(1) :: start,maxdim
         integer(hid_t) :: treal
         integer(hid_t) :: flplID, xferID, dcplID, memspaceID, aid
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
         integer :: ori, des, nvyp, stageid, mid, message, info
         integer, dimension(10) :: istat
         character(len=:), allocatable :: filename
         character(len=8) :: st

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'

         ierr = 0
         tpo = 0
         ldim(1) = 1
         call h5open_f(ierr)
         treal = detect_precision()
         tnpp = int(npp/dspl)
         nvyp = num_procs_loc()
         ori = id_proc() - nvyp 
         des = id_proc() + nvyp 
  
         if (ori >= 0) then
            call MPI_IRECV(message,1,p_dtype_int,ori,rtag,comm_world(),&
            &mid,ierr)
            call MPI_WAIT(mid,istat,ierr)
         endif

         call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,comm_loc(),ierr)
         if (tp == 0) then
            if (ori < 0) then
               call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
               call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
               info = MPI_INFO_NULL
               call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
               call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
               &access_prp=flplID) 
               call wrattr_file(file,file_id,xferID)
               call h5gopen_f(file_id, '/', rootID, ierr)
               call h5screate_simple_f(1, ldim, aspace_id, ierr)
               call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
               &aid, ierr )
               call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
               call h5aclose_f(aid, ierr)
               call h5sclose_f(aspace_id, ierr)
               call h5pclose_f(xferID, ierr)
               call h5pclose_f(flplID, ierr)
               call h5gclose_f(rootID, ierr)
               call h5fclose_f(file_id, ierr)
            endif
            if (des < num_procs()) then
               call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
               &id,ierr)
            else
               id = MPI_REQUEST_NULL         
            endif
            call h5close_f(ierr)
            return
         else 
            if (tnpp > 0) then 
               color = 1
            else
               color = MPI_UNDEFINED
            endif
            call MPI_COMM_SPLIT(comm_loc(), color, 0, pgrp, ierr )

            if (tnpp > 0) then
               call MPI_COMM_RANK(pgrp, pid, ierr)
               call MPI_COMM_SIZE(pgrp, pnvp, ierr)
               allocate(np(pnvp), dims(2,pnvp), stat = ierr)
               call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
               &pgrp, ierr)
               dims(1, 1) = 1
               dims(2, 1) = np(1) 
               do i = 2, pnvp
                  dims(1,i) = dims(2,i-1) + 1
                  dims(2,i) = dims(1,i) + np(i) - 1
               enddo
               allocate(buff(tnpp), stat = ierr)

               call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
               call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
               info = MPI_INFO_NULL
               call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
               call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    

               if (ori < 0) then
                  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
                  &access_prp=flplID) 
                  call wrattr_file(file,file_id,xferID)
                  call h5gopen_f(file_id, '/', rootID, ierr)
                  call h5screate_simple_f(1, ldim, aspace_id, ierr)
                  call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
                  &aid, ierr )
                  call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
                  call h5aclose_f(aid, ierr)
                  call h5sclose_f(aspace_id, ierr)
               else
                  call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr,&
                  &access_prp=flplID)
                  call h5gopen_f(file_id, '/', rootID, ierr)
                  call h5aopen_f(rootID, 'tp', aid, ierr)
                  call h5aread_f(aid, H5T_NATIVE_INTEGER, tpo, ldim, ierr)
                  tp = tp + tpo
                  call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
                  call h5aclose_f(aid, ierr)
               endif               

               do i = 1, 3
                  buff(1:tnpp) = part(i,1:((tnpp-1)*dspl+1):dspl)*delta(i)+origin(i)
      
                  if (ori >=0 .and. tpo /= 0) then
                     call h5dopen_f(rootID, 'x'//char(iachar('0')+i), dset_id, ierr)
                     ldim(1) = tp
                     call h5dextend_f(dset_id, ldim, ierr)
                     call h5screate_simple_f(1, ldim, dspace_id, ierr)
                     ldim(1) = tnpp
                     call h5screate_simple_f(1, ldim, memspaceID, ierr )
                     start = tpo + dims(1,pid+1) - 1
                     call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
                     &ldim,ierr)
                     call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                     &dspace_id, xfer_prp=xferID)
                  else
                     maxdim = (/H5S_UNLIMITED_F/)
                     ldim(1) = 1
                     call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
                     call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
                     ldim(1) = tp
                     call h5pset_chunk_f(dcplID, 1, ldim, ierr)
                     call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), treal,&
                     &dspace_id, dset_id, ierr, dcplID)
                     ldim(1) = tp
                     call h5dextend_f(dset_id, ldim, ierr)
                     call h5sclose_f(dspace_id, ierr)
                     call h5screate_simple_f(1, ldim, dspace_id, ierr)
                     ldim(1) = tnpp
                     call h5screate_simple_f(1, ldim, memspaceID, ierr )
                     start = tpo + dims(1,pid+1) - 1
                     call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                     &ldim, ierr)
                     call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                     &dspace_id, xfer_prp=xferID)
                     call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                     &name='x_'//char(iachar('0')+i))
                     call h5pclose_f(dcplID, ierr)
                  endif
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               do i = 1, 3
                  buff(1:tnpp) = part((i+3),1:((tnpp-1)*dspl+1):dspl) 
                  if (ori >= 0 .and. tpo /= 0) then
                     call h5dopen_f(rootID, 'p'//char(iachar('0')+i), dset_id, ierr)
                     ldim(1) = tp
                     call h5dextend_f(dset_id, ldim, ierr)
                     call h5screate_simple_f(1, ldim, dspace_id, ierr)
                     ldim(1) = tnpp
                     call h5screate_simple_f(1, ldim, memspaceID, ierr )
                     start = tpo + dims(1,pid+1) - 1
                     call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
                     &ldim,ierr)
                     call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                     &dspace_id, xfer_prp=xferID)
                  else
                     maxdim = (/H5S_UNLIMITED_F/)
                     ldim(1) = 1
                     call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
                     call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
                     ldim(1) = tp
                     call h5pset_chunk_f(dcplID, 1, ldim, ierr)
                     call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), treal,&
                     &dspace_id, dset_id, ierr, dcplID)
                     ldim(1) = tp
                     call h5dextend_f(dset_id, ldim, ierr)
                     call h5sclose_f(dspace_id, ierr)
                     call h5screate_simple_f(1, ldim, dspace_id, ierr)
                     ldim(1) = tnpp
                     call h5screate_simple_f(1, ldim, memspaceID, ierr )
                     start = tpo + dims(1,pid+1) - 1
                     call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                     &ldim, ierr)
                     call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                     &dspace_id, xfer_prp=xferID)
                     call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                     &name='p_'//char(iachar('0')+i))
                     call h5pclose_f(dcplID, ierr)
                  endif
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               buff(1:tnpp) = part(7,1:((tnpp-1)*dspl+1):dspl) 
               if (ori >= 0 .and. tpo /= 0) then
                  call h5dopen_f(rootID, 'q', dset_id, ierr)
                  ldim(1) = tp
                  call h5dextend_f(dset_id, ldim, ierr)
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr )
                  start = tpo + dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
                  &ldim,ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
               else
                  maxdim = (/H5S_UNLIMITED_F/)
                  ldim(1) = 1
                  call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
                  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
                  ldim(1) = tp
                  call h5pset_chunk_f(dcplID, 1, ldim, ierr)
                  call h5dcreate_f(rootID, 'q', treal,&
                  &dspace_id, dset_id, ierr, dcplID)
                  ldim(1) = tp
                  call h5dextend_f(dset_id, ldim, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr )
                  start = tpo + dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                  &ldim, ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
                  call wrattr_dataset(file,dset_id,unit='a.u.',&
                  &name='q')
                  call h5pclose_f(dcplID, ierr)
               endif
               call h5sclose_f(memspaceID, ierr)
               call h5sclose_f(dspace_id, ierr)
               call h5dclose_f(dset_id, ierr)


               call h5pclose_f(xferID, ierr)
               call h5pclose_f(flplID, ierr)
               call h5gclose_f(rootID, ierr)
               call h5fclose_f(file_id, ierr)
               deallocate(np,dims,buff)
            endif
         endif

         if (des < num_procs()) then
            call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
            &id,ierr)
         else
            id = MPI_REQUEST_NULL         
         endif

         if (pgrp /= MPI_COMM_NULL) then
            call MPI_COMM_FREE(pgrp, ierr )
         endif
         
         call h5close_f(ierr)
         
      end subroutine pwpart_3d_pipe_orig
!
subroutine pwpart_3d_pipe(file,x,p,q,npp,dspl,z0,rtag,stag,id,ierr,s, dr, dz)

 implicit none

 class(hdf5file), intent(in) :: file
 real, dimension(:,:), intent(in) :: x, p
 real, dimension(:,:), intent(in), optional :: s
 real, intent(in), optional :: dr, dz
 real, dimension(:), intent(in) :: q
 real, intent(in) :: z0
 integer(kind=LG), intent(in) :: npp
 integer, intent(in) :: dspl
 integer, intent(in) :: rtag, stag
 integer, intent(inout) :: id, ierr
! local data
 integer :: tnpp, tp, tpo, color, pgrp, pid, pnvp, i
 integer(hsize_t), dimension(1) :: ldim
 integer, dimension(:), pointer :: np => null()
 integer, dimension(:,:), pointer:: dims => null()
 real, dimension(:), pointer :: buff => null()
 integer(hsize_t), dimension(1) :: start,maxdim
 integer(hid_t) :: treal
 integer(hid_t) :: flplID, xferID, dcplID, memspaceID, aid
 integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
 integer :: ori, des, nvyp, mid, message, info
 integer, dimension(10) :: istat
 character(len=:), allocatable :: filename
 character(len=8) :: st
 logical :: has_spin

! OPENPMD hierachy
integer(hid_t) :: dataID, iterID, partID, specID, momID, posID, weightID
! other OPENPMD variables
real, dimension(2) :: local_gridspacing
real, dimension(2) :: local_GlobalOffset
real, dimension(2) :: local_position
real :: qsum, total_charge, q0
character(len=20) :: iter_str
logical gexist

 has_spin = .false.
 if (present(s)) has_spin = .true.


 allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
 write (st,'(I8.8)') file%n
 filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'


 qsum = 0
 total_charge = 0
 q0 = 0
 ierr = 0
 tpo = 0
 ldim(1) = 1
 call h5open_f(ierr)
 treal = detect_precision()
 tnpp = int(npp/dspl)
 nvyp = num_procs_loc()
 ori = id_proc() - nvyp
 des = id_proc() + nvyp

 write(iter_str,'(I8)') file%n
! write(*,*) ' in pwpart_3d_pipe',filename

 if (ori >= 0) then
    call MPI_IRECV(message,1,p_dtype_int,ori,rtag,comm_world(),&
    &mid,ierr)
    call MPI_WAIT(mid,istat,ierr)
 endif

 qsum = SUM(q(1:((tnpp-1)*dspl+1):dspl))
 call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,comm_loc(),ierr)
 call MPI_ALLREDUCE(qsum,total_charge,1,p_dtype_real,MPI_SUM,comm_loc(),ierr)
 ! print *,total_charge
 total_charge = total_charge * (2 * pi * dr**2 * dz) * (q_e * file%n0)/(file%kp)**3
 if (tp == 0) then
    if (ori < 0) then
       call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
       info = MPI_INFO_NULL
       call h5pset_fapl_mpio_f(flplID, comm_loc(), info, ierr)
       call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
       &access_prp=flplID)
       call wrattr_file(file,file_id,xferID, type=.false.)
       call h5gopen_f(file_id, '/', rootID, ierr)

       ! added code
       call h5gcreate_f(rootID,'data',dataID,ierr)
       call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
       call h5gcreate_f(iterID,trim(file%particlespath),partID,ierr)
       specID = partID
       call add_h5_atribute(partID, 'speciesType',trim(file%spectype))
       call add_h5_atribute(partID, 'numParticles',tp)
       call add_h5_atribute(partID, 'chargeUnitSI',1.0)
       call add_h5_atribute(partID, 'totalCharge',total_charge)

       call h5gcreate_f(specID,'momentum',momID,ierr)
       call h5gcreate_f(specID,'position',posID,ierr)

       call wrattr_dataset(file,momID,unit='m_ec',&
          &name='p', type=.false.)

       call wrattr_dataset(file,posID,unit='c/\omega_p',&
          &name='x', type=.false.)
       call add_h5_atribute(iterID,'dt',file%dt)
       call add_h5_atribute(iterID,'time',file%t)
       call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)
       call h5gclose_f(momID, ierr)
       call h5gclose_f(posID, ierr)
       call h5gclose_f(partID, ierr)
       call h5gclose_f(iterID, ierr)
       call h5gclose_f(dataID, ierr)
       !end added


       call h5screate_simple_f(1, ldim, aspace_id, ierr)
       call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
       &aid, ierr )
       call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
       call h5aclose_f(aid, ierr)
       call h5sclose_f(aspace_id, ierr)
       call h5pclose_f(xferID, ierr)
       call h5pclose_f(flplID, ierr)
       call h5gclose_f(rootID, ierr)
       call h5fclose_f(file_id, ierr)
    endif
    if (des < num_procs()) then
       call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
       &id,ierr)
    else
       id = MPI_REQUEST_NULL
    endif
    call h5close_f(ierr)
    return
 else
    if (tnpp > 0) then
       color = 1
    else
       color = MPI_UNDEFINED
    endif
    call MPI_COMM_SPLIT(comm_loc(), color, 0, pgrp, ierr )

    if (tnpp > 0) then
       call MPI_COMM_RANK(pgrp, pid, ierr)
       call MPI_COMM_SIZE(pgrp, pnvp, ierr)
       allocate(np(pnvp), dims(2,pnvp), stat = ierr)
       call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
       &pgrp, ierr)
       dims(1, 1) = 1
       dims(2, 1) = np(1)
       do i = 2, pnvp
          dims(1,i) = dims(2,i-1) + 1
          dims(2,i) = dims(1,i) + np(i) - 1
       enddo
       allocate(buff(tnpp), stat = ierr)

       call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
       info = MPI_INFO_NULL
       call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
       call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)

       if (ori < 0) then
          call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr,&
          &access_prp=flplID)
          call wrattr_file(file,file_id,xferID, type = .false.)
          call h5gopen_f(file_id, '/', rootID, ierr)
          call h5gcreate_f(rootID,'data',dataID,ierr)
          call h5gcreate_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
          call h5gcreate_f(iterID,trim(file%particlespath),partID,ierr)
          specID = partID
          call add_h5_atribute(partID, 'speciesType',trim(file%spectype))
          call add_h5_atribute(partID, 'numParticles',tp)
          call add_h5_atribute(partID, 'chargeUnitSI',1.0)
          call add_h5_atribute(partID, 'totalCharge',total_charge)
          call h5gcreate_f(specID,'momentum',momID,ierr)
          call h5gcreate_f(specID,'position',posID,ierr)

          call wrattr_dataset(file,momID,unit='m_ec',&
             &name='p', type=.false.)

          call wrattr_dataset(file,posID,unit='c/\omega_p',&
             &name='x', type=.false.)
          ! call add_h5_atribute(momID,'unitDimension',[ 1., 1., -1., 0., 0., 0., 0.])
          ! call add_h5_atribute(posID,'unitDimension',[ 1., 0., 0., 0., 0., 0., 0.])

          call add_h5_atribute(iterID,'dt',file%dt)
          call add_h5_atribute(iterID,'time',file%t)
          call add_h5_atribute(iterID,'timeUnitSI',file%timeunitsi)

          call h5screate_simple_f(1, ldim, aspace_id, ierr)
          call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
          &aid, ierr )
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
          call h5aclose_f(aid, ierr)
          call h5sclose_f(aspace_id, ierr)

       else
          call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, ierr,&
          &access_prp=flplID)
          call h5gopen_f(file_id, '/', rootID, ierr)
          call h5gopen_f(rootID,'data',dataID,ierr)
          call h5gopen_f(dataID,adjustl(trim(iter_str)),iterID,ierr)
          call h5gopen_f(iterID,trim(file%particlespath),partID,ierr)
          ! call h5gopen_f(partID,'species',specID,ierr)
          specID = partID
          call h5gopen_f(specID,'momentum',momID,ierr)
          call h5gopen_f(specID,'position',posID,ierr)

          call h5aopen_f(rootID, 'tp', aid, ierr)
          call h5aread_f(aid, H5T_NATIVE_INTEGER, tpo, ldim, ierr)
          tp = tp + tpo
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
          call h5aclose_f(aid, ierr)

          call h5aopen_f(partID, 'numParticles', aid, ierr)
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
          call h5aclose_f(aid, ierr)

          call h5aopen_f(partID, 'totalCharge', aid, ierr)
          call h5aread_f(aid, treal, q0, ldim, ierr)
          total_charge = total_charge + q0
          call h5awrite_f(aid, treal, total_charge, ldim, ierr)
          call h5aclose_f(aid, ierr)
       endif


       do i = 1, 3
          ! if (i == 1) then
          !    buff(1:tnpp) = part(1,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else if (i == 2) then
          !    buff(1:tnpp) = part(1,1:((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else
          !    buff(1:tnpp) = part(i,1:((tnpp-1)*dspl+1):dspl)+z0
          ! end if

          if (i == 3) then
            ! buff(1:tnpp) = x(i,1:((tnpp-1)*dspl+1):dspl)+z0
            buff(1:tnpp) = 0
          else
            buff(1:tnpp) = x(i,1:((tnpp-1)*dspl+1):dspl)
          endif

          if (ori >=0 .and. tpo /= 0) then
             call h5dopen_f(posID, file%partlabel(i), dset_id, ierr)
             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
             &ldim,ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
          else
             maxdim = (/H5S_UNLIMITED_F/)
             ldim(1) = 1
             call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
             ldim(1) = tp
             call h5pset_chunk_f(dcplID, 1, ldim, ierr)
             call h5dcreate_f(posID, file%partlabel(i), treal,&
             &dspace_id, dset_id, ierr, dcplID)
            

             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5sclose_f(dspace_id, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
             &ldim, ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
             call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
             &name='x_'//char(iachar('0')+i), type=.false.)
             call h5pclose_f(dcplID, ierr)
          endif
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       do i = 1, 3
          ! if (i == 1) then
          !    buff(1:tnpp) = part(4,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl)) - &
          !    &part(5,1:((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else if (i == 2) then
          !    buff(1:tnpp) = part(4,1:((tnpp-1)*dspl+1):dspl)*&
          !    &sin(part(2,1:((tnpp-1)*dspl+1):dspl)) + &
          !    &part(5,1:((tnpp-1)*dspl+1):dspl)*&
          !    &cos(part(2,1:((tnpp-1)*dspl+1):dspl))
          ! else
          !    buff(1:tnpp) = part(6,1:((tnpp-1)*dspl+1):dspl)
          ! end if

          buff(1:tnpp) = p(i,1:((tnpp-1)*dspl+1):dspl) * m_e * c**2/q_e
          if (ori >= 0 .and. tpo /= 0) then
             call h5dopen_f(momID, file%partlabel(i), dset_id, ierr)
             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
             &ldim,ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
          else
             maxdim = (/H5S_UNLIMITED_F/)
             ldim(1) = 1
             call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
             ldim(1) = tp
             call h5pset_chunk_f(dcplID, 1, ldim, ierr)
             call h5dcreate_f(momID, file%partlabel(i), treal,&
             &dspace_id, dset_id, ierr, dcplID)

             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5sclose_f(dspace_id, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
             &ldim, ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
             call wrattr_dataset(file,dset_id,unit='m_ec',&
             &name='p_'//char(iachar('0')+i), type=.false.)
             call h5pclose_f(dcplID, ierr)
          endif
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       if (has_spin) then

       do i = 1, 3

          buff(1:tnpp) = s(i,1:((tnpp-1)*dspl+1):dspl)

          if (ori >= 0 .and. tpo /= 0) then
             call h5dopen_f(partID, 's'//char(iachar('0')+i), dset_id, ierr)
             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
             &ldim,ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
          else
             maxdim = (/H5S_UNLIMITED_F/)
             ldim(1) = 1
             call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
             call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
             ldim(1) = tp
             call h5pset_chunk_f(dcplID, 1, ldim, ierr)
             call h5dcreate_f(partID, 's'//char(iachar('0')+i), treal,&
             &dspace_id, dset_id, ierr, dcplID)
             ldim(1) = tp
             call h5dextend_f(dset_id, ldim, ierr)
             call h5sclose_f(dspace_id, ierr)
             call h5screate_simple_f(1, ldim, dspace_id, ierr)
             ldim(1) = tnpp
             call h5screate_simple_f(1, ldim, memspaceID, ierr )
             start = tpo + dims(1,pid+1) - 1
             call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
             &ldim, ierr)
             call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
             &dspace_id, xfer_prp=xferID)
             call wrattr_dataset(file,dset_id,unit='a.u.',&
             &name='s_'//char(iachar('0')+i), type=.false.)
             call h5pclose_f(dcplID, ierr)
          endif
          call h5sclose_f(memspaceID, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5dclose_f(dset_id, ierr)
       enddo

       endif


       buff(1:tnpp) = (file%t + (x(3,1:((tnpp-1)*dspl+1):dspl)+z0))* file%timeunitsi
       if (ori >= 0 .and. tpo /= 0) then
          call h5dopen_f(specID, 'time', dset_id, ierr)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
          &ldim,ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
       else
          maxdim = (/H5S_UNLIMITED_F/)
          ldim(1) = 1
          call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
          call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
          ldim(1) = tp
          call h5pset_chunk_f(dcplID, 1, ldim, ierr)
          call h5dcreate_f(specID, 'time', treal,&
          &dspace_id, dset_id, ierr, dcplID)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='a.u.',&
          &name='t', type=.false.)
          call h5pclose_f(dcplID, ierr)
       endif

       call h5sclose_f(memspaceID, ierr)
       call h5sclose_f(dspace_id, ierr)
       call h5dclose_f(dset_id, ierr)


       buff(1:tnpp) = q(1:((tnpp-1)*dspl+1):dspl) * (2 * pi * dr**2 * dz) * (-q_e * file%n0)/(file%kp)**3
       if (ori >= 0 .and. tpo /= 0) then
          call h5dopen_f(specID, 'weight', dset_id, ierr)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
          &ldim,ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
       else
          maxdim = (/H5S_UNLIMITED_F/)
          ldim(1) = 1
          call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
          call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
          ldim(1) = tp
          call h5pset_chunk_f(dcplID, 1, ldim, ierr)
          call h5dcreate_f(specID, 'weight', treal,&
          &dspace_id, dset_id, ierr, dcplID)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='a.u.',&
          &name='weight', type=.false.)
          call h5pclose_f(dcplID, ierr)
       endif

       call h5sclose_f(memspaceID, ierr)
       call h5sclose_f(dspace_id, ierr)
       call h5dclose_f(dset_id, ierr)
       
       if (ori >= 0 .and. tpo /= 0) then
          call h5dopen_f(specID, 'mass', dset_id, ierr)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
          &ldim,ierr)
          call h5dwrite_f(dset_id, treal, buff/(-q_e) * m_e, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
       else
          maxdim = (/H5S_UNLIMITED_F/)
          ldim(1) = 1
          call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
          call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
          ldim(1) = tp
          call h5pset_chunk_f(dcplID, 1, ldim, ierr)
          call h5dcreate_f(specID, 'mass', treal,&
          &dspace_id, dset_id, ierr, dcplID)
          ldim(1) = tp
          call h5dextend_f(dset_id, ldim, ierr)
          call h5sclose_f(dspace_id, ierr)
          call h5screate_simple_f(1, ldim, dspace_id, ierr)
          ldim(1) = tnpp
          call h5screate_simple_f(1, ldim, memspaceID, ierr )
          start = tpo + dims(1,pid+1) - 1
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
          &ldim, ierr)
          call h5dwrite_f(dset_id, treal, buff/(-q_e) * m_e, ldim, ierr, memspaceID,&
          &dspace_id, xfer_prp=xferID)
          call wrattr_dataset(file,dset_id,unit='a.u.',&
          &name='mass', type=.false.)
          call h5pclose_f(dcplID, ierr)
       endif


       call h5sclose_f(memspaceID, ierr)
       call h5sclose_f(dspace_id, ierr)
       call h5dclose_f(dset_id, ierr)



       



       ! if (ori >= 0 .and. tpo /= 0) then
       !    call h5dopen_f(specID, 'weight', dset_id, ierr)
       !    ldim(1) = tp
       !    call h5dextend_f(dset_id, ldim, ierr)
       !    call h5screate_simple_f(1, ldim, dspace_id, ierr)
       !    ldim(1) = tnpp
       !    call h5screate_simple_f(1, ldim, memspaceID, ierr )
       !    start = tpo + dims(1,pid+1) - 1
       !    call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,start,&
       !    &ldim,ierr)
       !    call h5dwrite_f(dset_id, treal, buff/(-q_e), ldim, ierr, memspaceID,&
       !    &dspace_id, xfer_prp=xferID)
       ! else
       !    maxdim = (/H5S_UNLIMITED_F/)
       !    ldim(1) = 1
       !    call h5screate_simple_f(1, ldim, dspace_id, ierr, maxdim)
       !    call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
       !    ldim(1) = tp
       !    call h5pset_chunk_f(dcplID, 1, ldim, ierr)
       !    call h5dcreate_f(specID, 'weight', treal,&
       !    &dspace_id, dset_id, ierr, dcplID)
       !    ldim(1) = tp
       !    call h5dextend_f(dset_id, ldim, ierr)
       !    call h5sclose_f(dspace_id, ierr)
       !    call h5screate_simple_f(1, ldim, dspace_id, ierr)
       !    ldim(1) = tnpp
       !    call h5screate_simple_f(1, ldim, memspaceID, ierr )
       !    start = tpo + dims(1,pid+1) - 1
       !    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
       !    &ldim, ierr)
       !    call h5dwrite_f(dset_id, treal, buff/(-q_e), ldim, ierr, memspaceID,&
       !    &dspace_id, xfer_prp=xferID)
       !    call wrattr_dataset(file,dset_id,unit='a.u.',&
       !    &name='weight', type=.false.)
       !    call h5pclose_f(dcplID, ierr)
       ! endif

       ! call h5sclose_f(memspaceID, ierr)
       ! call h5sclose_f(dspace_id, ierr)
       ! call h5dclose_f(dset_id, ierr)

       

       

       
       

       call h5gclose_f(posID, ierr)
       call h5gclose_f(momID, ierr)
       call h5gclose_f(partID, ierr)
       call h5gclose_f(iterID, ierr)

       call h5gclose_f(dataID, ierr)
       call h5gclose_f(rootID, ierr)
       call h5pclose_f(xferID, ierr)
       call h5pclose_f(flplID, ierr)
       call h5fclose_f(file_id, ierr)
       deallocate(np,dims,buff)
    endif
 endif

 if (des < num_procs()) then
    call MPI_ISEND(message,1,p_dtype_int,des,stag,comm_world(),&
    &id,ierr)
 else
    id = MPI_REQUEST_NULL
 endif

 if (pgrp /= MPI_COMM_NULL) then
    call MPI_COMM_FREE(pgrp, ierr )
 endif

 call h5close_f(ierr)

end subroutine pwpart_3d_pipe
!
!
      subroutine wpart_orig(file,part,npp,dspl,ierr)
     
         implicit none

         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: part
         integer, intent(in) :: npp,dspl
         integer, intent(inout) :: ierr
! local data
         integer :: tp
         integer(hsize_t), dimension(1) :: ldim
         integer(hsize_t), dimension(2) :: dim
         integer(hid_t) :: treal
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
         integer(hid_t) :: memspaceID, aid
         character(len=:), allocatable :: filename
         character(len=8) :: st

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'

         ierr = 0
         ldim(1) = 1
         call h5open_f(ierr)
         treal = detect_precision()
         tp = int(npp/dspl)
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
         call wrattr_file(file,file_id,H5P_DEFAULT_F)
         call h5gopen_f(file_id, '/', rootID, ierr)
         call h5screate_simple_f(1, ldim, aspace_id, ierr)
         call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
         &aid, ierr )
         call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
         call h5aclose_f(aid, ierr)
         call h5sclose_f(aspace_id, ierr)

         if (tp == 0) then
            call h5gclose_f(rootID, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
         else
            dim(1) = size(part,1)
            dim(2) = tp 
            call h5screate_simple_f(2, dim, dspace_id, ierr)
            call h5screate_simple_f(2, dim, memspaceID, ierr)
            call h5dcreate_f(rootID, file%dataname, treal, dspace_id, dset_id,&
            &ierr)
            call wrattr_dataset(file,dset_id)
            call h5dwrite_f(dset_id, treal, part(:,1:(1+(tp-1)*dspl):dspl),&
            &dim, ierr, memspaceID, dspace_id)
            call h5sclose_f(memspaceID, ierr)
            call h5sclose_f(dspace_id, ierr)
            call h5gclose_f(rootID, ierr)
            call h5dclose_f(dset_id, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
         end if
         
      end subroutine wpart_orig
!
      subroutine rpart_orig(file,part,npp,ierr)
     
         implicit none

         ! class(parallel_pipe), intent(in), pointer :: pp
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(inout) :: part
         integer, intent(out) :: npp
         integer, intent(inout) :: ierr
! local data
         integer :: tp
         integer(hsize_t), dimension(1) :: ldim
         integer(hsize_t), dimension(2) :: dim
         integer(hid_t) :: treal
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
         integer(hid_t) :: memspaceID, aid
         character(len=:), allocatable :: filename
         character(len=8) :: st

                  
         allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
         write (st,'(I8.8)') file%n
         filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'

         ierr = 0
         ldim(1) = 1
         call h5open_f(ierr)
         treal = detect_precision()
         call h5fopen_f(filename,H5F_ACC_RDONLY_F, file_id, ierr)
         call h5gopen_f(file_id, '/', rootID, ierr)
         call h5aopen_f(rootID, 'tp', aid, ierr)
         call h5aread_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
         call h5aclose_f(aid, ierr)
         npp = tp

         if (tp == 0) then
            call h5gclose_f(rootID, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
         else
            dim(1) = size(part,1)
            dim(2) = tp 
            call h5screate_simple_f(2, dim, dspace_id, ierr)
            call h5screate_simple_f(2, dim, memspaceID, ierr)
            call h5dopen_f(rootID, file%dataname, dset_id, ierr)
            call h5dread_f(dset_id, treal, part,&
            &dim, ierr, memspaceID, dspace_id)
            call h5sclose_f(memspaceID, ierr)
            call h5sclose_f(dspace_id, ierr)
            call h5gclose_f(rootID, ierr)
            call h5dclose_f(dset_id, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
         end if
         
      end subroutine rpart_orig

subroutine wpart(file,x,p,q,npp,dspl,ierr,s)

 implicit none

 class(hdf5file), intent(in) :: file
 real, dimension(:,:), intent(in) :: x, p
 real, dimension(:,:), intent(in), optional :: s
 real, dimension(:), intent(in) :: q
 integer(kind=LG), intent(in) :: npp
 integer, intent(in) :: dspl
 integer, intent(inout) :: ierr
! local data
 integer :: tp
 integer(hsize_t), dimension(1) :: ldim
 integer(hsize_t), dimension(2) :: dim
 integer(hid_t) :: treal
 integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
 integer(hid_t) :: memspaceID, aid
 character(len=:), allocatable :: filename
 character(len=8) :: st
 logical :: has_spin

 has_spin = .false.
 if (present(s)) has_spin = .true.

 allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
 write (st,'(I8.8)') file%n
 filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'

 ierr = 0
 ldim(1) = 1
 call h5open_f(ierr)
 treal = detect_precision()
 tp = int(npp/dspl)
 call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
 call wrattr_file(file,file_id,H5P_DEFAULT_F)
 call h5gopen_f(file_id, '/', rootID, ierr)
 call h5screate_simple_f(1, ldim, aspace_id, ierr)
 call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
 &aid, ierr )
 call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
 call h5aclose_f(aid, ierr)
 call h5sclose_f(aspace_id, ierr)

 if (tp == 0) then
    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
 else
    ! dim(1) = size(part,1)
    ! dim(2) = tp
    ! call h5screate_simple_f(2, dim, dspace_id, ierr)
    ! call h5screate_simple_f(2, dim, memspaceID, ierr)
    ! call h5dcreate_f(rootID, file%dataname, treal, dspace_id, dset_id,&
    ! &ierr)
    ! call wrattr_dataset(file,dset_id)
    ! call h5dwrite_f(dset_id, treal, part(:,1:(1+(tp-1)*dspl):dspl),&
    ! &dim, ierr, memspaceID, dspace_id)
    ! call h5sclose_f(memspaceID, ierr)
    ! call h5sclose_f(dspace_id, ierr)
    ! call h5gclose_f(rootID, ierr)
    ! call h5dclose_f(dset_id, ierr)
    ! call h5fclose_f(file_id, ierr)
    ! call h5close_f(ierr)

    ! write position
    dim(1) = size(x,1)
    dim(2) = tp
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dcreate_f(rootID, 'x', treal, dspace_id, dset_id,&
    &ierr)
    call wrattr_dataset(file,dset_id)
    call h5dwrite_f(dset_id, treal, x(:,1:(1+(tp-1)*dspl):dspl),&
    &dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    ! write momentum
    dim(1) = size(p,1)
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dcreate_f(rootID, 'p', treal, dspace_id, dset_id,&
    &ierr)
    call wrattr_dataset(file,dset_id)
    call h5dwrite_f(dset_id, treal, p(:,1:(1+(tp-1)*dspl):dspl),&
    &dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    ! write spin
    if (has_spin) then

    dim(1) = size(s,1)
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dcreate_f(rootID, 's', treal, dspace_id, dset_id,&
    &ierr)
    call wrattr_dataset(file,dset_id)
    call h5dwrite_f(dset_id, treal, s(:,1:(1+(tp-1)*dspl):dspl),&
    &dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    endif

    ! write charge
    call h5screate_simple_f(1, dim(2), dspace_id, ierr)
    call h5screate_simple_f(1, dim(2), memspaceID, ierr)
    call h5dcreate_f(rootID, 'q', treal, dspace_id, dset_id,&
    &ierr)
    call wrattr_dataset(file,dset_id)
    call h5dwrite_f(dset_id, treal, q(1:(1+(tp-1)*dspl):dspl),&
    &dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
 end if

end subroutine wpart
!
subroutine rpart(file,x,p,q,npp,ierr,s)

 implicit none

 class(hdf5file), intent(in) :: file
 real, dimension(:,:), intent(inout) :: x, p
 real, dimension(:,:), intent(inout), optional :: s
 real, dimension(:), intent(inout) :: q
 integer(kind=LG), intent(out) :: npp
 integer, intent(inout) :: ierr
! local data
 integer :: tp
 integer(hsize_t), dimension(1) :: ldim
 integer(hsize_t), dimension(2) :: dim
 integer(hid_t) :: treal
 integer(hid_t) :: file_id, rootID, dset_id, dspace_id
 integer(hid_t) :: memspaceID, aid
 character(len=:), allocatable :: filename
 character(len=8) :: st
 logical :: has_spin

 has_spin = .false.
 if (present(s)) has_spin = .true.

 allocate(character(len(trim(file%filename))+len(trim(file%dataname))+11) :: filename)
 write (st,'(I8.8)') file%n
 filename = trim(file%filename)//trim(file%dataname)//'_'//st//'.h5'

 ierr = 0
 ldim(1) = 1
 call h5open_f(ierr)
 treal = detect_precision()
 call h5fopen_f(filename,H5F_ACC_RDONLY_F, file_id, ierr)
 call h5gopen_f(file_id, '/', rootID, ierr)
 call h5aopen_f(rootID, 'tp', aid, ierr)
 call h5aread_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
 call h5aclose_f(aid, ierr)
 npp = tp

 if (tp == 0) then
    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
 else
    ! dim(1) = size(part,1)
    ! dim(2) = tp
    ! call h5screate_simple_f(2, dim, dspace_id, ierr)
    ! call h5screate_simple_f(2, dim, memspaceID, ierr)
    ! call h5dopen_f(rootID, file%dataname, dset_id, ierr)
    ! call h5dread_f(dset_id, treal, part,&
    ! &dim, ierr, memspaceID, dspace_id)
    ! call h5sclose_f(memspaceID, ierr)
    ! call h5sclose_f(dspace_id, ierr)
    ! call h5gclose_f(rootID, ierr)
    ! call h5dclose_f(dset_id, ierr)
    ! call h5fclose_f(file_id, ierr)
    ! call h5close_f(ierr)

    ! read position
    dim(1) = size(x,1)
    dim(2) = tp
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dopen_f(rootID, 'x', dset_id, ierr)
    call h5dread_f(dset_id, treal, x, dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    ! read momentum
    dim(1) = size(p,1)
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dopen_f(rootID, 'p', dset_id, ierr)
    call h5dread_f(dset_id, treal, p, dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    ! read spin
    if (has_spin) then

    dim(1) = size(s,1)
    call h5screate_simple_f(2, dim, dspace_id, ierr)
    call h5screate_simple_f(2, dim, memspaceID, ierr)
    call h5dopen_f(rootID, 's', dset_id, ierr)
    call h5dread_f(dset_id, treal, s, dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    endif

    ! read charge
    call h5screate_simple_f(1, dim(2), dspace_id, ierr)
    call h5screate_simple_f(1, dim(2), memspaceID, ierr)
    call h5dopen_f(rootID, 'q', dset_id, ierr)
    call h5dread_f(dset_id, treal, q, dim, ierr, memspaceID, dspace_id)
    call h5sclose_f(memspaceID, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5gclose_f(rootID, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
 end if

end subroutine rpart
!

!
      function detect_precision()
         integer(hid_t) :: detect_precision
! local data
         real :: small

         small = 1.0e-12
         small = 1.0 + small
         if (small>1.0) then 
            detect_precision = H5T_NATIVE_DOUBLE 
         else
            detect_precision = H5T_NATIVE_REAL
         endif       

      end function detect_precision     
!      
      end module hdf5io_class

