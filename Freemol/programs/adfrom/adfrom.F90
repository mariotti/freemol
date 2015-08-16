!
!H
!H  
!H----------------------------------------------------------------------
!H----------------------------------------------------------------------
!H  PROGRAM adfrom
!H----------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H----------------------------------------------------------------------
!H
!H  Read data from ADF Tape21 and produce ASCII readable format
!H  
!H
!H----------------------------------------------------------------------
!H
!
program adfrom
  !
  !  Modules
  !
  use vartypes
  use adfrommodule
  use adfmodule
!  use FROMADF
!  use KF
!  use frstdio
!  use adfrom_general_data
  !
  implicit none
  !
  ! Some default parameters: see makefile definitions
  !
#ifdef PRGPRINT
  integer(KINT), parameter :: adfrom_print=1
#else
  integer(KINT), parameter :: adfrom_print=0
#endif
  !
  !  declare
  !

  !
  !  code
  !
  ! Initialize adfrom
  call adfrom_init()
  !
  ! Initialize adf and input tape21 file
  call adf99_init(unit_in,file_in)
  !
  !
  call read_general(unit_in)
  !
  call read_geometry(unit_in)
  !
  call write_general(unit_out)
  !
  call write_geometry(unit_out)
  !
  ! we try to save before errors at least the geometry
  call flush(unit_out)
  !
  call read_basis(unit_in)
  !
  call read_symmetry(unit_in)
  !
  call write_basis(unit_out,sto_fmt)
  !
  call write_symmetry(unit_out,mo_fmt)
  !
  !  call af_proc_cmdline(in_unit,FileInput,out_unit,FileOutput)
  !
  !  call init_adf99(in_unit,FileInput)
  !
  !  call rw_gen_section()
  !
  !  call rw_basis_section()
  !
  !  call rw_sym_section()
  !
  !  if(runtype(1:8).eq.'FREQUENC') then
  !     call rw_freq_section()
  !  endif
  !
  ! TERMINATE
  call adf99_deallocate()
  call adfrom_end()
  !
end program adfrom
