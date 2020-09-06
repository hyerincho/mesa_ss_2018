      program run
      use run_star_support, only: do_read_star_job
      use run_star_extras_astero, only: do_run_star_astero
      use utils_lib, only: mesa_error
      
      implicit none
      
      integer :: ierr
      
      ierr = 0
      call do_read_star_job('inlist', ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      call do_run_star_astero
      
      end program
