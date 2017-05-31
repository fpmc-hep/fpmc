c---------------------------------------------------------------
c initialise all the parameters
c---------------------------------------------------------------

      subroutine fpmc_var_ini(readcard)
      implicit none
      logical readcard
      external ffread

      include "ffcard.inc" ! get all common blocks in place

      call ffread

      end subroutine

c---------------------------------------------------------------
c dummy function for the end of herwig execution
c---------------------------------------------------------------

      subroutine hwaend
      end subroutine

