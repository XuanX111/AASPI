MODULE aaspi_horizon_util
!______________________________________________________
! a suite of programs to read various horizon formats
!
! Author: Kurt J. Marfurt
!         Tao Dash Zhao
! AASPI, University of Oklahoma
! Dec 12, 2014        
!______________________________________________________

!_______________________________________________________________________
! The following subtroutines are contained here:
!
! SUBROUTINE read_horizon
! SUBROUTINE read_interpolated_horizon
! SUBROUTINE read_gridded_horizon
! SUBROUTINE extract_gridded_horizon_value
! SUBROUTINE extract_interpolated_horizon_value
! SUBROUTINE eight_point_interpolation
! SUBROUTINE read_picked_horizon
!
!_______________________________________________________________________





IMPLICIT NONE
!______________________________________________________________________
! variables defining the dimensions of the input seismic data volume
!______________________________________________________________________
TYPE horizon_file  
   REAL, DIMENSION(:,:), POINTER    :: time =>  NULL()     	! time picks at line and cdp bins
   INTEGER 			    :: nx_grid		   	! nx value of the grid
   INTEGER  		    	    :: ny_grid		   	! ny value of the grid
   REAL    			    :: dx_grid		   	! dx value of the grid
   REAL     		    	    :: dy_grid		   	! dy value of the grid
   REAL    			    :: hmin   		   	! mininum horizon time or depth (t axis positive down)
   REAL     		    	    :: hmax   		   	! maximum horizon time or depth (t axis positive down)
   REAL    			    :: xmin   		   	! minimum horizon x value (origin of grid position 1)
   REAL     		    	    :: ymin   		   	! minimum horizon y value (origin of grid position 1)
   REAL    			    :: xmax   		   	! maximum horizon x value 
   REAL     		    	    :: ymax   		   	! maximum horizon y value 
   REAL     		    	    :: znull  		   	! dead pick flag           
   REAL     		    	    :: eps    		   	! round off value          
   REAL     		    	    :: scalet 		   	! scale factor applied to picks to convert from ms to s, etc.
   INTEGER                 	    :: first_line 	   	! first line of the horizon 
   INTEGER                          :: last_line  	   	! last line of the horizon
   INTEGER                          :: first_cdp  	   	! first cdp of the horizon
   INTEGER                          :: last_cdp   	   	! last cdp of the horizon 
   LOGICAL                          :: is_gridded_horizon  	! if == .TRUE. picks are defined as orthogonal (N-S, E-W) B-spline knots
   LOGICAL                          :: is_interpolated_horizon  ! if == .TRUE. picks have been interpolated to each cdp_no,line_no pair
   LOGICAL                          :: vertical_axis_down       ! if == .TRUE. the vertical axis is positive down
END TYPE

CONTAINS
!======================================================================
!======================================================================

SUBROUTINE read_horizon(time,x,y,cdp_no,line_no,live,ncdp,nline,horizon,        &
                        file_command_line_arg,lu_horizon,			&
                        input_error,stderr)
!______________________________________________________________________
! Read in a horizon interpolated to cdp_no and line_no (e.g. seisx format)
!______________________________________________________________________
!     Author
!     Kurt J. Marfurt
!     AASPI
!     The University of Oklahoma
!     July 21, 2016   
!_____________________________________________________________________
! modules to allow checking of subroutine argument types
!_____________________________________________________________________
USE aaspi_io                    ! aaspi io module
USE string_util                 ! aaspi string utility module
USE aaspi_file_util             ! aaspi file utilities
USE fftw3_util                  ! fftw version 3 utilities
!______________________________________________________________________
! set all variables to have a default type of 'NONE'.
! doing allows the compiler to find typos, mispellings, and other simple errors
!______________________________________________________________________
IMPLICIT NONE
!______________________________________________________________________
! variables passed to/from the calling routine
!______________________________________________________________________
TYPE(horizon_file)      :: horizon

INTEGER, INTENT(IN)     :: ncdp              		! number of cdps
INTEGER, INTENT(IN)     :: nline             		! number of lines
INTEGER, INTENT(IN)     :: cdp_no(ncdp,nline)		! cdp no of each trace
INTEGER, INTENT(IN)     :: line_no(ncdp,nline)		! line no of each trace
REAL,    INTENT(IN)     :: x(ncdp,nline)		! x-coordinate of each trace
REAL,    INTENT(IN)     :: y(ncdp,nline)		! y-coordinate of each trace
LOGICAL, INTENT(INOUT)  :: live(ncdp,nline)		! y-coordinate of each trace
REAL,    INTENT(OUT)    :: time(ncdp,nline)		! sample value of the horizon at each trace
CHARACTER (LEN=*), INTENT(IN)	:: file_command_line_arg! command line argument to call file name
INTEGER,           INTENT(IN)   :: lu_horizon		! horizon i/o unit number
INTEGER,           INTENT(INOUT):: input_error		! input error counter    
INTEGER,           INTENT(IN)   :: stderr		! standard error i/o unit
!______________________________________________________________________
! local variables 
!______________________________________________________________________
CHARACTER (LEN=256)     :: horizon_fn	        	! horizon file name
INTEGER			:: jline			! line index
INTEGER			:: jcdp  			! cdp index

!______________________________________________________________________
! read in the horizon in one of two formats
!______________________________________________________________________
CALL read_picked_horizon(horizon,horizon_fn,file_command_line_arg,lu_horizon, &
                         horizon%hmin,horizon%hmax,input_error,stderr)
IF(input_error > 0) RETURN
loop_over_line:  DO jline=1,nline
   loop_over_cdp: DO jcdp=1,ncdp
!___________________________________________________________________
!     extract the times values of the shallower and deeepr horizons
!     traces that fall outside the defined horizons are set to be dead!
!___________________________________________________________________
      IF(horizon%is_gridded_horizon) THEN
         CALL extract_gridded_horizon_value(time(jcdp,jline),x(jcdp,jline),y(jcdp,jline),live(jcdp,jline),horizon,horizon%znull)

      ELSE IF(horizon%is_interpolated_horizon) THEN
!___________________________________________________________________
!        check that the current trace falls within the horizon bounds, then extract
!___________________________________________________________________
         CALL extract_interpolated_horizon_value(time(jcdp,jline),cdp_no(jcdp,jline),line_no(jcdp,jline),live(jcdp,jline),horizon,horizon%znull)
      END IF
      IF(live(jcdp,jline)) THEN
!___________________________________________________________________
!        For all live traces, check again if there is a valid horizon pick. If not, set traces to be dead
!___________________________________________________________________
         IF(ABS(time(jcdp,jline)) < 1.0E-5) THEN
            live(jcdp,jline)=.FALSE.
            CYCLE loop_over_cdp
         END IF
      END IF
   END DO loop_over_cdp
END DO loop_over_line


END SUBROUTINE read_horizon

!======================================================================
!======================================================================


SUBROUTINE read_interpolated_horizon(interpolated_horizon,nskip,ncol,ZNULL,	&
                                     vertical_axis_down,scalet,line_col,cdp_col,time_col,lu_in,stderr)
!______________________________________________________________________
! Read in a horizon interpolated to cdp_no and line_no (e.g. seisx format)       
!______________________________________________________________________
!     Author
!     Tao Dash Zhao          
!     AASPI                             
!     The University of Oklahoma
!     October 1, 2014   
!_____________________________________________________________________
! modules to allow checking of subroutine argument types
!_____________________________________________________________________
USE aaspi_io                    ! aaspi io module
USE string_util                 ! aaspi string utility module
USE aaspi_file_util             ! aaspi file utilities
USE fftw3_util                  ! fftw version 3 utilities
!______________________________________________________________________
! set all variables to have a default type of 'NONE'.
! doing allows the compiler to find typos, mispellings, and other simple errors
!______________________________________________________________________
IMPLICIT NONE
!______________________________________________________________________
! I/O variables
!______________________________________________________________________
TYPE(horizon_file) 	:: interpolated_horizon 

INTEGER, INTENT(IN)     :: nskip 	! number of header lines to skip
INTEGER, INTENT(IN)     :: ncol  	! number of columns in the horizon file
INTEGER, INTENT(IN)     :: line_col  	! column number of line_no 
INTEGER, INTENT(IN)     :: cdp_col  	! column number of cdp_no 
INTEGER, INTENT(IN)     :: time_col  	! column number of time
INTEGER, INTENT(IN)     :: stderr  	! standard error unit 
INTEGER, INTENT(IN) 	:: lu_in        ! horizon i/o unit
REAL,    INTENT(IN)     :: ZNULL        ! a 'ZNULL' value indicating an undefined horizon value
REAL,    INTENT(IN)     :: scalet       ! scale factor to convert horizon picks to be the same units as the seismic data
LOGICAL, INTENT(IN)     :: vertical_axis_down		! if == .TRUE. vertical axis is positive down
!______________________________________________________________________
! local variables
!______________________________________________________________________
INTEGER      :: jchar1 			! index used in character strings
INTEGER      :: jchar2			! index used in character strings
INTEGER      :: cdp_no     		! current cdp number        
INTEGER      :: line_no    		! current line number        
CHARACTER(LEN=256)	:: horizon_name ! seisx horizon name 
CHARACTER(LEN=256)	:: grid_name    ! seisx grid_name          
CHARACTER(LEN=256)	:: horizon_fn   ! seisx-format horizon file name
REAL,DIMENSION(1:ncol)  :: cols         ! array for reading in all columns
REAL         :: eps           ! a fraction of ZNULL used to protect against roundoff error

!______________________________________________________________________
! miscellaneous loop indices
!______________________________________________________________________
INTEGER      :: io_error	 	! i/o error counter
INTEGER      :: iskip                   ! loop index for skipping header lines
!______________________________________________________________________
! variables used in defining the horizon picks
!______________________________________________________________________
INTEGER      :: idummy     		! dummy integer variable
REAL         :: xa         		! x coordinate of the pick
REAL         :: ya         		! y coordinate of the pick
REAL         :: t   	 		! time of the pick
!_____________________________________________________________________
! control variables                              
!_____________________________________________________________________
INTEGER            :: input_error      	! counter for command line and other data errors
INTEGER            :: allocation_status	! return code for Fortran90 allocate command
INTEGER            :: iostat            ! status of input/output attempt
LOGICAL		   :: verbose		! if == .TRUE. output verbose output
!_______________________________________________________________________
! special constants
!_______________________________________________________________________
!_______________________________________________________________________
! read all horizon times into memory
! if they do not exist, set the output seismic data trace to be dead
! if they fall outside the bounds of the seismic survey, throw them away
! initialize all horizon time values to be undefined
!_______________________________________________________________________
!_______________________________________________________________________
! skip all header lines to reach the first line of actual data
!_______________________________________________________________________
REWIND(lu_in)
iskip=0
eps=.0001*ABS(ZNULL)
loop_over_headers1: DO
   iskip=iskip+1
   IF(iskip.GT.nskip) THEN
      EXIT loop_over_headers1
   ELSE
      READ(lu_in,*)
   END IF 
END DO loop_over_headers1
!_______________________________________________________________________
! loop over all the input lines of data and determine the minimum and maximum values
!______________________________________________________________________
interpolated_horizon%first_cdp=+1000000
interpolated_horizon%last_cdp=-1000000
interpolated_horizon%first_line=+1000000
interpolated_horizon%last_line=-1000000
interpolated_horizon%hmin=+1.E+20
interpolated_horizon%hmax=-1.E+20

loop_over_input_cards1: DO
!_______________________________________________________________________
!  match the columns to corresponding dimensions
!_______________________________________________________________________
   READ(lu_in,*,IOSTAT=iostat) cols(1:ncol)
      line_no=INT(cols(line_col))
      cdp_no=INT(cols(cdp_col))
      t=cols(time_col)
!_______________________________________________________________________
!  check if there was an error or end of file encountered in reading this string
!  if so, we've hit the end of the ascii pick file
!_______________________________________________________________________
   IF(iostat /= 0) EXIT loop_over_input_cards1
   IF(ABS(t-ZNULL) > eps) THEN
      interpolated_horizon%first_cdp=MIN(interpolated_horizon%first_cdp,cdp_no)
      interpolated_horizon%last_cdp=MAX(interpolated_horizon%last_cdp,cdp_no)
      interpolated_horizon%first_line=MIN(interpolated_horizon%first_line,line_no)
      interpolated_horizon%last_line=MAX(interpolated_horizon%last_line,line_no)
   END IF
END DO loop_over_input_cards1
!_______________________________________________________________________
! the number of grid points has been determined and the memory has been allocated
!_______________________________________________________________________
ALLOCATE(interpolated_horizon%time(interpolated_horizon%first_cdp:interpolated_horizon%last_cdp,     &
                                   interpolated_horizon%first_line:interpolated_horizon%last_line),  & 
         STAT=allocation_status)               
write(0,*) 'allocate memory for time.  interpolated_horizon%first_cdp:interpolated_horizon%last_cdp =',interpolated_horizon%first_cdp,interpolated_horizon%last_cdp
write(0,*) 'allocate memory for time.  interpolated_horizon%first_line:interpolated_horizon%last_line =',interpolated_horizon%first_line,interpolated_horizon%last_line
!_______________________________________________________________________
! initialize all picks to be ZNULL
!_______________________________________________________________________
interpolated_horizon%time(:,:)=ZNULL  
!_______________________________________________________________________
! set the read pointer back to the first line of actual data, then start reading
!_______________________________________________________________________
REWIND(lu_in)
iskip=0
loop_over_headers2: DO
   iskip=iskip+1
   IF(iskip.GT.nskip) THEN
      EXIT loop_over_headers2
   ELSE
      READ(lu_in,*)
   END IF 
END DO loop_over_headers2            
loop_over_input_cards2: DO
!_______________________________________________________________________
!  match the columns to corresponding dimensions
!_______________________________________________________________________
   READ(lu_in,*,IOSTAT=iostat) cols(1:ncol)
   line_no=INT(cols(line_col))
   cdp_no=INT(cols(cdp_col))
   t=cols(time_col)
!_______________________________________________________________________
!  check if there was an error or end of file encountered in reading this string
!  if so, we've hit the end of the ascii pick file
!_______________________________________________________________________
   IF(iostat /= 0) EXIT loop_over_input_cards2
   IF(ABS(t-ZNULL) > eps) THEN
!_______________________________________________________________________
!     change vertical axis to be positive down if it was not
!_______________________________________________________________________
      IF(.NOT. vertical_axis_down) t=-t
      interpolated_horizon%time(cdp_no,line_no)=scalet*t
      interpolated_horizon%hmin=MIN(interpolated_horizon%time(cdp_no,line_no),interpolated_horizon%hmin)
      interpolated_horizon%hmax=MAX(interpolated_horizon%time(cdp_no,line_no),interpolated_horizon%hmax)    
   END IF            
END DO loop_over_input_cards2

END SUBROUTINE read_interpolated_horizon 

!======================================================================
!
!======================================================================
SUBROUTINE read_gridded_horizon(gridded_horizon,ZNULL,vertical_axis_down,scalet,lu_in,stderr)
!______________________________________________________________________
! Read in an earthvision surface                                  
!======================================================================
!     Author
!     Tao Dash Zhao          
!     AASPI                             
!     The University of Oklahoma
!     October 1, 2014   
!_____________________________________________________________________
! I use modules to allow checking of subroutine argument types
!_____________________________________________________________________
USE aaspi_io                    ! aaspi io module
USE string_util                 ! aaspi string utility module
USE aaspi_file_util             ! aaspi file utilities
USE fftw3_util                  ! fftw version 3 utilities
!______________________________________________________________________
! set all variables to have a default type of 'NONE'.
! doing allows the compiler to find typos, mispellings, and other simple errors
!______________________________________________________________________
IMPLICIT NONE
!______________________________________________________________________
! variables defining the dimensions of the input seismic data volume
!______________________________________________________________________
TYPE(horizon_file) :: gridded_horizon 
REAL,INTENT(IN)    :: ZNULL              	! a 'ZNULL' value indicating an undefined horizon value
REAL,   INTENT(IN) :: scalet       ! scale factor to convert horizon picks to be the same units as the seismic data
LOGICAL,INTENT(IN) :: vertical_axis_down	! if == .TRUE. vertical axis is positive down
INTEGER,INTENT(IN) :: lu_in        		! horizon i/o unit
INTEGER,INTENT(IN) :: stderr  			! standard error unit 
!______________________________________________________________________
! local variables               
!______________________________________________________________________
INTEGER      :: jchar1 			! index used in character strings
INTEGER      :: jchar2			! index used in character strings
INTEGER      :: jx                      ! current or closest horizon node in x
INTEGER      :: jy                      ! current or closest horizon node in y
REAL         :: dx_half                 ! dx/2
REAL         :: dy_half                 ! dy/2
INTEGER      :: io_error	 	! i/o error counter
!______________________________________________________________________
! variables used in defining the horizon picks
!______________________________________________________________________
INTEGER      :: node_count 		! the total number of horizon picks that fall within the survey
INTEGER      :: j          		! dummy loop index
REAL         :: xa         		! x coordinate of the pick
REAL         :: ya         		! y coordinate of the pick
REAL         :: time			! time of the pick
INTEGER      :: cdp_no     		! current cdp number        
INTEGER      :: line_no    		! current line number        
CHARACTER(LEN=256)  :: string           ! a character string used in reading ascii horizon pick file
REAL         :: eps           		! a fraction of ZNULL used to protect against roundoff error

!_____________________________________________________________________
! control variables                              
!_____________________________________________________________________
INTEGER            :: input_error      	! counter for command line and other data errors
INTEGER            :: allocation_status	! return code for Fortran90 allocate command
INTEGER            :: iostat            ! status of input/output attempt
LOGICAL		   :: verbose		! if == .TRUE. output verbose output
!_______________________________________________________________
!_______________________________________________________________
! initialize the AASPI software package    
!_______________________________________________________________
CALL aaspi_initialize()

eps=.0001*ABS(ZNULL)
node_count=0
!_______________________________________________________________________
! read all horizon times into memory
! if they do not exist, set the output seismic data trace to be dead
! if they fall outside the bounds of the seismic survey, throw them away
! initialize all horizon time values to be undefined
!_______________________________________________________________________
gridded_horizon%hmin=1.E+20
gridded_horizon%hmax=-1.E+20
loop_over_input_card_deck: DO 
   READ(lu_in,'(a80)',IOSTAT=iostat) string
!_______________________________________________________________________
!  check if there was an error or end of file encountered in reading this string
!  if so, we've hit the end of the ascii pick file
!_______________________________________________________________________
   IF(iostat /= 0) EXIT loop_over_input_card_deck 
   IF(string(1:1) == '#') THEN
     !WRITE(stderr,*) ' comment card ',j,' : ',string
!_______________________________________________________________________
!     look for the card containing the number of B-spline knots (nodes)
!_______________________________________________________________________
      jchar1=INDEX(string,'Grid_size:')
      IF(jchar1 > 0) THEN
         jchar1=jchar1+10
         jchar2=INDEX(string,'x') 
         READ(string(jchar1:jchar2-1),*) gridded_horizon%nx_grid
         READ(string(jchar2+1:LEN(string)),*) gridded_horizon%ny_grid
         WRITE(stderr,*) 'Prepare to allocate memory for horizon. Number of grid knots in x = ',gridded_horizon%nx_grid,' in y ',gridded_horizon%ny_grid
         ALLOCATE(gridded_horizon%time(1:gridded_horizon%nx_grid,1:gridded_horizon%ny_grid), STAT=allocation_status)               
         IF(allocation_status /= 0) STOP 'Could not allocate space to store horizon nodes!'
         gridded_horizon%time(:,:)=ZNULL
      END IF

!_______________________________________________________________________
!     look for the card containing the limits of 'surface'              
!_______________________________________________________________________
      jchar1=INDEX(string,'Grid_space:') 
      IF(jchar1 > 0) THEN
         jchar1=jchar1+11
         !jchar2=INDEX(string,',')
         READ(string(jchar1:LEN(string)),*) gridded_horizon%xmin,gridded_horizon%xmax,gridded_horizon%ymin,gridded_horizon%ymax
         WRITE(stderr,'(a25,f15.3)') 'gridded_horizon%xmin',gridded_horizon%xmin,'gridded_horizon%xmax',gridded_horizon%xmax,'gridded_horizon%ymin',gridded_horizon%ymin,'gridded_horizon%ymax',gridded_horizon%ymax
         gridded_horizon%dx_grid=(gridded_horizon%xmax-gridded_horizon%xmin)/(gridded_horizon%nx_grid-1)
         gridded_horizon%dy_grid=(gridded_horizon%ymax-gridded_horizon%ymin)/(gridded_horizon%ny_grid-1)
         WRITE(stderr,'(a25,f15.3)') 'gridded_horizon%dx_grid',gridded_horizon%dx_grid,'gridded_horizon%dy_grid',gridded_horizon%dy_grid
      END IF 
      CYCLE loop_over_input_card_deck
   END IF
   READ(string,*) xa,ya,time
   node_count=0
   jx=NINT((xa-gridded_horizon%xmin)/gridded_horizon%dx_grid)+1
   jy=NINT((ya-gridded_horizon%ymin)/gridded_horizon%dy_grid)+1
   IF(jy >=1 .AND. jy <= gridded_horizon%ny_grid) THEN
      IF(jx >=1 .AND. jx <= gridded_horizon%nx_grid) THEN
!___________________________________________________________________
!        grid nodes fall within the defined limits.
!        if not a ZNULL value, convert horizon time from ms to s
!_______________________________________________________________________
         IF(ABS(time-ZNULL) > eps) THEN
            IF(.NOT. vertical_axis_down) time=-time
            gridded_horizon%time(jx,jy)=scalet*time
            gridded_horizon%hmin=MIN(gridded_horizon%hmin,gridded_horizon%time(jx,jy))
            gridded_horizon%hmax=MAX(gridded_horizon%hmax,gridded_horizon%time(jx,jy))
         END IF
         node_count=node_count+1
         
      ELSE
         WRITE(stderr,*) 'node at (jx,jy) = (',jx,',',jy,') falls outside the defined limits!'
         STOP 666
      END IF
   ELSE
      WRITE(stderr,*) 'node at (jx,jy) = (',jx,',',jy,') falls outside the defined limits!'
      STOP 666
   END IF
END DO loop_over_input_card_deck
WRITE(stderr,'(a25,f15.3)') 'gridded_horizon%hmin',gridded_horizon%hmin,'gridded_horizon%hmax',gridded_horizon%hmax

END SUBROUTINE read_gridded_horizon

!=========================================================================================================
!
!=========================================================================================================
SUBROUTINE extract_interpolated_horizon_value(time,cdp_no,line_no,live,horizon,ZNULL)
!_______________________________________________________________________________________________________________________
! determine time value at a trace at location x,y for a gridded horzon (e.g. Earthvision format)
!_______________________________________________________________________________________________________________________
IMPLICIT NONE
TYPE(horizon_file) 	:: horizon 
INTEGER,        INTENT(IN)      ::  cdp_no              ! CDP No. of input trace
INTEGER,        INTENT(IN)      ::  line_no             ! Line No. of input trace
REAL,           INTENT(IN)      ::  ZNULL               ! ZNULL value indicating horizon value does not exist
REAL,           INTENT(OUT)     ::  time                ! two-way travel time of horizon at x,y
LOGICAL,        INTENT(INOUT)   ::  live                ! live trace flag

!___________________________________________________________________________________
! local arrays
!___________________________________________________________________________________
REAL         :: eps           ! a fraction of ZNULL used to protect against roundoff error

IF(live) THEN
   eps=.0001*ABS(ZNULL)
!___________________________________________________________________
!  check to see if trace falls within the bounds of the interpolated horizon
!___________________________________________________________________
   IF(cdp_no < horizon%first_cdp .OR. cdp_no > horizon%last_cdp   .OR. &
      line_no < horizon%first_line .OR. line_no > horizon%last_line  ) THEN
!___________________________________________________________________
!     current trace falls out of the range of the interpolated horizons
!___________________________________________________________________
      live=.FALSE.
      time=0.0
   ELSE
!___________________________________________________________________
!     current trace falls within the range of the interpolated horizon
!     check for znull!
!___________________________________________________________________
      time=horizon%time(cdp_no,line_no)
      IF( (ABS(time-znull) < eps) ) THEN
         live=.FALSE.               
      END IF
   END IF
ELSE
!___________________________________________________________________
!  trace was already dead. set time to 0.0
!___________________________________________________________________
   time=0.0
END IF

END SUBROUTINE extract_interpolated_horizon_value
!=========================================================================================================
!
!=========================================================================================================
SUBROUTINE extract_gridded_horizon_value(time,x,y,live,horizon,ZNULL)
!_______________________________________________________________________________________________________________________
! determine time value at a trace at location x,y for a gridded horzon (e.g. Earthvision format)
!_______________________________________________________________________________________________________________________
IMPLICIT NONE
TYPE(horizon_file) 	:: horizon 
REAL,           INTENT(IN)      ::  x,y                 ! trace location
REAL,           INTENT(IN)      ::  ZNULL               ! ZNULL value indicating horizon value does not exist
REAL,           INTENT(OUT)     ::  time                ! two-way travel time of horizon at x,y
LOGICAL,        INTENT(INOUT)   ::  live                ! live trace flag

!___________________________________________________________________________________
! local arrays
!___________________________________________________________________________________
INTEGER      :: jx1,jy1       ! node indices to the NE of the current trace
INTEGER      :: jx2,jy2       ! node indices to the NW of the current trace
INTEGER      :: jx3,jy3       ! node indices to the SW of the current trace
INTEGER      :: jx4,jy4       ! node indices to the SE of the current trace
REAL         :: xksi          ! local coordinate along x within a rectangular region
REAL         :: eta           ! local coordinate along y within a rectangular region
REAL         :: x3            ! x coordinate of node SW of current trace
REAL         :: y3            ! y coordinate of node SW of current trace
REAL         :: dx_half       ! dx/2
REAL         :: dy_half       ! dy/2
INTEGER      :: jx,jy         ! increment
REAL         :: eps           ! a fraction of ZNULL used to protect against roundoff error

IF(live) THEN
!___________________________________________________________________
!  input trace is currently live
!  proceed to see if it falls within the area of the defined horizon
!___________________________________________________________________

   eps=.0001*ABS(ZNULL)
   dx_half=horizon%dx_grid/2.0
   dy_half=horizon%dy_grid/2.0
!___________________________________________________________________
!  find the horizon node just to the SW of the current trace for horizon
!___________________________________________________________________
   jx3=((x-horizon%xmin)/horizon%dx_grid)+1
   jy3=((y-horizon%ymin)/horizon%dy_grid)+1
!___________________________________________________________________
!     define a square containing the point
!     number the nodes counterclockwise from the NE corner
!___________________________________________________________________
!
!       (x2,y2)---------(x1,y1)
!          |               |
!          |               |
!          |               |            ^
!          |               |            |
!          |  (x,y)        |            |
!          |               |            N
!       (x3,y3)---------(x4,y4)
!
!
!___________________________________________________________________
   jx1=jx3+1
   jy1=jy3+1
   jx2=jx3
   jy2=jy3+1
   jx4=jx3+1
   jy4=jy3
!___________________________________________________________________
!   make sure all the neighboring nodes fall within the limits of defined horizon
!___________________________________________________________________
   IF(jx3 < 1 .OR. jx1 > horizon%nx_grid .OR. jy3 <1 .OR. jy3 > horizon%ny_grid) THEN
!___________________________________________________________________
!     trace is beyond the limits of the defined horizon
!___________________________________________________________________
      time=0.0
!___________________________________________________________________
!     make sure all the neighboring nodes are defined
!___________________________________________________________________
   ELSE IF( (ABS(horizon%time(jx1,jy1)-ZNULL) > eps) .AND. ( ABS(horizon%time(jx2,jy2)-ZNULL) > eps) .AND. &
       (ABS(horizon%time(jx3,jy3)-ZNULL) > eps) .AND. ( ABS(horizon%time(jx4,jy4)-ZNULL) > eps) ) THEN
!___________________________________________________________________
!      interpolate using a bilinear (finite element) interpolator of the form t=axy+bx+cy+d
!      first step: compute the local coordinates in the rectangle (xksi,eta) where
!      -1.0 <= xksi <= +1.0 and  -1.0 <= eta <= +1.0
!___________________________________________________________________
      x3=horizon%xmin+(jx3-1)*horizon%dx_grid
      y3=horizon%ymin+(jy3-1)*horizon%dy_grid
      xksi=(x-(x3+dx_half))/dx_half
      eta=(y-(y3+dy_half))/dy_half
!___________________________________________________________________
!     second step: interpolate using the shape functions
!___________________________________________________________________
      time=0.25*( (1.0+xksi)*(1+eta)*horizon%time(jx1,jy1)   &
                 +(1.0-xksi)*(1+eta)*horizon%time(jx2,jy2)   &
                 +(1.0-xksi)*(1-eta)*horizon%time(jx3,jy3)   &
                 +(1.0+xksi)*(1-eta)*horizon%time(jx4,jy4))
   ELSE
!___________________________________________________________________
!     see if the nearest node is defined
!___________________________________________________________________
      jx=NINT((x-horizon%xmin)/horizon%dx_grid)+1
      jy=NINT((y-horizon%ymin)/horizon%dy_grid)+1
      IF( ABS(horizon%time(jx,jy)-ZNULL) > eps) THEN
         time=horizon%time(jx,jy)
if(time <= 0 .or. time >= 3) then
end if
      ELSE
!___________________________________________________________________
!        set this trace to be dead
!___________________________________________________________________
         live=.FALSE.
if(time <= 0 .or. time >= 3) then
end if
      END IF

   END IF
END IF
IF(.NOT. live) THEN
!___________________________________________________________________
!  input trace is dead. set value of time = 0.0 
!___________________________________________________________________
   time=0.0
END IF
   

END SUBROUTINE extract_gridded_horizon_value

!=========================================================================================================
!
!=========================================================================================================
SUBROUTINE eight_point_interpolation(u_interp,u,x,y,z,xmin,ymin,zmin,dx,dy,dz,          &
                                     live,left_cdp,right_cdp,line_index,nz,mx,my,mz)
! determine time value at a trace at location x,y for a gridded horzon (e.g. Earthvision format)
!_______________________________________________________________________________________________________________________
IMPLICIT NONE
REAL,           INTENT(IN)      ::  x,y,z               ! location of data to be interpolated
REAL,           INTENT(IN)      ::  dx,dy,dz            ! grid increments
REAL,           INTENT(IN)      ::  xmin,ymin,zmin      ! grid origin
INTEGER,        INTENT(IN)      ::  left_cdp,right_cdp  ! left and right cdp
INTEGER,        INTENT(IN)      ::  mx,my,mz            ! padding
INTEGER,        INTENT(IN)      ::  nz                  ! number of vertical samples
INTEGER,        INTENT(IN)      ::  line_index(-my:+my) ! indirect addressing to line index
REAL,           INTENT(IN)      ::  u(-mz:nz-1+mz,left_cdp-mx:right_cdp+mx,-my:+my)
REAL,           INTENT(OUT)     ::  u_interp            ! interpolated value
LOGICAL,        INTENT(INOUT)   ::  live                ! live trace flag
!___________________________________________________________________________________
! local arrays
!___________________________________________________________________________________
INTEGER      :: jx1,jy1,jz1       ! grid index to the upper right and above (x,y,z)
INTEGER      :: jx2,jy2,jz2       ! grid index to the upper left  and above (x,y,z)
INTEGER      :: jx3,jy3,jz3       ! grid index to the lower left  and above (x,y,z)
INTEGER      :: jx4,jy4,jz4       ! grid index to the lower right and above (x,y,z)
INTEGER      :: jx5,jy5,jz5       ! grid index to the upper right and below (x,y,z)
INTEGER      :: jx6,jy6,jz6       ! grid index to the upper left  and below (x,y,z)
INTEGER      :: jx7,jy7,jz7       ! grid index to the lower left  and below (x,y,z)
INTEGER      :: jx8,jy8,jz8       ! grid index to the lower right and below (x,y,z)
REAL         :: xksi          ! local coordinate along x within a rectangular region
REAL         :: eta           ! local coordinate along y within a rectangular region
REAL         :: zeta          ! local coordinate along z within a rectangular region
REAL         :: dx_half       ! dx/2
REAL         :: dy_half       ! dy/2
REAL         :: dz_half       ! dz/2
REAL         :: x3,y3,z3          ! coord of upper lower left corner

dx_half=dx/2.0
dy_half=dy/2.0
dz_half=dz/2.0
!___________________________________________________________________
!  find the horizon node just to the SW of the current trace for horizon 1
!___________________________________________________________________
jx3=((x-xmin)/dx)+1
jy3=((y-ymin)/dy)+1
jz3=((z-zmin)/dz)+1
!___________________________________________________________________
!     define a cube containing the point
!     number the nodes counterclockwise from the NE corner
!___________________________________________________________________
!
!       (x2,y2,z2)-----(x1,y1,z1)
!          |               |
!          |               |
!          |               |           TOP OF CUBE
!          |               |
!          |  (x,y)        |
!          |               |
!       (x3,y3,z3)----(x4,y4,z4)
!
!
!            (x,y)
!
!       (x6,y6,z6)-----(x5,y5,z5)
!          |               |
!          |               |
!          |               |          BOTTOM OF CUBE
!          |               |
!          |  (x,y)        |
!          |               |
!       (x7,y7,z7)----(x8,y8,z8)
!
!
!___________________________________________________________________
!___________________________________________________________________
jx1=jx3+1
jy1=jy3+1
jz1=jz3
jx2=jx3
jy2=jy3+1
jz2=jz3
jx4=jx3+1
jy4=jy3
jz4=jz3

jx5=jx1
jy5=jy1
jz5=jz1+1
jx6=jx2
jy6=jy2
jz6=jz2+1
jx7=jx3
jy7=jy3
jz7=jz3+1
jx8=jx4
jy8=jy4
jz8=jz4+1

!___________________________________________________________________
!        make sure all the neighboring nodes fall within the limits of defined horizon
!___________________________________________________________________
IF(jx3 < left_cdp .OR. jx1 > right_cdp) THEN
   live=.FALSE.
ELSE IF (jy3 < line_index(-my) .OR. jy2 > line_index(+my)) THEN
   live=.FALSE.
ELSE IF(jz3 < 0 .OR. jz5 > nz-1) THEN
   live=.FALSE.
ELSE
!___________________________________________________________________
!  interpolate using a bilinear (finite element) interpolator of the form t=axy+bx+cy+d
!   first step: compute the local coordinates in the rectangle (xksi,eta) where
!   -1.0 <= xksi <= +1.0 and  -1.0 <= eta <= +1.0 and -1.0 <= zeta <= +1.0
!___________________________________________________________________
   x3=xmin+(jx3-1)*dx
   y3=ymin+(jy3-1)*dy
   z3=zmin+(jz3-1)*dz
   xksi=(x-(x3+dx_half))/dx_half
   eta=(y-(y3+dy_half))/dy_half
   zeta=(z-(z3+dz_half))/dz_half
   live=.TRUE.
!___________________________________________________________________
!  second step: interpolate using the shape functions
!___________________________________________________________________
   u_interp=0.125*( (1.0+xksi)*(1+eta)*(1-zeta)*u(jx1,jy1,jz1)   &
                   +(1.0-xksi)*(1+eta)*(1-zeta)*u(jx2,jy2,jz2)   &
                   +(1.0-xksi)*(1-eta)*(1-zeta)*u(jx3,jy3,jz3)   &
                   +(1.0+xksi)*(1-eta)*(1-zeta)*u(jx4,jy4,jz4)   &
                   +(1.0+xksi)*(1+eta)*(1+zeta)*u(jx5,jy5,jz5)   &
                   +(1.0-xksi)*(1+eta)*(1+zeta)*u(jx6,jy6,jz6)   &
                   +(1.0-xksi)*(1-eta)*(1+zeta)*u(jx7,jy7,jz7)   &
                   +(1.0+xksi)*(1-eta)*(1+zeta)*u(jx8,jy8,jz8)   )
END IF

END SUBROUTINE eight_point_interpolation
!=========================================================================================================
!
!=========================================================================================================
SUBROUTINE read_picked_horizon(horizon,horizon_fn,file_command_line_arg,lu_horizon, &
                               hmin_picked,hmax_picked,input_error,stderr)
!________________________________________________________________________________________________________
! read in a horizon using command line arguments
! Authors: Kurt J. Marfurt and Tao Zhao
!________________________________________________________________________________________________________
!_____________________________________________________________________
! modules to allow checking of subroutine argument types
!_____________________________________________________________________
USE aaspi_io                    ! aaspi io module

TYPE(horizon_file)	:: horizon 
!________________________________________________________________________________________________________
! variables passed to/from a calling routine      
!________________________________________________________________________________________________________
CHARACTER (LEN=*), INTENT(IN)	:: file_command_line_arg! command line argument to call file name
INTEGER,           INTENT(IN)   :: lu_horizon		! horizon i/o unit number
CHARACTER (LEN=*), INTENT(OUT)  :: horizon_fn		! horizon file name
INTEGER,           INTENT(INOUT):: input_error		! input error counter    
REAL,              INTENT(OUT)  :: hmin_picked		! min value of picked horizon
REAL,              INTENT(OUT)  :: hmax_picked		! max value of picked horizon
INTEGER,           INTENT(IN)   :: stderr		! standard error i/o unit
!_____________________________________________________________________________________________
! local arrays
!_____________________________________________________________________________________________
REAL                            :: ZNULL              	! a 'ZNULL' value indicating an undefined horizon value
REAL                          	:: eps			! variable to prevent roundoff
INTEGER                         :: nskip                ! number of header lines to skip when reading a horizon
INTEGER      			:: ncol                 ! number of columns in the horizon file
INTEGER      			:: line_col             ! column number of line_no for interpolated horizon
INTEGER      			:: cdp_col              ! column number of cdp_no for interpolated horizon
INTEGER      			:: time_col             ! column number of time_no for interpolated horizon
CHARACTER(LEN=256)  		:: horizon_type         ! a character string used to read in the horizon type
!______________________________________________________________________
! control variables
!______________________________________________________________________
INTEGER                 :: allocation_status      	! return code for Fortran90 allocate command
INTEGER                 :: return_code            	! return code from sep calls $$$$
!REAL          		:: eps               		! fraction of znull to allow for roundoff in testing
LOGICAL       		:: file_exists      		! if == .TRUE. the file specified exists
LOGICAL                 :: test1,test2


!______________________________________________________________________
! read file descriptors from the command line
!______________________________________________________________________
horizon_fn			= aaspi_get_param(file_command_line_arg,'')
horizon_type           		= aaspi_get_param('horizon_type','gridded')
horizon%znull     		= aaspi_get_param('znull',-999999.0)
horizon%vertical_axis_down     	= (aaspi_get_param('vertical_axis_down','y')=='y')
horizon%scalet                 	= aaspi_get_param('vertical_horizon_units_scale_factor',1.000)

horizon%is_gridded_horizon	= (TRIM(horizon_type)=='gridded')
horizon%is_interpolated_horizon	= (TRIM(horizon_type)=='interpolated')
horizon%eps                    	= .0001*ABS(horizon%znull)

IF(horizon%is_interpolated_horizon) THEN
   nskip=aaspi_get_param('nskip',0)
   ncol=aaspi_get_param('ncol',3)
   IF(nskip < 0)THEN
      WRITE(stderr,*) 'Error! Number of header lines to skip in horizon files must be a non-negative integer!'
      input_error=input_error+1
      RETURN
   ENDIF
   IF(ncol < 3) THEN
      WRITE(stderr,*) 'error in routine read_picked_horizon! horizons are in a bad format! Not enough columns in the horizon files.'
      input_error=input_error+1
      RETURN
   END IF
   line_col=aaspi_get_param('line_col',1)
   cdp_col =aaspi_get_param('cdp_col',2)
   time_col=aaspi_get_param('time_col',3)
   IF((line_col > ncol).OR.(cdp_col > ncol).OR.(time_col > ncol)) THEN
      WRITE(stderr,*) 'error in routine read_picked_horizon! Column numbers of line, cdp and time picks should not be greater than the number of total columns in the horizon!'
      input_error=input_error+1
      RETURN
   ELSE IF((line_col <= 0).OR.(cdp_col <= 0).OR.(time_col <= 0)) THEN 
      WRITE(stderr,*) 'Column numbers of line, cdp and time picks should all be positive integers...'
      input_error=input_error+1
      RETURN
   END IF
END IF
WRITE(stderr,'(a25,L5)') 'interpolated_horizon',horizon%is_interpolated_horizon
WRITE(stderr,'(a25,L5)') 'gridded_horizon',horizon%is_gridded_horizon
WRITE(stderr,'(a25,L5)') 'vertical_axis_down',horizon%vertical_axis_down
WRITE(stderr,'(a25,f12.3)') 'znull',horizon%znull
!_______________________________________________________________
! check input horizons
!_______________________________________________________________
IF(horizon_fn == '') THEN
   WRITE(stderr,*) 'command line error in routine read_picked_horizon!  horizon file name must be entered after ',TRIM(file_command_line_arg)//'= on command line!'
   input_error=input_error+1
   RETURN
END IF
INQUIRE(file=horizon_fn,EXIST=file_exists)
IF( .NOT. file_exists) THEN
   WRITE(stderr,*) 'error! error in routine read_picked_horizon! horizon file = ',horizon_fn,' does not exist!'
   WRITE(stderr,*) 'check spelling and read permission'
   input_error=input_error+1
   RETURN
END IF

OPEN(lu_horizon,file=horizon_fn,status='UNKNOWN')
!_______________________________________________________________________
! read all horizon times into memory
! if they do not exist, set the output seismic data trace to be dead
! if they fall outside the bounds of the seismic survey, throw them away
! initialize all horizon time values to be undefined
!_________________________________________________________________
!_________________________________________________________________
! read in horizon 
!_________________________________________________________________
IF(horizon%is_gridded_horizon) THEN
!__________________________________________________________________
!  read in the two gridded horizons
!__________________________________________________________________
   CALL read_gridded_horizon(horizon,horizon%znull,horizon%vertical_axis_down,horizon%scalet,lu_horizon,stderr)
   hmin_picked=horizon%hmin
   hmax_picked=horizon%hmax
   WRITE(stderr,*) 'gridded horizon now in memory.'
ELSE IF(horizon%is_interpolated_horizon) THEN
!__________________________________________________________________
!  read in the two interpolated horizons
!__________________________________________________________________
   CALL read_interpolated_horizon(horizon,nskip,ncol,horizon%znull,  &
                               horizon%vertical_axis_down,horizon%scalet,line_col,cdp_col,time_col,lu_horizon,stderr)
   hmin_picked=horizon%hmin
   hmax_picked=horizon%hmax
   WRITE(stderr,*) 'interpolated horizon now in memory.'
ELSE
   WRITE(stderr,*) 'error! horizon type must be either gridded_horizon or interpolated_horizon!'
   input_error=input_error+1
   RETURN
END IF

RETURN
END SUBROUTINE read_picked_horizon

!=========================================================================================================
!
!=========================================================================================================
END MODULE aaspi_horizon_util
