PROGRAM EM
!______________________________________________________________________
! Gaussian Mixture Models with EM             
! written in order to compare to more useful som3d and gtm3d algorithms
!______________________________________________________________________
! USAGE:
! Compute the probability of given datapoint in certain cluster
!______________________________________________________________________
! Attribute-Assisted Seismic Processing and Interpretation (AASPI)
! The University of Oklahoma
! Norman, OK, USA
! Xuan Qi,
!______________________________________________________________________
! access the necessary modules from the library
!______________________________________________________________________
USE aaspi_io                    ! aaspi io module
USE mpi_send_receive_array_util ! aaspi mpi send and receive utility module
USE string_util                 ! aaspi string utility module
USE aaspi_file_util             ! aaspi file initialization and updating utilities
USE mutes                       ! mute application utilities
USE get_time                    ! wall time clock utilities
USE data_util
USE prog_util                   ! program initialization utilities (licensing data, mpi communication, etc.)
USE aaspi_horizon_util          ! aaspi read horizon module 
USE aaspi_cluster_util		! aaspi clustering utilities
!______________________________________________________________________
! initialize all variable types to be NONE (undefined)
!______________________________________________________________________
IMPLICIT NONE
!______________________________________________________________________
! parameters describing the input data volume
!______________________________________________________________________
INTEGER                 :: nt_in,ncdp_in,nline_in                     ! number of samples along each axes
REAL                    :: o1_in,o2_in,o3_in                          ! first sample of each axis
REAL                    :: dt,d2,d3                                   ! sample increment of each axis
CHARACTER (LEN=80)      :: unit1,unit2,unit3                          ! units along each axis
CHARACTER (LEN=80)      :: label1,label2,label3                       ! label of each axis
REAL                    :: dcdp,dline                                 ! cdp and line separation in m or ft (measured for a given single change in index no.)
INTEGER                 :: esize                                      ! size in bytes of each data word
REAL                    :: inline_azimuth                             ! azimuth in degrees from North of increasing CDP indices
REAL                    :: crossline_azimuth                          ! azimuth in degrees from North of increasing Line_no indices
INTEGER                 :: ndimension                                 ! number of dimensions on the input data
INTEGER                 :: fixed_scale                                ! if == 1, fixed scale, if == 0, auto scale the data when plotting
CHARACTER (LEN=256)     :: unique_project_name                        ! project name used to facilitate links to subsequent processes
CHARACTER (LEN=256)     :: suffix                                     ! suffix typically attached to a file name to allow one to differentiate run parameters
CHARACTER (LEN=256)     :: colorbar_fn                                ! color bar file name
CHARACTER (LEN=256)     :: title                                      ! title of the input data volume
INTEGER                 :: n_keys                                     ! number of header keys

INTEGER      :: nattr	       ! number of input attribute volumes
INTEGER      :: nclusters       ! number of output clusters
INTEGER      :: ntrace         ! number of traces
INTEGER      :: ncdp_out         ! number of output seismic CDPs
INTEGER      :: nline_out       ! number of output seismic lines
INTEGER      :: first_line_in   ! first line in input data  
INTEGER      :: last_line_in    ! last line in input data  
INTEGER      :: first_line_out  ! first line to be processed
INTEGER      :: last_line_out   ! last line to be processed
INTEGER      :: begin_line_out  ! index of first line to be processed (increments of 1, starting at 1)
INTEGER      :: end_line_out    ! index of last line to be processed (increments of 1, starting at 1)
INTEGER      :: first_cdp_in    ! first cdp in input data  
INTEGER      :: last_cdp_in     ! last cdp in input data  
INTEGER      :: first_cdp_out   ! first cdp to be processed
INTEGER      :: last_cdp_out    ! last cdp to be processed
INTEGER      :: begin_cdp_out   ! index of first CDP to be processed (increments of 1, starting at 1)
INTEGER      :: end_cdp_out     ! index of last CDP to be processed (increments of 1, starting at 1)
INTEGER      :: cdp_order       ! = +1 if d2 > 0, = -1 if d2 < 0.
INTEGER      :: line_order      ! = +1 if d3 > 0, = -1 if d3 < 0.
INTEGER      :: line_decimation	! skip line increment for training
INTEGER      :: cdp_decimation 	! skip cdp increment for training
INTEGER      :: vertical_sample_decimation	! skip sample decimation to extract training vectors
INTEGER      :: jt_start			! first data sample in memory
INTEGER      :: jt_end  			! first data sample in memory
INTEGER      :: jt				! loop index

REAL         :: o1_out,o2_out,o3_out ! first values of o1,o2,o3 axes of output data 
REAL         :: d2_abs,d3_abs    ! ABS(d2), ABS(d3)
REAL         :: t_start        ! start time of analysis
REAL         :: t_end 	       	 ! end time of analysis    
REAL         :: t_min          ! mimimum time of input data                  
REAL         :: t_max 	         ! maximum time of input data
INTEGER      :: first_sample_in   ! first sample of input data
INTEGER      :: last_sample_in   ! last sample of input data 
INTEGER      :: nt_out        ! number of output samples 
INTEGER      :: j_training_vector 	! training vector index              
INTEGER(KIND=8)      :: n_training_vectors	! actual number of training vectors
INTEGER      :: max_training_vectors	! maximum number of training vectors

INTEGER      :: nrow           ! order of the rcov matrix
REAL	     :: nstd           ! the number of standard deviations to span in initialization
INTEGER(KIND=8)      :: n_data_vectors    	! Number of live data vectors
INTEGER      :: n_data_training_iterations 	! Number of iterations to train the attribute volume
INTEGER      :: jiter          ! current iteration no.
INTEGER      :: ierror	       ! i/o error counter
INTEGER      :: count	       ! counter of number of traces
INTEGER      :: temp           ! temporary integer number
INTEGER      :: nlive          ! number of live traces
INTEGER      :: jattr,kattr     ! attribute loop indices 
INTEGER      :: jcluster        ! cluster loop indices 
INTEGER      :: jline           ! line index
INTEGER      :: jcdp            ! cdp index
INTEGER      :: jp              ! prototype vector index
INTEGER      :: jsamp		! sample index
REAL         :: fraction, center! values in intializing cluster centers
REAL (KIND=8):: total_dist2     ! sum of all distances, which should decrease with iterations
REAL (KIND=8):: total_dist2_previous     ! sum of all distances from previous iteration
REAL (KIND=8):: convergence     ! convergence fraction from previous iteration
CHARACTER (LEN=256) :: command_line_arg 		!for different Input files
CHARACTER (LEN=256) :: command_line_arg1 		!for different Input files
CHARACTER (LEN=3)   :: string3					! 
LOGICAL		    :: verbose				! if == .TRUE., generate verbose output
REAL                :: cluster_number_dummy		! place holder in calling routine
!______________________________________________________________________
! control variables 
!______________________________________________________________________
INTEGER               :: allocation_status  	! return code for Fortran90 allocate command   
INTEGER               :: input_error        	! input error counter   
INTEGER               :: return_code        	! return code from sep calls $$$$
INTEGER,   PARAMETER  :: stderr=0	  	! standard error unit 
INTEGER               :: info                   ! return code on LAPACK subroutines
CHARACTER (LEN=1), PARAMETER :: qu='"'      	! quotation mark
LOGICAL       :: is_gridded_horizon           ! if == .TRUE. horizons are formatted by x/y coordinates
LOGICAL       :: is_interpolated_horizon      ! if == .TRUE. horizons are formatted by line/cdp numbers
REAL          :: znull             ! value indicating an undefined horizon value
REAL          :: eps               ! fraction of znull to allow for roundoff in testing
LOGICAL       :: file_exists      ! if == .TRUE. the file specified exists
!______________________________________________________________________
! allocatable arrays
!______________________________________________________________________
REAL,          ALLOCATABLE, DIMENSION(:,:)   :: data_vector_slave	! nattr input data vector on slave   
REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: data_vector_master	! nattr input data vectors on master   
REAL,          ALLOCATABLE, DIMENSION(:)     :: cluster_number_slave	! cluster number on slave   
REAL,          ALLOCATABLE, DIMENSION(:,:)   :: cluster_number_master	! cluster number on master   
REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: d_vector          	! an input line of data vectors
REAL,          ALLOCATABLE, DIMENSION(:,:)   :: d_training_vector    	! decimated vector used in training
LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: live			! if == .TRUE., this is a live trace
REAL,          ALLOCATABLE, DIMENSION(:)     :: d_mean       		! the mean value of each attribute
REAL,          ALLOCATABLE, DIMENSION(:,:)   :: rcov         		! covariance matrix amongst the attributes
REAL,          ALLOCATABLE, DIMENSION(:)     :: sigma        		! SQRT  of diagonal of covariance matrix
REAL,          ALLOCATABLE, DIMENSION(:,:)   :: rhs          		! right hand side in matrix inversion (identity matrix on input, inverse of rcov on output)
REAL (KIND=4), ALLOCATABLE, DIMENSION(:,:)   :: cov_inv      		! inverse of the covariance matrix used in computing Mahalanobis distance
REAL (KIND=4), ALLOCATABLE, DIMENSION(:,:)   :: mu           		! centroids of each cluster
INTEGER      , ALLOCATABLE, DIMENSION(:)     :: nvector         	! number of vectors assigned to each cluster
REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: attribute_cluster_sum   ! sum of attributes of each cluster in attribute space
INTEGER,        ALLOCATABLE, DIMENSION(:,:)  :: top_sample     		! first sample to be analyzed
INTEGER,        ALLOCATABLE, DIMENSION(:,:)  :: bottom_sample   	! last sample to be analyzed
CHARACTER (LEN=256), ALLOCATABLE, DIMENSION(:)  :: input_title          ! input file titles

!___________________________________________________________________
! AASPI-format file names
!___________________________________________________________________
CHARACTER (LEN=256), ALLOCATABLE, DIMENSION(:)  :: input_fn     ! input file name
CHARACTER (LEN=256) :: d_mean_fn                 		! attribute mean file name
CHARACTER (LEN=256) :: d_std_fn                 		! attribute standard file name
CHARACTER (LEN=256) :: covariance_fn           			! covariance matrix file name
CHARACTER (LEN=256) :: cluster_number_fn         		! kmeans cluster number file name
CHARACTER (LEN=256) :: error_fn                  		! output error file name  
!_________________________________________________________________
! trace header variables
!______________________________________________________________________
INTEGER,ALLOCATABLE :: header(:,:)           ! input trace header block
INTEGER             :: live_in_key           ! live input trace key
INTEGER             :: sx_in_key             ! shot x coordinate key
INTEGER             :: sy_in_key             ! shot y coordinate key
INTEGER             :: gx_in_key             ! group x coordinate key
INTEGER             :: gy_in_key             ! group y coordinate key
INTEGER             :: scalco_in_key         ! scale xy coordinate key
INTEGER             :: live_out_key          ! live output trace key
INTEGER             :: sx_out_key            ! shot x coordinate key
INTEGER             :: sy_out_key            ! shot y coordinate key
INTEGER             :: gx_out_key            ! group x coordinate key
INTEGER             :: gy_out_key            ! group y coordinate key
INTEGER             :: scalco_out_key        ! scale xy coordinate key
INTEGER             :: next_header_in        ! location of next header to be read in 
INTEGER             :: next_header_out       ! location of next header to be written out  
INTEGER             :: line_key,cdp_key      ! line and cdp no. keys
INTEGER             :: naxis                 ! number of header data axes (needs to be = 2!)
INTEGER             :: cdp_x_key             ! cdp x coordinate key
INTEGER             :: cdp_y_key             ! cdp y coordinate key
INTEGER             :: line_no_key           ! line number key         
INTEGER             :: cdp_no_key            ! cdp number key 
!_____________________________________________________________________
! variables associated with MPI
!_____________________________________________________________________
INTEGER,PARAMETER  :: max_process=999   ! maximum number of processes (used in constructing tags)
INTEGER            :: n_process         ! total number of processes
INTEGER            :: n_slave           ! number of slave processes
INTEGER            :: this_process      ! the rank of this process
INTEGER            :: mpi_rc            ! mpi return code
INTEGER            :: j_process         ! loop index
INTEGER            :: slave_cnt         ! process number of slave(s)
INTEGER            :: slave_process     ! process number of slave(s)
INTEGER            :: master_process=0  ! master process no. (always 0)
INTEGER            :: msg_tag
INTEGER            :: ibuf
INTEGER            :: nwords            ! number of words to be transferred
LOGICAL            :: master            ! if .TRUE., this is the master process (this_process == 0)
LOGICAL            :: slave             ! if .TRUE., this is a slave process (a master is also its own slave if n_process=1)
LOGICAL            :: use_mpi           ! if = .TRUE., run under MPI
INTEGER            :: more_traces       ! more_traces = 0 finished, = 1 more traces coming 
INTEGER            :: num_sent          ! number of traces sent
INTEGER            :: mpi_comm=MPI_COMM_WORLD  ! set up default MPI communicator as world
CHARACTER (LEN=256):: process_name      ! name of each process
INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: process_status
INTEGER            :: more_traces_tag = 1
INTEGER            :: ready_tag = 2
INTEGER            :: cdp_tag = 3
INTEGER            :: data_in_tag = 10
INTEGER            :: cluster_number_tag = 11
INTEGER            :: live_flag
!______________________________________________________________________
! variables defining the grid of the input seismic data volume
!______________________________________________________________________
INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: cdp_no               ! cdp number
INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: line_no              ! current line number
REAL,    ALLOCATABLE, DIMENSION(:,:)    :: xcoord               ! current trace x
REAL,    ALLOCATABLE, DIMENSION(:,:)    :: ycoord               ! current trace y
!______________________________________________________________________
! variables defining bounding horizons
! horizon will have the following two possible formats:
! 1) an interpolated horizon with time picks given for every cdp_no and line_no (e.g. seisx format)
! 2) a gridded horizon with time picks on an uninterpolated NS by EW rectilinear grid (e.g. from Petrel: earthvision format)
!______________________________________________________________________
TYPE(horizon_file) :: horizon_upper             ! upper horizon type
TYPE(horizon_file) :: horizon_lower             ! lower horizon type
!______________________________________________________________________
! variables used in defining the horizon picks
!______________________________________________________________________
REAL         :: scalet                  ! scale factor to convert horizon units to those of the seismic data      
REAL         :: w1         		! weight for interpolation1
REAL         :: w2         		! weight for interpolation2
COMPLEX	     :: arg  	   		! complex argument for horizon (0.0,+omega*t)
INTEGER      :: nskip      		! number of header lines to skip when reading a horizon
INTEGER      :: ncol  	   		! number of columns in the horizon file
INTEGER      :: line_col   		! column number of line_no for interpolated horizon 
INTEGER      :: cdp_col    		! column number of cdp_no for interpolated horizon 
INTEGER      :: time_col   		! column number of time_no for interpolated horizon
CHARACTER(LEN=256)  :: horizon_type     ! a character string used to read in the horizon type
REAL         :: tstart     		! start time with respect to the picked horizon (positive = below)
REAL         :: tend       		! end time with respect to the picked horizon (positive = below)
LOGICAL      :: vertical_axis_down      ! if == .TRUE. the vertical axis is positive down
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: t_horizon_upper  ! two-way travel times of the upper horizon for inlines and crosslines
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: t_horizon_lower  ! two-way travel times of the lower horizon for inlines and crosslines
INTEGER      :: t_horizon_count         ! counter for samples within operation window
LOGICAL      :: use_horizon             ! use horizon selection
CHARACTER(LEN=1)   :: use_horizon_char  ! the use_horizon in character 
INTEGER, PARAMETER :: lu_horizon_upper=98  ! horizon_upper file logical unit
INTEGER, PARAMETER :: lu_horizon_lower=97  ! horizon_lower file logical unit
CHARACTER (LEN=256):: horizon_upper_fn	   ! horizon_upper file name
CHARACTER (LEN=256):: horizon_lower_fn	   ! horizon_lower file name
!_____________________________________________________________________
! timing and op count arrays. make these double precision to avoid truncation errors.
!_____________________________________________________________________
REAL(KIND=8) :: last_time(0:20)        			! last time clock was called for a given operation.
REAL(KIND=8) :: total_time(0:20)=1.e-3 			! sum of time for a given operation.
REAL(KIND=8) :: ops(0:20)=0.           			! op count.
INTEGER, PARAMETER	:: NOT_DEFINED=-99999		! an integer znull value
!_______________________________________________________________
! initialize aaspi_io
!_______________________________________________________________
CALL aaspi_initialize()
!_______________________________________________________________
! read in command line arguments
! if MPI is NOT going to be used, then there are no other processors to read it!
!_______________________________________________________________
verbose=(aaspi_get_param('verbose','n') == 'y')
use_mpi=(aaspi_get_param('use_mpi','n') == 'y')
!_______________________________________________________________
! determine if we will run under MPI
! if so, determine the number of processors,
! this processor number, and whether it is a master or slave
!_______________________________________________________________
CALL invoke_mpi(n_process,max_process,n_slave,this_process,slave_process,         &
                master,slave,use_mpi,stderr)
IF(use_mpi) THEN
   CALL MPI_Bcast(verbose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_rc)
END IF
IF(master) THEN
!________________________________________________________________________________________________
!  check whether the license is valid. (times are hardcoded in the parameter statements above!)
!________________________________________________________________________________________________
   input_error=license_expiration_prog('kmeans3d',release_date,expiration_month,expiration_year,stderr)
   error_fn=aaspi_get_param('error_fn','kmeans3d_error_file')
END IF

IF(master) THEN
   input_error=0
!___________________________________________________________
!   read output auxilliary data set names from the command line
!_______________________________________________________________
   nattr=aaspi_get_param('nattr', NOT_DEFINED)   
   IF(nattr == NOT_DEFINED) THEN
      WRITE(stderr,*) 'error! number of input attribute volumes must follow nattr= on the command line!'
      STOP 666
   END IF
!___________________________________________________________
!  read in the number of clusters                                 
!_______________________________________________________________
   nclusters=aaspi_get_param('nclusters',3)
   IF(nclusters <= 1)  THEN
      WRITE(stderr,*) 'error! must have at least 2 clusters following nclusters= on the command line!'
      STOP 666
   END IF
   WRITE(stderr,*) 'nclusters = ',nclusters
!_______________________________________________________________
!  allocate memory for input attribute file names
!_______________________________________________________________
   ALLOCATE(input_fn(nattr),input_title(nattr),STAT=allocation_status)
   IF(allocation_status /= 0) THEN
      WRITE(stderr,*) 'error! could not allocate memory for ',nattr,' input attribute file names!'
      STOP 666
   END IF

!_______________________________________________________________
!  read in the input file names
!_______________________________________________________________
   CALL read_multiple_input_file_names(input_fn,nattr,input_error,stderr)
!_______________________________________________________________
!  read in output file names
!_______________________________________________________________
   cluster_number_fn=aaspi_get_param('cluster_number_fn','')
   IF(cluster_number_fn == '') THEN
     WRITE(stderr,*) 'command line error! valid output file name must follow cluster_number_fn'
     input_error=input_error+1
   END IF
   covariance_fn=aaspi_get_param('covariance_fn','covariance_kmeans.txt')
   d_mean_fn=aaspi_get_param('d_mean_fn','attribute_mean_kmeans.txt')
   d_std_fn=aaspi_get_param('d_std_fn','attribute_std_kmeans.txt')
END IF
IF(use_mpi) THEN
   CALL check_for_errors_found_only_on_master(this_process,input_error,use_mpi,master,               &
                  'one or more input files do not exist or cannot be read!',stderr,verbose,error_fn)
ELSE IF(input_error /= 0) THEN
   WRITE(stderr,*) 'Program terminated due to ',input_error,' input errors!'
   STOP 666
END IF
   
IF(master) THEN

   CALL read_params
!_______________________________________________________________
!  read in the titles of the input files
!_______________________________________________________________
   DO jattr=1,nattr
      input_title(jattr)=aaspi_get_hist(input_fn(jattr),'title','UNKNOWN',return_code)
   END DO

!_______________________________________________________________
!  compute decimated data vectors for training
!_______________________________________________________________
   write(0,*) 'compute max_training_vectors. line_decimation,cdp_decimation,vertical_sample_decimation ',line_decimation,cdp_decimation,vertical_sample_decimation
   max_training_vectors=((nline_out-1)/line_decimation+1)*((ncdp_out-1)/cdp_decimation+1)*((nt_out-1)/vertical_sample_decimation+1)
   write(0,*) 'max_training_vectors = ',max_training_vectors
   write(0,*) 'maximum number of training vectors ', max_training_vectors
END IF
IF(use_mpi) THEN
   CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
!_______________________________________________________________
!  broadcast relevant input parameters from the master to all the slaves.
!_______________________________________________________________
   CALL MPI_Bcast(nattr,                1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(nclusters,            1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(jt_start,             1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(jt_end,               1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(max_training_vectors, 1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(begin_line_out,       1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(end_line_out,         1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(nattr,                1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(nt_out,               1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(ncdp_in,              1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(nline_in,             1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(ncdp_out,             1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(nline_out,            1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(begin_cdp_out,        1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(end_cdp_out,          1, MPI_INTEGER , 0, MPI_COMM_WORLD, mpi_rc)
END IF

IF(master) THEN
!_______________________________________________________________
!  read in the number of AASPI header keys from the hff file   		
!  currently, we only need the live trace (trid) trace header value
!_______________________________________________________________
   n_keys=aaspi_get_number_keys(input_fn(1),return_code)
   IF(return_code /= 0) STOP "error! cannot read in number of keys from the hff (.H@@) file!"
   WRITE(stderr,*) 'n_keys read in = ',n_keys

   CALL aaspi_get_key_index(input_fn(1),"trid",live_in_key,return_code)
   IF(return_code /= 0) STOP "trouble obtaining trid key index"
   
   CALL aaspi_get_key_index(input_fn(1),"gx",gx_in_key,return_code)      
   IF(return_code /= 0) THEN
        gx_in_key=0
   END IF

   CALL aaspi_get_key_index(input_fn(1),"gy",gy_in_key,return_code)      
   IF(return_code /= 0) THEN
        gy_in_key=0
   END IF

   CALL aaspi_get_key_index(input_fn(1),"sx",sx_in_key,return_code)      
   IF(return_code /= 0) THEN
        sx_in_key=0
   END IF

   CALL aaspi_get_key_index(input_fn(1),"sy",sy_in_key,return_code)      
   IF(return_code /= 0) THEN
        sy_in_key=0
   END IF

   CALL aaspi_get_key_index(input_fn(1),"scalco",scalco_in_key,return_code)      
   IF(return_code /= 0) THEN
        scalco_in_key=0
   END IF

   CALL aaspi_get_key_index(input_fn(1),'cdp_x',cdp_x_key,return_code)
   IF(return_code /= 0) THEN
      WRITE(stderr,*) 'error! cdp_x key index does not exist on input file!'
      input_error=input_error+1
   END IF

   CALL aaspi_get_key_index(input_fn(1),'cdp_y',cdp_y_key,return_code)
   IF(return_code /= 0) THEN
      WRITE(stderr,*) 'error! cdp_y key index does not exist on input file!'
      input_error=input_error+1
   END IF

   CALL aaspi_get_key_index(input_fn(1),'cdp_no',cdp_no_key,return_code)
   IF(return_code /= 0) THEN
      WRITE(stderr,*) 'error! cdp_no key index does not exist on input file!'
      input_error=input_error+1
   END IF

   CALL aaspi_get_key_index(input_fn(1),'line_no',line_no_key,return_code)
   IF(return_code /= 0) THEN
      WRITE(stderr,*) 'error! line_no key index does not exist on input file!'
      input_error=input_error+1
   END IF

   WRITE(stderr,*) 'live_in_key',live_in_key
   WRITE(stderr,*) 'cdp_x_key',cdp_x_key
   WRITE(stderr,*) 'cdp_y_key',cdp_y_key
   WRITE(stderr,*) 'cdp_no_key',cdp_no_key
   WRITE(stderr,*) 'line_no_key',line_no_key
!_______________________________________________________________
! write out history for projected means
!_______________________________________________________________
write(0,*) 'init cluster_number_fn ',cluster_number_fn
   title='k-means cluster number'
   colorbar_fn='rainbow.sep'
   fixed_scale=1
   CALL init_3d_history_file_2015(nt_out,ncdp_out,nline_out,                    &
                                  dt,d2,d3,                                     &
                                  t_start,o2_out,o3_out,                        &
                                  label1,label2,label3,                         &
                                  unit1,unit2,unit3,                            &
                                  dcdp,dline,esize,                             &
                                  colorbar_fn,fixed_scale,unique_project_name,  &
                                  suffix,title,input_fn(1),cluster_number_fn,stderr)
   CALL aaspi_put_hist(cluster_number_fn,'min_amplitude',0.0,return_code)
   CALL aaspi_put_hist(cluster_number_fn,'max_amplitude',REAL(nclusters),return_code)


END IF

write(0,*) 'after intialization'
IF(master) THEN
!______________________________________________________________________
!  allocate arrays needed only on master
!______________________________________________________________________
   ALLOCATE(data_vector_master(nattr,jt_start:jt_end,begin_cdp_out:end_cdp_out),&
            d_vector(first_sample_in:last_sample_in,1:ncdp_in,nattr),		&
            cluster_number_master(jt_start:jt_end,1:ncdp_in),			&
            live(1:ncdp_in,1:nline_in),						&
    	    d_mean(nattr),							&
            rcov(nattr,nattr),                     				&
            rhs(nattr,nattr),                     				&
            t_horizon_upper(1:ncdp_in,1:nline_in), 				&
            t_horizon_lower(1:ncdp_in,1:nline_in), 				&
            header(1:n_keys,1:ncdp_in),						&
            line_no(1:ncdp_in,1:nline_in),                                      &
            cdp_no(1:ncdp_in,1:nline_in),                                       &
            xcoord(1:ncdp_in,1:nline_in),                                       &
            ycoord(1:ncdp_in,1:nline_in),                                       &
            STAT=allocation_status)
   IF(allocation_status /= 0) THEN
      WRITE(stderr,*) this_process,': error! could not allocate arrays on master!'
      STOP 666
   END IF
   WRITE(0,*) 'after allocate on master'
END IF

IF(slave) THEN
!______________________________________________________________________
!  allocate arrays needed only on slaves
!______________________________________________________________________
   ALLOCATE(data_vector_slave(nattr,jt_start:jt_end),		&
	    cluster_number_slave(jt_start:jt_end),		&
            STAT=allocation_status)
   IF(allocation_status /= 0) THEN
      WRITE(stderr,*) this_process,': error! could not allocatedata_vector_slave on slaves! '
      STOP 666
   END IF
END IF
IF(use_mpi) THEN
   CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
END IF
!______________________________________________________________________
!  allocate arrays needed on both master and slaves
!______________________________________________________________________
ALLOCATE(top_sample(1:ncdp_in,1:nline_in),				&
         bottom_sample(1:ncdp_in,1:nline_in),				&
         d_training_vector(nattr,max_training_vectors),			&
         cov_inv(nattr,nattr),						&
         sigma(nattr),							&
         mu(nattr,nclusters),						&
         attribute_cluster_sum(nattr,nclusters),			&
         nvector(nclusters),						&
         STAT=allocation_status)
IF(allocation_status /= 0) THEN
   WRITE(stderr,*) this_process,' : could not allocate required memory on both master and slaves'
   STOP 'program aborted!'
ELSE
   WRITE(stderr,*) 'arrays successfully allocated'
END IF
IF(use_mpi) THEN
   CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
END IF
IF(master) THEN
   write(0,*) 'after allocate on master AND slaves'
END IF

!==========================================================================================================
! BEGIN DEFININTION OF HORIZONS
!==========================================================================================================
!_____________________________________________________________________________________________
! read in the 3d geometry as defined in the seismic trace headers
!_____________________________________________________________________________________________
IF(master) THEN
   CALL read_3d_geometry(input_fn(1),xcoord,ycoord,line_no,cdp_no,live,ncdp_in,nline_in,input_error,stderr)
END IF
CALL check_for_errors_found_only_on_master(this_process,input_error,use_mpi,master,               &
                  'errors in reading geometry!',stderr,verbose,error_fn)
IF(master) WRITE(stderr,*) 'geometry successfully read in'

IF(master) THEN
   use_horizon=(aaspi_get_param('use_horizon','n') == 'y')
   read_in_horizon_time: IF(use_horizon) THEN
!________________________________________________________________________________________________
!     read in user-defined horizon files in either gridded or interpolated format.
!     Note that his routine uses command line arguments of the form:
!
!     upper_horizon_fn=
!     horizon_type=
!     znull=
!     vertical_axis_down=
!     scalet=
!     ncol=
!     line_col=
!     cdp_col=
!     time_col=
!
!     in general, these values as well as use_horizon= will be written by a GUI calling subroutine HorizonTabItem.cpp 
!________________________________________________________________________________________________
      CALL read_horizon(t_horizon_upper,xcoord,ycoord,cdp_no,line_no,live,ncdp_in,nline_in,horizon_upper,&
                        'upper_horizon_fn',lu_horizon_upper,                       				&
                        input_error,stderr)
!________________________________________________________________________________________________
!     test that the upper horizon picks fall within the operation window range
!________________________________________________________________________________________________
      IF(input_error == 0 .AND. horizon_upper%hmin < t_start .OR. horizon_upper%hmax > t_end) THEN
         WRITE(stderr,*) 'error! upper horizon falls beyond the limits of the operation window!'
         WRITE(stderr,'(a30,f12.4)') 't_start',t_start,'t_end',t_end,&
                                     'upper horizon top',horizon_upper%hmin,'shallower horizon bottom',horizon_upper%hmax
         input_error=input_error+1
      END IF

      CALL read_horizon(t_horizon_lower,xcoord,ycoord,cdp_no,line_no,live,ncdp_in,nline_in,horizon_lower,&
                        'lower_horizon_fn',lu_horizon_lower,                       				&
                        input_error,stderr)
!________________________________________________________________________________________________
!     test that the lower horizon picks fall within the operation window range
!________________________________________________________________________________________________
      IF(input_error == 0 .AND. horizon_lower%hmin < t_start .OR. horizon_lower%hmax > t_end) THEN
         WRITE(stderr,*) 'error! lower horizon falls beyond the limits of the operation window!'
         WRITE(stderr,'(a30,f12.4)') 't_start',t_start,'t_end',t_end,&
                                     'lower horizon top',horizon_lower%hmin,'shallower horizon bottom',horizon_lower%hmax
         input_error=input_error+1
      END IF
      IF(input_error == 0) THEN
!________________________________________________________________________________________________
!        convert time to samples
!________________________________________________________________________________________________
         DO jline=begin_line_out,end_line_out
            top_sample(begin_cdp_out:end_cdp_out,jline)   =NINT(t_horizon_upper(begin_cdp_out:end_cdp_out,jline)/dt)
            bottom_sample(begin_cdp_out:end_cdp_out,jline)=NINT(t_horizon_lower(begin_cdp_out:end_cdp_out,jline)/dt)
         END DO
      END IF
   ELSE
!________________________________________________________________________________________________
!     no horizons read in. simply crop the data volume between time slices at t_start and t_end
!________________________________________________________________________________________________
      write(0,*) 'set top_sample and bottom_sample to be fixed!. jt_start,jt_end ',jt_start,jt_end
      top_sample(:,:)=jt_start
      bottom_sample(:,:)=jt_end
   END IF read_in_horizon_time
END IF

CALL check_for_errors_found_only_on_master(this_process,input_error,use_mpi,master,               &
                  'errors in reading horizons!',stderr,verbose,error_fn)
IF(master) WRITE(stderr,*) 'horizons successfully read in!'
!==========================================================================================================
! END DEFININTION OF HORIZONS
!==========================================================================================================


IF(use_mpi) THEN
   CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
END IF
IF(use_mpi) THEN
!_______________________________________________________________
!  broadcast relevant input parameters from the master to all the slaves.
!_______________________________________________________________
   CALL MPI_Bcast(top_sample,          ncdp_in*nline_in, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_rc)
   CALL MPI_Bcast(bottom_sample,       ncdp_in*nline_in, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_rc)
END IF   

compute_cov_inv: IF(master) THEN
!___________________________________________________________________
!    read the  data volume into memory.
!     read in the headers in the same manner.
!     header         = an array of all headers present
!     n__keys        = the number of headers in the .hff input data file.
!     next_header_in = the sequential trace position in the trace header file.
!___________________________________________________________________
   WRITE(stderr,*) 'Read in the data volume and calculate the mean of each attribute'
!___________________________________________________________________
! initialize
!___________________________________________________________________
   CALL init_time(last_time(0))
   CALL init_time(last_time(1)) 
!___________________________________________________________________
!  read through the data and compute the mean for each attribute
!___________________________________________________________________
   CALL compute_mean_values(top_sample,bottom_sample,live,input_fn,d_mean,      &
                            first_sample_in,last_sample_in,                     &
                            ncdp_in,begin_cdp_out,end_cdp_out,begin_line_out,end_line_out,      &
                            nattr,n_data_vectors,verbose,stderr)
   CALL add_time(total_time(1),last_time(1))

   WRITE(stderr,'(2a15)') 'attribute no.','mean'
   WRITE(stderr,'(i15,e15.3)') (jattr,d_mean(jattr),jattr=1,nattr)
   CALL init_time(last_time(2))
!____________________________________________________________________
!  compute the covariance matrix
!____________________________________________________________________
   CALL init_time(last_time(2))
   CALL compute_covariance_matrix(top_sample,bottom_sample,live,input_fn,d_mean,rcov,   &
                                  first_sample_in,last_sample_in,                       &
                                  ncdp_in,begin_cdp_out,end_cdp_out,begin_line_out,end_line_out,        &
                                  nattr,verbose,stderr)
!____________________________________________________________________
!  Extract the diagonal elements
!____________________________________________________________________
   DO jattr=1,nattr
      sigma(jattr)=SQRT(rcov(jattr,jattr))
   END DO
!____________________________________________________________________
!  divide each row of the rcov matrix by the diagonal to normalize the data
!____________________________________________________________________
   DO jattr=1,nattr
      DO kattr=1,nattr
         rcov(jattr,kattr)=rcov(jattr,kattr)/(sigma(kattr)*sigma(jattr))
      END DO
   END DO
   write(0,*) 'after  dividing by sigma'
   WRITE(stderr,'(a10,8i12)') 'row',(kattr,kattr=1,nattr)
   DO jattr=1,nattr
      WRITE(stderr,'(i10,8e12.3)') jattr,(rcov(jattr,kattr),kattr=1,nattr)
   END DO

   CALL print_statistics_for_excel_plotting(d_mean_fn,d_std_fn,covariance_fn,'','',                          &
                                            input_title,d_mean,sigma,rcov,rcov,rcov,nattr,stderr)

   CALL add_time(total_time(2),last_time(2))
!_____________________________________________________________________
! compute the inverse of the covariance matrix
!_____________________________________________________________________
   CALL init_time(last_time(3))
   CALL compute_covariance_matrix_inverse(rcov,cov_inv,nattr,verbose,stderr)
   CALL add_time(total_time(3),last_time(3))
END IF compute_cov_inv

construct_training_data: IF(master) THEN
!____________________________________________________________________
!  Decimate the data for training
!  Reset file positions to the beginning of each attribute files
!____________________________________________________________________
   DO jattr=1,nattr
      CALL aaspi_seek(input_fn(jattr),0,SEEK_SET,return_code)
   END DO
!____________________________________________________________________
!  read in the data
!  no need to remove mean or scale at this point
!  decimate input and store result in d_vector
!____________________________________________________________________
   write(0,*) 'line_decimation,cdp_decimation,vertical_sample_decimation ',line_decimation,cdp_decimation,vertical_sample_decimation
   n_training_vectors=0  
   decimation_loop_over_lines: DO jline=1,end_line_out
      IF(MOD(jline,10) == 0) WRITE(stderr,'(a,3i10)') 'decimate data for training: jline,end_line_out',jline,end_line_out
     
      DO jattr=1,nattr 
        CALL aaspi_read(input_fn(jattr),d_vector(:,:,jattr),return_code)
      END DO
      IF( MOD( (jline-begin_line_out),line_decimation ) == 0) THEN
         decimation_loop_over_cdps: DO jcdp=begin_cdp_out,end_cdp_out,cdp_decimation
           IF(live(jcdp,jline)) THEN   
              DO jt=top_sample(jcdp,jline),bottom_sample(jcdp,jline),vertical_sample_decimation
                 n_training_vectors=n_training_vectors+1
                 DO jattr=1,nattr
                    d_training_vector(jattr,n_training_vectors)=d_vector(jt,jcdp,jattr)
                 END DO
              END DO    
            END IF
          END DO decimation_loop_over_cdps
       END IF
   END DO decimation_loop_over_lines
END IF construct_training_data
CALL add_time(total_time(3),last_time(3))
training_if_block: IF(master) THEN
!____________________________________________________________________
!    initialize cluster centers
!____________________________________________________________________  
write(0,*) 'initialize mu'
write(0,'(4a12)') 'jcluster','fraction','jattr','mu'
   center=REAL(nclusters-1)/2.0
   DO jcluster=1,nclusters
      fraction=(REAL(jcluster-1)-center)/center
      DO jattr=1,nattr
         mu(jattr,jcluster)=d_mean(jattr)+fraction*sigma(jattr) 
write(0,'(i12,e12.5,i12,e12.5)') jcluster,fraction,jattr,mu(jattr,jcluster)
      END DO
   END DO


   training_loop: DO jiter=1,n_data_training_iterations
      total_dist2=0.0
!____________________________________________________________________
!     loop over training data and find the cluster means
!____________________________________________________________________  
      attribute_cluster_sum(:,:)=0.0
      nvector(:)=0
      loop_over_training_vectors: DO j_training_vector=1,n_training_vectors
         CALL go_cluster(cov_inv,d_training_vector(:,j_training_vector),mu,attribute_cluster_sum,	&
                         nattr,nclusters,cluster_number_dummy,nvector,total_dist2)
      END DO loop_over_training_vectors
!____________________________________________________________________  
!      compute new centers
!____________________________________________________________________  
       DO jcluster=1,nclusters
          IF(nvector(jcluster) /= 0) THEN
             mu(:,jcluster)=attribute_cluster_sum(:,jcluster)/nvector(jcluster)
         ELSE
            WRITE(stderr,*) ' cluster ',jcluster,' not updated. iteration = ',jiter
         END IF
         WRITE(stderr,*) 'cluster = ',jcluster,' nvector = ',nvector(jcluster)
         WRITE(stderr,'(8e12.5)') (mu(jattr,jcluster),jattr=1,nattr)
      END DO
      IF(jiter > 1) THEN
         convergence=(total_dist2_previous-total_dist2)/total_dist2
         total_dist2_previous=total_dist2
      ELSE
         convergence=0.0
         total_dist2_previous=total_dist2
      END IF
      WRITE(stderr,'(a,i5,a,e12.5,a,e12.5)') 'iteration = ',jiter,' total_dist2 = ',total_dist2,' convergence = ',convergence
      IF(jiter > 1 .AND. convergence < .001) EXIT training_loop
   END DO training_loop
END IF training_if_block
IF(use_mpi) THEN
   CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
END IF
!___________________________________________________________________
! broadcast the inverse covariance matrix and cluster centers to all processors
!___________________________________________________________________
CALL MPI_Bcast(cov_inv,      nattr*nattr,     MPI_REAL,   0, MPI_COMM_WORLD, mpi_rc)
CALL MPI_Bcast(sigma,        nattr,           MPI_REAL,   0, MPI_COMM_WORLD, mpi_rc)
CALL MPI_Bcast(mu,           nattr*nclusters, MPI_REAL,   0, MPI_COMM_WORLD, mpi_rc)
CALL MPI_Bcast(top_sample,   ncdp_in*nline_in,MPI_INTEGER,0, MPI_COMM_WORLD, mpi_rc)
CALL MPI_Bcast(bottom_sample,ncdp_in*nline_in,MPI_INTEGER,0, MPI_COMM_WORLD, mpi_rc)

IF(master) WRITE(stderr,*) 'Start Processing. end_line_out =',end_line_out
IF(master) THEN
!____________________________________________________________________
!     Reset file positions for each of the attribute files
!____________________________________________________________________
      DO jattr=1,nattr
          CALL aaspi_seek(input_fn(jattr),0,SEEK_SET,return_code)
      END DO
END IF

!______________________________________________________________________
! Do the calculations
!______________________________________________________________________
loop_over_lines:  DO jline=1,end_line_out
   IF(master) THEN 
!______________________________________________________________________
!     read in the attributes, one cdp at a time 
!     scale data. window into target zone of interest. 
!______________________________________________________________________
      IF(MOD(jline,1) == 0) WRITE(stderr,'(a,3i10)') 'calculate posterior mean projections: jline,end_line_out',jline,end_line_out               
      next_header_in=(jline-1)*ncdp_in+1
!______________________________________________________________________
!     read in a line of headers. these will be output! 
!______________________________________________________________________
      CALL aaspi_get_val_headers(input_fn(1),next_header_in,ncdp_in,header(:,1:ncdp_in),return_code)  
      DO jattr=1,nattr
         CALL aaspi_read(input_fn(jattr),d_vector(:,:,jattr),return_code)
      END DO
      DO jcdp=begin_cdp_out,end_cdp_out  
         IF(live(jcdp,jline)) THEN          
            DO jattr=1,nattr
               DO jt=top_sample(jcdp,jline),bottom_sample(jcdp,jline)
                  data_vector_master(jattr,jt,jcdp)=d_vector(jt,jcdp,jattr)
               END DO  
            END DO
         ELSE
            data_vector_master(:,:,jcdp)=0.0
         END IF
      END DO
   END IF

   IF (jline < begin_line_out)  CYCLE loop_over_lines
   CALL init_time(last_time(13))
!______________________________________________________________________
!   Master Compute Block
!______________________________________________________________________
    mpi_block: IF (use_mpi) THEN
      CALL MPI_Barrier(mpi_comm,mpi_rc)
      IF (master) THEN
         num_sent=0
         jcdp=begin_cdp_out
!______________________________________________________________________
!        initiate by sending each slave a live trace
!        if there are fewer live traces, tell the slave there are no more traces.
!______________________________________________________________________
         startup_slaves: DO slave_process=1,n_slave  !Distribute work to slaves
           DO WHILE (jcdp <= end_cdp_out)
             IF (live(jcdp,jline)) EXIT		! a live trace has been found. go to next step.
!______________________________________________________________________
!            skip dead traces
!______________________________________________________________________
             if(verbose) write(0,*) 'dead trace! jcdp = ',jcdp,' jline = ', jline
             jcdp=jcdp+1
           END DO
           IF ( jcdp <= end_cdp_out ) THEN
             more_traces=.TRUE.
             CALL MPI_Send(more_traces,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
             CALL send_data_to_slave(slave_process,jcdp)
             live_flag=live(jcdp,jline)
             CALL MPI_Send(live_flag,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
             jcdp=jcdp+1
             num_sent=num_sent+1
           ELSE             ! Otherwise unassigned processors wait
             more_traces=.FALSE.
             if(verbose) write(0,*) 'tell slave_process ',slave_process,' that there are no more traces on this jline = ', jline
             CALL MPI_Send(more_traces,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
           END IF
         END DO startup_slaves 

         receive_results_and_send_more_work: DO WHILE (jcdp <= end_cdp_out)
            CALL MPI_Recv(slave_process,1,MPI_INTEGER,MPI_ANY_SOURCE,ready_tag,MPI_COMM_WORLD,process_status,mpi_rc)
            CALL receive_data_on_master(slave_process,jcdp)
            more_traces=.TRUE.
            IF(verbose) write(0,*) this_process,':  prepare to send more_traces. more_traces_tag =  ',more_traces_tag,' slave_process = ',slave_process
            CALL MPI_Send(more_traces,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
            CALL send_data_to_slave(slave_process,jcdp)
            live_flag=live(jcdp,jline)
            CALL MPI_Send(live_flag,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
            jcdp=jcdp+1
            DO WHILE (jcdp <= end_cdp_out)
               IF (live(jcdp,jline)) EXIT	! a live trace has been found. continue to next step.
               if(verbose) write(0,*) 'dead trace! jcdp = ',jcdp,' jline = ', jline
               jcdp=jcdp+1
            END DO
         END DO receive_results_and_send_more_work
         collect_work_on_master: DO slave_cnt=1,num_sent ! Stop slaves
            CALL MPI_Recv(slave_process,1,MPI_INTEGER,MPI_ANY_SOURCE,ready_tag,MPI_COMM_WORLD,process_status,mpi_rc)
            CALL receive_data_on_master(slave_process,jcdp)
            more_traces=.FALSE.
            CALL MPI_Send(more_traces,1,MPI_LOGICAL,slave_process,more_traces_tag,MPI_COMM_WORLD,mpi_rc)
         END DO collect_work_on_master
      ELSE   
!_______________________________________________________________________________________
!        slaves. receive data. perform the work. send back the results.
!_______________________________________________________________________________________
         perform_work_on_slaves: DO
            IF(verbose) write(0,*) this_process,':  prepare to receive more_traces. more_traces_tag =  ',more_traces_tag
            CALL MPI_Recv(more_traces,1,MPI_LOGICAL,0,more_traces_tag,MPI_COMM_WORLD,process_status,mpi_rc)
            IF (verbose) WRITE(stderr,*) this_process,'more_traces=',more_traces
            IF (.NOT. more_traces) THEN
                IF(verbose) WRITE(stderr,*) this_process,' go to barrier. no more traces on jline = ',jline
                EXIT perform_work_on_slaves
            END IF
!            IF(verbose) WRITE(stderr,*) this_process,' waiting for data from master. jcdp = ',jcdp
            CALL receive_data_on_slave(jcdp)
            CALL MPI_Recv(live_flag,1,MPI_LOGICAL,0,more_traces_tag,MPI_COMM_WORLD,process_status,mpi_rc)
            IF(live_flag) THEN
               cluster_number_slave(:)=0.0
               DO jt=top_sample(jcdp,jline),bottom_sample(jcdp,jline)
                  CALL go_cluster(cov_inv,data_vector_slave(:,jt),mu,attribute_cluster_sum,	&
                               nattr,nclusters,cluster_number_slave(jt),nvector,total_dist2)
               END DO
            ELSE
               cluster_number_slave(:)=0.0
            END IF
            CALL MPI_Send(this_process,1,MPI_INTEGER,0,ready_tag,MPI_COMM_WORLD,mpi_rc)
            CALL send_data_to_master(jcdp)
         END DO perform_work_on_slaves
      END IF
      IF (verbose) WRITE(stderr,*) this_process,'at barrier. jline = ',jline
      CALL MPI_Barrier(MPI_COMM_WORLD,mpi_rc)
   ELSE mpi_block 
!____________________________________________________________________
!     not running on MPI
!____________________________________________________________________
      over_cdps2: DO jcdp=begin_cdp_out,end_cdp_out
         IF (.NOT. live(jcdp,jline)) CYCLE over_cdps2
         DO jt=top_sample(jcdp,jline),bottom_sample(jcdp,jline)
            CALL go_cluster(cov_inv,data_vector_master(:,jt,jcdp),mu,attribute_cluster_sum,	&
                            nattr,nclusters,cluster_number_master(jt,jcdp),nvector,total_dist2)
         END DO
      END DO over_cdps2
   END IF mpi_block
   IF(verbose) WRITE(stderr,*) this_process,': after barrier.'

   CALL add_time(total_time(13),last_time(13))

   IF(master) THEN
!____________________________________________________________________
!     write out the final Mean projections into AASPI file format
!____________________________________________________________________
      CALL init_time(last_time(8))
      next_header_out=(jline-begin_line_out)*ncdp_out+1
      CALL aaspi_put_val_headers(cluster_number_fn,next_header_out,ncdp_out,header(:,begin_cdp_out:end_cdp_out),return_code)
      DO jcdp=1,end_cdp_out
         IF (.NOT. live(jcdp,jline)) THEN
            cluster_number_master(:,jcdp)=0.0
         END iF
      END DO
      CALL aaspi_write(cluster_number_fn,cluster_number_master(:,:))
      CALL add_time(total_time(8),last_time(8))
   END IF
END DO loop_over_lines


CALL add_time(total_time(0),last_time(0))
WRITE(stderr,'(a9,a40,2a16)') 'process','task','time (hr)','time/trace (s)'

IF(master) THEN
   total_time(4)=total_time(4)-total_time(11) ! remove dist2
   total_time(6)=total_time(6)-total_time(12) ! remove dist2

   WRITE(stderr,'(i4,a1,a40,2f16.3)')                                                         	   &
     this_process,':',' read  attribute data',     total_time(1)/3600. ,total_time(1)/n_data_vectors, &
     this_process,':',' compute covariance matrix',total_time(2)/3600. ,total_time(2)/n_data_vectors, &
     this_process,':',' compute inverse covariance matrix',    total_time(3)/3600. ,total_time(3)/n_data_vectors, &
     this_process,':',' compute dist2 over data',  total_time(11)/3600. ,total_time(11)/n_data_vectors, &
     this_process,':',' other train over data',    total_time(4)/3600. ,total_time(4)/n_data_vectors, &
     this_process,':',' write results to disk',    total_time(8)/3600. ,total_time(8)/n_data_vectors, &
     this_process,':',' do the calculations',      total_time(13)/3600. ,total_time(13)/n_data_vectors, &
     this_process,':',' total time',               total_time(0)/3600. ,total_time(0)/n_data_vectors
END IF


WRITE(stderr,*) this_process,': Normal completion of program kmeans3d'
WRITE(stderr,*) this_process,': before call to mpi_barrier. use_mpi = ',use_mpi
IF (use_mpi) THEN
   IF(verbose) WRITE(stderr,*) this_process,': call MPI_Barrier'
   CALL MPI_Barrier(MPI_COMM_WORLD, mpi_rc)
   IF(verbose) WRITE(stderr,*) this_process,': call MPI_Finalize'
   CALL MPI_Finalize(mpi_rc)
END IF


CONTAINS
!====================================================================
!
!====================================================================
	
SUBROUTINE send_data_to_slave(slave_process,jcdp)
IMPLICIT NONE
INTEGER,INTENT(IN) :: slave_process
INTEGER,INTENT(IN) :: jcdp

IF (verbose) THEN
   WRITE(stderr,'(i3,a,i3,a,i8,a,i12,a,i12)') this_process,': sending data to slave',slave_process,' for jcdp=',jcdp,' cdp_tag=',cdp_tag,' data_in_tag=',data_in_tag
END IF
CALL MPI_Send(jcdp,1,MPI_INTEGER,slave_process,cdp_tag,MPI_COMM_WORLD,mpi_rc)
CALL MPI_Send(data_vector_master(1:nattr,jt_start:jt_end,jcdp),nattr*nt_out,MPI_REAL, slave_process,data_in_tag,MPI_COMM_WORLD,mpi_rc)
END SUBROUTINE send_data_to_slave
!====================================================================
!
!====================================================================


SUBROUTINE receive_data_on_slave(jcdp)
IMPLICIT NONE
INTEGER,INTENT(OUT) :: jcdp

IF (verbose) THEN
   WRITE(stderr,'(i3,a,i3,a,i8,a,i12,a,i12)') this_process,': receiving data from master',slave_process,' for jcdp=',jcdp,' cdp_tag=',cdp_tag,' data_in_tag=',data_in_tag
END IF
CALL MPI_Recv(jcdp,1,MPI_INTEGER,0,cdp_tag,MPI_COMM_WORLD,process_status,mpi_rc)
CALL MPI_Recv(data_vector_slave(:,:),nattr*nt_out,MPI_REAL,0,data_in_tag,MPI_COMM_WORLD,process_status,mpi_rc)
IF (verbose) THEN
   WRITE(stderr,'(i3,a,i8)') this_process,': received data for jcdp=',jcdp
END IF
END SUBROUTINE receive_data_on_slave
!====================================================================
!
!====================================================================
SUBROUTINE send_data_to_master(jcdp)
IMPLICIT NONE
INTEGER,INTENT(IN) :: jcdp

IF (verbose) THEN
   WRITE(stderr,'(i3,a,i8)') this_process,': sending data to master for jcdp=',jcdp
END IF
CALL MPI_Send(jcdp,1,MPI_INTEGER,0,cdp_tag,MPI_COMM_WORLD,mpi_rc)
CALL MPI_Send(cluster_number_slave(:),nt_out,MPI_REAL,0,cluster_number_tag,MPI_COMM_WORLD,mpi_rc)

END SUBROUTINE send_data_to_master

!====================================================================
!
!====================================================================


SUBROUTINE receive_data_on_master(slave_process,jcdp)
IMPLICIT NONE

INTEGER, INTENT(IN)	:: slave_process
INTEGER, INTENT(IN)	:: jcdp           
INTEGER  :: jcdp_received       !local variable in order to not contaminate global

IF (verbose) THEN
   WRITE(stderr,'(i3,a,i8)') this_process,': receiving data from slave',slave_process
END IF
CALL MPI_Recv(jcdp_received,1,MPI_REAL,slave_process,cdp_tag,MPI_COMM_WORLD,process_status,mpi_rc)
CALL MPI_Recv(cluster_number_master(:,jcdp_received),nt_out,MPI_REAL,slave_process,cluster_number_tag,MPI_COMM_WORLD,process_status,mpi_rc)

IF (verbose) THEN
   WRITE(stderr,'(i3,a,i3,a,i8,a,i8)') this_process,': got data from slave ',slave_process,' for jcdp_received=',jcdp_received,' jcdp=',jcdp
END IF

END SUBROUTINE receive_data_on_master

!====================================================================
!
!====================================================================


SUBROUTINE read_params
!____________________________________________________________________
!  read input parameters from command line and input files
!____________________________________________________________________
!_______________________________________________________________
!  read input data values from the history file.
!_______________________________________________________________
CALL read_3d_history_file(input_fn(1),nt_in,ncdp_in,nline_in,dt,d2,d3,o1_in,o2_in,o3_in,      &
                             unit1,unit2,unit3,label1,label2,label3,                    &
                             dcdp,dline,inline_azimuth,crossline_azimuth,               &
                             unique_project_name,suffix,fixed_scale,                    &
                             colorbar_fn,title,esize,stderr)
!_______________________________________________________________
! read in output time range
!_______________________________________________________________
t_min=o1_in
t_max=t_min+(nt_in-1)*dt
t_start=aaspi_get_param('t_start',t_min)   
t_end=aaspi_get_param('t_end',t_max)

first_sample_in=NINT(t_min/dt)
last_sample_in=first_sample_in+nt_in-1
jt_start=NINT(t_start/dt)
jt_end=NINT(t_end/dt)
nt_out=jt_end-jt_start+1

!_______________________________________________________________
! read in output cdp range 
!_______________________________________________________________

first_cdp_in=NINT(o2_in)
first_line_in=NINT(o3_in)
d2_abs=NINT(ABS(d2))
cdp_order=INT(SIGN(1.,d2))
last_cdp_in=first_cdp_in+(ncdp_in-1)*d2_abs*cdp_order
first_cdp_out=aaspi_get_param('first_cdp_out',first_cdp_in)       
last_cdp_out=aaspi_get_param('last_cdp_out',last_cdp_in) 
first_cdp_out=MAX(first_cdp_in,first_cdp_out)
last_cdp_out=MIN(last_cdp_in,last_cdp_out)
begin_cdp_out=(first_cdp_out-first_cdp_in)/(d2_abs*cdp_order)+1
end_cdp_out=(last_cdp_out-first_cdp_in)/(d2_abs*cdp_order)+1
ncdp_out=end_cdp_out-begin_cdp_out+1
o2_out=first_cdp_out

!_______________________________________________________________
! read in output line range 
!_______________________________________________________________

d3_abs=NINT(ABS(d3))
line_order=NINT(SIGN(1.,d3))
first_line_in=NINT(o3_in)
last_line_in=first_line_in+(nline_in-1)*d3_abs*line_order
first_line_out=aaspi_get_param('first_line_out',first_line_in)       
last_line_out=aaspi_get_param('last_line_out',last_line_in)       
first_line_out=MAX(first_line_in,first_line_out)
last_line_out=MIN(last_line_in,last_line_out)
begin_line_out=(first_line_out-first_line_in)/(d3_abs*line_order)+1
end_line_out=(last_line_out-first_line_in)/(d3_abs*line_order)+1
nline_out=end_line_out-begin_line_out+1
o3_out=first_line_out

line_decimation=aaspi_get_param('line_decimation',10)
cdp_decimation=aaspi_get_param('cdp_decimation',10)
vertical_sample_decimation=aaspi_get_param('vertical_sample_decimation',10)
n_data_training_iterations=aaspi_get_param('n_data_training_iterations',25)

WRITE(stderr,'(a20,i10)') 'nt_in',nt_in,'first_sample_in',first_sample_in,         &
          'last_sample_in',last_sample_in,'jt_start',jt_start,                    &
          'jt_end',jt_end,'nt_out',nt_out,'vertical_sample_decimation',vertical_sample_decimation,                    &
          'first_line_in',first_line_in,'nline_in',nline_in,'first_line_out',first_line_out,       &
          'begin_line_out',begin_line_out,'last_line_in',last_line_in,'line_decimation',line_decimation,          &
          'last_line_out',last_line_out,'end_line_out',end_line_out,'nline_out',nline_out,            &
          'line_order',line_order,'d3_abs',NINT(d3_abs),                              &
          'first_cdp_in',first_cdp_in,'first_cdp_out',first_cdp_out,            &
          'begin_cdp_out',begin_cdp_out,'last_cdp_in',last_cdp_in,              &
          'last_cdp_out',last_cdp_out,'end_cdp_out',end_cdp_out,'cdp_decimation',cdp_decimation,                &
          'ncdp_in',ncdp_in,'ncdp_out',ncdp_out,'cdp_order',cdp_order,'d2_abs',NINT(d2_abs)
!_______________________________________________________________
!  check validity of input parameters
!_______________________________________________________________
   input_error=0
   IF(t_start < t_min) THEN
      WRITE(stderr,*) 'input error! t_start = ',t_start,' < t_min = ',t_min
      input_error=input_error+1
   END IF
   IF(t_end > t_max) THEN
      WRITE(stderr,*) 'input error! t_end = ',t_end,' < t_max = ',t_max
      input_error=input_error+1
   END IF
   IF(input_error /= 0) THEN
      WRITE(stderr,*) 'program aborted due to ',input_error,' input errors'
      STOP 'check input parameters and resubmit!'
   END IF

END SUBROUTINE read_params

SUBROUTINE go_cluster(cov_inv,a,mu,attribute_cluster_sum,		&
                      nattr,nclusters,cluster_number,nvector,total_dist2)
!_______________________________________________________________
! assign attribute vectors to the current cluster centers and
! update centers for the next iteration
!_______________________________________________________________
IMPLICIT NONE
!_______________________________________________________________
!  variables passed to/from calling routine
!_______________________________________________________________
INTEGER,INTENT(IN)	:: nattr               		! the number of attributes         
INTEGER,INTENT(IN)	:: nclusters            		! the number of desired clusters   
REAL,   INTENT(IN)	:: cov_inv(nattr,nattr)		! the inverse covariance matrix
REAL,   INTENT(IN)	:: a(nattr)      		! an attribute vector            
REAL,   INTENT(IN)	:: mu(nattr,nclusters)  		! the current cluster centers    
REAL (KIND=8),   INTENT(INOUT)   :: attribute_cluster_sum(nattr,nclusters)	! sum of attributes vectors about current cluster center
REAL (KIND=8),   INTENT(INOUT)   :: total_dist2                          	! sum of all distances, which should decrease with iterations
INTEGER,INTENT(INOUT)   :: nvector(nclusters)					! number of data vectors per cluster                
REAL,   INTENT(OUT)     :: cluster_number					! cluster number assigned to current data vector
!_______________________________________________________________
!  local variables 
!_______________________________________________________________
REAL                    :: sum(nattr)			! intermediate sum
REAL                    :: dist(nattr)          	! the distance along each attribute axis between attribute vector and each cluster center
REAL                    :: mah_dist2(nclusters)	        ! square of the Mahalanobis distance to current cluster centers
INTEGER                 :: jcluster			! loop index
INTEGER                 :: jattr,kattr			! loop indices
REAL                    :: dist2_min                    ! minimum squared distance
INTEGER                 :: nearest_cluster              ! nearest cluster center
!_______________________________________________________________
! compute the Mahalanobis distance
!_______________________________________________________________
DO jcluster=1,nclusters
   mah_dist2(jcluster)=0.0
!_______________________________________________________________
! compute Euclidian distance along each axis to each cluster and its transpose
!_______________________________________________________________
   dist(:)=a(:)-mu(:,jcluster)
   DO kattr=1,nattr
      sum(kattr)=0.0
      DO jattr=1,nattr
         sum(kattr)=sum(kattr)+cov_inv(jattr,kattr)*dist(jattr)
      END DO
      mah_dist2(jcluster)=mah_dist2(jcluster)+dist(kattr)*sum(kattr)
   END DO
END DO
!_______________________________________________________________
! determine which cluster is smallest
!_______________________________________________________________
nearest_cluster=1
dist2_min=mah_dist2(nearest_cluster)
DO jcluster=2,nclusters
   IF(mah_dist2(jcluster) < dist2_min) THEN 
      dist2_min=mah_dist2(jcluster)
      nearest_cluster=jcluster
   END IF
END DO
!_______________________________________________________________
! accumulate for future update cluster centers
!_______________________________________________________________
attribute_cluster_sum(:,nearest_cluster)=attribute_cluster_sum(:,nearest_cluster)+a(:)
nvector(nearest_cluster)=nvector(nearest_cluster)+1
total_dist2=total_dist2+mah_dist2(nearest_cluster)
cluster_number=nearest_cluster

END SUBROUTINE go_cluster


END PROGRAM kmeans3d

