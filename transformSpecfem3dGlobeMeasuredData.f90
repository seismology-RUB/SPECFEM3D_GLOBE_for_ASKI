!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
program transformSpecfem3dGlobeMeasuredData
  use specfem3dForASKI_mod
  use inversionBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use discreteFourierTransform
  use componentTransformation
  use complexKernelFrequency
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage
  use mathConstants

  implicit none

  ! command line
  type (argument_parser) :: ap
  character(len=max_length_string) :: parfile,str,forward_method_consistent_with_data
  character(len=max_length_string), dimension(:), pointer :: str_vec

  ! basics
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=35) :: myname = 'transformSpecfem3dGlobeMeasuredData'
  type (inversion_basics) :: invbasics

  ! component transformation
  integer :: ncomp,icomp
  character(len=character_length_component), dimension(:), allocatable :: comp
  double precision, dimension(:,:), pointer :: trans_coef
  real, dimension(:,:,:), allocatable :: trans_coef_all_transpose

  ! specfem seismograms
  character(len=2) :: band_instrument_code
  character(len=100) :: seisfile_extension
  integer :: NSTEP
  real :: DT
  real, dimension(:,:), pointer :: traces
  real, dimension(:), pointer :: stf

  ! fourier transformation
  type (discrete_fourier_transform) :: DFT
  complex, dimension(:,:), allocatable :: spectra,spectrum_rotated
  complex, dimension(:), allocatable :: stf_spectrum
  real :: df
  integer :: ifreq,nfreq
  integer, dimension(:), pointer :: jf
  complex, dimension(:), allocatable :: f

  double precision :: unit_factor,inverse_unit_factor

  ! filtering
  logical :: apply_filter_request,apply_filter,apply_event_filter,apply_station_filter
  character(len=400) :: path_event_filter,path_station_filter
  complex, dimension(:), pointer :: event_filter,station_comp_filter
  complex, dimension(:), allocatable :: filter

  ! other stuff
  integer :: istat,lu
  logical :: one_event_only,deconvolve_stf,print_usage_and_stop,&
       diff_time_series,scale_time_series,use_complex_freq
  real :: ts_scale_factor
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_specfem_seismograms,path_measured_data,file_measured_data
  double complex :: two_pi_i

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(str_vec,trans_coef,traces,stf,jf,event_filter,station_comp_filter)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Transform standard SPECFEM3D GLOBE 7.0.0 output to ASKI 1.0-1.2 spectral data in measured-data format")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-bicode",.true.,"(mandatory) bandcode and instrument code: the first two characters before "//&
       "the component in seismogram filename, e.g. 'LH' if your filenames look like 'network.staname.LH*.semd'",&
       "sval","")
  call addOption(ap,"-dt",.true.,"(mandatory) time step of the seismograms (as in SPECFEM3D Par_file)","rval","0.0")
  call addOption(ap,"-nstep",.true.,"(mandatory) number of samples NSTEP as in SPECFEM3D Par_file","ival","0")
  call addOption(ap,"-ocomp",.true.,"(mandatory) receiver components for which measured data is produced; valid "//&
       "components: '"//trim(all_valid_components)//"')","svec","")
  call addOption(ap,"-filter",.false.,"(optional) if set, the respective event filters and station (component) "//&
       "filters (as defined in main parfile) will be applied to the spectra")
  call addOption(ap,'-evid',.true.,"(optional) indicates a single event for which measured data is produced, "//&
       "otherwise measured data is produced for all events (as defined in ASKI FILE_EVENT_LIST)",'sval','')
  call addOption(ap,"-cfreq",.true.,"(optional) use complex frequencies producing data consistent with the given "//&
       "forward code (e.g. giving 'GEMINI' here, will produce Gemini-consistent spectral data at complex frequencies). "//&
       "Any filter values are also assumed to be given at those frequencies!","sval","")
  call addOption(ap,"-dconv",.false.,"(optional) if set, the source time function will be deconvolved from "//&
       "SPECFEM seismograms; consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI")
  call addOption(ap,"-diffts",.false.,"(optional) if set, additionally the time series will be differentiated "//&
       "(in the frequency domain after Fourier transform)")
  call addOption(ap,"-scale",.true.,"(optional) factor (different from 0) by which the time series are scaled "//&
       "before further processing","rval","1.0")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) goto 3
!
  str = ap.sval."main_parfile"
  parfile = str
  if (.level.(.errmsg.ap) == 2) goto 3
!
  print_usage_and_stop = .false.
!
  ! -bicode
  if(.not.(ap.optset."-bicode")) then
     write(*,*) "ERROR: please indicate -bicode"
     print_usage_and_stop = .true.
  else
     str = ap.sval."-bicode"
     band_instrument_code = str
  end if
!
  ! -evid
  one_event_only = (ap.optset."-evid")
  if(one_event_only) then
     str = ap.sval."-evid"
     evid_one_event_only = str
  end if
!
  ! -dt
  if(.not.(ap.optset."-dt")) then
     write(*,*) "ERROR: please indicate -dt"
     print_usage_and_stop = .true.
  else
     DT = ap.rval."-dt"
  end if  
!
  ! -nstep
  if(.not.(ap.optset."-nstep")) then
     write(*,*) "ERROR: please indicate -nstep"
     print_usage_and_stop = .true.
  else
     NSTEP = ap.ival."-nstep"
  end if  
!
  ! -filter
  apply_filter_request = (ap.optset."-filter")
!
  ! -ocomp
  if(.not.(ap.optset."-ocomp")) then
     print *, "ERROR: please indicate -ocomp"
     print_usage_and_stop = .true.
  else
     str_vec => ap.svec.'-ocomp'
     if (.level.(.errmsg.ap) == 2) goto 3
     if(.not.associated(str_vec)) then
        write(*,*) "ERROR: for some reason, there is no list of station components returned by argument parser, "//&
             "even though there was no error parsing argument -ocomp. This is strange..."
        write(*,*) ""
        goto 3
     end if
     ncomp = size(str_vec)
     allocate(comp(ncomp))
     do icomp=1,ncomp
        comp(icomp) = str_vec(icomp)
     end do
     deallocate(str_vec)
     if(.not.allValidComponents(comp,i_invalid=icomp)) then
        write(*,*) "ERROR: ",icomp,"'th output component '"//trim(comp(icomp))//"' not valid. Valid components are '"//&
             all_valid_components//"'"
        goto 3
     end if
  end if
!
  ! -dconv
  deconvolve_stf = (ap.optset."-dconv")
!
  ! -diffts
  diff_time_series = (ap.optset."-diffts")
!
  ! -scale
  scale_time_series = (ap.optset."-scale")
  if(scale_time_series) then
     ts_scale_factor = ap.rval."-scale"
     if (.level.(.errmsg.ap) == 2) goto 3
     if(ts_scale_factor == 0.0) then
        write(*,*) "ERROR: scaling factor for time series ('-scale' option) must not be zero!"
        print_usage_and_stop = .true.
     end if
  end if ! scale_time_series
!
  ! -cfreq
  use_complex_freq = (ap.optset."-cfreq")
  if(use_complex_freq) then
     forward_method_consistent_with_data = ap.sval."-cfreq"
  end if
!
  print_usage_and_stop = print_usage_and_stop .or. .level.(.errmsg.ap) == 2
!
  if(print_usage_and_stop) goto 3
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(parfile),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  preliminary processing
!
  if(apply_filter_request) then
     apply_event_filter = lval(.inpar.invbasics,'APPLY_EVENT_FILTER')
     apply_station_filter = lval(.inpar.invbasics,'APPLY_STATION_FILTER')
     apply_filter = apply_event_filter .or. apply_station_filter
  end if
!
  ! in case of one_event_only, check if evid_one_event_only is valid
  if(one_event_only) then
     errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid_one_event_only)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg == 2) then
        write(*,*) "ERROR: eventID '"//trim(evid_one_event_only)//"' given by option -evid is not contained in event list"
        goto 3
     end if
  end if
!
  ! get all transformation matrices here for transformations from specfem seismograms orientation
  ! get all transformation matrices here for transformations from specfem seismograms orientation
  ! NEZ to the requested output components contained in array comp
  ! THIS IS ACTUALLY QUITE INEFFICIENT, SINCE FOR MOST APPLICATIONS, ONLY THE NEZ SEISMOGRAMS ARE NEEDED!
  ! READING ALL NEZ SEISMOGRAMS AND DO "TRANSFORMATIONS" TO NEZ (or S,N,DOWN etc.) IS REDUNDANT COMPUTATION!
  ! KEEP IT LIKE THIS ANYWAY, IN ORDER BE ABLE TO USE THE SIMPLE ROUTINE readTraces FROM specfem3dForASKI_mod
  ! (MIGHT IMPROVE THIS IN THE FUTURE)
  allocate(trans_coef_all_transpose(3,ncomp,.nstat.(.statlist.invbasics)))
  istat = 0
  do while (nextStationSeismicNetwork(.statlist.invbasics,station))
     istat = istat + 1
     ! transpose trans_coef here by shwitching NEZ comp (i.e. coef_in = comp and coef_out = N,E,UP),
     ! as we need the transpose in matmul operation when actually transforming later on
     trans_coef => transform(.comptrans.invbasics,comp,(/'N ','E ','UP'/),.staname.station)
     if(.not.associated(trans_coef)) then
        write(*,*) "ERROR: no transformation coefficients for ",istat,"'th station '"//trim(.staname.station)//"'"
        goto 1
     end if
     trans_coef_all_transpose(:,:,istat) = trans_coef
     deallocate(trans_coef)
  end do ! while next station
!
  unit_factor = .ufmdata.invbasics
  inverse_unit_factor = 1.d0/unit_factor
!
  ! compute fourier transformation factors efactors here
  ! taper parameters like in SPECFEM3D-for-ASKI codes
  df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
  nfreq = ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')
  jf => ivecp(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',nfreq)
  if(.not.associated(jf)) then
     write(*,*) "ERROR: could not read ",nfreq," frequency indices from vector 'MEASURED_DATA_INDEX_OF_FREQ' in main parfile '"&
          //trim(parfile)//"'"
     goto 1
  end if
  two_pi_i = mc_two_pid * mc_cid
  if(apply_filter) allocate(filter(nfreq))
  allocate(f(nfreq))
  call new(errmsg,myname)
  if(use_complex_freq) then
     do ifreq = 1,nfreq
        ! get the frequency by module complexKernelFrequency (will return correct real-valued 
        ! frequency in complex variable if the forward method has no complex frequencies)
        f(ifreq) = getComplexKernelFrequency(forward_method_consistent_with_data,df,jf(ifreq))
     end do ! ifreq
     call initiateForwardDFT(DFT,DT,0,NSTEP-1,f,errmsg,hanning_taper=0.05)
  else
     f = jf*df
     call initiateForwardDFT(DFT,DT,0,NSTEP-1,f,errmsg,hanning_taper=0.05)
  end if
  if(.level.errmsg /=0) call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
  allocate(spectra(nfreq,3*.nstat.(.statlist.invbasics)),spectrum_rotated(nfreq,ncomp))
  if(deconvolve_stf) allocate(stf_spectrum(nfreq))
!
  seisfile_extension = '.sem.ascii' ! put hard-coded seismogram file extension here. must be adapted, if SPECFEM3D_GLOBE standard changes
!
!------------------------------------------------------------------------
!  write some info about this run now
!
  if(one_event_only) then
     write(*,*) "creating ASKI measured data from SPECFEM3D seismograms for one event and ",&
          .nstat.(.statlist.invbasics)," stations, "
  else
     write(*,*) "creating ASKI measured data from SPECFEM3D seismograms for ",.nev.(.evlist.invbasics)," events and ",&
          .nstat.(.statlist.invbasics)," stations, "
  end if
  write(*,*) "as of main parameter file '"//trim(parfile)//"'"
  write(*,*) ""
  write(*,*) "input SPECFEM3D seismograms: "
  write(*,*) "   NSTEP =  ",NSTEP
  write(*,*) "   DT =  ",DT
  write(*,*) "   band and instrument code = ",band_instrument_code
  write(*,*) "   seismogram orientation = NEZ"
  write(*,*) ""
  write(*,*) "output spectra: "
  write(*,*) "   number of frequencies = ",nfreq
  write(*,*) "   frequency step df = ",df
  write(*,*) "   frequency indices jf = ",jf
  if(use_complex_freq) then
     if(methodHasComplexKernelFrequency(forward_method_consistent_with_data)) then
        write(*,*) "   assuming complex frequencies consistent with forward method '",&
             trim(forward_method_consistent_with_data),&
             "', for which this frequency discretization corresponds to frequencies [Hz]:"
        write(*,*) "   ",f
     else ! methodHasComplexKernelFrequency
        write(*,*) "   WARNING: even though requested by flag -cfreq, method '",trim(forward_method_consistent_with_data),&
             "' does not use complex frequencies (or is unknown): assuming real-valued frequencies f = jf*df [Hz]:"
        write(*,*) "   ",real(f)
     end if ! methodHasComplexKernelFrequency
  else
     write(*,*) "   assuming real-valued frequencies f = jf*df [Hz]:"
     write(*,*) "   ",real(f)
  end if
  if(scale_time_series) then
     write(*,*) "   time series will be scaled by factor ",ts_scale_factor," before further processing"
  else
     write(*,*) "   time series will NOT be scaled by a factor"
  end if
  if(deconvolve_stf) then
     write(*,*) "   normalized differentiated source time function will be deconvolved from velocity "//&
          "seismograms (producing displacement spectra w.r.t. dirac)"
  else
     write(*,*) "   NO deconvolution of source time function will be applied. Reading in displacements ",&
          ", using seismogram file extension '"//trim(seisfile_extension)//"'"
  end if
  if(diff_time_series) then
     write(*,*) "   ADDITIONALLY, time series will be differentiated (in the frequency domain after Fourier transform)"
  else
     write(*,*) "   time series will NOT be differentiated additionally"
  end if
  if(apply_filter_request) then
     if(apply_filter) then
        if(apply_event_filter) then
           path_event_filter = (.inpar.invbasics).sval.'PATH_EVENT_FILTER'
           write(*,*) "   event filters will be applied, as defined by filter files in path '"//trim(path_event_filter)//"'"
        end if
        if(apply_station_filter) then
           path_station_filter = (.inpar.invbasics).sval.'PATH_STATION_FILTER'
           write(*,*) "   station filters will be applied, as defined by filter files in path '"//trim(path_station_filter)//"'"
        end if
     else ! apply_filter
        write(*,*) "   ALTHOUGH THERE IS FILTERING REQUESTED, the ASKI main parfile switches event and station filters OFF, "//&
             "so NO filtering will be applied here!!"
     end if ! apply_filter
  else ! apply_filter_request
     write(*,*) "   no filtering requested, the specfem time series will simply be transformed to frequency domain"
  end if ! apply_filter_request
  write(*,*) "   the displacement output spectra are computed according to the unit factor ",unit_factor
  write(*,*) ""

  path_measured_data = (.inpar.invbasics).sval.'PATH_MEASURED_DATA'
!
!
!------------------------------------------------------------------------
!  loop on all events (even if one_event_only, for simplicity of coding) 
!
  lu = get(fuh)
!
  do while (nextEventSeismicEventList(.evlist.invbasics,event))
!
     if(one_event_only) then
        evid = evid_one_event_only
     else
        evid = .evid.event
     end if
!
     ! in case of filtering, read in event specific filter values
     if(apply_event_filter) then
        if(associated(event_filter)) deallocate(event_filter)
        errmsg = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu,event_filter,ndata=nfreq)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg == 2) then
           write(*,*) "ERROR: could not read source filter from ascii file '"//trim(path_event_filter)//&
                "filter_"//trim(evid)
           goto 1
        end if
     end if
!
     path_specfem_seismograms = trim(path_measured_data)//'data_'//trim(evid)//'_OUTPUT_FILES/'
!
     write(*,*) "read all traces of event '",trim(evid),"' from path '",trim(path_specfem_seismograms),"'"
     ! read in all traces
     if(associated(traces)) deallocate(traces)
     call readTraces(traces,NSTEP,.statlist.invbasics,path_specfem_seismograms,band_instrument_code,&
          .true.,seisfile_extension,.false.,lu)
     if(.not.associated(traces)) then
        write(*,*) "no spectra produced for this event"
        goto 2
     end if
     if(scale_time_series) then
        traces = traces*ts_scale_factor
     end if
     ! if(diff_time_series) then ! differentiate by central differences
     !    do itrace = 1,size(traces,2)
     !       ! treat first sample separately (right finite difference)
     !       tmp_sample = traces(1,itrace) ! remember current sample for differentiation of the next sample
     !       traces(1,itrace) = real(  ( dble(traces(2,itrace)) - dble(traces(1,itrace)) )*one_over_DT  )
     !       do isamp = 2,size(traces,1)-1
     !          rtmp = real(  ( dble(traces(isamp+1,itrace)) - dble(tmp_sample) )*one_over_2DT  )
     !          tmp_sample = traces(isamp,itrace)
     !          traces(isamp,itrace) = rtmp
     !       end do ! isamp
     !       ! treat last sample separately (left finite difference)
     !       traces(size(traces,1),itrace) = real(  (dble(traces(size(traces,1),itrace))-dble(tmp_sample))*one_over_DT  )
     !    end do ! itrace
     ! end if ! diff_time_series
!
     if(deconvolve_stf) then
        write(*,*) "read, normalize and differentiate source time function of event '",trim(evid),"' from file '",&
             trim(path_specfem_seismograms),"plot_source_time_function.txt'"
        if(associated(stf)) deallocate(stf)
        call getDiffStf(stf,NSTEP,DT,trim(path_specfem_seismograms)//'plot_source_time_function.txt',lu,normalize=.true.)
        if(.not.associated(stf)) then
           write(*,*) "no spectra produced for this event"
           goto 2
        end if
        ! call diffTraces(traces,DT) ! FS FS INSTEAD OF DIFFERENTIATING THE TRACES IN THE TIME DOMAIN,
        !                                    IT IS BETTER TO MULTIPLY BY  2*pi*i*f  BELOW IN THE DECONVOLUTION!
     end if ! deconvolve_stf
!
     write(*,*) "compute spectra from traces"
     ! fourier transform to frequency domain of all traces at once
     call new(errmsg,myname)
     call transformForwardDFT(DFT,traces,spectra,errmsg)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg==2) goto 1
     call dealloc(errmsg)
     if(deconvolve_stf) then
        write(*,*) "compute spectrum of derivative of source time function and deconvolve synthetics"
        ! fourier transform to frequency domain of the single trace of stf'
        call new(errmsg,myname)
        call transformForwardDFT(DFT,stf,stf_spectrum,errmsg)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg==2) goto 1
        call dealloc(errmsg)
        ! deconvolve spectra
        !   THE ONLY DIFFERENCE TO THE CODE IN specfem3d_for_ASKI.f90 IS:
        !   HERE, spectra AND stf_spectrum ARE ALWAYS SINGLE PRECISION!
        !   IN specfem3d_for_ASKI.f90 , stf_spectrum IS ALWAYS DOULBLE PRECISION, AND BY OPTIONAL FLAG spectra MAY BE DOUBLE PRECISION, TOO
        do ifreq = 1,nfreq
           ! FS FS INSTEAD OF DIFFERENTIATING THE TRACES IN THE TIME DOMAIN (above), MULTIPLY HERE BY 2*pi*i*f
           ! (this seems to be more stable from a precisioin point of view, since the DFT is done in double precision above
           !  but time-domain differentiation would be done in single precision)
           spectra(ifreq,:) = two_pi_i * f(ifreq) * spectra(ifreq,:) / stf_spectrum(ifreq)
        end do ! ifreq
     end if ! deconvolve_stf
!
     if(diff_time_series) then ! differentiate the time series by multiplying their spectra by two_pi_i * f
        write(*,*) "at last, differentiate the traces by multiplying their spectra by two_pi_i * f"
        spectra(ifreq,:) = two_pi_i * f(ifreq) * spectra(ifreq,:)
     end if ! diff_time_series
!
     ! write spectra to files
     write(*,*) "write measured data files to path '",trim(path_measured_data),"'"
     istat = 0
     do while (nextStationSeismicNetwork(.statlist.invbasics,station))
        istat = istat + 1
!
        ! Rotate to requested components and apply the conversion factor to the requested unit.
        ! Note that the standard SPECFEM3D_GLOBE seismograms are in SI units (displacements in meters). Thus, multiply by inverse_unit_factor here.
        spectrum_rotated(:,:) = inverse_unit_factor * &
             matmul(spectra(:,(istat-1)*3+1:(istat-1)*3+3) , trans_coef_all_transpose(:,:,istat))
!
        do icomp = 1,ncomp
!
           ! in case of filtering, read in station_comp_filter for this path and component
           if(apply_station_filter) then
              if(associated(station_comp_filter)) deallocate(station_comp_filter)
              errmsg = readAsciiData(trim(path_station_filter)//"filter_"//trim(.staname.station)//"_"//trim(comp(icomp)),&
                   lu,station_comp_filter,ndata=nfreq)
              if(.level.errmsg /=0) call print(errmsg)
              if(.level.errmsg == 2) then
                 write(*,*) "ERROR: could not read station component filter from ascii file '"//trim(path_station_filter)//&
                      "filter_"//trim(.staname.station)//"_"//trim(comp(icomp))
                 goto 1
              end if
           end if
!
           ! finally set filter (if any)
           if(apply_filter) then
              filter = (1.,0.)
              if(apply_event_filter) filter = filter*event_filter
              if(apply_station_filter) filter = filter*station_comp_filter
           end if
!
           ! define filename of output file
           file_measured_data = "data_"//trim(evid)//"_"//trim(.staname.station)//"_"//trim(comp(icomp))
!
           write(*,*) "writing measured data file '",trim(file_measured_data),"'"
           if(apply_filter) then
              errmsg = writeAsciiData(trim(path_measured_data)//file_measured_data,lu,&
                   filter*spectrum_rotated(:,icomp))
           else
              errmsg = writeAsciiData(trim(path_measured_data)//file_measured_data,lu,spectrum_rotated(:,icomp))
           end if
           if(.level.errmsg/=0) call print(errmsg)
           if(.level.errmsg==2) goto 1
           call dealloc(errmsg)
        end do ! icomp
!
     end do ! while next station
!
2    write(*,*) ""
     if(one_event_only) exit
  end do ! while next event
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "good bye"
!
1 call dealloc(invbasics)
  call add(fuh,lu); call dealloc(fuh)
  call dealloc(ap)
  call dealloc(DFT)
  call dealloc(errmsg)
  if(allocated(comp)) deallocate(comp)
  if(associated(trans_coef)) deallocate(trans_coef)
  if(allocated(trans_coef_all_transpose)) deallocate(trans_coef_all_transpose)
  if(associated(traces)) deallocate(traces)
  if(associated(stf)) deallocate(stf)
  if(allocated(spectra)) deallocate(spectra)
  if(allocated(spectrum_rotated)) deallocate(spectrum_rotated)
  if(allocated(stf_spectrum)) deallocate(stf_spectrum)
  if(associated(jf)) deallocate(jf)
  if(allocated(f)) deallocate(f)
  if(associated(event_filter)) deallocate(event_filter)
  if(associated(station_comp_filter)) deallocate(station_comp_filter)
  if(allocated(filter)) deallocate(filter)
!
  stop
!
3 if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
  call usage(ap)
  goto 1
end program transformSpecfem3dGlobeMeasuredData
