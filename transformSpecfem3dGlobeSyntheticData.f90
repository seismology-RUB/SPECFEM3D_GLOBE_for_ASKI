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
program transformSpecfem3dGlobeSyntheticData
  use specfem3dForASKI_mod
  use inversionBasics
  use iterationStepBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use discreteFourierTransform
  use componentTransformation
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage
  use mathConstants

  implicit none

  ! command line
  type (argument_parser) :: ap
  character(len=max_length_string) :: parfile,str
  character(len=max_length_string), dimension(:), pointer :: str_vec

  ! basics
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=36) :: myname = 'transformSpecfem3dGlobeSyntheticData'
  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

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
  real, dimension(:), allocatable :: f

  double precision :: unit_factor,inverse_unit_factor

  ! other stuff
  integer :: istat,lu
  logical :: one_event_only,deconvolve_stf,print_usage_and_stop
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_specfem_seismograms,path_synthetic_data,file_synthetic_data
  double complex :: two_pi_i

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(str_vec,trans_coef,traces,stf,jf)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Transform standard SPECFEM3D GLOBE 7.0.0 output to ASKI 1.0-1.2 spectral data in synthetic-data format")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-bicode",.true.,"(mandatory) bandcode and instrument code: the first two characters before "//&
       "the component in seismogram filename, e.g. 'LH' if your filenames look like 'network.staname.LH*.sem'",&
       "sval","")
  call addOption(ap,"-dt",.true.,"(mandatory) time step of the seismograms (e.g. displayed in output_solver.txt)","rval","0.0")
  call addOption(ap,"-nstep",.true.,"(mandatory) number of samples NSTEP that should be accounted for (can be "//&
       "less than actual number which is e.g. displayed in output_solver.txt)","ival","0")
  call addOption(ap,"-ocomp",.true.,"(mandatory) receiver components for which synthetic data is produced; valid "//&
       "components: '"//trim(all_valid_components)//"')","svec","")
  call addOption(ap,'-evid',.true.,"(optional) indicates a single event for which synthetic data is produced, "//&
       "otherwise synthetic data is produced for all events (as defined in ASKI FILE_EVENT_LIST)",'sval','')
  call addOption(ap,"-dconv",.false.,"(optional) if set, the normalized and differentiated source time function "//&
       "(error function) will be deconvolved from the differentiated synthetics; consistend with "//&
       "'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI")
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
  ! -evid
  one_event_only = (ap.optset."-evid")
  if(one_event_only) then
     str = ap.sval."-evid"
     evid_one_event_only = str
  end if
!
  ! -dconv
  deconvolve_stf = (ap.optset."-dconv")
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
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  preliminary processing
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
  ! NEZ to the requested output components contained in array comp
  ! THIS IS ACTUALLY QUITE INEFFICIENT, SINCE FOR MOST APPLICATIONS, ONLY THE NEZ SEISMOGRAMS ARE NEEDED!
  ! READING ALL NEZ SEISMOGRAMS AND DO "TRANSFORMATIONS" TO NEZ (or S,N,DOWN etc.) IS REDUNDANT COMPUTATION!
  ! KEEP IT LIKE THIS ANYWAY, IN ORDER BE ABLE TO USE THE SIMPLE ROUTINE readTraces FROM specfem3dForASKI_mod
  ! (MIGHT IMPROVE THIS IN THE FUTURE)
  allocate(trans_coef_all_transpose(3,ncomp,.nstat.(.statlist.invbasics)))
  istat = 0
  do while (nextStationSeismicNetwork(.statlist.invbasics,station))
     istat = istat + 1
     ! transpose trans_coef here by switching comp with NEZ (i.e. coef_in = comp and coef_out = N,E,UP),
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
  nfreq = ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')
  jf => ivecp(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',nfreq)
  if(.not.associated(jf)) then
     write(*,*) "ERROR: could not read ",nfreq," frequency indices from vector 'ITERATION_STEP_INDEX_OF_FREQ' "//&
          "in iteration step parfile"
     goto 1
  end if
  two_pi_i = mc_two_pid * mc_cid
  allocate(f(nfreq))
  f = jf*df
  call new(errmsg,myname)
  call initiateForwardDFT(DFT,DT,0,NSTEP-1,f,errmsg,hanning_taper=0.05)
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
     write(*,*) "creating ASKI synthetic data from SPECFEM3D seismograms for one event and ",&
          .nstat.(.statlist.invbasics)," stations, "
  else
    write(*,*) "creating ASKI synthetic data from SPECFEM3D seismograms for ",.nev.(.evlist.invbasics)," events and ",&
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
  write(*,*) "   frequency step = ",df
  write(*,*) "   frequency indices = ",jf
  write(*,*) "   output spectra will be produced for the ",ncomp," receiver components: ",comp//", "
  if(deconvolve_stf) then
     write(*,*) "   normalized differentiated source time function will be deconvolved from velocity "//&
          "seismograms (producing displacement spectra w.r.t. dirac)"
  else
     write(*,*) "   NO deconvolution of source time function will be applied. Reading in displacements ",&
          ", using seismogram file extension '"//trim(seisfile_extension)//"'"
  end if
  write(*,*) "   the displacement output spectra are computed according to the unit factor ",unit_factor
  write(*,*) ""
!
  path_synthetic_data = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SYNTHETIC_DATA')
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
     path_specfem_seismograms = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//&
          'kernel_displ_'//trim(evid)//'_OUTPUT_FILES/'
!
     write(*,*) "read all traces of event '",trim(evid),"' from path '",trim(path_specfem_seismograms),"'"
!
     ! read in all traces
     if(associated(traces)) deallocate(traces)
     call readTraces(traces,NSTEP,.statlist.invbasics,path_specfem_seismograms,band_instrument_code,&
          .true.,seisfile_extension,.false.,lu)
     if(.not.associated(traces)) then
        write(*,*) "no spectra produced for this event"
        goto 2
     end if
!
     if(deconvolve_stf) then
        ! in case that the source time function should be deconvolved, it is assumed here that it was a
        ! steep error function.
        ! --> use velocity seismograms (i.e. differentiate the displacement seismograms) and deconvolve the 
        ! differentiated source time function (for reasons of numerical stability at very small frequencies)
        !                                                                                                                                                                                                                                    
        ! ALSO IMPORTANT TO NOTE HERE:
        ! for stability reasons differentiate the source-time function (error function) in the time domain (get Gaussian)
        ! before it is Fourier transformed and deconvolved (Gaussian spectrum is close to amplitude 1 everywhere).
        ! Apply the very same DFT object here, which was initiated INCLUDING hanning taper. THIS SHOULD ONLY BE DONE
        ! WHEN DIFFERENTIATING THE ERROR FUNCTION. If a source time function is Fourier transformed which is NOT
        ! zero on the tail, you should initiate a separate DFT object WITHOUT taper!
!
        if(associated(stf)) deallocate(stf)
        write(*,*) "read, normalize and differentiate source time function of event '",trim(evid),"' from file '",&
             trim(path_specfem_seismograms),"plot_source_time_function.txt'"
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
        ! fourier transform to frequency domain of all traces at once
        call new(errmsg,myname)
        call transformForwardDFT(DFT,stf,stf_spectrum,errmsg)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg==2) goto 1
        call dealloc(errmsg)
        ! deconvolve spectra
        !   THE ONLY DIFFERENCE TO THE CODE IN specfem3d_par_ASKI.f90 IS:
        !   HERE, spectra AND stf_spectrum ARE ALWAYS SINGLE PRECISION!
        !   IN specfem3d_par_ASKI.f90 , stf_spectrum IS ALWAYS DOULBLE PRECISION, AND BY OPTIONAL FLAG spectra MAY BE DOUBLE PRECISION, TOO
        do ifreq = 1,nfreq
           ! FS FS INSTEAD OF DIFFERENTIATING THE TRACES IN THE TIME DOMAIN (above), MULTIPLY HERE BY 2*pi*i*f
           ! (this seems to be more stable from a precisioin point of view, since the DFT is done in double precision above
           !  but time-domain differentiation would be done in single precision)
           spectra(ifreq,:) = two_pi_i * f(ifreq) * spectra(ifreq,:) / stf_spectrum(ifreq)
        end do ! ifreq
     end if ! deconvolve_stf
!
     ! write spectra to files
     write(*,*) "write synthetic data files to path '",trim(path_synthetic_data),"'"
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
           ! define filename of output file
           file_synthetic_data = "synthetics_"//trim(evid)//"_"//trim(.staname.station)//"_"//trim(comp(icomp))

           write(*,*) "writing synthetic data file '",trim(file_synthetic_data),"'"
           errmsg = writeAsciiData(trim(path_synthetic_data)//file_synthetic_data,lu,spectrum_rotated(:,icomp))
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
1 call dealloc(invbasics); call dealloc(iterbasics)
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
!
  stop
!
3 if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
  call usage(ap)
  goto 1
end program transformSpecfem3dGlobeSyntheticData
