# SPECFEM3D_GLOBE for ASKI

### Extension package for [ASKI](https://github.com/seismology-RUB/ASKI): use [SPECFEM3D_GLOBE](https://github.com/geodynamics/specfem3d_globe) for solving the seismic forward problem

SPECFEM3D_GLOBE for ASKI, as well as ASKI and some of its components, 
documentation and examples are available under terms of the 
[GNU General Public License](https://github.com/seismology-RUB/ASKI/blob/master/LICENSE)
(version 2 or higher) via [github](https://github.com/seismology-RUB). 
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/aski in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the respective repository on gitHub (e.g.
[here for the ASKI main package](https://github.com/seismology-RUB/ASKI/issues))

The main author is Florian Schumacher, Ruhr-University Bochum, Germany. 


## Usage, Documentation

Read the [manual](doc/SPECFEM3D_GLOBE_for_ASKI_manual.pdf) for information on 
how to set parameters correctly and how to use SPECFEM3D_GLOBE as a forward 
solver for ASKI. 


## Requirements

1. You require an installation of the [ASKI main package](https://github.com/seismology-RUB/ASKI):
   ```
   git clone --depth 1 --branch master https://github.com/seismology-RUB/ASKI
   ```
   
   The directory created by the git clone command will be referred to below as `ASKI/`

2. You need a functioning installation of the SPECFEM3D_GLOBE code, including 
   modifications for usage with ASKI:
   * You can either use the modified slim version of SPECFEM3D_GLOBE_V7.0.0 \[2015-07-10\]
     that comes with this package in subdirectory [specfem3d_globe/](specfem3d_globe/)
   * or use your running installation of SPECFEM3D_GLOBE and extend it for usage
     with ASKI (see "Extend regular SPECFEM3D" below).
     
   Also refer to the [manual](doc/SPECFEM3D_GLOBE_for_ASKI_manual.pdf), sections 1.3, 1.4.
   
   *In both cases, you still must install this package* (item "Installation" below).
       
3. You need basic experience in using the regular SPECFEM3D_GLOBE software!


## Installation

You should clone the latest version of the master branch of the 
[gitHub repository](https://github.com/seismology-RUB/SPECFEM3D_GLOBE_for_ASKI) 
to *the same* directory where you have cloned the ASKI main package to (in the 
[ASKI documentation](https://github.com/seismology-RUB/ASKI/blob/master/doc/ASKI_manual.pdf)
exemplarily called `/your/programs/`). That is, make sure that the git clone command
```
git clone --depth 1 --branch master https://github.com/seismology-RUB/SPECFEM3D_GLOBE_for_ASKI
```

creates the directory `/your/programs/SPECFEM3D_GLOBE_for_ASKI` (also referred to 
below simply as `SPECFEM3D_GLOBE_for_ASKI/`) and ASKI was installed to directory
`/your/programs/ASKI` (also referred to below simply as `ASKI/`) .

You need to compile few more ASKI binaries:
* in [SPECFEM3D_GLOBE_for_ASKI/Makefile](Makefile), set `COMPILER` appropriately, 
  adjust `FFLAGS` if required and set the variables `BLAS`, `LAPACK`, just as you did 
  for installing the ASKI main package
* Execute command
  ```
  make all
  ```
  
  from path `SPECFEM3D_GLOBE_for_ASKI/`.
  
After that, `ASKI/bin/` should contain the new binaries.


## Extend regular SPECFEM3D

If you have a regular SPECFEM3D_GLOBE installation which has not significantly
different functionality compared with SPECFEM3D_GLOBE release version 7.0.0 
\[2015-07-10\], you can extend it for ASKI by the following steps:

1. install SPECFEM3D_GLOBE on your system and make it run, gain 
   experience in using it (below, the installation path is refered to as 
   `specfem3d_globe/`)
2. append content of file [SPECFEM3D_GLOBE_for_ASKI/specfem3D_par_ASKI.f90](specfem3D_par_ASKI.f90) to file
   `specfem3d_globe/src/specfem3D/specfem3D_par.F90`
3. in `specfem3d_globe/src/specfem3D/prepare_timerun.F90` in subroutine `prepare_timerun`:<br>
   add the following line at the beginning of the subroutine, after the `use ...` statements:
   ```
   use specfem_for_ASKI_par
   ```
   
   add the following line close to the end of the subroutine, before `synchronize_all()` is called:
   ```
   call prepare_timerun_ASKI()
   ```
   
4. in `specfem3d_globe/src/specfem3D/iterate_time.F90` in subroutine `iterate_time`:<br>
   add the following line at the beginning of the subroutine, after the `use ...` statements:
   ```
   use specfem_for_ASKI_par
   ```
   
   add the following line just before the `enddo` of the main time loop:
   ```
   call write_ASKI_output()
   ```
   
5. append content of file [SPECFEM3D_GLOBE_for_ASKI/ASKI_external_model.f90](ASKI_external_model.f90) to file
   `specfem3d_globe/src/meshfem3D/meshfem3D_par.f90`
6. in `specfem3d_globe/src/meshfem3D/setup_model.f90` in subroutine `setup_model`:<br>
   add the following line at the beginning of the subroutine, after the `use ...` statements:
   ```
   use ASKI_external_model
   ```
   
   add the following line just before info output is written to `IMAIN`, after the 3D models are broadcasted:
   ```
   call broadcast_ASKI_external_model(myrank)
   ```
   
7. in `specfem3d_globe/src/meshfem3D/get_model.F90` in `subroutine get_model`:<br>
   add the following line at the beginning of the subroutine, after the `use ...` statements:
   ```
   use ASKI_external_model
   ```
   
   add the following lines just before define elastic parameters in the model (i.e. setting all 
   arrays `rhostore`, `kappavstore`, `muvstore`, ...) , just after all other `get model` routines:
   ```
   call values_ASKI_external_model(iregion_code,xmesh,ymesh,zmesh,r, &
            vpv,vph,vsv,vsh,rho,Qmu,Qkappa,eta_aniso,dvp, &
            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
            c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
   ```

8. append content of file [SPECFEM3D_GLOBE_for_ASKI/parallel_ASKI.f90](parallel_ASKI.f90)
   to file `specfem3d_globe/src/shared/parallel.f90`
9. recompile the relevant SPECFEM3D binaries by executing
   ```
   make xmeshfem3D xspecfem3D
   ```
   
   in directory `specfem3d_globe/`
10. in order to produce ASKI output in SPECFEM3D simulations, copy file
    [SPECFEM3D_GLOBE_for_ASKI/Par_file_ASKI](Par_file_ASKI) to your respective `DATA/` path
    (which is e.g. `specfem3d_globe/EXAMPLES/my_example/DATA/`, or `specfem3d_globe/DATA/`). This 
    file must be adjusted for any specific simulation (just as all other parameter files), 
    refer to the [manual](doc/SPECFEM3D_GLOBE_for_ASKI_manual.pdf) on how to use it.

Additionally, you may refer to the [manual](doc/SPECFEM3D_GLOBE_for_ASKI_manual.pdf)
section 1.4 for details on extending a regular SPECFEM3D code copy for use with ASKI.

