# A Fortran Sample code for particle packing in a rectangle box.
User guide: 
1) main codes: particle_packing.f, common.f(not necessary for these code sample)</br>
2) input files: hopp3d.inp, restart.dat </br>
3) output files </br>
   3-1) data generated during runing the job: </br>
   ppor.dat for bed porosity related information</br>
   preflow.dat for particle information data, need for calculation from existing data</br>
   *.vtu data for particle/container wall geometry visulation in software paraview </br>
   3-2) data generated during post-treatment: averageforce.dat, orenXY.dat, coordination.dat, etc.</br>
4) other related files: location-1-015.dat (needed only for ellipsoidal particles)</br>
5) compile command in linux (Massive3): "gfortran particle_packing.f"</br>
6) run the case: use "sbatch case" to submit a job from job script "case". </br>

For DEM packing theory, please refer to paper:</br> 
[1] J. Gan, A. Yu, Z. Zhou, DEM simulation on the packing of fine ellipsoids, Chemical Engineering Science, 156 (2016) 64-76.
