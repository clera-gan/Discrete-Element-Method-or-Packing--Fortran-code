# A Fortran Sample code for particle packing in a rectangle box.
User guide: 
1) main codes: particle_packing.f, common.f(not necessary for these code sample)</br>
2) input files: </br>
   hopp3d.inp, to set general input information </br>
   restart.dat, set the first variable to 1 if start runing from non-existing "preflow.dat" data file, otherwise set it to 0 </br>
4) output files </br>
   3-1) data generated during runing the job: </br>
   ppor.dat for bed porosity related information</br>
   preflow.dat for particle information data, need for calculation from existing data</br>
   *.vtu data for particle/container wall geometry visulation in software paraview </br>
   3-2) data generated during post-treatment: averageforce.dat, orenXY.dat, coordination.dat, etc.</br>
4) other related files: location-1-015.dat (needed only for ellipsoidal particles)</br>
5) compile command in linux (Massive3): "gfortran particle_packing.f"</br>
6) to run the case: use "sbatch case" to submit a job from job script "case". Please refer to Massive user guide at https://docs.massive.org.au/M3/slurm/simple-batch-jobs.html#running-simple-batch-jobs for more details on runing jobs and job status check. </br>

For DEM packing theory, please refer to paper:</br> 
[1] J. Gan, A. Yu, Z. Zhou, DEM simulation on the packing of fine ellipsoids, Chemical Engineering Science, 156 (2016) 64-76.

Contact person:
Dr. Jieqing Gan: jieqing.gan@monash.edu
