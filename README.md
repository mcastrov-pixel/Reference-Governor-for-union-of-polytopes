# Reference Governor for union of non overlaping polytopes
Reference governor implementation for constraint sets represented as a union of polytopes

The approach is suited for cases where a non overlaping cover of the free space is available. Safe transition zones are created automatically using the concept of (weak) extensions and restrictions.

This code (and example) are taken from our work [1]. Please make sure to properly cite us if you found this useful :)



Two examples are available:

a drone navigating in an urban environment (see figure 1)
  main file is **drone_main.m**


an on-orbit proximity operation (see figure 2)
  main file is **CWH_main.m**



 <figure align="center">
  <img src="https://github.com/user-attachments/files/29059246/TRAJ_v1drone3D_2.pdf" width="75%">
  <figcaption><b>Figure 1:</b> Several 3D trajectories for a drone navigating in an urban environment.</figcaption>
</figure>


  
 <figure align="center">
  <img src="https://github.com/user-attachments/files/29059240/ISS_CWH3D_1_.pdf" width="75%">
  <figcaption><b>Figure 2:</b> Free space decomposition using hyperrectangles around a large spacecraft </figcaption>
</figure>
