MaximumSize_on_TopOpt_in_PETSc
=================================================================================
A 3D large-scale & maximum size-constrained topology optimization code using PETSc
=================================================================================

This code solves the compliance minimization problem subject to volume and 
maximum size restrictions. The problem is formulated using the robust design 
approach based on eroded, intermediate and dilated designs. A quarter of the 
3D MBB beam is modeled and the density filter weights are adjusted to simulate 
an extension of the design domain. The theoretical details are presented in 
the manuscript: 

E. Fernandez, K. Yang, S. Koppen, P. Alarcon, S. Bauduin, P. Duysinx (2020). 
Imposing minimum and maximum member size, minimum cavity size, and minimum 
separation distance between solid members in topology optimization. 
(arXiv preprint arXiv:2003.00263).

This code is based on the version 2017 of the TopOpt_in_PETSc code. The
original code can be found here: https://github.com/topopt/TopOpt_in_PETSc 
(N. Aage, E. Andreassen, B. S. Lazarov (2014), Topology optimization using PETSc: 
An easy-to-use, fully parallel, open source topology optimization framework, 
Structural and Multidisciplinary Optimization, DOI 10.1007/s00158-014-1157-0)

Modifications made by Eduardo Fernandez & Kaike Yang, University of Liege, Belgium.
 
The following files were modified:

- main(cc):  Continuation over the beta and penal parameters / Scaling the upper 
             bound of the volume constraint / Export additional data (MPIIO output).
			
- TopOpt(cc/h): Projection parameters (beta, thresholds, minimum size) / Maximum size 
                dimensions and parameters / Declaration and Initialization of arrays.
				
- Filter(cc/h): Solid Passive elements / Change of the filter Matrix H / Numerical 
                treatment on H wrt the boundaries / eroded, intermediate and dilated 
		projections / Maximum size Constraint.

- LinearElasticity(cc/h): Robust formulation (on eroded) / 3D-MBB boundary conditions;

- MPIIO(cc): Export additional data.							


NOTE: The code works with PETSc version 3.7.4
