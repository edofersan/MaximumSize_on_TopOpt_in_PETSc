MaximumSize_on_TopOpt_in_PETSc
=================================================================================
A 3D large-scale & maximum size-constrained topology optimization code using PETSc
=================================================================================

This code solves the compliance minimization problem subject to volume and 
maximum size restrictions, as shown in **(II)**. The problem is formulated using the 
robust design approach based on eroded, intermediate and dilated designs. A 
quarter of the 3D MBB beam is modeled and the density filter weights are adjusted 
to simulate an extension of the design domain.

<p align="center">
  <img width="1000" src="https://i.ibb.co/gzPvWpk/Logo-1.png" />
</p>

The theoretical details are presented in the manuscript: 

E. Fernandez, K. Yang, S. Koppen, P. Alarcon, S. Bauduin, P. Duysinx (2020). 
Imposing minimum and maximum member size, minimum cavity size, and minimum 
separation distance between solid members in topology optimization. 
(arXiv preprint arXiv:2003.00263).

## The Source

This code is based on the version 2017 of the TopOpt_in_PETSc code. The
original code can be found here: https://github.com/topopt/TopOpt_in_PETSc 
(N. Aage, E. Andreassen, B. S. Lazarov (2014), Topology optimization using PETSc: 
An easy-to-use, fully parallel, open source topology optimization framework, 
Structural and Multidisciplinary Optimization, DOI 10.1007/s00158-014-1157-0)

> **NOTE**: The code works with **PETSc version 3.7.4**

## Implementation

Implementations made by Eduardo Fernandez & Kaike Yang, University of Liege, Belgium.

The following files were modified from the original code:

- main(cc):  Continuation over the beta and penal parameters / Scaling the upper 
             bound of the volume constraint / Export additional data (MPIIO output).
			
- TopOpt(cc/h): Projection parameters (beta, thresholds, minimum size) / Maximum size 
                dimensions and parameters / Declaration and Initialization of arrays.
				
- Filter(cc/h): Solid Passive elements / Change of the filter Matrix H / Numerical 
                treatment on H wrt the boundaries / eroded, intermediate and dilated 
		projections / Maximum size Constraint.

- LinearElasticity(cc/h): Robust formulation (on eroded) / 3D-MBB boundary conditions;

- MPIIO(cc): Export additional data.

> **NOTE**: By default, the code solves the problem **(II)**. To solve the problem 
**(I)**, define m=1 in TopOpt.cc and comment out the maximum size restriction as 
indicated in Filter.cc (line 322). 

## Compilation and Execution

See requirements and instructions at https://github.com/topopt/TopOpt_in_PETSc

- To compile, e.g: make topopt
- To execute, e.g: mpiexec -np 6 ./topopt
- Visualize using ParaView.

The default problem takes about 1 hour to be solved on a laptop. The default 
solution looks like this:  

<p align="center">
  <img width="600" src="https://i.ibb.co/sKH8Wy7/Logo-2.png" />
</p>

To recover the full MBB beam, the solution must be reflected with respect to the
two planes of symmetry. The full optimized solution looks as follows:

<p align="center">
  <img width="600" src="https://i.ibb.co/6r99ng5/Logo-3.png" />
</p>
