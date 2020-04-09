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

This code is based on the version 2017 of the TopOpt_in_PETSc code -> (N. Aage, 
E. Andreassen, B. S. Lazarov (2014), Topology optimization using PETSc: 
An easy-to-use, fully parallel, open source topology optimization framework, 
Structural and Multidisciplinary Optimization, DOI 10.1007/s00158-014-1157-0)

Modifications made by Eduardo Fernandez & Kaike Yang, at the University of Liege, Belgium.

NOTE: The code works with PETSc version 3.7.4
