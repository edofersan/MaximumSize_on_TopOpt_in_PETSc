#include <petsc.h>
#include <TopOpt.h>
#include <LinearElasticity.h>
#include <MMA.h>
#include <Filter.h>
#include <MPIIO.h>
#include <mpi.h>


/*
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

Disclaimer:                                                              
The authors reserves all rights but does not guaranty that the code is   
free from errors. Furthermore, we shall not be liable in any event     
caused by the use of the program.                                     
 */

/* 
==========================================================
Modified by: Eduardo Fernandez, Kaike Yang, December 2019.
==========================================================
*/ 
 
 
static char help[] = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";

int main(int argc, char *argv[]){

	// Error code for debugging
	PetscErrorCode ierr;
	
	// Initialize PETSc / MPI and pass input arguments to PETSc
	PetscInitialize(&argc,&argv,PETSC_NULL,help);

	// STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
	TopOpt *opt = new TopOpt();

	// STEP 2: THE PHYSICS
	LinearElasticity *physics = new LinearElasticity(opt);

	// STEP 3: THE FILTERING
	Filter *filter = new Filter(opt);
	
	// STEP 4: VISUALIZATION USING VTK
	MPIIO *output = new MPIIO(opt->da_nodes,3,"ux, uy, uz",9,"x, xPhys, TestFilter, Xs, Xv, Ero, Dil, gms, gmsDIL");
	
	// STEP 5: THE OPTIMIZER MMA
	MMA *mma;
	PetscInt itr=0;
	opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !

	// STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
	ierr = filter->FilterProject(opt); CHKERRQ(ierr);
	
	// STEP 7: OPTIMIZATION LOOP   
	PetscScalar ch = 1.0;
	double t1,t2;
	while (itr < opt->maxItr && ch > 0.001){
		// Update iteration counter
		itr++;
		
		// UPDATING THE CONTINUATION PARAMETERS
		if (itr % opt->IterProj == 0){
			opt->penal   = PetscMin(opt->penal+0.25,3); // penal = penal : 0.25 : 3.0
			// move limits: initial 0.6, final 0.05
			opt->movlim  = (opt->movlimEnd-opt->movlimIni)/(3.0-1.0)*(opt->penal-1.0)+opt->movlimIni; 
			opt->Beta    = PetscMin((1.50*opt->Beta),38.0); // beta = beta*1.5 
			PetscPrintf(PETSC_COMM_WORLD,"===================================================\n");
			PetscPrintf(PETSC_COMM_WORLD,"It.: %i, Beta: %3.2f, Penal: %2.2f   movlim: %f\n", itr,opt->Beta,opt->penal,opt->movlim);
			PetscPrintf(PETSC_COMM_WORLD,"===================================================\n");
		}
		
		// Scaling the upper bound volume constraint of the Dilated Design
		if (itr==3 || itr%10 == 3){ 
			opt->volfrac = opt->VolDil_VolInt*opt->volfracREF;
			PetscPrintf(PETSC_COMM_WORLD,"=====================\n");
			PetscPrintf(PETSC_COMM_WORLD,"volfrac: %3.2f\n",opt->volfrac);
			PetscPrintf(PETSC_COMM_WORLD,"=====================\n");
		}
		
		// start timer
		t1 = MPI_Wtime();

		// Compute (a) obj+const, (b) sens, (c) obj+const+sens 
		ierr = physics->ComputeObjectiveConstraintsSensitivities(opt); CHKERRQ(ierr);
		
		// Compute objective scale
		if (itr==1 || (itr % opt->IterProj == 0)){ 
		  PetscPrintf(PETSC_COMM_WORLD,"It.: %i, obj_0: %f\n",itr,opt->fx);
			opt->fscale = 30.0/opt->fx; 
		}
		// Scale objective and sens
		opt->fx = opt->fx*opt->fscale;
		VecScale(opt->dfdx,opt->fscale);

		// Filter sensitivities (chainrule)
		ierr = filter->Gradients(opt); CHKERRQ(ierr);

		// Sets outer movelimits on design variables
		ierr = mma->SetOuterMovelimit(opt->Xmin,opt->Xmax,opt->movlim,opt->x,opt->xmin,opt->xmax); CHKERRQ(ierr);

		// Update design by MMA
		ierr = mma->Update(opt->x,opt->dfdx,opt->gx,opt->dgdx,opt->xmin,opt->xmax); CHKERRQ(ierr);

		// Inf norm on the design change
		ch = mma->DesignChange(opt->x,opt->xold);
		
		// Filter design field
		ierr = filter->FilterProject(opt); CHKERRQ(ierr);

		// stop timer
		t2 = MPI_Wtime();

		// Print to screen
		PetscPrintf(PETSC_COMM_WORLD,"It.: %i, obj.: %f, g[0]: %f, g[1]: %f, g[2]: %f, ch.: %f, time: %f\n",
				itr,(opt->fx/opt->fscale),opt->gx[0],opt->gx[1],opt->gx[2],ch,t2-t1);

		// Write field data: first 10 iterations and then every 20th
		if (itr<2 || itr%20==0){
			output->WriteVTK(opt->da_nodes,physics->GetStateField(),opt, itr);
		}

		// Dump data needed for restarting code at termination
		if (itr%3==0)	{
			opt->WriteRestartFiles(&itr, mma);
			physics->WriteRestartFiles();
		}
	}
	// Write restart WriteRestartFiles
	opt->WriteRestartFiles(&itr, mma);  
	physics->WriteRestartFiles();

	// Dump final design
	output->WriteVTK(opt->da_nodes,physics->GetStateField(),opt, itr+1);

	// STEP 7: CLEAN UP AFTER YOURSELF
	delete mma;
	delete output;
	delete filter;
	delete opt;  
	delete physics;

	// Finalize PETSc / MPI
	PetscFinalize();
	return 0;
}
