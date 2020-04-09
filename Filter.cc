#include <Filter.h>


/* -----------------------------------------------------------------------------
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013 
Copyright (C) 2013-2014,

This Filter implementation is licensed under Version 2.1 of the GNU
Lesser General Public License.  

This MMA implementation is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This Module is distributed in the hope that it will be useful,implementation 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
-------------------------------------------------------------------------- */

/* 
==========================================================
Modified by: Eduardo Fernandez, Kaike Yang, December 2019.
==========================================================
*/ 

Filter::Filter(TopOpt *opt){
	// Set all pointers to NULL
	H=NULL;  
	Hs=NULL;
	da_elem=NULL;
	pdef=NULL;
	
	// ==================
	// MAXIMUM SIZE =====
	DII=NULL;
	Avcorr=NULL; 			
	Vv=NULL; 					
	dVvdxPhys=NULL; 		
	Av=NULL; 					
	smax=NULL; 								
	dAggdx=NULL;
	// ROBUST:
	DII_DIL=NULL;
	// ==================
	// ==================
	
  // Call the setup method
	SetUp(opt);
}

Filter::~Filter(){
	// Deallocate data
	if (Hs!=NULL){ VecDestroy(&Hs); }
	if (H!=NULL){ MatDestroy(&H); }
	if (da_elem!=NULL){ DMDestroy(&da_elem); }
	if (pdef!=NULL){delete pdef; }
	
	// ========================================
	// MAXIMUM SIZE ===========================
	if (DII!=NULL){ MatDestroy(&DII); }
	if (Avcorr!=NULL){ VecDestroy(&Avcorr); }
	if (Vv!=NULL){ VecDestroy(&Vv); }
	if (dVvdxPhys!=NULL){ VecDestroy(&dVvdxPhys); }
	if (Av!=NULL){ VecDestroy(&Av); }
	if (smax!=NULL){ VecDestroy(&smax); }
	if (dAggdx!=NULL){ VecDestroy(&dAggdx); }
	// ROBUST:
	if (DII_DIL!=NULL){ MatDestroy(&DII_DIL); }
	if (Avcorr_DIL!=NULL){ VecDestroy(&Avcorr_DIL); }
	// ========================================
	// ========================================

}

// Filter design variables
PetscErrorCode Filter::FilterProject(TopOpt *opt){
	PetscErrorCode ierr;

	// Filter the design variables or copy to xPhys
	// STANDARD FILTER
	if (opt->filter == 1){
		// Filter the densitities
		ierr = MatMult(H,opt->x,opt->xPhys); CHKERRQ(ierr);
		VecPointwiseDivide(opt->xPhys,opt->xPhys,Hs);
	}
	// PDE FILTER
	else if (opt->filter == 2 ){
		ierr = pdef->FilterProject(opt->x, opt->xPhys); CHKERRQ(ierr);
		// Check for bound violation: simple, but cheap check!
		PetscScalar *xp;
		PetscInt locsiz;
		VecGetArray(opt->xPhys,&xp);
		VecGetLocalSize(opt->xPhys,&locsiz);
		for (PetscInt i=0;i<locsiz;i++){
			if (xp[i] < 0.0){
				if (PetscAbsReal(xp[i]) > 1.0e-4){
					PetscPrintf(PETSC_COMM_WORLD,"BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH RESOLUTION: xPhys = %f\n",xp[i]);
				}
				xp[i]= 0.0;
			}
			if (xp[i] > 1.0){
				if (PetscAbsReal(xp[i]-1.0) > 1.0e-4){
					PetscPrintf(PETSC_COMM_WORLD,"BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH RESOLUTION: xPhys = %f\n",xp[i]);
				}
				xp[i]=1.0;
			}

		}
		VecRestoreArray(opt->xPhys,&xp);
	}
	// Projection
	else if (opt->filter == 3){

	  // ====================================================================
		// ============== PASSIVE ELEMENTS SOLID ==============================
		if (opt->penal>=1.00){
		VecPointwiseMult(opt->x,opt->x,opt->xSolid0); // x = x*S0
		VecAXPY(opt->x,1.0,opt->xSolid1);             // x = x*S0 + 1.0*S1
		}
		// ============== END PASSIVE ELEMENTS ================================
		// ====================================================================
	
		// Filter the densitities
		ierr = MatMult(H,opt->x,opt->xFilter); CHKERRQ(ierr);
		
		/*
		VecPointwiseDivide(opt->xFilter,opt->xFilter,Hs);
		*/
		
		// ====================================================================
		// ============== PASSIVE ELEMENTS SOLID ==============================
		if (opt->penal>=1.00){
		VecPointwiseMult(opt->xFilter,opt->xFilter,opt->xSolid0);
		VecAXPY(opt->xFilter,1.0,opt->xSolid1);
		}
		// ============== END PASSIVE ELEMENTS ================================
		// ====================================================================
		
		
		// For Robust Formulation =============================================
		// ====================================================================
		PetscInt rstart,rend,k,i; 					// For inner loop
		PetscScalar *xF,*xDil,*xEro,*xInt; 			// For getting fields
		PetscScalar *dDil,*dEro,*dInt; 				// For getting derivatives
		PetscScalar Beta,etaDil,etaEro,etaInt; 		// For getting parameters
		etaDil = opt->etaDil; 						
		etaEro = opt->etaEro; 						
		etaInt = opt->eta;
		Beta   = opt->Beta;
		VecGetOwnershipRange(opt->xFilter,&rstart,&rend);
		
		VecGetArray(opt->xFilter,&xF);
		VecGetArray(opt->xPhysDil,&xDil);
		VecGetArray(opt->xPhysEro,&xEro);
		VecGetArray(opt->xPhys,&xInt);
		VecGetArray(opt->dDildFilt,&dDil);
		VecGetArray(opt->dErodFilt,&dEro);
		VecGetArray(opt->dPhysdFilt,&dInt);
		
		k = 0;
		for (i=rstart; i<rend; i++) {
			xDil[k] = Projection(xF[k],Beta,etaDil);
			xEro[k] = Projection(xF[k],Beta,etaEro);
			xInt[k] = Projection(xF[k],Beta,etaInt);
			dDil[k] = ProjectionSensitivity(xF[k],Beta,etaDil);
			dEro[k] = ProjectionSensitivity(xF[k],Beta,etaEro);
			dInt[k] = ProjectionSensitivity(xF[k],Beta,etaInt);
			k++;
		}
		
		VecRestoreArray(opt->xFilter,&xF);
		VecRestoreArray(opt->xPhysDil,&xDil);
		VecRestoreArray(opt->xPhysEro,&xEro);
		VecRestoreArray(opt->xPhys,&xInt);
		VecRestoreArray(opt->dDildFilt,&dDil);
		VecRestoreArray(opt->dErodFilt,&dEro);
		VecRestoreArray(opt->dPhysdFilt,&dInt);
		
		// ...........PASSIVE ELEMENTS SOLID 
		if (opt->penal>=1.00){
		VecPointwiseMult(opt->xPhys,opt->xPhys,opt->xSolid0);
		VecAXPY(opt->xPhys,1.0,opt->xSolid1);
		VecPointwiseMult(opt->xPhysEro,opt->xPhysEro,opt->xSolid0);
		VecAXPY(opt->xPhysEro,1.0,opt->xSolid1);
		VecPointwiseMult(opt->xPhysDil,opt->xPhysDil,opt->xSolid0);
		VecAXPY(opt->xPhysDil,1.0,opt->xSolid1);
		}
		// ......... END PASSIVE ELEMENTS
		
		// ROBUST FORMULATION VOLUME ADAPTATION ............
		VolDil = 0.0;
		VolInt = 0.0;
		VecSum(opt->xPhysDil,&VolDil);
		VecSum(opt->xPhys,&VolInt);
		opt->VolDil_VolInt = VolDil/VolInt;
		// .................................................
		
		// End Robust Formulation =============================================
		// ====================================================================
		
		
	} 
	// COPY IN CASE OF SENSITIVITY FILTER
	else {	ierr = VecCopy(opt->x,opt->xPhys); CHKERRQ(ierr); }

	return ierr;
}


// NEW FOR PROJECTION ================================================================================
// ===================================================================================================
PetscScalar Filter::Projection(PetscScalar xFiltered, PetscScalar Beta, PetscScalar eta) {
	PetscScalar tbn, tbrn, dtbn;
	tbn  = tanh(Beta*eta);
	tbrn = tanh(Beta*(xFiltered-eta));
	dtbn = tbn + tanh(Beta*(1.0-eta));
	return (tbn+tbrn)/(dtbn);
}

PetscScalar Filter::ProjectionSensitivity(PetscScalar xFiltered, PetscScalar Beta, PetscScalar eta) {
	PetscScalar tbn, dtbn;
	tbn  = 1.0/cosh(Beta*(xFiltered-eta));
	tbn  = Beta*tbn*tbn;
	dtbn = tanh(Beta*eta)+tanh(Beta*(1.0-eta));
	return tbn/dtbn;
}
// END FOR PROJECTION ================================================================================
// ===================================================================================================



// Filter the sensitivities
PetscErrorCode Filter::Gradients(TopOpt *opt){

	PetscErrorCode ierr;
	// Chainrule/Filter for the sensitivities
	if (opt->filter == 0)
		// Filter the sensitivities, df,dg
	{
		Vec xtmp;
		ierr = VecDuplicate(opt->x,&xtmp);  CHKERRQ(ierr);
		VecPointwiseMult(xtmp,opt->dfdx,opt->x);
		MatMult(H,xtmp,opt->dfdx);
		VecPointwiseDivide(xtmp,opt->dfdx,Hs);
		VecPointwiseDivide(opt->dfdx,xtmp,opt->x);
		VecDestroy(&xtmp);
	}
	else if (opt->filter == 1) {
		// Filter the densities, df,dg: STANDARD FILTER
		Vec xtmp;
		ierr = VecDuplicate(opt->x,&xtmp);  CHKERRQ(ierr);
		// dfdx
		VecPointwiseDivide(xtmp,opt->dfdx,Hs);
		MatMult(H,xtmp,opt->dfdx);
		// dgdx
		VecPointwiseDivide(xtmp,opt->dgdx[0],Hs);
		MatMult(H,xtmp,opt->dgdx[0]);
		// tidy up
		VecDestroy(&xtmp);
	}
	else if (opt->filter == 2){
		// Filter the densities, df,dg: PDE FILTER  
		ierr = pdef->Gradients(opt->dfdx,opt->dfdx); CHKERRQ(ierr);
		ierr = pdef->Gradients(opt->dgdx[0],opt->dgdx[0]); CHKERRQ(ierr);
	}
		// Sensitivities of the projection and filter
	else if (opt->filter == 3) {
		Vec xtmp;
		ierr = VecDuplicate(opt->x,&xtmp);  CHKERRQ(ierr);

		/* (for robust formulation)
		ierr = VecPointwiseMult(xtmp,opt->dPhysdFilt,opt->dfdx); CHKERRQ(ierr);		
		ierr = VecCopy(xtmp,opt->dfdx); CHKERRQ(ierr);			

		ierr = VecPointwiseMult(xtmp,opt->dPhysdFilt,opt->dgdx[0]); CHKERRQ(ierr);		
		ierr = VecCopy(xtmp,opt->dgdx[0]); CHKERRQ(ierr);
		*/
		
		// Robust formulation ===========================
		VecPointwiseMult(xtmp,opt->dErodFilt,opt->dfdx);
		VecCopy(xtmp,opt->dfdx);
		// VOLUME ADAPTATION ..............................
		VecPointwiseMult(xtmp,opt->dDildFilt,opt->dgdx[0]);		
		VecCopy(xtmp,opt->dgdx[0]);
		// ==============================================
		
		/* (for change of H matrix)
		// Filter the densities, df,dg: STANDARD FILTER
		// dfdx
		VecPointwiseDivide(xtmp,opt->dfdx,Hs);
		MatMult(H,xtmp,opt->dfdx);
		// dgdx
		VecPointwiseDivide(xtmp,opt->dgdx[0],Hs);
		MatMult(H,xtmp,opt->dgdx[0]);
		*/
		
		
		// ==== MODIFICATION OF H MATRIX ============
		// =========( Hs is integrated )=============
		// dfdx
		VecCopy(opt->dfdx,xtmp);
		MatMultTranspose(H,xtmp,opt->dfdx);// dfdx = H'*dfdx
		// dgdx
		VecCopy(opt->dgdx[0],xtmp);
		MatMultTranspose(H,xtmp,opt->dgdx[0]);// dgdx = H'*dgdx		
		// ==========================================
		// ==========================================
		
		
		
		
		// =======================================================================
		// ================ MAXIMUM SIZE CONSTRAINT ==============================
		// =======================================================================
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// To deactivate the MaxSize restrictions, comment out this entire block.
    // and define m=1 in the TopOpt.cc file.		
		
		// INTERMEDIATE ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		PetscScalar qq; 
		qq = PetscMin(opt->penal + 1.0, 3.0) ;
		
		// voids vector with penalization ----------------------------------------
		VecCopy(opt->xPhys,Vv);  						// Vv =  xPhys
		VecScale(Vv,-1); 								// Vv = -xPhys
		VecShift(Vv,1);  								// Vv = 1 - xPhys
		VecPow(Vv,qq); 									// Vv = (1 - xPhys)^qq
		
		// and its derivative
		VecCopy(Vv,dVvdxPhys); 							// dVvdxPhys = Vv = (1 - xPhys)^penal
		VecPow(dVvdxPhys,1.0-1.0/qq); 					// dVvdxPhys = Vv^(qq-1)
		VecScale(dVvdxPhys,-1.0*qq);  					// dVvdxPhys = -qq*(1 - xPhys)^(qq-1)
		
		// Amount of voids inside the test regions
		MatMult(DII,Vv,Av); 							// Av = DII*Vv	
		VecAXPY(Av,1,Avcorr); 							// Av = DII*Vv + Avcorr
		
		// local constraints
		VecScale(Av,-1); 								// Av = -Av
		VecShift(Av,opt->epsiMS);  						// Av = epsi-Av
		VecCopy(Av,opt->gms); 							// To plot
		
		// Feature for the Aggregation		
		VecCopy(Av,smax); 								// smax = epsi-(DII*Vv + Avcorr)
		VecShift(smax,(1.0-opt->epsiMS));  				// smax = smax + 1 - epsi = SI
				
		// Aggregation
		VecPow(smax,opt->pagg);							// smax = smax^(r)
		VecSum(smax,&Agg); 								// Agg = sum(smax^r)
		Agg = Agg/opt->Nagg;							// Agg = sum(smax^r)/N
		Agg = PetscPowScalar(Agg,1/opt->pagg); 			// Agg = (sum(smax^r)/N)^(1/r) 
		opt->gx[1] = Agg-1.0+opt->epsiMS; 				// Gms = Agg - 1 + epsi
		
		// Sensitivity------------------------------------------------------------
		VecPow(smax,(1.0-1.0/opt->pagg));				// smax = SI^(r*(1-1/r)) = SI^(r-1)
		Agg = PetscPowScalar(Agg,(1.0-opt->pagg));  	// Agg  = (sum(SI^r)/N)^(1/r-1) 
		Agg = (-1.0/opt->Nagg)*Agg; 					// Agg  = -1/N*(sum(SI^r)/N)^(1/r-1)
		VecScale(smax,Agg); 					    	// smax = SI^(r-1)*(-1/N)*(sum(SI^r)/N)^(1/r-1)
		
		MatMultTranspose(DII,smax,Vv); 			    	// Vv     = DII'*smax
		VecPointwiseMult(dAggdx,dVvdxPhys,Vv);	    	// dAggdx = dVvdxPhys.*(DII'*smax)		

		VecPointwiseMult(xtmp,opt->dPhysdFilt ,dAggdx); // xtmp   = dPhysdFilt.*dAggdx
		VecCopy(xtmp,dAggdx); 							// dAggdx = dAggdFilt
		
		// MODIFICATION OF H MATRIX ============
		MatMultTranspose(H,dAggdx,opt->dgdx[1]); 		// dgdx[1] = H'*dAggdFilt
		// =====================================
		
		
		
		
		// DILATE ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		// voids vector with penalization ----------------------------------------
		VecCopy(opt->xPhysDil,Vv);  					// Vv =  xPhys
		VecScale(Vv,-1); 								// Vv = -xPhys
		VecShift(Vv,1);  								// Vv = 1 - xPhys
		VecPow(Vv,qq); 									// Vv = (1 - xPhys)^qq
		
		// and its derivative
		VecCopy(Vv,dVvdxPhys); 							// dVvdxPhys = Vv = (1 - xPhys)^penal
		VecPow(dVvdxPhys,1.0-1.0/qq); 					// dVvdxPhys = Vv^(qq-1)
		VecScale(dVvdxPhys,-1.0*qq);  					// dVvdxPhys = -qq*(1 - xPhys)^(qq-1)
		
		// Amount of voids inside the test regions
		MatMult(DII_DIL,Vv,Av); 						// Av = DII*Vv	
		VecAXPY(Av,1,Avcorr_DIL); 						// Av = DII*Vv + Avcorr
		
		// local constraints
		VecScale(Av,-1); 								// Av = -Av
		VecShift(Av,opt->epsiMS);  						// Av = epsi-Av
		VecCopy(Av,opt->gmsDIL); 						// To plot
		
		// Feature for the Aggregation		
		VecCopy(Av,smax); 								// smax = epsi-(DII*Vv + Avcorr)
		VecShift(smax,(1.0-opt->epsiMS));  				// smax = smax + 1 - epsi = SI
				
		// Aggregation
		VecPow(smax,opt->pagg);							// smax = smax^(r)
		VecSum(smax,&Agg); 								// Agg = sum(smax^r)
		Agg = Agg/opt->Nagg;							// Agg = sum(smax^r)/N
		Agg = PetscPowScalar(Agg,1/opt->pagg); 			// Agg = (sum(smax^r)/N)^(1/r) 
		opt->gx[2] = Agg-1.0+opt->epsiMS; 				// Gms = Agg - 1 + epsi
		
		// Sensitivity------------------------------------------------------------
		VecPow(smax,(1.0-1.0/opt->pagg));				// smax = SI^(r*(1-1/r)) = SI^(r-1)
		Agg = PetscPowScalar(Agg,(1.0-opt->pagg));  	// Agg  = (sum(SI^r)/N)^(1/r-1) 
		Agg = (-1.0/opt->Nagg)*Agg; 					// Agg  = -1/N*(sum(SI^r)/N)^(1/r-1)
		VecScale(smax,Agg); 					    	// smax = SI^(r-1)*(-1/N)*(sum(SI^r)/N)^(1/r-1)
		
		MatMultTranspose(DII_DIL,smax,Vv); 	    		// Vv     = DII'*smax
		VecPointwiseMult(dAggdx,dVvdxPhys,Vv);	    	// dAggdx = dVvdxPhys.*(DII'*smax)		

		VecPointwiseMult(xtmp,opt->dDildFilt,dAggdx);   // xtmp   = dPhysdFilt.*dAggdx
		VecCopy(xtmp,dAggdx); 							// dAggdx = dAggdFilt
		
		// MODIFICATION OF H MATRIX ============
		MatMultTranspose(H,dAggdx,opt->dgdx[2]); 		// dgdx[1] = H'*dAggdFilt
		// =====================================
		
		// =======================================================================
		// =======================================================================
		// ============== END MAXIMUM SIZE CONSTRAINT ============================
		// =======================================================================
		
		
		// tidy up
		VecDestroy(&xtmp);
	}
	return ierr;
}


PetscErrorCode Filter::SetUp(TopOpt *opt){

	PetscErrorCode ierr;

	if (opt->filter==0 || opt->filter==1 || opt->filter==3){

		// Extract information from the nodal mesh
		PetscInt M,N,P,md,nd,pd; 
		DMBoundaryType bx, by, bz;
		DMDAStencilType stype;
		ierr = DMDAGetInfo(opt->da_nodes,NULL,&M,&N,&P,&md,&nd,&pd,NULL,NULL,&bx,&by,&bz,&stype); CHKERRQ(ierr);

		// Find the element size
		Vec lcoor;
		DMGetCoordinatesLocal(opt->da_nodes,&lcoor);
		PetscScalar *lcoorp;
		VecGetArray(lcoor,&lcoorp);

		PetscInt nel, nen;
		const PetscInt *necon;
		DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon);

		PetscScalar dx,dy,dz;
		// Use the first element to compute the dx, dy, dz
		dx = lcoorp[3*necon[0*nen + 1]+0]-lcoorp[3*necon[0*nen + 0]+0];
		dy = lcoorp[3*necon[0*nen + 2]+1]-lcoorp[3*necon[0*nen + 1]+1];
		dz = lcoorp[3*necon[0*nen + 4]+2]-lcoorp[3*necon[0*nen + 0]+2];
		VecRestoreArray(lcoor,&lcoorp);

		// -----------------------------------------------------------------------------------------------------
		// -------------------------- From here the code is slightly changed -----------------------------------

		// Create the minimum element connectivity shit
		PetscInt ElemConn;
		PetscInt ElemConnFilt;
		
		// Check dx,dy,dz and find max conn for a given rmin
		ElemConnFilt = (PetscInt)PetscMax(ceil(opt->rmin/dx)-1,PetscMax(ceil(opt->rmin/dy)-1,ceil(opt->rmin/dz)-1));
		
		// FOR 3D USE THE FOLLOWING LINE
		ElemConnFilt = PetscMin(ElemConnFilt,PetscMin((M-1)/2,PetscMin((N-1)/2,(P-1)/2)));

		// The following is needed due to roundoff errors 
		PetscInt tmp;
		MPI_Allreduce(&ElemConnFilt, &tmp, 1,MPIU_INT, MPI_MAX,PETSC_COMM_WORLD );
		ElemConn = tmp;

		// =====================================================================================================
		// ======== FOR MAXIMUM SIZE, THE PREVIOUS LINES ARE REPLACED BY: ======================================
		// =====================================================================================================
		PetscScalar Rref = PetscMax(opt->Rmax,opt->rmin);
		Rref 		 = PetscMax(opt->RmaxDil,Rref);
		ElemConn = (PetscInt)PetscMax(ceil(Rref/dx)-1,PetscMax(ceil(Rref/dy)-1,ceil(Rref/dz)-1));
		ElemConn = PetscMin(ElemConn,PetscMin((M-1)/2,PetscMin((N-1)/2,(P-1)/2)));
		MPI_Allreduce(&ElemConn, &tmp, 1,MPIU_INT, MPI_MAX,PETSC_COMM_WORLD );
		ElemConn = tmp;
		// =====================================================================================================
		// =====================================================================================================
		// =====================================================================================================

		// Print to screen: mesh overlap!
		PetscPrintf(PETSC_COMM_WORLD,"#  Filter radius rmin  = %f results in a stencil of %i elements \n",opt->rmin,ElemConnFilt);
		PetscPrintf(PETSC_COMM_WORLD,"# MaxSize radius Rmax  = %f results in a stencil of %i elements \n",opt->Rmax,ElemConn);
		PetscPrintf(PETSC_COMM_WORLD,"# Reference radius Rref= %f results in a stencil of %i elements \n",Rref,ElemConn);
		
		// -------------------------- Up to here the code was slightly changed ---------------------------------
		// -----------------------------------------------------------------------------------------------------
		
		
		PetscPrintf(PETSC_COMM_WORLD,"# M,N,P = %i,%i,%i \n",M,N,P);

		// Find the geometric partitioning of the nodal mesh, so the element mesh will coincide 
		PetscInt *Lx=new PetscInt[md];
		PetscInt *Ly=new PetscInt[nd];
		PetscInt *Lz=new PetscInt[pd];

		// get number of nodes for each partition
		const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
		DMDAGetOwnershipRanges(opt->da_nodes, &LxCorrect, &LyCorrect, &LzCorrect); 

		// subtract one from the lower left corner.
		for (int i=0; i<md; i++){
			Lx[i] = LxCorrect[i];
			if (i==0){Lx[i] = Lx[i]-1;}
		}
		for (int i=0; i<nd; i++){
			Ly[i] = LyCorrect[i];
			if (i==0){Ly[i] = Ly[i]-1;}
		}
		for (int i=0; i<pd; i++){
			Lz[i] = LzCorrect[i];
			if (i==0){Lz[i] = Lz[i]-1;}
		}

		// Create the element grid:
		DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,M-1,N-1,P-1,md,nd,pd,
				1,ElemConn,Lx,Ly,Lz,&da_elem);

		// Set the coordinates: from 0+dx/2 to xmax-dx/2 and so on
		PetscScalar xmax = (M-1)*dx;
		PetscScalar ymax = (N-1)*dy;
		PetscScalar zmax = (P-1)*dz;
		DMDASetUniformCoordinates(da_elem , dx/2.0,xmax-dx/2.0, dy/2.0,ymax-dy/2.0, dz/2.0,zmax-dz/2.0);

		// Allocate and assemble
		DMCreateMatrix(da_elem,&H);
		DMCreateGlobalVector(da_elem,&Hs);

		// Set the filter matrix and vector
		DMGetCoordinatesLocal(da_elem,&lcoor);
		VecGetArray(lcoor,&lcoorp);
		DMDALocalInfo info;
		DMDAGetLocalInfo(da_elem,&info);
		
		
		// =============================================================
		// FOR PASSIVE ELEMENTS ========================================
		DMCreateGlobalVector(da_elem,&opt->xSolid0);
		VecSet(opt->xSolid0,1.0);
		DMCreateGlobalVector(da_elem,&opt->xSolid1);
		VecSet(opt->xSolid1,0.0);
		PetscScalar valCero = 0.0;
		PetscScalar valOne = 1.0;
		PetscScalar RMIN   = opt->rmin/2.0; // half of the filter radius
		// =============================================================
		// =============================================================
		
		
		// NEW FOR SYMMETRY ===================
		// ====================================
		PetscScalar xjj = 1.0;
		PetscScalar yjj = 1.0;
		PetscScalar zjj = 1.0;
		PetscScalar distijj = 1.0;
		// ====================================
		// ====================================
		
		
		// ===================================
		// ====== NEW FOR MAXIMUM SIZE =======
		PetscScalar oneval = 1.0;
		DMCreateMatrix(da_elem,&DII);
		// ROBUST:
		DMCreateMatrix(da_elem,&DII_DIL);
		// ====== END FOR MAXIMUM SIZE =======
		// ===================================
		
		
		
		// The variables from info that are used are described below:
		// -------------------------------------------------------------------------
		// sw = Stencil width
		// mx, my, mz = Global number of "elements" in each direction 
		// xs, ys, zs = Starting point of this processor, excluding ghosts
		// xm, ym, zm = Number of grid points on this processor, excluding ghosts
		// gxs, gys, gzs = Starting point of this processor, including ghosts
		// gxm, gym, gzm = Number of grid points on this processor, including ghosts
		// -------------------------------------------------------------------------
		
		// Outer loop is local part = find row
		// What is done here, is:
		// 
		// 1. Run through all elements in the mesh - should not include ghosts
		for (PetscInt k=info.zs; k<info.zs+info.zm; k++) {
			for (PetscInt j=info.ys; j<info.ys+info.ym; j++) {
				for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
					// The row number of the element we are considering:
					PetscInt row = (i-info.gxs) + (j-info.gys)*(info.gxm) + (k-info.gzs)*(info.gxm)*(info.gym);
					//
					// 2. Loop over nodes (including ghosts) within a cubic domain with center at (i,j,k)
					//    For each element, run through all elements in a box of size stencilWidth * stencilWidth * stencilWidth 
					//    Remark, we want to make sure we are not running "out of the domain", 
					//    therefore k2 etc. are limited to the max global index (info.mz-1 etc.)
					for (PetscInt k2=PetscMax(k-info.sw,0);k2<=PetscMin(k+info.sw,info.mz-1);k2++){
						for (PetscInt j2=PetscMax(j-info.sw,0);j2<=PetscMin(j+info.sw,info.my-1);j2++){
							for (PetscInt i2=PetscMax(i-info.sw,0);i2<=PetscMin(i+info.sw,info.mx-1);i2++){
								PetscInt col = (i2-info.gxs) + (j2-info.gys)*(info.gxm) + (k2-info.gzs)*(info.gxm)*(info.gym);
								PetscScalar dist = 0.0;
								// Compute the distance from the "col"-element to the "row"-element
								for(PetscInt kk=0; kk<3; kk++){
									dist = dist + PetscPowScalar(lcoorp[3*row+kk]-lcoorp[3*col+kk],2.0);
								}
								dist = PetscSqrtScalar(dist);
								
								// =========================================================================
								// =========== NEW FOR PASSIVE =============================================
								// 3D CANT BEAM
								if (dist<(dx*0.1) && dist<(dy*0.1) && dist<(dz*0.1)){
									VecSetValueLocal(opt->xSolid0,row,valOne, INSERT_VALUES);
									VecSetValueLocal(opt->xSolid1,row,valCero, INSERT_VALUES);																		
									// 2 line of elements around the applied force									
									if ( ((lcoorp[3*row+0])-(opt->xc[0]))<(RMIN) && 
									     ((opt->xc[5])-(lcoorp[3*row+2]))<(RMIN*2.0) ){
										VecSetValueLocal(opt->xSolid0,row,valCero, INSERT_VALUES);
									  VecSetValueLocal(opt->xSolid1,row,valOne, INSERT_VALUES);
									}
									// 2 line of elements in support 
									if ( (-(lcoorp[3*row+0])+(opt->xc[1]))<(RMIN*2.0) && 
									     (-(opt->xc[4])+(lcoorp[3*row+2]))<(RMIN*2.0) ){
									  VecSetValueLocal(opt->xSolid0,row,valCero, INSERT_VALUES);
									  VecSetValueLocal(opt->xSolid1,row,valOne,  INSERT_VALUES);
									}
								}
								// =========== END PASSIVE =================================================
								// =========================================================================
								
								
								
								// ======================================================================
								// ============ NEW FOR MAXIMUM SIZE ====================================
								// INTERMEDIATE :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
								if (dist<=opt->Rmax){
									oneval = 0.0;
									// ########### Symmetry in Y-plane ####################
									if ( (lcoorp[3*row+1] + opt->Rmax) > opt->xc[3]) {
										xjj = lcoorp[3*col];
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->Rmax) && (distijj>=opt->Rint)){
											oneval = oneval + 1.0;
										}
									}
									// ########### Symmetry in X-plane ####################
									if ((lcoorp[3*row+0] - opt->Rmax) < opt->xc[0]) {
										xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1];
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->Rmax) && (distijj>=opt->Rint)){
											oneval = oneval + 1.0;
										}
									}
									// ########### Symmetry in XY-plane ####################
									if (((lcoorp[3*row+0] - opt->Rmax) < opt->xc[0]) && ( (lcoorp[3*row+1] + opt->Rmax) > opt->xc[3]) ){
										xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->Rmax) && (distijj>=opt->Rint)){
											oneval = oneval + 1.0;
										}
									}
									// Normal condition withouth symmetry 
									if (dist>=opt->Rint){
										oneval = oneval + 1.0;
									}
									// Final Weight
									MatSetValuesLocal(DII, 1, &row, 1, &col, &oneval, INSERT_VALUES);
								}
								
								// DILATE ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
								if (dist<=opt->RmaxDil){
									oneval = 0.0;
									// ########### Symmetry in Y-plane ####################
									if ( (lcoorp[3*row+1] + opt->RmaxDil) > opt->xc[3]) {
										xjj = lcoorp[3*col];
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->RmaxDil) && (distijj>=opt->RintDil)){
											oneval = oneval + 1.0;
										}
									}
									// ########### Symmetry in X-plane ####################
									if ((lcoorp[3*row+0] - opt->RmaxDil) < opt->xc[0]) {
										xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1];
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->RmaxDil) && (distijj>=opt->RintDil)){
											oneval = oneval + 1.0;
										}
									}
									// ########### Symmetry in XY-plane ####################
									if (((lcoorp[3*row+0] - opt->RmaxDil) < opt->xc[0]) && ( (lcoorp[3*row+1] + opt->RmaxDil) > opt->xc[3]) ){
										xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if ((distijj<=opt->RmaxDil) && (distijj>=opt->RintDil)){
											oneval = oneval + 1.0;
										}
									}
									// Normal condition withouth symmetry 
									if (dist>=opt->RintDil){
										oneval = oneval + 1.0;
									}
									// Final Weight
									MatSetValuesLocal(DII_DIL, 1, &row, 1, &col, &oneval, INSERT_VALUES);
								}
								
								
								// ================ END FOR MAXIMUM SIZE ================================
								// ======================================================================
								
								
								
								
								if (dist<opt->rmin){						
									
									// Longer distances should have less weight
									dist = 1.0 - dist/opt->rmin;
									
									// ============= SYMETRY CONDITION OF THE FILTER =======================
									// =====================================================================
									// row element is close to the plane Ymax and the col element is mirrored
									if ((lcoorp[3*row+1] + opt->rmin) > opt->xc[3]) 	{
									  xjj = lcoorp[3*col];
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										
										if (distijj<opt->rmin){
											dist = dist + (1.0-distijj/opt->rmin);
									    MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
										else {
											MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
									}
									
									// row element is close to the plane Xmin and the col element is mirrored
									if ((lcoorp[3*row+0] - opt->rmin) < opt->xc[0]) 	{
									  xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										zjj = lcoorp[3*col+2];
										yjj = lcoorp[3*col+1];
										
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										
										if (distijj<opt->rmin){
											dist = dist + (1.0-distijj/opt->rmin);
									    MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
										else {
											MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
									}
									
									if (((lcoorp[3*row+1] + opt->rmin) > opt->xc[3]) && ((lcoorp[3*row+0] - opt->rmin) < opt->xc[0])) {
										xjj = lcoorp[3*col] + 2.0*(opt->xc[0]-lcoorp[3*col]);
										yjj = lcoorp[3*col+1] + 2.0*(opt->xc[3]-lcoorp[3*col+1]);
										zjj = lcoorp[3*col+2];
										distijj = 0.0;
										distijj = distijj + PetscPowScalar(lcoorp[3*row]-xjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+1]-yjj,2.0);
										distijj = distijj + PetscPowScalar(lcoorp[3*row+2]-zjj,2.0);
										distijj = PetscSqrtScalar(distijj);
										if (distijj<opt->rmin){
											dist = dist + (1.0-distijj/opt->rmin);
									    MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
										else {
											MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
										}
									}
									
									// ============= END SYMMETRY OF FILTER ================================
									// =====================================================================
									
									else {
									MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
									}									
								}
							}
						}
					}
				}
			}
		}
		// Assemble H:
		MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
		
		// Compute the Hs, i.e. sum the rows
		Vec dummy;
		VecDuplicate(Hs,&dummy);
		VecSet(dummy,1.0);
		MatMult(H,dummy,Hs);
		
		// This is to take into acount the extension of the filter
	 	PetscScalar maxHs=1;	
		VecMax(Hs,NULL,&maxHs); 
		VecSet(Hs,maxHs);
		
		// ================================================
		// ============== MODIFICATION OF H MATRIX ========
		VecSet(Hs,1.0/maxHs);
		VecSet(dummy,1.0);
		MatDiagonalScale(H,Hs,dummy); // H = diag(1/Hs)*H*diag(1)
				
		// TO TEST the extension of the filter:
		// if the weights are well defined, the filtered field
		// at symmetry planes should have same values than in 
		// the interior of the design. Filtered field at other 
		// boundaries should have lower values. Check this
		// on paraview by plotting the TestFilter field.
		DMCreateGlobalVector(da_elem,&opt->TestFilter);
		VecSet(opt->TestFilter,0);
		VecSet(dummy,1.0);
		MatMult(H,dummy,opt->TestFilter);
		//=================================================
		//=================================================
		
		
		// ================================================
		// =========== NEW FOR PASSIVE ====================
		VecAssemblyBegin(opt->xSolid0);
		VecAssemblyEnd(opt->xSolid0);
		VecAssemblyBegin(opt->xSolid1);
		VecAssemblyEnd(opt->xSolid1);
		// ================================================
		// ================================================		
		
		
		
		// ===============================================================
		// ===================== MAXIMUM SIZE MATRIX =====================
		// ===============================================================
		
		// Assemble DII:
		MatAssemblyBegin(DII, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(DII, MAT_FINAL_ASSEMBLY);
		
		MatAssemblyBegin(DII_DIL, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(DII_DIL, MAT_FINAL_ASSEMBLY);

		// Compute the sum of the neighborhoods ----------------------------------
		Vec DIIsumrow;								// Vector for the neighborhoods
		VecDuplicate(Hs,&DIIsumrow);				// defining the size of the vector
		VecSet(dummy,1.0);							// a vectror of 1
		MatMult(DII,dummy,DIIsumrow);				// to sum the rows
		VecMax(DIIsumrow,NULL,&maxHs);  			// Taking the maximum for extending the domain
		maxHs = 1/maxHs; 							// factor to scale each weight wij
		VecSet(DIIsumrow,maxHs);					// all the test regions with the same size
		// Final matrix to compute the amount of voids
		VecSet(dummy,1.0);
		MatDiagonalScale(DII,DIIsumrow,dummy); 		// DII = diag(a)*DII*diag(b)
		
		VecSet(dummy,1.0);							// a vectror of 1
		MatMult(DII_DIL,dummy,DIIsumrow);			// to sum the rows
		VecMax(DIIsumrow,NULL,&maxHs);  			// Taking the maximum for extending the domain
		maxHs = 1/maxHs; 							// factor to scale each weight wij
		VecSet(DIIsumrow,maxHs);					// all the test regions with the same size
		// Final matrix to compute the amount of voids
		VecSet(dummy,1.0);
		MatDiagonalScale(DII_DIL,DIIsumrow,dummy); 	// DII = diag(a)*DII*diag(b)

		
		
		// Computing the correction of the voids counting in the boundaries
		DMCreateGlobalVector(da_elem,&Avcorr);
		VecDuplicate(Hs,&Avcorr);
		MatMult(DII,dummy,Avcorr); 						// Avcorr = DII*1
		VecAYPX(Avcorr,-1,dummy); 						// Avcorr = 1 - DII*1
		
		DMCreateGlobalVector(da_elem,&Avcorr_DIL);
		VecDuplicate(Hs,&Avcorr_DIL);
		MatMult(DII_DIL,dummy,Avcorr_DIL); 			// Avcorr = DII*1
		VecAYPX(Avcorr_DIL,-1,dummy); 				// Avcorr = 1 - DII*1
		
		
		// This is to initialize the vectors for the derivatives
		VecDuplicate(Hs,&dVvdxPhys); 				// dVvdxPhys
		VecSet(dVvdxPhys,1.0); 			
		VecDuplicate(Hs,&smax); 					// smax
		VecSet(smax,1.0); 			 
		VecDuplicate(Hs,&Vv); 						// Vv
		VecSet(Vv,1);
		VecDuplicate(Hs,&Av); 						// Av
		VecSet(Av,1); 		
		VecDuplicate(Hs,&dAggdx); 					// dAggdx
		VecSet(dAggdx,1); 	

		// destroying
		VecDestroy(&DIIsumrow);
		// ===============================================================
		// ================== END MAXIMUM SIZE MATRIX ====================
		// ===============================================================
		
		
		
		// Clean up
		VecRestoreArray(lcoor,&lcoorp);
		VecDestroy(&dummy);
		delete [] Lx;
		delete [] Ly;
		delete [] Lz;

	} 
	else if (opt->filter==2){
		// ALLOCATE AND SETUP THE PDE FILTER CLASS
		pdef = new PDEFilt(opt);
	}

	return ierr;

}


PetscErrorCode Filter::DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]) {
	PetscErrorCode ierr;
	DM_DA          *da = (DM_DA*)dm->data;
	PetscInt       i,xs,xe,Xs,Xe;
	PetscInt       j,ys,ye,Ys,Ye;
	PetscInt       k,zs,ze,Zs,Ze;
	PetscInt       cnt=0, cell[8], ns=1, nn=8;
	PetscInt       c; 
	if (!da->e) {
		if (da->elementtype == DMDA_ELEMENT_Q1) {ns=1; nn=8;}
		ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xe,&ye,&ze);
		CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(dm,&Xs,&Ys,&Zs,&Xe,&Ye,&Ze);
		CHKERRQ(ierr);
		xe    += xs; Xe += Xs; if (xs != Xs) xs -= 1;
		ye    += ys; Ye += Ys; if (ys != Ys) ys -= 1;
		ze    += zs; Ze += Zs; if (zs != Zs) zs -= 1;
		da->ne = ns*(xe - xs - 1)*(ye - ys - 1)*(ze - zs - 1);
		PetscMalloc((1 + nn*da->ne)*sizeof(PetscInt),&da->e);
		for (k=zs; k<ze-1; k++) {
			for (j=ys; j<ye-1; j++) {
				for (i=xs; i<xe-1; i++) {
					cell[0] = (i-Xs  ) + (j-Ys  )*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[1] = (i-Xs+1) + (j-Ys  )*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[2] = (i-Xs+1) + (j-Ys+1)*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[3] = (i-Xs  ) + (j-Ys+1)*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[4] = (i-Xs  ) + (j-Ys  )*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[5] = (i-Xs+1) + (j-Ys  )*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[6] = (i-Xs+1) + (j-Ys+1)*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[7] = (i-Xs  ) + (j-Ys+1)*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					if (da->elementtype == DMDA_ELEMENT_Q1) {
						for (c=0; c<ns*nn; c++) da->e[cnt++] = cell[c];
					}
				}
			}
		}
	}
	*nel = da->ne;
	*nen = nn;
	*e   = da->e;
	return(0);
}
