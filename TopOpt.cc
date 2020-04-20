#include <TopOpt.h>
#include <cmath>

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

TopOpt::TopOpt(){
      
  x          = NULL;
  xPhys      = NULL;
  xFilter    = NULL;
  dfdx       = NULL;
  dgdx       = NULL;
  dPhysdFilt = NULL;
  gx         = NULL;
  da_nodes   = NULL;
  da_elem    = NULL;
  xo1        = NULL;
  xo2        = NULL;
  U          = NULL;
  L          = NULL;
	
  // ==== PASSIVE elements
  xSolid0    = NULL;
  xSolid1    = NULL;
  // =====================

  // ==== Robust Formulation
  xPhysEro   = NULL;
  xPhysDil   = NULL;
  dErodFilt  = NULL;
  dDildFilt  = NULL;
  // =======================  
  
  // === MAX SIZE ====
  gms=NULL;
  // Robust:
  gmsDIL=NULL;
  // =================
	
	
  TestFilter = NULL;
	
  SetUp();
}

TopOpt::~TopOpt(){
  
  // Delete vectors
  if (x!=NULL){ VecDestroy(&x); }
  if (dfdx!=NULL){ VecDestroy(&dfdx); }
  if (dgdx!=NULL){ VecDestroyVecs(m,&dgdx); }
  // Densities
  if (xPhys!=NULL){ VecDestroy(&xPhys); }
  if (xold!=NULL){ VecDestroy(&xold); }
  if (xmin!=NULL){ VecDestroy(&xmin); }
  if (xmax!=NULL){ VecDestroy(&xmax); }
    
  if (da_nodes!=NULL){ DMDestroy(&(da_nodes)); }
  if (da_elem!=NULL){ DMDestroy(&(da_elem)); }
    
  // Delete constraints
  if (gx!=NULL){ delete [] gx; }
    
  // mma restart method    		
  if (xo1!=NULL){ VecDestroy(&xo1); }
  if (xo2!=NULL){ VecDestroy(&xo2); }
  if (L!=NULL){ VecDestroy(&L); }
  if (U!=NULL){ VecDestroy(&U);  }

  // PROJECTION ====================================
  if (xFilter!=NULL){ VecDestroy(&xFilter); }
  if (dPhysdFilt!=NULL){ VecDestroy(&dPhysdFilt); }
  // ===============================================
    
	// ====================== PASSIVE elements
	if (xSolid0!=NULL){ VecDestroy(&xSolid0); }
	if (xSolid1!=NULL){ VecDestroy(&xSolid1); }
	// =======================================
	
	// ========================== ROBUST FORMULATION 
	if (xPhysEro!=NULL){ VecDestroy(&xPhysEro); }
	if (xPhysDil!=NULL){ VecDestroy(&xPhysDil); }
	if (dErodFilt!=NULL){ VecDestroy(&dErodFilt); }
	if (dDildFilt!=NULL){ VecDestroy(&dDildFilt); }
	// ================================================
	
	// ======================= Max Size
	if (gms!=NULL){ VecDestroy(&gms); }
	if (gmsDIL!=NULL){ VecDestroy(&gmsDIL); }
	// ================================
	
	if (TestFilter!=NULL){ VecDestroy(&TestFilter); }
}


// NO METHODS !
//PetscErrorCode TopOpt::SetUp(Vec CRAPPY_VEC){
PetscErrorCode TopOpt::SetUp(){
	PetscErrorCode ierr;
  
	// SET DEFAULTS for FE mesh and levels for MG solver
        nxyz[0] = 144+1;
        nxyz[1] = 24+1;
        nxyz[2] = 48+1;
        xc[0]   = 0.0;
        xc[1]   = 3.0;
        xc[2]   = 0.0;
        xc[3]   = 0.5;
        xc[4]   = 0.0;
        xc[5]   = 1.0;
		nu      = 0.3;
        nlvls   = 4;

	// SET DEFAULTS for optimization problems
	volfrac= 0.25       ;  // Volume fraction on Dilated (it's scaled on the main file)
	volfracREF= volfrac ;  // Volume fraction in Intermediate
	maxItr = 449        ;  // Maximum number of iterations (check IterProj above: 40*9=360)
	rmin   = (3.0/(144.0)*1.01)*2.0 ; // Filter radius (it is not the minimum size radius)
	Emin   = 1.0e-6     ; // 
	Emax   = 1.0        ;
	filter = 3          ; // Option 3 = density + Heaviside
	Xmin   = 0.0        ;
	Xmax   = 1.0        ;	
		
	// SET PARAMETERS OF THE PROJECTION ========
	Beta    = 1.5; // Initial value of Beta
	IterProj= 50 ; // Iterations by step
	penal   = 1.0; // Initial value for penal
	// -----------------------------------------
	
	// Robust formulation =======================
	etaEro  = 0.75 ; // erosion threshold
	eta     = 0.50 ; // Intermediate threshold
	etaDil  = 0.25 ; // dilation threshold
	tero    = rmin/2.0*0.58 ; // offset distance
	tdil    = rmin/2.0*0.58 ; // offset distance
	// ==========================================
	
	// MAXIMUM SIZE PARAMETERS ==================
	Rmax     = (3.0/(144.0)*3.01);   // Maximum size radius (intermediate)
	Rint     = (3.0/(144.0)*1.01);   // Minimum size radius (inner ring)
	RmaxDil  = Rmax + tdil;          // Maximum size radius in dilated design
	RintDil  = Rint + tdil;          // Minimum size radius in dilated design
	epsiMS   = 0.05;                 // Epsilon parameter
	pagg     = 100.0;                // Aggregation exponent
	Nagg     = 144.0*24.0*48.0;      // number of data (for p-mean)
	m        = 3;                    // volume and MaxSize (interm. and dilated)
	movlim   = 0.5;                  // move limits MMA (interpolated in the main file)
	movlimIni= 0.5;                  // Initial move limits
	movlimEnd= 0.05;                 // final move limits
	// ==========================================
	
	restart = PETSC_TRUE;
	ierr = SetUpMESH(); CHKERRQ(ierr);
	ierr = SetUpOPT(); CHKERRQ(ierr);
	return(ierr);
}


PetscErrorCode TopOpt::SetUpMESH(){
	
	PetscErrorCode ierr;
	
	// Read input from arguments
	PetscBool flg;
	
	// Physics parameters
	PetscOptionsGetInt(NULL,NULL,"-nx",&(nxyz[0]),&flg);
	PetscOptionsGetInt(NULL,NULL,"-ny",&(nxyz[1]),&flg);
	PetscOptionsGetInt(NULL,NULL,"-nz",&(nxyz[2]),&flg);
	PetscOptionsGetReal(NULL,NULL,"-xcmin",&(xc[0]),&flg);	
	PetscOptionsGetReal(NULL,NULL,"-xcmax",&(xc[1]),&flg);
	PetscOptionsGetReal(NULL,NULL,"-ycmin",&(xc[2]),&flg);
	PetscOptionsGetReal(NULL,NULL,"-ycmax",&(xc[3]),&flg);
	PetscOptionsGetReal(NULL,NULL,"-zcmin",&(xc[4]),&flg);
	PetscOptionsGetReal(NULL,NULL,"-zcmax",&(xc[5]),&flg);
   PetscOptionsGetReal(NULL,NULL,"-penal",&penal,&flg);
	PetscOptionsGetInt(NULL,NULL,"-nlvls",&nlvls,&flg);

	
	// Write parameters for the physics _ OWNED BY TOPOPT
	PetscPrintf(PETSC_COMM_WORLD,"########################################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"############################ FEM settings ##############################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Number of nodes: (-nx,-ny,-nz):        (%i,%i,%i) \n",nxyz[0],nxyz[1],nxyz[2]);
        PetscPrintf(PETSC_COMM_WORLD,"# Number of degree of freedom:           %i \n",3*nxyz[0]*nxyz[1]*nxyz[2]);
	PetscPrintf(PETSC_COMM_WORLD,"# Number of elements:                    (%i,%i,%i) \n",nxyz[0]-1,nxyz[1]-1,nxyz[2]-1);
	PetscPrintf(PETSC_COMM_WORLD,"# Dimensions: (-xcmin,-xcmax,..,-zcmax): (%f,%f,%f)\n",xc[1]-xc[0],xc[3]-xc[2],xc[5]-xc[4]);
	PetscPrintf(PETSC_COMM_WORLD,"# -nlvls: %i\n",nlvls);
       	PetscPrintf(PETSC_COMM_WORLD,"########################################################################\n");

	// Check if the mesh supports the chosen number of MG levels
	PetscScalar divisor = PetscPowScalar(2.0,(PetscScalar)nlvls-1.0);
	// x - dir
	if ( std::floor((PetscScalar)(nxyz[0]-1)/divisor) != (nxyz[0]-1.0)/((PetscInt)divisor) ) {
		PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"X - number of nodes %i is cannot be halfened %i times\n",nxyz[0],nlvls-1);
		exit(0);
	}	
	// y - dir
        if ( std::floor((PetscScalar)(nxyz[1]-1)/divisor) != (nxyz[1]-1.0)/((PetscInt)divisor) ) {
                PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"Y - number of nodes %i is cannot be halfened %i times\n",nxyz[1],nlvls-1);
		exit(0);
        }
	// z - dir
        if ( std::floor((PetscScalar)(nxyz[2]-1)/divisor) != (nxyz[2]-1.0)/((PetscInt)divisor) ) {
                PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"Z - number of nodes %i is cannot be halfened %i times\n",nxyz[2],nlvls-1);
		exit(0);
        }


	// Start setting up the FE problem
	// Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_PERIODIC
	DMBoundaryType bx = DM_BOUNDARY_NONE;
	DMBoundaryType by = DM_BOUNDARY_NONE;
	DMBoundaryType bz = DM_BOUNDARY_NONE;

	// Stencil type - box since this is closest to FEM (i.e. STAR is FV/FD)
	DMDAStencilType  stype = DMDA_STENCIL_BOX;

	// Discretization: nodes:
	// For standard FE - number must be odd
	// FOr periodic: Number must be even
	PetscInt nx = nxyz[0];
	PetscInt ny = nxyz[1];
	PetscInt nz = nxyz[2];

	// number of nodal dofs
	PetscInt numnodaldof = 3;

	// Stencil width: each node connects to a box around it - linear elements
	PetscInt stencilwidth = 1;

	// Coordinates and element sizes: note that dx,dy,dz are half the element size
	PetscReal xmin=xc[0], xmax=xc[1], ymin=xc[2], ymax=xc[3], zmin=xc[4], zmax=xc[5];
	dx = (xc[1]-xc[0])/(PetscScalar(nxyz[0]-1));
	dy = (xc[3]-xc[2])/(PetscScalar(nxyz[1]-1));
	dz = (xc[5]-xc[4])/(PetscScalar(nxyz[2]-1));

	// Create the nodal mesh
	ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
			numnodaldof,stencilwidth,0,0,0,&(da_nodes));
	CHKERRQ(ierr);

	// Set the coordinates
	ierr = DMDASetUniformCoordinates(da_nodes, xmin,xmax, ymin,ymax, zmin,zmax);
	CHKERRQ(ierr);
	
	// Set the element type to Q1: Otherwise calls to GetElements will change to P1 !
	// STILL DOESN*T WORK !!!!
	ierr = DMDASetElementType(da_nodes, DMDA_ELEMENT_Q1);
	CHKERRQ(ierr);
	
	// Create the element mesh: NOTE THIS DOES NOT INCLUDE THE FILTER !!!
	// find the geometric partitioning of the nodal mesh, so the element mesh will coincide 
	// with the nodal mesh
	PetscInt md,nd,pd; 
	ierr = DMDAGetInfo(da_nodes,NULL,NULL,NULL,NULL,&md,&nd,&pd,NULL,NULL,NULL,NULL,NULL,NULL);
	CHKERRQ(ierr);
	
	// vectors with partitioning information
	PetscInt *Lx=new PetscInt[md];
	PetscInt *Ly=new PetscInt[nd];
	PetscInt *Lz=new PetscInt[pd];

	// get number of nodes for each partition
	const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
	ierr = DMDAGetOwnershipRanges(da_nodes, &LxCorrect, &LyCorrect, &LzCorrect); 
	CHKERRQ(ierr);
	
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

	// Create the element grid: NOTE CONNECTIVITY IS 0
	PetscInt conn = 0;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx-1,ny-1,nz-1,md,nd,pd,
			1,conn,Lx,Ly,Lz,&(da_elem));
	CHKERRQ(ierr);
	
	// Set element center coordinates
	ierr = DMDASetUniformCoordinates(da_elem , xmin+dx/2.0,xmax-dx/2.0, ymin+dy/2.0,ymax-dy/2.0, zmin+dz/2.0,zmax-dz/2.0);
	CHKERRQ(ierr);

	// Clean up
	delete [] Lx;
	delete [] Ly;
	delete [] Lz;
  
  
  	return(ierr);
}

PetscErrorCode TopOpt::SetUpOPT(){
  
	PetscErrorCode ierr;
  
	//ierr = VecDuplicate(CRAPPY_VEC,&xPhys); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_elem,&xPhys);  CHKERRQ(ierr);

	// FOR PROJECTION ============================================
	ierr = DMCreateGlobalVector(da_elem,&xFilter);  CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_elem,&dPhysdFilt);  CHKERRQ(ierr);
	// ===========================================================
	
	// FOR ROBUST FORMULATION ===========
	DMCreateGlobalVector(da_elem,&xPhysEro);
  DMCreateGlobalVector(da_elem,&xPhysDil);
	// ==================================
	
	// MAXIMUM SIZE ==================
	DMCreateGlobalVector(da_elem,&gms);
	DMCreateGlobalVector(da_elem,&gmsDIL);
	// ===============================

	// Total number of design variables
	VecGetSize(xPhys,&n);

	PetscBool flg;
	
	// Optimization paramteres
	PetscOptionsGetReal(NULL,NULL,"-Emin",&Emin,&flg);
	PetscOptionsGetReal(NULL,NULL,"-Emax",&Emax,&flg);
	PetscOptionsGetReal(NULL,NULL,"-volfrac",&volfrac,&flg);
   PetscOptionsGetReal(NULL,NULL,"-penal",&penal,&flg);
	PetscOptionsGetReal(NULL,NULL,"-rmin",&rmin,&flg);
	PetscOptionsGetInt(NULL,NULL,"-maxItr",&maxItr,&flg);
	PetscOptionsGetInt(NULL,NULL,"-filter",&filter,&flg);
	PetscOptionsGetReal(NULL,NULL,"-Xmin",&Xmin,&flg);
        PetscOptionsGetReal(NULL,NULL,"-Xmax",&Xmax,&flg);
	PetscOptionsGetReal(NULL,NULL,"-movlim",&movlim,&flg);
        
	PetscPrintf(PETSC_COMM_WORLD,"################### Optimization settings ####################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Problem size: n= %i, m= %i\n",n,m);
	PetscPrintf(PETSC_COMM_WORLD,"# -filter: %i  (0=sens., 1=dens, 2=PDE)\n",filter);
	PetscPrintf(PETSC_COMM_WORLD,"# -rmin: %f\n",rmin);
	PetscPrintf(PETSC_COMM_WORLD,"# -volfrac: %f\n",volfrac);
        PetscPrintf(PETSC_COMM_WORLD,"# -penal: %f\n",penal);
	PetscPrintf(PETSC_COMM_WORLD,"# -Emin/-Emax: %e - %e \n",Emin,Emax);
	PetscPrintf(PETSC_COMM_WORLD,"# -maxItr: %i\n",maxItr);
	PetscPrintf(PETSC_COMM_WORLD,"# -movlim: %f\n",movlim);
       	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");

        // Allocate after input
        gx = new PetscScalar[m];
	if (filter==0){
		Xmin = 0.001; // Prevent division by zero in filter
	}
	
	// Allocate the optimization vectors
	ierr = VecDuplicate(xPhys,&x); CHKERRQ(ierr);
	VecSet(x,volfrac); // Initialize to volfrac !
	VecSet(xPhys,volfrac); // Initialize to volfrac !

	// PROJECTION =============================
	VecSet(xFilter,volfrac);
	VecDuplicate(x,&dPhysdFilt); 
	// ========================================

	// ROBUST FORMULATION =====================
	VecSet(xPhysEro,volfrac);
	VecSet(xPhysDil,volfrac);
	VecDuplicate(x,&dPhysdFilt);
	VecDuplicate(x,&dErodFilt);
	VecDuplicate(x,&dDildFilt);
	// ========================================
	
	// MAXIMUM SIZE =======
	VecSet(gms,volfrac);
	VecSet(gmsDIL,volfrac);
	// ====================
	
	// Sensitivity vectors
	ierr = VecDuplicate(x,&dfdx); CHKERRQ(ierr);
	ierr = VecDuplicateVecs(x,m, &dgdx); CHKERRQ(ierr);

	// Bounds and 
	VecDuplicate(x,&xmin);
	VecDuplicate(x,&xmax);
	VecDuplicate(x,&xold);	
	VecSet(xold,volfrac);
  
  	return(ierr);
}

PetscErrorCode TopOpt::AllocateMMAwithRestart(PetscInt *itr, MMA **mma)  {

	PetscErrorCode ierr = 0;

	// Set MMA parameters (for multiple load cases)
	PetscScalar aMMA[m];
	PetscScalar cMMA[m];
	PetscScalar dMMA[m];
	for (PetscInt i=0;i<m;i++){
	    aMMA[i]=0.0;
	    dMMA[i]=0.0;
	    cMMA[i]=1000.0;
	}

	// Check if restart is desired
	restart = PETSC_TRUE; // DEFAULT USES RESTART
	flip = PETSC_TRUE;     // BOOL to ensure that two dump streams are kept
	PetscBool onlyLoadDesign = PETSC_FALSE; // Default restarts everything

	// Get inputs
	PetscBool flg;
	char filenameChar[PETSC_MAX_PATH_LEN];
	PetscOptionsGetBool(NULL,NULL,"-restart",&restart,&flg);
	PetscOptionsGetBool(NULL,NULL,"-onlyLoadDesign",&onlyLoadDesign,&flg);

	if (restart) {
	  ierr = VecDuplicate(x,&xo1); CHKERRQ(ierr);
	  ierr = VecDuplicate(x,&xo2); CHKERRQ(ierr);
	  ierr = VecDuplicate(x,&U); CHKERRQ(ierr);
	  ierr = VecDuplicate(x,&L); CHKERRQ(ierr);
	}
	
	// Determine the right place to write the new restart files
	std::string filenameWorkdir = "./";
	PetscOptionsGetString(NULL,NULL,"-workdir",filenameChar,sizeof(filenameChar),&flg);
	if (flg){
		filenameWorkdir = "";
		filenameWorkdir.append(filenameChar);
	}
	filename00 = filenameWorkdir;
	filename00Itr = filenameWorkdir;
	filename01 = filenameWorkdir;
	filename01Itr = filenameWorkdir;

	filename00.append("/Restart00.dat");
	filename00Itr.append("/Restart00_itr_f0.dat");
	filename01.append("/Restart01.dat");
	filename01Itr.append("/Restart01_itr_f0.dat");

	// Where to read the restart point from
	std::string restartFileVec = ""; // NO RESTART FILE !!!!!
	std::string restartFileItr = ""; // NO RESTART FILE !!!!!

	PetscOptionsGetString(NULL,NULL,"-restartFileVec",filenameChar,sizeof(filenameChar),&flg);
	if (flg) {
	   restartFileVec.append(filenameChar);
	}
	PetscOptionsGetString(NULL,NULL,"-restartFileItr",filenameChar,sizeof(filenameChar),&flg);
	if (flg) {
		restartFileItr.append(filenameChar);
	}

	// Which solution to use for restarting
	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Continue from previous iteration (-restart): %i \n",restart);
	PetscPrintf(PETSC_COMM_WORLD,"# Restart file (-restartFileVec): %s \n",restartFileVec.c_str());
	PetscPrintf(PETSC_COMM_WORLD,"# Restart file (-restartFileItr): %s \n",restartFileItr.c_str());
	PetscPrintf(PETSC_COMM_WORLD,"# New restart files are written to (-workdir): %s (Restart0x.dat and Restart0x_itr_f0.dat) \n",filenameWorkdir.c_str());

	// Check if files exist:
	PetscBool vecFile = fexists(restartFileVec);
	if (!vecFile) { PetscPrintf(PETSC_COMM_WORLD,"File: %s NOT FOUND \n",restartFileVec.c_str()); }
	PetscBool itrFile = fexists(restartFileItr);
	if (!itrFile) { PetscPrintf(PETSC_COMM_WORLD,"File: %s NOT FOUND \n",restartFileItr.c_str()); }
	
	// Read from restart point
	
	PetscInt nGlobalDesignVar;
	VecGetSize(x,&nGlobalDesignVar); // ASSUMES THAT SIZE IS ALWAYS MATCHED TO CURRENT MESH
	if (restart && vecFile && itrFile){
		
		PetscViewer view;
		// Open the data files 
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,restartFileVec.c_str(),FILE_MODE_READ,&view);	
				
		VecLoad(x,view);
		VecLoad(xPhys,view);
		VecLoad(xo1,view);
		VecLoad(xo2,view);
		VecLoad(U,view);
		VecLoad(L,view);
		PetscViewerDestroy(&view);
		
		// Read iteration and fscale
		std::fstream itrfile(restartFileItr.c_str(), std::ios_base::in);
		itrfile >> itr[0];
		itrfile >> fscale;
		
		
		// Choose if restart is full or just an initial design guess
		if (onlyLoadDesign){
			PetscPrintf(PETSC_COMM_WORLD,"# Loading design from file: %s \n",restartFileVec.c_str());
			*mma = new MMA(nGlobalDesignVar,m,x, aMMA, cMMA, dMMA);
		}
		else {
			PetscPrintf(PETSC_COMM_WORLD,"# Continue optimization from file: %s \n",restartFileVec.c_str());
			*mma = new MMA(nGlobalDesignVar,m,*itr,xo1,xo2,U,L,aMMA,cMMA,dMMA);
		}

		PetscPrintf(PETSC_COMM_WORLD,"# Successful restart from file: %s and %s \n",restartFileVec.c_str(),restartFileItr.c_str());
	}
	else {
		*mma = new MMA(nGlobalDesignVar,m,x,aMMA,cMMA,dMMA);
	}  

	return ierr;
} 


PetscErrorCode TopOpt::WriteRestartFiles(PetscInt *itr, MMA *mma) {

	PetscErrorCode ierr=0;
	// Only dump data if correct allocater has been used
	if (!restart){
		return -1;
	}

	// Get restart vectors
	mma->Restart(xo1,xo2,U,L);
	
	// Choose previous set of restart files
	if (flip){ flip = PETSC_FALSE; 	}	
	else {     flip = PETSC_TRUE; 	}

	// Write file with iteration number of f0 scaling
	// and a file with the MMA-required vectors, in the following order:
	// : x,xPhys,xold1,xold2,U,L
	PetscViewer view; // vectors
	PetscViewer restartItrF0; // scalars
	
	PetscViewerCreate(PETSC_COMM_WORLD, &restartItrF0);
	PetscViewerSetType(restartItrF0, PETSCVIEWERASCII);
	PetscViewerFileSetMode(restartItrF0, FILE_MODE_WRITE);
	
	// Open viewers for writing
	if (!flip){
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename00.c_str(),FILE_MODE_WRITE,&view);
		PetscViewerFileSetName(restartItrF0, filename00Itr.c_str());
	}
	else if (flip){
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename01.c_str(),FILE_MODE_WRITE,&view);
		PetscViewerFileSetName(restartItrF0, filename01Itr.c_str());
	}

	// Write iteration and fscale
	PetscViewerASCIIPrintf(restartItrF0, "%d ", itr[0]);
	PetscViewerASCIIPrintf(restartItrF0," %e",fscale);
	PetscViewerASCIIPrintf(restartItrF0,"\n");

	// Write vectors
	VecView(x,view); // the design variables
	VecView(xPhys,view);
	VecView(xo1,view);
	VecView(xo2,view);
	VecView(U,view);	
	VecView(L,view);
	
	// Clean up
	PetscViewerDestroy(&view);
	PetscViewerDestroy(&restartItrF0);

	//PetscPrintf(PETSC_COMM_WORLD,"DONE WRITING DATA\n");
	return ierr;
}
