/**********************************************************************
* This example aims to demonstrate how to use p-FEM using the NeoPZ   *
* library. It solves the equation div(-k grad(p)) = f for the varia-  *
* ble p using a Discontinuous Galerkin formulation.                   *
* The domain is a unit square ranging from (0,0) to (1,1), k=-1 and   *
* dirichlet homogeneous conditions are imposed in all boundaries.     *
**********************************************************************/

#include <TPZGmshReader.h>
#include <pzgmesh.h>
#include <pzanalysis.h>
#include <TPZMatLaplacian.h>
#include <pzbndcond.h>
#include <pzstepsolver.h>
#ifdef USING_MKL
#include <StrMatrix/TPZSSpStructMatrix.h>
#endif
#include <Pre/pzgengrid.h>
#include <pzintel.h>

#include <string>


enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2
};

/**
 * @brief Routine from creating the geometric mesh of the domain. It is a square with nodes
 * (0,0), (1,0), (1,1), (0,1).
 * @param dim dimension of the problem
 * @param nelx number of divisions on the x direction
 * @param nely number of divisions on the y direction
 * @param meshType defines with elements will be created (triangular, square, trapezoidal)
 * @param matIds store the material identifiers
 */
TPZGeoMesh *CreateGeoMesh(const int dim, int nelx, int nely, EElementType meshType, TPZVec<int> &matIds);


/**
* Generates a computational mesh that implements the problem to be solved
*/
static TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TPZVec<int> &matIds, const int initialPOrder);

int main(int argc, char **argv)
{
    constexpr int numthreads{8};
    #ifdef USING_MKL
    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
    mkl_set_num_threads(numthreads);
    #endif

    constexpr int dim{2};
    constexpr int nDiv{50};
    TPZVec<int> matIdVec;
    TPZGeoMesh *gMesh = CreateGeoMesh(dim,nDiv,nDiv,ETriangular,matIdVec);


    constexpr int initialPOrder{0};
    TPZCompMesh *cMesh = CreateCompMesh(gMesh,matIdVec,initialPOrder);

    //Setting up the analysis object
    constexpr bool optimizeBandwidth{true};
    TPZAnalysis an(cMesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
#ifdef USING_MKL
    TPZSymetricSpStructMatrix matskl(cMesh);
#else
    TPZSkylineStructMatrix matskl(cMesh); // asymmetric case ***
#endif
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    //setting reference solution
    auto exactSolution = [] (const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &solDx){
        const auto x = pt[0], y = pt[1];
        sol[0] = std::pow(10,-std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x);
    };
    an.SetExact(exactSolution);
    an.SetThreadsForError(numthreads);

    //setting variables for post processing
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("POrder");
    std::string plotfile= "solution.vtk";
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    constexpr int postProcessResolution{0};

    constexpr int nPRefinements{3};
    //global errror values in energy, L2 and H1 (semi) norms
    TPZVec<REAL> errorVec(3,0);
    const int nElems = cMesh->NElements();
    cMesh->ElementSolution().Resize(nElems,3);
    constexpr int nRefEls{20};
    for(auto it = 0 ; it < nPRefinements; it++){
        std::cout<<"============================"<<nPRefinements<<std::endl;
        std::cout<<"\tIteration "<<it+1<<" out of "<<nPRefinements<<std::endl;
        std::cout<<"\tAssemble matrix with NDoF = "<<cMesh->NEquations()<<"."<<std::endl;
        an.Assemble(); //Assembles the global stiffness matrix (and load vector)
        std::cout<<"\tAssemble finished."<<std::endl;
        an.Solve();
        std::cout<<"\tPost processing..."<<std::endl;
        an.PostProcess(postProcessResolution);
        std::cout<<"\tPost processing finished."<<std::endl;
        if(it < nPRefinements - 1){
            an.PostProcessError(errorVec,true);
            auto &allElementErrors = cMesh->ElementSolution();
            TPZManVector<std::pair<int,STATE>,nRefEls> biggestElementErrors(nRefEls,std::make_pair(-1,-1));
            for(int iel = 0 ; iel < nElems; iel++){
                int pos = 0;
                while( biggestElementErrors[pos].second > allElementErrors(iel,0) ){
                    pos++;
                }
                for(int iPos = pos; iPos < nRefEls - 1; iPos++){
                    biggestElementErrors[iPos + 1] = biggestElementErrors[iPos];
                }
                auto id = cMesh->Element(iel)->Index();
                auto elError = allElementErrors(iel,0);
                biggestElementErrors[pos] = std::make_pair(id,elError);
            }
            for(int iel = 0 ; iel < nRefEls; iel++){
                const auto elId = biggestElementErrors[iel].first;
                auto refEl =  dynamic_cast<TPZCompElDisc *> (cMesh->Element(elId));
                const int currentPorder = refEl->GetgOrder();
                refEl->PRefine(currentPorder+1);
            }
        }
    }

    delete cMesh;
    delete gMesh;
    return 0;
};


TPZGeoMesh *CreateGeoMesh(const int dim, int nelx, int nely, EElementType meshType, TPZVec<int> &matIds) {
    //Creating geometric mesh, nodes and elements.
    //Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(dim);

    //Auxiliary vector to store coordinates:
    TPZVec<REAL> coord1(3, 0.);
    TPZVec<REAL> coord2(3, 0.);
    coord1[0] = 0;coord1[1] = 0;coord1[2] = 0;
    coord2[0] = 1;coord2[1] = 1;coord2[2] = 0;

    TPZManVector<int> nelem(2, 1);
    nelem[0] = nelx;
    nelem[1] = nely;

    TPZGenGrid gengrid(nelem, coord1, coord2);

    switch (meshType) {
        case ESquare:
            gengrid.SetElementType(EQuadrilateral);
            break;
        case ETriangular:
            gengrid.SetElementType(ETriangle);
            break;
        case ETrapezoidal:
            gengrid.SetDistortion(0.25);
            break;
    }
    constexpr int matIdDomain = 1, matIdBoundary = 2;
    gengrid.Read(gmesh, matIdDomain);
    gengrid.SetBC(gmesh, 4, matIdBoundary);
    gengrid.SetBC(gmesh, 5, matIdBoundary);
    gengrid.SetBC(gmesh, 6, matIdBoundary);
    gengrid.SetBC(gmesh, 7, matIdBoundary);

    gmesh->BuildConnectivity();

    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    matIds.Resize(2);
    matIds[0] = matIdDomain;
    matIds[1] = matIdBoundary;
    //Printing geometric mesh:

    //ofstream bf("before.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
}

TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TPZVec<int> &matIds, const int initialPOrder){
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    //Definition of the approximation space
    const int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(initialPOrder);
    cmesh->SetDimModel(dim);


    const int matId = matIds[0];
    constexpr REAL perm{-1};

    //Inserting material
    TPZMatLaplacian * mat = new TPZMatLaplacian(matId, dim);
    mat->SetPermeability(perm);
    mat->SetSymmetric();


    auto forcingFunction = [](const TPZVec<REAL>& pt, TPZVec<STATE> &result){
        REAL x = pt[0];
        REAL y = pt[1];
        result[0] =  -9*std::pow(2,1 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*
                    std::pow(5,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x)*std::log(2) -
                    9*std::pow(2,1 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*
                    std::pow(5,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x)*std::log(5) -
                    9*std::pow(2,1 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*
                    std::pow(5,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x)*std::log(10) -
                    9*std::pow(10,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x)*
                    std::log(10) + 9*std::pow(10,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*
                    std::pow(-2*M_PI + 15*x,3)*std::pow(std::log(10),2) +
                    9*std::pow(10,2 - std::pow(-2*M_PI + 15*x,2) - std::pow(-2*M_PI + 15*y,2))*(-2*M_PI + 15*x)*
                    std::pow(-2*M_PI + 15*y,2)*std::pow(std::log(10),2);
    };
    constexpr int pOrderForcingFunction{10};
    mat->SetForcingFunction(forcingFunction,pOrderForcingFunction);

    //Inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat);

    //Boundary conditions
    constexpr int dirichlet = 0;
    constexpr int neumann = 1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);


    const int &matIdBc1 = matIds[1];
    val2(0,0)=0.0;
    auto bc1 = mat->CreateBC(mat, matIdBc1, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc1);

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    cmesh->ApproxSpace().CreateInterfaces(*cmesh);
    return cmesh;
}