#include <Mesh/pzgmesh.h>
#include <Mesh/pzcmesh.h>
#include <Pre/pzgengrid.h>
#include <Material/pzpoisson3d.h>
#include <Analysis/pzanalysis.h>
#include <StrMatrix/pzskylstrmatrix.h>
#ifdef USING_MKL
#include <StrMatrix/TPZSSpStructMatrix.h>
#endif
#include <Matrix/pzstepsolver.h>

//#define PZDEBUG

#ifdef PZDEBUG
#include <Post/TPZVTKGeoMesh.h>
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.discontinuous"));
#endif

//------------------Elasticity Problem------------------------

enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2
};

/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidim5ensional formada por nel elementos de tamanho elsize
 * @param dim dimension of the problem
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(const int dim, int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

/**
 * @brief Routine for creating a computational mesh for the Poisson problem
 * @note It creates the appropriate approximation space
 * @param gmesh geometric mesh
 * @param pOrder polynomial order of the approximation space
 * @param dim dimension of the problem
 * @param matId identifier for the material
 */
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int dim, const int matId);

void Error(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);

//Variáveis globais do problema:

const int gDim = 2; // Dimension of the problem
const int gMatID = 1; // Material of the volumetric element
const int gMatLagrange = -10; // Material of the Lagrange multipliers
const int gMatBCbott = -1, gMatBCtop = -2, gMatBCleft = -3, gMatBCright = -4; // Materials of the boundary conditions
const int gDirichlet = 0, gNeumann = 1, gMixed = 2, gDirichletvar = 4, gPointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
//    EConfig conf = EThiago;
    int initial_p = 1;
    int final_p = 1;
    int initial_h = 1;
    int final_h = 9;
    bool plotting = false;
    EElementType elementType = ESquare;
    int numthreads = 0;

    switch (argc) {
        case 9:
            numthreads = atoi(argv[8]);
        case 8:
            elementType = EElementType(atoi(argv[7]));
        case 7:
            plotting = atoi(argv[6]);
        case 6:
            final_h = atoi(argv[5]);
        case 5:
            initial_h = atoi(argv[4]);
        case 4:
            final_p = atoi(argv[3]);
        case 3:
            initial_p = atoi(argv[2]);
        case 2:
            break;
//            conf = EConfig(atoi(argv[1]));
    };
    int n_ref_p = final_p - initial_p + 1;
    int n_ref_h = final_h - initial_h + 1;

#ifdef USING_MKL
    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
    mkl_set_num_threads(numthreads);
#endif

    std::string rootname;
    double hx = 2, hy = 2; //Dimensões em x e y do domínio
    double x0 = -1;
    double y0 = -1;

//    //Problem data:
//    switch (conf) {
//        case EThiago:
//        case EThiagoPlus:
//        case EThiagoPlusPlus:
//            TElasticityExample1::fProblemType = TElasticityExample1::EThiago;
//            TElasticityExample1::fStressState = TElasticityExample1::EPlaneStrain;
//            TElasticityExample1::fElast = 206.8150271873455;
//            TElasticityExample1::fNu = 0.3040039545229857;
//            hx = 1;
//            hy = 1;
//            x0 = 0;
//            y0 = 0;
//
//            rootname = ConfigRootname[conf] + "_Thiago";
//            break;
//        case EAxiSymmetric:
//        case EAxiSymmetricPlus:
//            TElasticityExample1::fProblemType = TElasticityExample1::Etest1;
//            TElasticityExample1::fStressState = TElasticityExample1::EAxiSymmetric;
//            TElasticityExample1::fElast = 100.;
//            TElasticityExample1::fNu = 0.;
//            hx = 2;
//            hy = 2;
//            x0 = 1;
//            y0 = -1;
//            rootname = ConfigRootname[conf] + "_Test1";
//            break;
//        default:
//            DebugStop();
//            break;
//    }

    for (unsigned int pref = initial_p - 1; pref < final_p; ++pref) {
        for (unsigned int href = initial_h; href <= final_h; ++href) {
            unsigned int h_level = 1 << href;
            unsigned int nelx = h_level, nely = h_level; //Number of elements in x and y directions
            std::cout << "********* " << "Number of h refinements: " << href << " (" << nelx << "x" << nely << " elements). p order: " << pref + 1 << ". *********" << std::endl;
            unsigned int nx = nelx + 1, ny = nely + 1; //Number of nodes in x and y directions
            unsigned int stressPOrder = pref + 1; //Polynomial order of the approximation

            TPZGeoMesh *gmesh = CreateGMesh(gDim, nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh

#ifdef PZDEBUG

            std::string filename("../gmesh.pz");
            {
                std::string meshName("testMesh");
                gmesh->SetName(meshName);
                TPZPersistenceManager::OpenWrite(filename);
                TPZPersistenceManager::WriteToFile(gmesh);
                TPZPersistenceManager::CloseWrite();
                delete gmesh;
                TPZPersistenceManager::OpenRead(filename);
                gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
                std::ofstream fileg("GeoMesh_1.txt"); //Prints the geometric mesh in txt format
                std::ofstream filegvtk("GeoMesh_1.vtk"); //Prints the geometric mesh in vtk format
                gmesh->Print(fileg);
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
            }
#endif
            //Creating computational mesh:
            TPZCompMesh *cmesh = CreateCMesh(gmesh, stressPOrder, gDim, gMatID);

#ifdef PZDEBUG
            {
                std::ofstream filec("CompMesh.txt"); //Prints the computational mesh in txt format
                cmesh->Print(filec);
                std::ofstream fileg1("GeoMesh_2.txt");
                gmesh->Print(fileg1); //Prints the geometric mesh in txt format


                filename = "../cmesh.pz";
                TPZPersistenceManager::OpenWrite(filename);
                TPZPersistenceManager::WriteToFile(cmesh);
                TPZPersistenceManager::WriteToFile(gmesh);//just to show off
                TPZPersistenceManager::CloseWrite();
                delete cmesh;
                delete gmesh;

                TPZPersistenceManager::OpenRead(filename);
                cmesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
                gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
            }
#endif



            //Solving the system:
            bool optimizeBandwidth = true;
            TPZAnalysis an(cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
#ifdef USING_MKL
            TPZSymetricSpStructMatrix matskl(cmesh);
#else
            TPZSkylineStructMatrix matskl(cmesh); // asymmetric case ***
#endif
            matskl.SetNumThreads(numthreads);
            an.SetStructuralMatrix(matskl);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);

            std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << "." << std::endl;
            an.Assemble(); //Assembles the global stiffness matrix (and load vector)
            std::cout << "Assemble finished." << std::endl;

//            TPZManVector<REAL, 6> Errors(cmesh->FindMaterial(gMatID)->NEvalErrors());
//            TElasticityExample1 example;
//            an.SetExact(example.Exact());
            //            an.PostProcessError(Errors,std::cout);

#ifdef PZDEBUG
            //Imprimir Matriz de rigidez Global:
            if (false) {
                std::ofstream filestiff("stiffness.nb");
                an.Solver().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);

                std::ofstream filerhs("rhs.nb");
                an.Rhs().Print("R = ", filerhs, EMathematicaInput);
            }
#endif

            std::cout << "Solving." << std::endl;
            an.Solve();
            std::cout << "Solved." << std::endl;


            {
                TPZStepSolver<STATE> solver;
                an.SetSolver(solver);
            }
#ifdef PZDEBUG
            if (0) {
                std::ofstream file("file.txt");
                an.Solution().Print("sol=", file, EMathematicaInput);

            }
#endif

            if (plotting) {
                std::string plotfile;
                {
                    std::stringstream sout;
                    sout << rootname << ".vtk";
                    plotfile = sout.str();
                }
                TPZStack<std::string> scalnames, vecnames;
                scalnames.Push("Pressure");
                vecnames.Push("Flux");
                int count = href * n_ref_p + pref - (initial_p - 1);
                an.SetStep(count);
                an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
                an.PostProcess(2);
            }

#ifdef PZDEBUG
            //Imprimindo vetor solução:
            {
                TPZFMatrix<STATE> solucao = cmesh->Solution(); //Pegando o vetor de solução, alphaj
                std::ofstream solout("sol.nb");
                solucao.Print("Sol", solout, EMathematicaInput); //Imprime na formatação do Mathematica

                std::ofstream fileAlpha("alpha.nb");
                an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
            }
#endif

            //   matids.clear();
            //   matids.insert(-1);
            //   TPZManVector<STATE,3> result;
            //  result = cmesh_m_HDiv->Integrate("state",matids);
            //  std::cout << "Sigma Y"  << result << std::endl;


            //    //Calculo do erro
            //    std::cout << "Computing Error " << std::endl;

            std::stringstream sout;
            sout << rootname;
            switch (elementType) {
                case ETriangular:
                    sout << "_tria";
                    break;
                case ESquare:
                    sout << "_quad";
                    break;
                case ETrapezoidal:
                    sout << "_trap";
                    break;
            }
//            sout << "_" << stressPOrder << "_Error.nb";
//            std::ofstream ErroOut(sout.str(), std::ios::app);
//            ErroOut << "(* Type of simulation " << rootname << " *)\n";
//            ErroOut << "(* Number of elements " << h_level << " *)" << std::endl;
//            ErroOut << "(* Type of Element ";
//            switch (elementType) {
//                case ETriangular:
//                    ErroOut << "triangular ";
//                    break;
//                case ESquare:
//                    ErroOut << "square ";
//                    break;
//                case ETrapezoidal:
//                    ErroOut << "trapezoidal ";
//                    break;
//            }
//            ErroOut << " *)\n";
//            ErroOut << "(* Number of Condensed equations " << cmesh_m_Hybrid->NEquations() << " *)" << std::endl;
//            ErroOut << "(* Number of equations before condensation " << cmesh_m_Hybrid->Solution().Rows() << " *)" << std::endl;
//            ErroOut << "(*\n";
//            an.SetExact(example.Exact());
//            an.SetThreadsForError(numthreads);
//            bool store_errors = true;
//            cmesh_m_Hybrid->ElementSolution().Redim(cmesh_m_Hybrid->NElements(), Errors.size());
//            std::cout << "Computing errors." << std::endl;
//            an.PostProcessError(Errors, store_errors, ErroOut);
//            std::cout << "Computed errors." << std::endl;
//            ErroOut << "nelx ribporder internalporder n_condensed - n_total - error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as\n";
//            ErroOut << "*)\n";
//            TPZManVector<STATE, 10> output(Errors.size() + 5, 0);
//            output[0] = h_level;
//            output[1] = stressPOrder;
//            output[2] = stressInternalPOrder;
//            output[3] = cmesh_m_Hybrid->NEquations();
//            output[4] = cmesh_m_Hybrid->Solution().Rows();
//            for (int i = 0; i < Errors.size(); i++) {
//                output[5 + i] = Errors[i];
//            }
//            ErroOut << "Error[[" << href + 1 << "," << pref + 1 << "]] = {" << output << "};\n";
//
//            std::cout << "Errors = " << Errors << std::endl;


            an.CleanUp();
            delete cmesh;
            delete gmesh;
        }
    }
    std::cout << "FINISHED!" << std::endl;

    return 0;
}


TPZGeoMesh *CreateGMesh(const int dim, int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType) {
    //Creating geometric mesh, nodes and elements.
    //Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(dim);

    //Auxiliary vector to store coordinates:
    //TPZVec <REAL> coord(3, 0.);
    TPZVec<REAL> gcoord1(3, 0.);
    TPZVec<REAL> gcoord2(3, 0.);
    gcoord1[0] = x0;
    gcoord1[1] = y0;
    gcoord1[2] = 0;
    gcoord2[0] = x0 + hx;
    gcoord2[1] = y0 + hy;
    gcoord2[2] = 0;
    //Inicialização dos nós:

    TPZManVector<int> nelem(2, 1);
    nelem[0] = nelx;
    nelem[1] = nely;

    TPZGenGrid gengrid(nelem, gcoord1, gcoord2);

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

    gengrid.Read(gmesh, gMatID);
    gengrid.SetBC(gmesh, 4, gMatBCbott);
    gengrid.SetBC(gmesh, 5, gMatBCright);
    gengrid.SetBC(gmesh, 6, gMatBCtop);
    gengrid.SetBC(gmesh, 7, gMatBCleft);

    gmesh->BuildConnectivity();

    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }

    //Printing geometric mesh:

    //ofstream bf("before.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

}

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int dim, const int matId) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

//    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
//    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZMaterial *material = new TPZMatPoisson3d(matId, dim); //criando material que implementa a formulacao fraca do problema modelo
    cmesh->InsertMaterialObject(material); //Insere material na malha


    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    cmesh->ApproxSpace().CreateInterfaces(*cmesh);
    return cmesh;

}
