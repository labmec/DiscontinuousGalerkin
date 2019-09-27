#include <Mesh/pzgmesh.h>
#include <Mesh/pzcmesh.h>
#include <Pre/pzgengrid.h>
#include <Material/pzpoisson3d.h>

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
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

/**
 * @brief Funcao para criar a malha computacional
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder);

void Error(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);

//Variáveis globais do problema:

const int dim = 2; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; // Materials of the boundary conditions
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType) {
    //Creating geometric mesh, nodes and elements.
    //Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

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
        case ETriangular:
            gengrid.SetElementType(ETriangle);
            break;
        case ETrapezoidal:
            gengrid.SetDistortion(0.25);
            break;
    }

    gengrid.Read(gmesh, matID);
    gengrid.SetBC(gmesh, 4, matBCbott);
    gengrid.SetBC(gmesh, 5, matBCright);
    gengrid.SetBC(gmesh, 6, matBCtop);
    gengrid.SetBC(gmesh, 7, matBCleft);

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

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZMaterial *material = new TPZMatPoisson3d(matID, dim); //criando material que implementa a formulacao fraca do problema modelo

    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}


int main(int argc, char *argv[]) {
//    TPZMaterial::gBigNumber = 1.e16;
//
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
//    EConfig conf = EThiago;
//    int initial_p = 1;
//    int final_p = 1;
//    int initial_h = 1;
//    int final_h = 9;
//    bool plotting = false;
//    EElementType elementType = ESquare;
//    int numthreads = 8;
//
//    switch (argc) {
//        case 9:
//            numthreads = atoi(argv[8]);
//        case 8:
//            elementType = EElementType(atoi(argv[7]));
//        case 7:
//            plotting = atoi(argv[6]);
//        case 6:
//            final_h = atoi(argv[5]);
//        case 5:
//            initial_h = atoi(argv[4]);
//        case 4:
//            final_p = atoi(argv[3]);
//        case 3:
//            initial_p = atoi(argv[2]);
//        case 2:
//            conf = EConfig(atoi(argv[1]));
//    };
//    int n_ref_p = final_p - initial_p + 1;
//    int n_ref_h = final_h - initial_h + 1;
//
//#ifdef USING_MKL
//    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
//    mkl_set_num_threads(numthreads);
//#endif
//
//    std::string rootname;
//    double hx = 2, hy = 2; //Dimensões em x e y do domínio
//    double x0 = -1;
//    double y0 = -1;
//
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
//    //    TPZManVector<STATE, 2> force(2);
//    //    TPZFNMatrix<4, STATE> sigma(2, 2);
//    //    TPZManVector<REAL, 3> x(3, 0.);
//    //    x[0] = x0 + hx / 2.;
//    //    x[1] = y0 + hy / 2.;
//    //    TElasticityExample1::Sigma(x, sigma);
//    //    TElasticityExample1::Force(x, force);
//
//    for (unsigned int pref = initial_p - 1; pref < final_p; ++pref) {
//        for (unsigned int href = initial_h; href <= final_h; ++href) {
//            unsigned int h_level = 1 << href;
//            unsigned int nelx = h_level, nely = h_level; //Number of elements in x and y directions
//            std::cout << "********* " << "Number of h refinements: " << href << " (" << nelx << "x" << nely << " elements). p order: " << pref + 1 << ". *********" << std::endl;
//            unsigned int nx = nelx + 1, ny = nely + 1; //Number of nodes in x and y directions
//            unsigned int stressPOrder = pref + 1; //Polynomial order of the approximation
//            int stressInternalPOrder = stressPOrder; //k
//            if (conf == EThiagoPlus || conf == EAxiSymmetricPlus) {
//                stressInternalPOrder += 1; //k+1
//            }
//            if (conf == EThiagoPlusPlus) {
//                stressInternalPOrder += 2; //k+2
//            }
//            int displacementPOrder = elementType == ETriangular ? stressInternalPOrder - 1 : stressInternalPOrder;
//            int rotationPOrder = displacementPOrder;
//            TPZGeoMesh *gmesh = CreateGMesh(nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh
//
//#ifdef PZDEBUG
//            std::ofstream fileg("MalhaGeo.txt"); //Prints the geometric mesh in txt format
//            std::ofstream filegvtk("MalhaGeo.vtk"); //Prints the geometric mesh in vtk format
//            gmesh->Print(fileg);
//            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
//#endif
//            //Creating computational mesh:
//            TPZCompMesh *cmesh_S_HDiv = CMesh_S(gmesh, stressPOrder); //Creates the computational mesh for the stress field
//            ChangeInternalOrder(cmesh_S_HDiv, stressInternalPOrder);
//            TPZCompMesh *cmesh_U_HDiv = CMesh_U(gmesh, displacementPOrder); //Creates the computational mesh for the displacement field
//            TPZCompMesh *cmesh_P_HDiv = CMesh_P(gmesh, rotationPOrder, hx / nelx); //Creates the computational mesh for the rotation field
//
//
//            //TPZCompMesh *cmesh_m_HDiv = CMesh_Girk(gmesh, RibpOrder); //Creates the multi-physics computational mesh
//            TPZCompMesh *cmesh_m_HDiv = CMesh_m(gmesh, stressInternalPOrder);
//            //TPZCompMesh *cmesh_m_HDiv = CMesh_AxiS(gmesh, InternalpOrder,  Example);
//#ifdef PZDEBUG
//            {
//                std::ofstream filecS("MalhaC_S.txt"); //Prints the stress computational mesh in txt format
//                std::ofstream filecU("MalhaC_U.txt"); //Prints the displacement computational mesh in txt format
//                std::ofstream filecP("MalhaC_P.txt"); //Prints the rotation computational mesh in txt format
//                cmesh_S_HDiv->Print(filecS);
//                cmesh_U_HDiv->Print(filecU);
//                cmesh_P_HDiv->Print(filecP);
//            }
//#endif
//
//            TPZManVector<TPZCompMesh*, 3> meshvector_HDiv(3);
//            meshvector_HDiv[0] = cmesh_S_HDiv;
//            meshvector_HDiv[1] = cmesh_U_HDiv;
//            meshvector_HDiv[2] = cmesh_P_HDiv;
//            TPZBuildMultiphysicsMesh::AddElements(meshvector_HDiv, cmesh_m_HDiv);
//            TPZBuildMultiphysicsMesh::AddConnects(meshvector_HDiv, cmesh_m_HDiv);
//            TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector_HDiv, cmesh_m_HDiv);
//            cmesh_m_HDiv->LoadReferences();
//
//            //            AddMultiphysicsInterfaces(*cmesh_m_HDiv);
//
//            //CreateCondensedElements(cmesh_m_HDiv);
//
//#ifdef PZDEBUG
//            std::ofstream fileg1("MalhaGeo2.txt");
//            gmesh->Print(fileg1); //Prints the geometric mesh in txt format
//
//            std::ofstream filecm("MalhaC_m.txt");
//            cmesh_m_HDiv->Print(filecm); //Prints the multi-physics computational mesh in txt format
//#endif
//
//            //Solving the system:
//            bool optimizeBandwidth = true;
//            cmesh_m_HDiv->InitializeBlock();
//
//            TPZCompMesh * cmesh_m_Hybrid;
//            TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
//            TPZHybridizeHDiv hybridizer;
//            tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmesh_m_HDiv, meshvector_HDiv, true, -1.);
//            cmesh_m_Hybrid->InitializeBlock();
//
//            TPZAnalysis an(cmesh_m_Hybrid, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
//#ifdef USING_MKL
//            TPZSymetricSpStructMatrix matskl(cmesh_m_Hybrid);
//#else
//            TPZSkylineStructMatrix matskl(cmesh_m_Hybrid); // asymmetric case ***
//#endif
//            matskl.SetNumThreads(numthreads);
//            an.SetStructuralMatrix(matskl);
//            TPZStepSolver<STATE> step;
//            step.SetDirect(ELDLt);
//            an.SetSolver(step);
//
//            std::cout << "Assemble matrix with NDoF = " << cmesh_m_Hybrid->NEquations() << "." << std::endl;
//            an.Assemble(); //Assembles the global stiffness matrix (and load vector)
//            std::cout << "Assemble finished." << std::endl;
//
//            TPZManVector<REAL, 6> Errors(cmesh_m_HDiv->FindMaterial(matID)->NEvalErrors());
//            TElasticityExample1 example;
//            an.SetExact(example.Exact());
//            //            an.PostProcessError(Errors,std::cout);
//
//#ifdef PZDEBUG
//            //Imprimir Matriz de rigidez Global:
//            if (false) {
//                std::ofstream filestiff("stiffness.nb");
//                an.Solver().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);
//
//                std::ofstream filerhs("rhs.nb");
//                an.Rhs().Print("R = ", filerhs, EMathematicaInput);
//            }
//#endif
//
//            std::cout << "Solving." << std::endl;
//            an.Solve();
//            std::cout << "Solved." << std::endl;
//
//
//            {
//                TPZStepSolver<STATE> solver;
//                an.SetSolver(solver);
//            }
//#ifdef PZDEBUG
//            if (0) {
//                std::ofstream file("file.txt");
//                an.Solution().Print("sol=", file, EMathematicaInput);
//
//            }
//#endif
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector_Hybrid, cmesh_m_Hybrid);
//
//            if (plotting) {
//                std::string plotfile;
//                {
//                    std::stringstream sout;
//                    sout << rootname << ".vtk";
//                    plotfile = sout.str();
//                }
//                TPZStack<std::string> scalnames, vecnames;
//                scalnames.Push("SigmaX");
//                scalnames.Push("SigmaY");
//                scalnames.Push("TauXY");
//                vecnames.Push("Flux");
//                vecnames.Push("displacement");
//                vecnames.Push("Stress");
//                int count = href * n_ref_p + pref - (initial_p - 1);
//                an.SetStep(count);
//                an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
//                an.PostProcess(2);
//            }
//
//#ifdef PZDEBUG
//            //Imprimindo vetor solução:
//            {
//                TPZFMatrix<STATE> solucao = cmesh_m_Hybrid->Solution(); //Pegando o vetor de solução, alphaj
//                std::ofstream solout("sol.nb");
//                solucao.Print("Sol", solout, EMathematicaInput); //Imprime na formatação do Mathematica
//
//                std::ofstream fileAlpha("alpha.nb");
//                an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
//            }
//#endif
//
//            //   matids.clear();
//            //   matids.insert(-1);
//            //   TPZManVector<STATE,3> result;
//            //  result = cmesh_m_HDiv->Integrate("state",matids);
//            //  std::cout << "Sigma Y"  << result << std::endl;
//
//
//            //    //Calculo do erro
//            //    std::cout << "Computing Error " << std::endl;
//
//            std::stringstream sout;
//            sout << rootname;
//            switch (elementType) {
//                case ETriangular:
//                    sout << "_tria";
//                    break;
//                case ESquare:
//                    sout << "_quad";
//                    break;
//                case ETrapezoidal:
//                    sout << "_trap";
//                    break;
//            }
//            sout << "_" << stressPOrder << "_Error.nb";
//            ofstream ErroOut(sout.str(), std::ios::app);
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
//
//
//            an.CleanUp();
//
//            delete cmesh_m_Hybrid;
//            for (int i = meshvector_Hybrid.size() - 1; i >= 0; i--) {
//                meshvector_Hybrid[i]->CleanUp();
//                delete meshvector_Hybrid[i];
//            }
//
//            delete cmesh_m_HDiv;
//            for (int i = meshvector_HDiv.size() - 1; i >= 0; i--) {
//                meshvector_HDiv[i]->CleanUp();
//                delete meshvector_HDiv[i];
//            }
//            delete gmesh;
//
//        }
//    }
//
//    //
//    //
//    //
//    //    //Pós-processamento (paraview):
//    //    std::cout << "Post Processing " << std::endl;
//    //    std::string plotfile("ElasticityTest.vtk");
//    //    TPZStack<std::string> scalnames, vecnames;
//    //    vecnames.Push("Displacement");
//    //    vecnames.Push("Stress");
//    //    vecnames.Push("Rotation");
//    ////    vecnames.Push("V_exact");
//    ////    vecnames.Push("P_exact");
//    //    //        vecnames.Push("V_exactBC");
//    //
//    //
//    //    int postProcessResolution = 3; //  keep low as possible
//    //
//    //    int dim = gmesh->Dimension();
//    //    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//    //    an.PostProcess(postProcessResolution,dim);
//
//    std::cout << "FINISHED!" << std::endl;

    return 0;
}
