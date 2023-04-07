#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"


class FieldTools
{
private:

    TPZManVector<TPZGeoMesh*,3> fgmesh;
    TPZManVector<TPZCompMesh*,3> fcmesh;
    REAL flx;//total slope heigth
    REAL fly;//water heigth
    REAL fM;//total heigth
    REAL ftype;
    int fporder;

public:
    FieldTools();
    FieldTools(FieldTools & copy);
    FieldTools(TPZManVector<TPZGeoMesh*,3> gmesh,REAL lx, REAL ly,int M,int type,int porder);
    ~FieldTools();


    void ComputeField (  );
    TPZManVector<TPZCompMesh*,3>  CreateCompMeshKL ( );
    TPZFMatrix<REAL> CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata );
    TPZManVector<TPZCompMesh *,3> SettingCreateFilds (  );
    void PrintMat ( std::string out,TPZFMatrix<REAL> mat );
    void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );
    template <class T>
    std::vector<T> str_vec ( std::vector<std::string> &vs );


};

FieldTools::FieldTools()
{

}

FieldTools::FieldTools(FieldTools & copy): fgmesh(copy.fgmesh),fcmesh(copy.fcmesh),flx(copy.flx),
fly(copy.fly),fM(copy.fM),ftype(copy.ftype),fporder(copy.fporder)
{

}
FieldTools::FieldTools(TPZManVector<TPZGeoMesh*,3> gmesh,REAL lx, REAL ly,int M,int type,int porder)
{
    fgmesh =gmesh;
    flx=lx;
    fly=ly;
    fM=M;
    ftype=type;
    fporder=porder;
    fcmesh = CreateCompMeshKL ( );
}

FieldTools::~FieldTools()
{

}

TPZManVector<TPZCompMesh*,3>  FieldTools::CreateCompMeshKL ( )
{

    TPZManVector<TPZCompMesh*,3>  veccmesh(3);
    int id=1;
    int dim = fgmesh[0]->Dimension();
    TPZCompMesh * cmesh1 = new TPZCompMesh ( fgmesh[0] );
    TPZCompMesh * cmesh2 = new TPZCompMesh ( fgmesh[1] );
    TPZCompMesh * cmesh3 = new TPZCompMesh ( fgmesh[2] );

    REAL lz=1.;

    KLMaterial * klmat1 = new KLMaterial ( id,20.,2.,lz,dim,ftype,fM );
    KLMaterial * klmat2 = new KLMaterial ( id,20.,2.,lz,dim,ftype,fM );
    KLMaterial * klmat3 = new KLMaterial ( id,1.,1.,lz,dim,ftype,fM );

    cmesh1->SetDefaultOrder ( fporder );
    cmesh2->SetDefaultOrder ( fporder );
    cmesh3->SetDefaultOrder ( fporder );
    cmesh1->SetDimModel ( dim );
    cmesh2->SetDimModel ( dim );
    cmesh3->SetDimModel ( dim );
    cmesh1->InsertMaterialObject ( klmat1 );
    cmesh2->InsertMaterialObject ( klmat2 );
    cmesh3->InsertMaterialObject ( klmat3 );
    cmesh1->SetAllCreateFunctionsContinuous();
    cmesh2->SetAllCreateFunctionsContinuous();
    cmesh3->SetAllCreateFunctionsContinuous();
    cmesh1->AutoBuild();
    cmesh2->AutoBuild();
    cmesh3->AutoBuild();
    veccmesh[0]=cmesh1;
    veccmesh[1]=cmesh2;
    veccmesh[2]=cmesh3;
    return veccmesh;
}

void FieldTools::ComputeField (  )
{


    KLAnalysis * klanal = new KLAnalysis ( fcmesh[0] );

    KLMaterial *mat = dynamic_cast<KLMaterial*> ( fcmesh[0]->MaterialVec() [1] );

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();



    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();
    REAL mean=10.;
    REAL cov=0.3;
    int samples=1000;
    string outdata = "coesao.txt";

    TPZFMatrix<REAL> field = CreateLogNormalRandomField ( eigenfunctions, mean, cov, samples, outdata );


    mean=30.*M_PI/180.;
    cov=0.2;
    string outdataatrito = "atrito.txt";
    TPZFMatrix<REAL> fieldphi = CreateLogNormalRandomField ( eigenfunctions, mean, cov, samples, outdataatrito );


    mean=1.;
    cov=0.2;
    string outdataperm = "permeability.txt";
    TPZFMatrix<REAL> fieldperm = CreateLogNormalRandomField ( eigenfunctions, mean, cov, samples, outdataperm );


//     string file1 ="coesao3.vtk";
//
//     klanal->Post ( file1,2,0 );

}

TPZFMatrix<REAL> FieldTools::CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata )
{
    TPZFMatrix<REAL>  PHIt;

    PHI.Transpose ( &PHIt );

    int M = PHI.Cols();

    cout << "M = " << M <<endl;

    std::normal_distribution<double> distribution ( 0., 1. );

    TPZFMatrix<REAL>  THETA ( M, samples, 0. );
    for ( int isample = 0; isample < samples; isample++ )
    {
        for ( int irdvar = 0; irdvar < M; irdvar++ )
        {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA ( irdvar,isample ) = xic;
        }
    }


    TPZFMatrix<REAL>  hhat;

    PHI.Multiply ( THETA, hhat );

    REAL sdev = cov * mean;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
    REAL lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhat.Rows(); i++ )
    {
        for ( int j = 0; j < hhat.Cols(); j++ )
        {
            REAL temp =  hhat ( i,j );
            hhat ( i,j ) = exp ( lambda + xi * temp );
        }
    }

    //string out="testefield.txt";

    PrintMat ( outdata,hhat );

    return hhat;
}

void FieldTools::PrintMat ( std::string out,TPZFMatrix<REAL> mat )
{
    std::ofstream print ( out );
    int row=mat.Rows();
    int cols = mat.Cols();

    for ( int i=0; i<row; i++ )
    {
        for ( int j=0; j<cols; j++ )
        {
            print << mat ( i,j ) << " ";
        }
        print<< std::endl;
    }


}


TPZManVector<TPZCompMesh *,3> FieldTools::SettingCreateFilds ( )
{
    string outco="/home/diogo/projects/pz-random-build-release/Projects/randomdarcy/coesao.txt";
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    fcmesh[0]->LoadSolution ( readco );

    TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis ( fcmesh[0] );
    analysis1x->LoadSolution ( readco );


    string outphi="/home/diogo/projects/pz-random-build-release/Projects/randomdarcy/atrito.txt";
    TPZFMatrix<REAL> readphi;
    ReadFile ( outphi,readphi );
    fcmesh[1]->LoadSolution ( readphi );

    TPZElastoPlasticAnalysis * analysis2x  = new TPZElastoPlasticAnalysis ( fcmesh[1] );
    analysis2x->LoadSolution ( readphi );

    string outperm="/home/diogo/projects/pz-random-build-release/Projects/randomdarcy/permeability.txt";
    TPZFMatrix<REAL> readperm;
    ReadFile ( outperm,readperm );
    fcmesh[2]->LoadSolution ( readperm );

    TPZElastoPlasticAnalysis * analysis3x  = new TPZElastoPlasticAnalysis ( fcmesh[2] );
    analysis3x->LoadSolution ( readperm );


    return fcmesh;
}

void  FieldTools::ReadFile ( std::string file,TPZFMatrix<REAL> &out )
{
    string line,line2, temp;

    ifstream myfile ( file );

    std::vector<std::vector<double>> coords;
    int counter=0;
    while ( getline ( myfile, line ) )
    {
        std::vector<string> tokens;
        istringstream iss ( line );
        while ( iss >> temp ) tokens.push_back ( temp );
        std::vector<double> input_doub_temp= str_vec<double> ( tokens );
        coords.push_back ( input_doub_temp );
        counter++;
    }

    out.Resize ( coords.size(),coords[1].size() );
    for ( int iel=0; iel<coords.size(); iel++ )
    {
        for ( int elnode=0; elnode<coords[iel].size(); elnode++ )
        {
            out ( iel,elnode ) = coords[iel][elnode];
        }
    }
}

template <class T>
std::vector<T> FieldTools::str_vec ( std::vector<std::string> &vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it )
    {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

