#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"


class FieldTools
{
private:

    TPZGeoMesh * fgmesh1;
    TPZGeoMesh * fgmesh2;
    TPZCompMesh * fcmesh;
    REAL flx;//total slope heigth
    REAL fly;//water heigth
    REAL fM;//total heigth
    REAL ftype;
    int fporder;

public:
    FieldTools();
    FieldTools(FieldTools & copy);
    FieldTools(TPZGeoMesh * gmesh1,TPZGeoMesh * gmesh2,REAL lx, REAL ly,int M,int type,int porder);
    ~FieldTools();


    TPZCompMesh * ComputeField (  );
    TPZCompMesh * CreateCompMeshKL ( TPZGeoMesh* gmesh );
    TPZFMatrix<REAL> CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata );
    TPZManVector<TPZCompMesh *,2> SettingCreateFilds (  );
    void PrintMat ( std::string out,TPZFMatrix<REAL> mat );
    void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );
    template <class T>
    std::vector<T> str_vec ( std::vector<std::string> &vs );


};

FieldTools::FieldTools()
{

}

FieldTools::FieldTools(FieldTools & copy): fgmesh1(copy.fgmesh1),fgmesh2(copy.fgmesh2),fcmesh(copy.fcmesh),flx(copy.flx),
fly(copy.fly),fM(copy.fM),ftype(copy.ftype),fporder(copy.fporder)
{

}
FieldTools::FieldTools(TPZGeoMesh * gmesh1,TPZGeoMesh * gmesh2,REAL lx, REAL ly,int M,int type,int porder)
{
    fgmesh1=gmesh1;
    fgmesh2=gmesh2;
    flx=lx;
    fly=ly;
    fM=M;
    ftype=type;
    fporder=porder;
}

FieldTools::~FieldTools()
{

}

TPZCompMesh * FieldTools::CreateCompMeshKL ( TPZGeoMesh* gmesh)
{

    int id=1;
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    REAL lz=1.;

    KLMaterial * klmat = new KLMaterial ( id,flx,fly,lz,dim,ftype,fM );

    cmesh->SetDefaultOrder ( fporder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    return cmesh;


}

TPZCompMesh * FieldTools::ComputeField (  )
{

    TPZCompMesh * cmesh = CreateCompMeshKL (fgmesh1);

    KLAnalysis * klanal = new KLAnalysis ( cmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> ( cmesh->MaterialVec() [1] );

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

    cmesh->LoadSolution ( field );

    string file1 ="coesao3.vtk";

    klanal->Post ( file1,2,0 );

    return cmesh;
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

TPZManVector<TPZCompMesh *,2> FieldTools::SettingCreateFilds ( )
{

    TPZManVector<TPZCompMesh *,2>  vecmesh ( 2 );
    TPZCompMesh * cmeshfield;
    TPZCompMesh*  cmeshfield2;
    cmeshfield =  CreateCompMeshKL ( fgmesh1);
    cmeshfield2 =  CreateCompMeshKL (fgmesh2 );
    string outco="/home/diogo/projects/pz-random-build-release/Projects/randomdarcy/coesao.txt";
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    cmeshfield->LoadSolution ( readco );

    TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis ( cmeshfield );
    analysis1x->LoadSolution ( readco );


    string outphi="/home/diogo/projects/pz-random-build-release/Projects/randomdarcy/atrito.txt";
    TPZFMatrix<REAL> readphi;
    ReadFile ( outphi,readphi );
    cmeshfield2->LoadSolution ( readphi );

    TPZElastoPlasticAnalysis * analysis2x  = new TPZElastoPlasticAnalysis ( cmeshfield2 );
    analysis2x->LoadSolution ( readphi );

    vecmesh[0]=cmeshfield;
    vecmesh[1]=cmeshfield2;

    return vecmesh;
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

