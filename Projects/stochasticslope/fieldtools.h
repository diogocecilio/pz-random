#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"


class FieldTools
{
private:

    TPZGeoMesh* fgmesh;
    TPZCompMesh* fcmesh;
    REAL flx;//total slope heigth
    REAL fly;//water heigth
    REAL fM;//total heigth
    REAL ftype;
    int fporder;
    string ffile;

public:
    FieldTools();
    FieldTools(FieldTools & copy);
    FieldTools(TPZGeoMesh* gmesh,REAL lx, REAL ly,int M,int type,int porder,string file);
    ~FieldTools();

    TPZCompMesh* GetCMesh()
    {
        return fcmesh;
    }

    void ComputeField (REAL mean,REAL cov,int samples );
    TPZCompMesh* CreateCompMeshKL ( );
    TPZFMatrix<REAL> CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata );
    TPZCompMesh * SettingCreateFild (  );
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
FieldTools::FieldTools(TPZGeoMesh*gmesh,REAL lx, REAL ly,int M,int type,int porder,string file)
{
    ffile=file;
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

TPZCompMesh* FieldTools::CreateCompMeshKL ( )
{

    int id=1;
    int dim = fgmesh->Dimension();
    TPZCompMesh * cmesh1 = new TPZCompMesh ( fgmesh );

    REAL lz=1.;

    KLMaterial * klmat1 = new KLMaterial ( id,flx,fly,lz,dim,ftype,fM );
    cmesh1->SetDefaultOrder ( fporder );
    cmesh1->SetDimModel ( dim );
    cmesh1->InsertMaterialObject ( klmat1 );
    cmesh1->SetAllCreateFunctionsContinuous();
    cmesh1->AutoBuild();

    return cmesh1;
}

void FieldTools::ComputeField ( REAL mean,REAL cov,int samples)
{

    KLAnalysis * klanal = new KLAnalysis ( fcmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> ( fcmesh->MaterialVec() [1] );

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();

    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();


    TPZFMatrix<REAL> field = CreateLogNormalRandomField ( eigenfunctions, mean, cov, samples, ffile );


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


TPZCompMesh * FieldTools::SettingCreateFild ( )
{
    cout << "outco = ?";
    TPZCompMesh * cmesh = CreateCompMeshKL ( );
    string outco="/home/diogo/projects/pz-random-build-release/Projects/stochasticslope/";
    outco+=ffile;
    cout << "outco = " << outco;
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    cmesh->LoadSolution ( readco );

    TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis ( cmesh );
    analysis1x->LoadSolution ( readco );

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

