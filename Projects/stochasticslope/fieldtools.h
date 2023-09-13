#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"

string filelocation2 = "/home/diogo/projects/pz-random-build-release/Projects/stochasticslope/";

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
    FieldTools(TPZGeoMesh* gmesh,REAL lx, REAL ly,int M,int type,int porder);
    ~FieldTools();

    TPZCompMesh* GetCMesh()
    {
        return fcmesh;
    }

    void ComputeField (REAL mean,REAL cov,int samples );

    void ComputeFields (REAL mean1,REAL mean2,REAL cov1, REAL cov2,int samples,REAL crossfac,string file1,string file2);

    void ComputeFields (TPZVec<REAL> mean,TPZVec<REAL> cov,TPZVec<string> file,int samples);

    TPZCompMesh * SettingCreateFild (string file );

    TPZCompMesh* CreateCompMeshKL ( );
    TPZFMatrix<REAL> CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata );
    TPZVec<TPZFMatrix<REAL>>  GenerateNonGaussinRandomField (TPZFMatrix<REAL> PHI, REAL meanc,REAL meanphi, REAL covc,REAL covphi,int samples,REAL crossfac);

    void  GenerateNonGaussinRandomField (TPZFMatrix<REAL> PHI, TPZVec<REAL> mean,TPZVec<REAL> cov,TPZVec<string> file,int samples);

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

FieldTools::FieldTools(TPZGeoMesh*gmesh,REAL lx, REAL ly,int M,int type,int porder)
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
    //delete fcmesh;
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

void FieldTools::ComputeFields ( REAL mean1,REAL mean2,REAL cov1, REAL cov2,int samples,REAL crossfac,string file1,string file2)
{

    KLAnalysis * klanal = new KLAnalysis ( fcmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> ( fcmesh->MaterialVec() [1] );

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();

    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();


    TPZVec<TPZFMatrix<REAL>>  fields = GenerateNonGaussinRandomField (eigenfunctions, mean1, mean2,  cov1, cov2, samples, crossfac);

    PrintMat ( file1,fields[0] );

    PrintMat ( file2,fields[1] );
}


void FieldTools::ComputeFields (TPZVec<REAL> mean,TPZVec<REAL> cov,TPZVec<string> file,int samples)
{
    KLAnalysis * klanal = new KLAnalysis ( fcmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> ( fcmesh->MaterialVec() [1] );

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();

    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();


    GenerateNonGaussinRandomField (eigenfunctions,   mean, cov,  file,  samples);


}

TPZFMatrix<REAL> FieldTools::CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata )
{
    TPZFMatrix<REAL>  PHIt;

   //PHI.Transpose ( &PHIt );

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


void  FieldTools::GenerateNonGaussinRandomField (TPZFMatrix<REAL> PHI, TPZVec<REAL> mean,TPZVec<REAL> cov,TPZVec<string> file,int samples)
{
      TPZFMatrix<REAL>  PHIt;

    int nsamples = samples;

    int M = PHI.Cols();

    int nfields=mean.size();


    for(int ifield=0;ifield<nfields;ifield++)
    {
        std::normal_distribution<REAL> distribution ( 0., 1. );
        TPZFMatrix<REAL> THETA (M,samples );

        for ( int n = 0; n < samples; n++ ) {
            for ( int iexp = 0; iexp < M; iexp++ ) {
                std::random_device rd{};
                std::mt19937 generator{ rd() };
                REAL xic = distribution ( generator );
                THETA(iexp,n) = xic;
            }
    }

        TPZFMatrix<REAL> hhat;
        PHI.Multiply ( THETA, hhat );


        REAL sdev = cov[ifield] * mean[ifield];
        REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean[ifield] ),2 ) ) );
        REAL lambda = log ( mean[ifield] ) - xi * xi / 2.;

        for ( int i = 0; i < hhat.Rows(); i++ ) {
            for ( int j = 0; j < hhat.Cols(); j++ ) {
                hhat(i,j) = exp ( lambda + xi * hhat(i,j) );
            }

        }

        PrintMat ( file[ifield],hhat );

    }

    std::cout << "\n Exiting  generalized eigenvalue prolem" << endl;

}

TPZVec<TPZFMatrix<REAL>>  FieldTools::GenerateNonGaussinRandomField (TPZFMatrix<REAL> PHI, REAL meanc,REAL meanphi, REAL covc,REAL covphi,int samples,REAL crossfac)
{
    TPZFMatrix<REAL>  PHIt;

    //PHI.Transpose ( &PHIt );

    int M = PHI.Cols();


    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine generator(seed);



    std::normal_distribution<REAL> distribution ( 0., 1. );

    int nsamples = samples;
    TPZFMatrix<REAL> THETA (M,samples ), THETA2( M, samples);
    REAL correlation = crossfac;
    for ( int n = 0; n < samples; n++ ) {
        for ( int iexp = 0; iexp < M; iexp++ ) {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA(iexp,n) = xic;
            THETA2(iexp,n) = xic * correlation + xiphi * sqrt ( 1 - correlation * correlation );
        }
    }

    TPZFMatrix<REAL> hhatphi, hhatcoes;

    PHI.Multiply ( THETA, hhatcoes );
    PHI.Multiply ( THETA2, hhatphi );


    REAL sdev = covc * meanc;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / meanc ),2 ) ) );
    REAL lambda = log ( meanc ) - xi * xi / 2.;

    for ( int i = 0; i < hhatcoes.Rows(); i++ ) {
        for ( int j = 0; j < hhatcoes.Cols(); j++ ) {
                hhatcoes(i,j) = exp ( lambda + xi * hhatcoes(i,j) );
        }

    }

    sdev = covphi * meanphi;
    xi = sqrt ( log ( 1 + pow ( ( sdev / meanphi ), 2 ) ) );
    lambda = log ( meanphi ) - xi * xi / 2.;
    for ( int i = 0; i < hhatphi.Rows(); i++ ) {
        for ( int j = 0; j < hhatphi.Cols(); j++ ) {
                hhatphi(i,j) = exp ( lambda + xi * hhatphi(i,j) );
        }

    }

    TPZVec<TPZFMatrix<REAL>> fields(2);
    fields[0] = hhatcoes;
    fields[1] = hhatphi;



    std::cout << "\n Exiting  generalized eigenvalue prolem" << endl;

    return fields;
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
    string outco=filelocation2;
    outco+=ffile;
    cout << "outco = " << outco;
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    cmesh->LoadSolution ( readco );

    TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis ( cmesh );
    analysis1x->LoadSolution ( readco );

    delete analysis1x;
    return fcmesh;
}

TPZCompMesh * FieldTools::SettingCreateFild (string file )
{
    TPZCompMesh * cmesh = CreateCompMeshKL ( );
    string outco=filelocation2;
    outco+=file;
    cout << "outco = " << outco;
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    cmesh->LoadSolution ( readco );

    TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis ( cmesh );
    analysis1x->LoadSolution ( readco );

    delete analysis1x;
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

