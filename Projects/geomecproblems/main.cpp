#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"
#include "pzgeoelbc.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"


#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplasticmem.h"
#include "pzelastoplastic2D.h"
#include "TPZPlasticStepPV.h"
#include "pzelastoplasticanalysis.h"
#include "readgidmesh.h"
#include "slope.h"
#include "slopecphi.h"
#include "footing.h"

int main() {

   SlopeCphi slopecphi;

   slopecphi.run();

 //    Slope slope;
//
//     slope.run();

//   Footing foot;

//  foot.runfoot();


    return 0;
}
