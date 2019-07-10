//
// Created by Joel Eaves on 2019-07-10.
//

#include <cstdlib>
#include "fix_example.h"

using namespace LAMMPS_NS;
using namespace FixConst;


Fix_Example::Fix_Example(LAMMPS *lmp, int narg, char **arg) :Fix(lmp,narg,arg){
    if(narg < 4)
        error -> all(FLERR,"Illegal fix print command.");
    nevery = atoi(arg[3]);
    if (nevery <= 0 )
        error->all(FLERR,"Illegal fix print command, you fuckup!");
}

int Fix_Example::setmask() {
    int mask = 0;
    mask |= FixConst::END_OF_STEP;
    return mask;
}

void Fix_Example::end_of_step() {
    int i;
    double** v = atom->v;
    int nlocal = atom->nlocal;
    std::unique_ptr < double [] > localAverageVelocity (new double [4]);
    std::unique_ptr < double [] > globalAverageVelocity (new double [4]);
    for(i=0;i<4;i++){
        localAverageVelocity[i] = 0;
        globalAverageVelocity[i] = 0;
    }
    for(i=0;i<nlocal;i++)
        MathExtra::add3( localAverageVelocity.get(),v[i],localAverageVelocity.get() );
    localAverageVelocity[3] = nlocal;
    /*
    double localAverageVelocity[4];
    memset(localAverageVelocity,0,4*sizeof(double));
    for(int i=0;i<nlocal;i++)
        MathExtra::add3(localAverageVelocity,v[i],localAverageVelocity);
    localAverageVelocity[3] = nlocal;
    double globalAverageVelocity[4];
    memset(globalAverageVelocity,0,4*sizeof(double));
    */
    MPI_Allreduce(localAverageVelocity.get(), globalAverageVelocity.get(), 4, MPI_DOUBLE, MPI_SUM, world);
    MathExtra::scale3( 1.0/globalAverageVelocity[3], globalAverageVelocity.get() );
    if(comm->me == 0) {
        std::cout << globalAverageVelocity[0] << ", " << globalAverageVelocity[1] << ", " << globalAverageVelocity[2]
                  << std::endl;
    }
}