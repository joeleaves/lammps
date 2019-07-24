//
// Created by Joel Eaves on 2019-07-20.
//

#include "compute_fskt.h"

using namespace std;
using namespace LAMMPS_NS

ComputeFskt::ComputeFskt(LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg),
        id_fix(NULL)
{
    if (narg < 3) error->all(FLERR,"Illegal compute vacf command");

    vector_flag = 1;
    size_vector = 4;
    extvector = 0;
    create_attribute = 1;

    // create a new fix STORE style
    // id = compute-ID + COMPUTE_STORE, fix group = compute group

    int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
    id_fix = new char[n];
    strcpy(id_fix,id);
    strcat(id_fix,"_COMPUTE_STORE");

    char **newarg = new char*[6];
    newarg[0] = id_fix;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "1";
    newarg[5] = (char *) "3";
    modify->add_fix(6,newarg);
    fix = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;

    // store current velocities in fix store array
    // skip if reset from restart file

    if (fix->restart_reset) fix->restart_reset = 0;
    else {
        double **roriginal = fix->astore;

        double **r = atom->r;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) {
                roriginal[i][0] = r[i][0];
                roriginal[i][1] = r[i][1];
                roriginal[i][2] = r[i][2];
            } else roriginal[i][0] = roriginal[i][1] = roriginal[i][2] = 0.0;
    }

    // displacement vector

    vector = new double[size_vector];
}
