//
// Created by Joel Eaves on 2019-07-20.
//


#ifdef COMPUTE_CLASS

ComputeStyle(Fskt,ComputeFskt)

#else

#ifndef LAMMPS_COMPUTE_FSKT_H
#define LAMMPS_COMPUTE_FSKT_H

#include "compute.h"
#include <cmath>

namespace LAMMPS_NS {

    class ComputeFskt : public compute {
    public:
        ComputeFskt(class LAMMPS *, int, char **);
        ~ComputeFskt();
        void init();
        virtual void compute_vector();
        void set_arrays(int);
    protected:
        bigint nfskt;
        char *id_fix;
        class FixStore *fix;
    };

}

#endif //LAMMPS_COMPUTE_FSKT_H
#endif