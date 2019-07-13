//
// Created by Joel Eaves on 2019-07-10.
//

#ifdef FIX_CLASS

FixStyle(heun,FixHeun)

#else

#ifndef LAMMPS_FIX_HEUN_H
#define LAMMPS_FIX_HEUN_H
#include "fix.h"
#include "force.h"
#include "random_mars.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "integrate.h"
#include "pair.h"
#include "neighbor.h"
#include "modify.h"
#include "timer.h"
#include "domain.h"
#include "kspace.h"
#include <memory>
#include <iostream>

namespace LAMMPS_NS {

    class FixHeun :public Fix{
    public:
        FixHeun(class LAMMPS*, int, char** arg);
        FixHeun();
        ~FixHeun();
        int setmask();
        void setup(int);
        virtual void initial_integrate(int);
        virtual void final_integrate();
        void grow_arrays(int);
        void copy_arrays(int, int, int);
        int pack_exchange(int, double *);
        int unpack_exchange(int, double *);
        void set_arrays(int i);
        double memory_usage();
        void force_clear();
        void compute_forces();
    private:
        class RanMars* random;
        int seed;
        double** f_previous;
        double** x_previous;
        class Compute* temperature;
        double t_target;//Target temperature
    };

}

#endif //LAMMPS_FIX_HEUN_H
#endif