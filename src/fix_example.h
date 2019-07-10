//
// Created by Joel Eaves on 2019-07-10.
//
#ifdef FIX_CLASS

FixStyle(example,Fix_Example)

#else

#ifndef LAMMPS_FIX_EXAMPLE_H
#define LAMMPS_FIX_EXAMPLE_H

#include "fix.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include <iostream>
#include <memory>

namespace LAMMPS_NS {

    class Fix_Example : public Fix {
    public:
        Fix_Example(LAMMPS *lmp, int narg, char **arg);
        ~Fix_Example() = default;
        int setmask() override;
        void end_of_step() override;
    };

}

#endif
#endif //LAMMPS_FIX_EXAMPLE_H
