//
// Created by Joel Eaves on 2019-07-10.
//

#include <cmath>
#include "fix_heun.h"

using namespace LAMMPS_NS;
using namespace FixConst;

//Use of the fix would look like
//fix   1 all heun seed temperature friction_coefficient
FixHeun::FixHeun(LAMMPS* lmp,int narg,char** arg):Fix(lmp,narg,arg),random(nullptr),temperature(nullptr),
    f_previous(nullptr),x_previous(nullptr){
    if( narg<5 ) error->all(FLERR,"Stop fucking up the Heun algorithm, idiot.");
    seed = force->inumeric(FLERR,arg[3]);
    t_target = force->numeric(FLERR, arg[4] );
    random = new RanMars(lmp, seed + comm->me);
    atom->add_callback(0);//This is required because we're passing extra array information around.
    grow_arrays(atom->nmax);
}

FixHeun::~FixHeun(){
    memory->destroy(f_previous);
    memory->destroy(x_previous);
    atom->delete_callback(id,0);
    if (random == nullptr)
        return;
    delete random;
}

int FixHeun::setmask(){
        int mask = 0;
        mask |= INITIAL_INTEGRATE;
        mask |= FINAL_INTEGRATE;
        return mask;
}

void FixHeun::setup(int){
    compute_forces();
}

void FixHeun::initial_integrate(int){
    int i,k;
    int* mask = atom->mask;
    int nlocal = atom->nlocal;
    double** x = atom->x;
    double** f = atom->f;
    double dt = update->dt;
    double gamma = sqrt(2*12*t_target*dt);//This contains a factor of 12 to give it the right std.
    double dw;
    for(i=0;i<nlocal;i++){
        if( mask[i] & groupbit ){
            for(k=0;k<3;k++) {
                dw = gamma*( random->uniform() - 0.5 );
                x_previous[i][k]=x[i][k];
                f_previous[i][k]=f[i][k] + dw/dt;
                x[i][k] += f_previous[i][k]*dt;
            }
        }
    }
}

void FixHeun::final_integrate(){
    /*Test to make sure that we can calculate the forces correctly.*/
    int i,k;
    int* mask = atom->mask;
    int nlocal = atom->nlocal;
    double** x = atom->x;
    double** f = atom->f;
    double** v = atom->v;
    double dt = update->dt;
    double gamma = sqrt(2*12*t_target*dt);
    double dw;
    for(i=0;i<nlocal;i++){
        if( mask[i] & groupbit ){
            for(k=0;k<3;k++) {
                dw = gamma*(random->uniform() - 0.5);
                x[i][k]= x_previous[i][k] + 0.5*dt*f_previous[i][k] +0.5*dt*f[i][k] + 0.5*dw;
                v[i][k] = ( x[i][k] - x_previous[i][k] )/dt;
            }
        }
    }
    compute_forces();//You need this because the positions were advanced during the previuos timestep.
}

void FixHeun::grow_arrays(int nmax)
{
    memory->grow(this->x_previous,nmax,3,"Heun: previous position.");
    memory->grow(this->f_previous,nmax,3,"Henu: previous force.");
    //std::cout << "Grew arrays no problem" << std::endl;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixHeun::copy_arrays(int i, int j, int /*delflag*/)
{
    //std::cout << "Copying arrays" << std::endl;
    memcpy(this->f_previous[j],this->f_previous[i],3*sizeof(double));
    memcpy(this->x_previous[j],this->x_previous[i],3*sizeof(double));
    //std::cout << "Copied arrays" << std::endl;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixHeun::pack_exchange(int i, double *buf)
{
    //std::cout << "Packing " << "\t" << std::endl;
    int m = 0;
    buf[m++] = x_previous[i][0];
    buf[m++] = x_previous[i][1];
    buf[m++] = x_previous[i][2];
    buf[m++] = f_previous[i][0];
    buf[m++] = f_previous[i][1];
    buf[m++] = f_previous[i][2];
    //std::cout << "Packed " << std::endl;
    return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixHeun::unpack_exchange(int nlocal, double *buf)
{
    //std::cout  << "Unpacking " << std::endl;
    int m = 0;
    x_previous[nlocal][0] = buf[m++];
    x_previous[nlocal][1] = buf[m++];
    x_previous[nlocal][2] = buf[m++];
    f_previous[nlocal][0] = buf[m++];
    f_previous[nlocal][1] = buf[m++];
    f_previous[nlocal][2] = buf[m++];
    //std::cout << "Unpacked " << std::endl;
    return m;
}

double FixHeun::memory_usage() {
    //std::cout << "Computing memory usage" << std::endl;
    int nmax = atom->nmax;
    double bytes = 0;
    bytes = 2*3*nmax*sizeof(double);//We're sending 2 arrays of size nmax x 3.
    //std::cout << "Memory size " << "\t" << bytes << std::endl;
    return bytes;
}

void FixHeun::set_arrays(int i) {
    //std::cout << "Setting" << std::endl;
    memset(this->x_previous[i],0,sizeof(double)*3);
    memset(this->f_previous[i],0,sizeof(double)*3);
    //std::cout << "Set " << std::endl;
}

void FixHeun::force_clear()//Adapted from verlet.cpp.
    {
        size_t nbytes;

        // clear force on all particles
        // if either newton flag is set, also include ghosts
        // when using threads always clear all forces.

        int nlocal = atom->nlocal;

        if (neighbor->includegroup == 0) {
            nbytes = sizeof(double) * nlocal;
            if (force->newton) nbytes += sizeof(double) * atom->nghost;

            if (nbytes) {
                memset(&atom->f[0][0],0,3*nbytes);
                //if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
                //if (extraflag) atom->avec->force_clear(0,nbytes);
            }

            // neighbor includegroup flag is set
            // clear force only on initial nfirst particles
            // if either newton flag is set, also include ghosts

        } else {
            nbytes = sizeof(double) * atom->nfirst;

            if (nbytes) {
                memset(&atom->f[0][0],0,3*nbytes);
                //if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
                //if (extraflag) atom->avec->force_clear(0,nbytes);
            }

            if (force->newton) {
                nbytes = sizeof(double) * atom->nghost;

                if (nbytes) {
                    memset(&atom->f[nlocal][0],0,3*nbytes);
                    //if (torqueflag) memset(&atom->torque[nlocal][0],0,3*nbytes);
                    //if (extraflag) atom->avec->force_clear(nlocal,nbytes);
                }
            }
        }
    }

void FixHeun::compute_forces() {
    int nflag = neighbor->decide();
    if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
    } else {
        if (modify->n_pre_exchange) {
            timer->stamp();
            modify->pre_exchange();
            timer->stamp(Timer::MODIFY);
        }
        //if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
            domain->reset_box();
            comm->setup();
            if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 &&
            update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        //if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (modify->n_pre_neighbor) {
            modify->min_pre_neighbor();
            timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_post_neighbor) {
            modify->min_post_neighbor();
            timer->stamp(Timer::MODIFY);
        }
    }
    int eflag = 0;//3;
    int vflag = 0;//5;
    //ev_set(update->ntimestep);
    force_clear();

    timer->stamp();

    if (modify->n_pre_force) {
        modify->pre_force(vflag);
        timer->stamp(Timer::MODIFY);
    }

    //if (pair_compute_flag) {
        force->pair->compute(eflag,vflag);
        timer->stamp(Timer::PAIR);
    //}

    if (atom->molecular) {
        /*
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
         */
        timer->stamp(Timer::BOND);
    }

    //if (kspace_compute_flag) {
    //force->kspace->compute(eflag,vflag);
    //    timer->stamp(Timer::KSPACE);
    //}

    if (modify->n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
    }

    if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
    }
} // Adapted from min.cpp


