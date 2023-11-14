#include "./pair_test.h"

#include "lammps.h"
#include "comm.h"
#include "utils.h"
#include "error.h"
#include "memory.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"


namespace LAMMPS_NS {

PairTestmlff::PairTestmlff(LAMMPS* lmp) : Pair(lmp) {
    this->me = this->comm->me;
    this->writedata = 1;
    this->int_tensor_options = c10::TensorOptions().dtype(torch::kInt64);
    this->float_tensor_options = c10::TensorOptions().dtype(torch::kFloat64);
    
    // Step . Log
    utils::logmesg(this->lmp, "*** PairTestmlff::PairTestmlff()\n");
}

PairTestmlff::~PairTestmlff() {
    this->memory->destroy(this->setflag);
    this->memory->destroy(this->cutsq);
}

void PairTestmlff::settings(int argc, char** argv) {
    if (argc != 3) 
        this->error->all(FLERR, "Illegal command pair_style!!!\n");
    
    this->pt_file = argv[0];
    this->cutgb = utils::numeric(FLERR, argv[1], false, this->lmp);
    this->cutgbs = utils::numeric(FLERR, argv[2], false, this->lmp);

    // Step . Log
    utils::logmesg(this->lmp, "*** PairTestmlff::settings()\n");
    std::cout << "\t + 1. this->pt_file = " << this->pt_file << std::endl;
    std::cout << "\t + 2. this->cutgb = " << this->cutgb << std::endl;
    std::cout << "\t + 3. this->cutgbs = " << this->cutgbs << std::endl;
}


void PairTestmlff::allocate() {
    this->memory->create<int>(this->setflag, this->atom->ntypes+1, this->atom->ntypes+1, "PairTestmlff::setflag");
    this->memory->create<double>(this->cutsq, this->atom->ntypes+1, this->atom->ntypes+1, "PairTestmlff::cutsq");

    // Step . Log
    utils::logmesg(this->lmp, "*** PairTestmlff::allocate()\n");
}

void PairTestmlff::coeff(int argc, char** argv) {
    if (!this->allocated)
        this->allocate();

    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, argv[0], 1, this->atom->ntypes, ilo, ihi, this->error);
    utils::bounds(FLERR, argv[1], 1, this->atom->ntypes, jlo, jhi, this->error);
    
    // Step 2. Assign `this->setflag`
    for (int ii=1; ii<this->atom->ntypes+1; ii++) {
        for (int jj=1; jj<this->atom->ntypes+1; jj++){
            this->setflag[ii][jj] = 1;
            this->cutsq[ii][jj] = this->cutgb * this->cutgb;
        }
    }

    // Step 3. 
    utils::logmesg(this->lmp, "*** PairTestmlff::coeff()\n");
    printf("\t + 1. ilo = %3d, ihi = %3d\n", ilo, ihi);
    printf("\t + 2. jlo = %3d, jhi = %3d\n", jlo, jhi);
}


void PairTestmlff::init_style() {
    this->neighbor->add_request(this, NeighConst::REQ_FULL);
    
    // Step . 
}


double PairTestmlff::init_one(int i, int j) {

}


void PairTestmlff::compute(int eflag, int vflag) {

}


}   // namespace : LAMMPS_NS