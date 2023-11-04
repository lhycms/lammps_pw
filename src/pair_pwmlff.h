#ifdef PAIR_CLASS
PairStyle(pwmlff, PairPwmlff);
#else
#ifndef LMP_PAIR_PWMLFF_H
#include "lammps.h"
#include "pair.h"
#include <string>
#include <iostream>
#include <stdio.h>
#include <torch/torch.h>
#include <torch/script.h>

namespace LAMMPS_NS {
class PairPwmlff : public Pair {
public:
    PairPwmlff(LAMMPS* lmp);

    void settings(int argc, char** argv) override;

    void coeff(int argc, char** argv) override;

    void compute(int eflag, int vflag) override;

    double init_one(int i, int j);

    void init_style();

protected:
    void allocate();

private:
    int me;
    std::string pt_file;
    double cutgb;
    double cutgbs;
    int max_num_neigh_atoms;

    torch::jit::script::Module module;

    at::Tensor davg;
    at::Tensor dstd;
    at::Tensor natoms_image_tensor; // [[12, 4, 8]]
    at::Tensor atom_types_tensor;
};  // class : PwmlffPair
}   // namespace : LAMMPS_NS
#endif
#endif