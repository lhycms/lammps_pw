#ifdef PAIR_CLASS
PairStyle(testmlff, PairTestmlff)
#else
#ifndef LMP_PAIR_TEST_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <torch/torch.h>
#include <torch/script.h>
#include "lammps.h"
#include "pair.h"

namespace LAMMPS_NS {
class PairTestmlff : public Pair {
public:
    PairTestmlff(LAMMPS* lmp);
    
    ~PairTestmlff();

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

    c10::TensorOptions int_tensor_options;
    c10::TensorOptions float_tensor_options;
};
}   // namespace : LAMMPS_NS

#endif
#endif