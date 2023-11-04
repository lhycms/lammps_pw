#include "./pair_pwmlff.h"

#include "lammps.h"
#include "pair.h"
#include "comm.h"
#include "utils.h"
#include "error.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"

#include <torch/torch.h>
#include <torch/jit.h>
#include <torch/script.h>
#include <stdlib.h>
#include <cmath>
#include "matersdk.h"


namespace LAMMPS_NS {

PairPwmlff::PairPwmlff(LAMMPS* lmp) : Pair(lmp) {
    this->me = this->comm->me;
    this->writedata = 1;

    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::Pwmlff() ***\n");
}


void PairPwmlff::settings(int argc, char** argv) {
    if (argc <= 3)
        this->error->all(FLERR, "Illegal command pair_style!!!\n");
    // Step 1. 读取 torch script module 文件路径
    this->pt_file = argv[0];
    this->cutgbs = utils::numeric(FLERR, argv[1], false, this->lmp);
    this->cutgb = utils::numeric(FLERR, argv[2], false, this->lmp);
    c10::TensorOptions int_tensor_options = c10::TensorOptions().dtype(torch::kInt32).device(c10::kCPU);
    this->atom_types_tensor = at::zeros({1, argc-3}, int_tensor_options);
    for (int ii=0; ii<argc-3; ii++) {
        int tmp_value = utils::inumeric(FLERR, argv[ii+3], false, this->lmp);
        this->atom_types_tensor[0][ii] = tmp_value;
    }

    // Step 2. 加载 torch script module 文件
    try {
        this->module = torch::jit::load(this->pt_file, c10::kCPU);
    } catch (c10::Error& e) {
        this->error->all(FLERR, "Error when loading module!!!\n");
    }

    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::settings()\n");
    std::cout << "\t + pt_file : " << this->pt_file << std::endl;
    std::cout << "\t + ntypes : " << this->module.attr("ntypes") << std::endl;
    std::cout << "\t + atom_type : " << this->module.attr("atom_type") << std::endl;
    std::cout << "\t + maxNeighborNum : " << this->module.attr("maxNeighborNum") << std::endl;
    std::cout << "\t + this->cutgbs : " << this->cutgbs << std::endl;
    std::cout << "\t + this->cutgb : " << this->cutgb << std::endl;

    this->max_num_neigh_atoms = module.attr("maxNeighborNum").toInt();
    //this->max_num_neigh_atoms = 20;
    this->davg = torch::index_select(this->module.attr("davg").toTensor(), 1, torch::arange(0, 4));
    this->dstd = torch::index_select(this->module.attr("dstd").toTensor(), 1, torch::arange(0, 4));
    auto atom_type_list = this->module.attr("atom_type").toList();


}


void PairPwmlff::allocate() {
    this->allocated = 1;

    // Step 1. Allocate memory for some `variables`
    this->memory->create<int>(this->setflag, this->atom->ntypes+1, this->atom->ntypes+1, "PairPwmlff:setflag");
    this->memory->create<double>(this->cutsq, this->atom->ntypes+1, this->atom->ntypes+1, "PairPwmlff:cutsq");
    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::allocate() ***\n");
}


void PairPwmlff::coeff(int argc, char** argv) {
    //if (argc != 2)
    //    this->error->all(FLERR, "Illegal pair_style command!!!\n")
    if (!this->allocated)
        this->allocate();
    

    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, argv[0], 1, this->atom->ntypes, ilo, ihi, this->error);
    utils::bounds(FLERR, argv[1], 1, this->atom->ntypes, jlo, jhi, this->error);

    // Step 2. Assign `this->setflag`
    for (int ii=1; ii<this->atom->ntypes+1; ii++) {
        for (int jj=MAX(ii, jlo); jj<this->atom->ntypes+1; jj++)
            this->setflag[ii][jj] = 1;
    }


    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::coeff() ***\n");
    printf("\t + ilo = %3d, ihi = %3d\n", ilo, ihi);
    printf("\t + jlo = %3d, jhi = %3d\n", jlo, jhi);
}


void PairPwmlff::init_style() {
    this->neighbor->add_request(this, NeighConst::REQ_FULL);

    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::init_style() ***\n");
}


double PairPwmlff::init_one(int i, int j) {
    if (this->setflag[i][j] == 0)
        this->error->all(FLERR, "Coefficients of pair i,j has not set!!!\n");
    
    // Step 1. Copy (i,j) to its mirror (j,i)
    this->setflag[j][i] = this->setflag[i][j];

    // Step . Log
    utils::logmesg(this->lmp, "*** PairPwmlff::init_one() ***\n");
    printf("\t + pair = [%3d, %3d]\n", i, j);
    return this->cutgb;
}


void PairPwmlff::compute(int eflag, int vflag) {
    if (eflag || vflag) 
        this->ev_setup(eflag, vflag);
    
    int local_inum = this->list->inum;
    int* local_ilist = (int*)malloc(sizeof(int) * local_inum);
    int* local_numneigh = (int*)malloc(sizeof(int) * local_inum);
    memset(local_numneigh, 0, sizeof(int) * local_inum);
    int* local_firstneigh = (int*)malloc(sizeof(int) * local_inum * this->max_num_neigh_atoms);
    memset(local_firstneigh, -1, sizeof(int) * local_inum * this->max_num_neigh_atoms);

    int tmp_center_idx;
    int tmp_neigh_idx;
    double* tmp_center_coord = (double*)malloc(sizeof(double) * 3);
    double* tmp_neigh_coord = (double*)malloc(sizeof(double) * 3);
    double* tmp_diff_coord = (double*)malloc(sizeof(double) * 3);
    double tmp_cutsq;

    // Step 1. Build neighbor list for `local atoms` which is input of `matersdk::deepPotSE::Se4pwOp`
    for (int ii=0; ii<this->list->inum; ii++) {
        tmp_center_idx = this->list->ilist[ii];
        tmp_center_coord[0] = this->atom->x[tmp_center_idx][0];
        tmp_center_coord[1] = this->atom->x[tmp_center_idx][1];
        tmp_center_coord[2] = this->atom->x[tmp_center_idx][2];

        local_ilist[ii] = this->list->ilist[ii];
        for (int jj=0; jj<this->list->numneigh[ii]; jj++) {
            tmp_neigh_idx = this->list->firstneigh[ii][jj];
            tmp_neigh_coord[0] = this->atom->x[tmp_neigh_idx][0];
            tmp_neigh_coord[1] = this->atom->x[tmp_neigh_idx][1];
            tmp_neigh_coord[2] = this->atom->x[tmp_neigh_idx][2];
            tmp_diff_coord[0] = tmp_neigh_coord[0] - tmp_center_coord[0];
            tmp_diff_coord[1] = tmp_neigh_coord[1] - tmp_center_coord[1];
            tmp_diff_coord[2] = tmp_neigh_coord[2] - tmp_center_coord[2];
            tmp_cutsq = (
                    tmp_diff_coord[0] * tmp_diff_coord[0] + 
                    tmp_diff_coord[1] * tmp_diff_coord[1] + 
                    tmp_diff_coord[2] * tmp_diff_coord[2]);
            if (tmp_cutsq < (this->cutgb*this->cutgb)) {
                local_firstneigh[ii*this->max_num_neigh_atoms+local_numneigh[ii]] = tmp_neigh_idx;
                local_numneigh[ii]++;
            }
        }
    }

    for (int ii=0; ii<local_inum; ii++) {
        if (local_numneigh[ii] > (this->max_num_neigh_atoms/this->atom->ntypes))
            this->error->all(FLERR, "The number of neighbor atoms is more than you set!!!");
    }

    // Step 2. Calculate tilde_r, tilde_r_deriv, relative_coord
    // Step 2.1. 
    c10::TensorOptions int_tensor_options = c10::TensorOptions().dtype(torch::kInt32).device(c10::kCPU);
    c10::TensorOptions float_tensor_options = c10::TensorOptions().dtype(torch::kFloat64).device(c10::kCPU);

    at::Tensor ilist_tensor = torch::from_blob(local_ilist, {1, local_inum}, int_tensor_options);
    at::Tensor numneigh_tensor = torch::from_blob(local_numneigh, {1, local_inum}, int_tensor_options);
    at::Tensor firstneigh_tensor = torch::from_blob(local_firstneigh, {1, local_inum, this->max_num_neigh_atoms}, int_tensor_options);

    at::Tensor types_tensor = torch::from_blob(this->atom->type, {1, this->atom->nlocal+this->atom->nghost}, int_tensor_options);
    types_tensor = types_tensor - 1;

    int* num_neigh_atoms_lst = (int*)malloc(sizeof(int) * this->atom->ntypes);
    for (int ii=0; ii<this->atom->ntypes; ii++) 
        num_neigh_atoms_lst[ii] = this->max_num_neigh_atoms / this->atom->ntypes;
    at::Tensor num_neigh_atoms_lst_tensor = torch::from_blob(num_neigh_atoms_lst, {1, this->atom->ntypes}, int_tensor_options);

    double* x_1d = (double*)malloc(sizeof(double) * (this->atom->nghost+this->atom->nlocal) * 3);
    memset(x_1d, 0, sizeof(double) * (this->atom->nghost+this->atom->nlocal) * 3);
    for (int ii=0; ii<(this->atom->nghost+this->atom->nlocal); ii++) {
        x_1d[ii*3 + 0] = this->atom->x[ii][0];
        x_1d[ii*3 + 1] = this->atom->x[ii][1];
        x_1d[ii*3 + 2] = this->atom->x[ii][2];
    }
    at::Tensor x_tensor = torch::from_blob(x_1d, {1, (this->atom->nghost+this->atom->nlocal), 3}, float_tensor_options);

    
    // Step 2.2. get tilde_r, tilde_r_deriv, relative_coord
    torch::autograd::variable_list features = matersdk::deepPotSE::Se4pwOp::forward(
            1,
            local_inum,
            ilist_tensor,
            numneigh_tensor,
            firstneigh_tensor,
            x_tensor,
            types_tensor,
            this->atom->ntypes,
            num_neigh_atoms_lst_tensor,
            this->cutgb,
            this->cutgbs);

    at::Tensor tilde_r_tensor = features[0];
    at::Tensor tilde_r_deriv_tensor = features[1];
    at::Tensor relative_coords_tensor = features[2];

    //std::cout << tilde_r_tensor << std::endl;
    //std::cout << tilde_r_deriv_tensor.sizes() << std::endl;
    //std::cout << relative_coords_tensor.sizes() << std::endl;
    
    // Step 2.3. Normalize `tilde_r`, `tilde_r_deriv` with `this->davg`, `this->dstd`
    at::Tensor prim_types_tensor = torch::index_select(types_tensor, 1, torch::arange(0, local_inum));
    for (int ii=0; ii<this->atom->ntypes; ii++) {
        auto mask = (prim_types_tensor.squeeze(0) == ii);
        at::Tensor selected_indices = mask.nonzero().squeeze(1).to(torch::kInt64);

        // Step 2.3.1. tilde_r
        at::Tensor selected_tilde_r_tensor = torch::index_select(tilde_r_tensor, 1, selected_indices);
        // selected_tilde_r: (1, N_c, N_b, 4); this->davg : (4,)
        selected_tilde_r_tensor = (selected_tilde_r_tensor - this->davg[ii].unsqueeze(0)) / this->dstd[ii].unsqueeze(0);
        torch::index_put_(tilde_r_tensor, {torch::tensor(0).to(torch::kInt64), selected_indices}, selected_tilde_r_tensor[0]);

        // Ste 2.3.2. tilde_r_deriv
        at::Tensor selected_tilde_r_deriv_tensor = torch::index_select(tilde_r_deriv_tensor, 1, selected_indices);
        // selected_tilde_r: (1, N_c, N_b, 4, 3); this->davg : (4,)
        selected_tilde_r_deriv_tensor = selected_tilde_r_deriv_tensor / this->dstd[ii].unsqueeze(0).unsqueeze(-1);
        torch::index_put_(tilde_r_deriv_tensor, {torch::tensor(0).to(torch::kInt64), selected_indices}, selected_tilde_r_deriv_tensor[0]);
    }

    // Step 2.4. Get `prim_indices_tensor`
    int* nstart_idxs = (int*)malloc(sizeof(int) * this->atom->ntypes);
    memset(nstart_idxs, 0, sizeof(int) * this->atom->ntypes);
    for (int ii=0; ii<this->atom->ntypes; ii++)
        for (int jj=0; jj<ii; jj++)
            nstart_idxs[ii] += num_neigh_atoms_lst[jj];
    int* nloop_idxs = (int*)malloc(sizeof(int) * this->atom->ntypes);
    memset(nloop_idxs, 0, sizeof(int) * this->atom->ntypes);
    int* prim_indices = (int*)malloc(sizeof(int) * local_inum * this->max_num_neigh_atoms);
    memset(prim_indices, 0, sizeof(int) * local_inum * this->max_num_neigh_atoms);  // global index in lammps starts from 1


    for (int ii=0; ii<local_inum; ii++) {
        tmp_center_idx = local_ilist[ii];
        for (int jj=0; jj<this->atom->ntypes; jj++)
            nloop_idxs[jj] = 0;
        
        for (int jj=0; jj<local_numneigh[ii]; jj++) {
            tmp_neigh_idx = local_firstneigh[ii*this->max_num_neigh_atoms + jj];

            int kk = this->atom->type[tmp_neigh_idx] - 1;   // start from 0/1

            prim_indices[(nstart_idxs[kk] + nloop_idxs[kk]) + ii*this->max_num_neigh_atoms] = this->atom->tag[tmp_neigh_idx];
            nloop_idxs[kk]++;
        }
    }
    at::Tensor prim_indices_tensor = torch::from_blob(prim_indices, {1, local_inum, this->max_num_neigh_atoms}, int_tensor_options);


    // Step 2.5. Get `natoms_image_tensor`
    at::Tensor natoms_image_tensor = at::zeros({1, 1+this->atom->ntypes}, int_tensor_options);
    natoms_image_tensor[0][0] = local_inum;
    for (int ii=0; ii<this->atom->ntypes; ii++) {
        for (int jj=0; jj<local_inum; jj++) {
            if (this->atom->type[local_ilist[jj]] == (ii+1)) {
                natoms_image_tensor[0][ii+1] = natoms_image_tensor[0][ii+1] + 1;
            }
        }
    }

    // Step 3. Forward
    tilde_r_tensor.requires_grad_(true);
    tilde_r_deriv_tensor.requires_grad_(true);
    relative_coords_tensor.requires_grad_(true);
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(tilde_r_tensor);
    inputs.push_back(tilde_r_deriv_tensor);
    inputs.push_back(prim_indices_tensor);
    inputs.push_back(natoms_image_tensor);
    inputs.push_back(atom_types_tensor);
    inputs.push_back(relative_coords_tensor);
    //std::cout << tilde_r_tensor.sizes() << std::endl;
    //std::cout << tilde_r_deriv_tensor.sizes() << std::endl;
    //std::cout << prim_indices_tensor.sizes() << std::endl;
    //std::cout << natoms_image_tensor.sizes() << std::endl;
    //std::cout << atom_types_tensor.sizes() << std::endl;
    //std::cout << relative_coords_tensor.sizes() << std::endl;
    auto outputs = this->module.forward(inputs);

    // Step 4. Return to lammps
    // Step 4.1. 
    at::Tensor etot_tensor = outputs.toTuple()->elements()[0].toTensor();   // [1, 1]
    at::Tensor eatom_tensor = outputs.toTuple()->elements()[1].toTensor();  // [1, Nc]
    at::Tensor fatom_tensor = outputs.toTuple()->elements()[2].toTensor();  // [1, Nc, 3]
    bool calc_virial_from_mlff = false;
    if (outputs.toTuple()->elements()[3].isTensor()) {
        calc_virial_from_mlff = true;
        at::Tensor vatom_tensor = outputs.toTuple()->elements()[3].toTensor();
    } else 
        auto vatom_tensor = outputs.toTuple()->elements()[3];   // may be `None` (IValue)

    // Step 4.2. 
    auto eatom_ptr = etot_tensor.accessor<double, 2>();
    auto fatom_ptr = fatom_tensor.accessor<double, 3>();
    
    double **f = atom->f;
    for (int ii=0; ii<local_inum; ii++) {
        f[ii][0] = fatom_ptr[0][ii][0];
        f[ii][1] = fatom_ptr[0][ii][1];
        f[ii][2] = fatom_ptr[0][ii][2];
    }

    if (eflag)
        eng_vdwl = etot_tensor[0][0].item<double>();
    if (eflag_atom) {
        for (int ii=0; ii<local_inum; ii++)
            eatom[ii] = eatom_ptr[0][ii];
    }

    if (vflag_fdotr) 
        virial_fdotr_compute();

    // Step . Free memory
    free(local_ilist);
    free(local_numneigh);
    free(local_firstneigh);
    free(tmp_center_coord);
    free(tmp_neigh_coord);
    free(tmp_diff_coord);
    free(nstart_idxs);
    free(nloop_idxs);
}


}   // namespace : LAMMPS_NS