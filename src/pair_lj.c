#include <data.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
//#define DEBUG_THIS_FILE
#include <log.h>

void pair_lj_setup(esmd_t *md, areal *cutoff, ireal *epsilon, ireal *sigma, ireal *mass, int ntypes){
  assert(ntypes < MAX_TYPES);
  md->pair_conf.cutoff = 0;
  for (int i = 0; i < ntypes; i ++){
    md->pair_conf.rmass[i] = 1.0 / mass[i];
    md->pair_conf.mass[i] = mass[i];
    for (int j = 0; j < ntypes; j ++){
      ireal sigma2 = sigma[i * ntypes + j] * sigma[i * ntypes + j];
      ireal sigma6 = sigma2 * sigma2 * sigma2;
      ireal sigma12 = sigma6 * sigma6;
      md->pair_conf.lj_param.c6[i][j] = 4.0 * 6.0 * epsilon[i * ntypes + j] * sigma6;
      md->pair_conf.lj_param.c12[i][j] = 4.0 * 12.0 * epsilon[i * ntypes + j] * sigma12;
      md->pair_conf.lj_param.ec6[i][j] = 4.0 * epsilon[i * ntypes + j] * sigma6;
      md->pair_conf.lj_param.ec12[i][j] = 4.0 * epsilon[i * ntypes + j] * sigma12;
      md->pair_conf.lj_param.cutoff2[i][j] = cutoff[i * ntypes + j] * cutoff[i * ntypes + j];
      if (cutoff[i * ntypes + j] > md->pair_conf.cutoff){
        md->pair_conf.cutoff = cutoff[i * ntypes + j];
      }
    }
  }
}

void pair_lj_force(esmd_t *md) {
  box_t *box = &(md->box);
  lj_param_t *lj_param = &(md->pair_conf.lj_param);
  areal evdwl = 0;
  areal virial = 0;
  int total_atoms = 0;
  for (int kk = -NCELL_CUT; kk < box->nlocal[2] + NCELL_CUT; kk ++){
    for (int jj = -NCELL_CUT; jj < box->nlocal[1] + NCELL_CUT; jj ++){
      for (int ii = -NCELL_CUT; ii < box->nlocal[0] + NCELL_CUT; ii ++){
        int self_off = get_cell_off(box, ii, jj, kk);
	cell_t *cell_self = box->cells + self_off;
	celldata_t *data_self = box->celldata + self_off;
	areal (*fi)[3] = data_self->f;
	for (int i = 0; i < cell_self->natoms; i ++) {
	  fi[i][0] = 0;
	  fi[i][1] = 0;
	  fi[i][2] = 0;
	}
      }
    }
  }

  double maxf = 0;
  int total_int = 0;
  for (int kk = 0; kk < box->nlocal[2]; kk ++){
    for (int jj = 0; jj < box->nlocal[1]; jj ++){
      for (int ii = 0; ii < box->nlocal[0]; ii ++){   
        int self_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell_self = box->cells + self_off;
        celldata_t *data_self = box->celldata + self_off;
        areal (*fi)[3] = data_self->f;
        areal (*xi)[3] = data_self->x;
        int   *ti = data_self->type;
        total_atoms += cell_self->natoms;
        for (int dx = -NCELL_CUT; dx <= 0; dx ++) {
          int dytop = (dx == 0) ? 0 : NCELL_CUT;
          for (int dy = -NCELL_CUT; dy <= dytop; dy ++) {
            int dztop = (dx == 0 && dy == 0) ? 0 : NCELL_CUT;
            for (int dz = -NCELL_CUT; dz <= dztop; dz ++) {
	      /* for (int dx = -NCELL_CUT; dx <= NCELL_CUT; dx ++) { */
	      /*   for (int dy = -NCELL_CUT; dy <= NCELL_CUT; dy ++) { */
	      /*     for (int dz = -NCELL_CUT; dz <= NCELL_CUT; dz ++) { */
              int neigh_x = dx + ii;
              int neigh_y = dy + jj;
              int neigh_z = dz + kk;
              int neigh_off = get_cell_off(box, ii + dx, jj + dy, kk + dz);
              
              cell_t *cell_neigh = box->cells + neigh_off;
              celldata_t *data_neigh = box->celldata + neigh_off;

              areal (*fj)[3] = data_neigh->f;
              areal (*xj)[3] = data_neigh->x;
              int   *tj = data_neigh->type;

              int self_interaction = 0;
              if (dx == 0 && dy == 0 && dz == 0) {
                self_interaction = 1;
              }

              for (int i = 0; i < cell_self->natoms; i ++) {
                int jtop = self_interaction ? i : cell_neigh->natoms;
                ireal *cutoff2_i = lj_param->cutoff2[ti[i]];
                ireal *c6_i = lj_param->c6[ti[i]];
                ireal *c12_i = lj_param->c12[ti[i]];
                ireal *ec6_i = lj_param->ec6[ti[i]];
                ireal *ec12_i = lj_param->ec12[ti[i]];
                for (int j = 0; j < jtop; j ++) {
                  int jtype = tj[j];
                  ireal delx = xj[j][0] - xi[i][0];
                  ireal dely = xj[j][1] - xi[i][1];
                  ireal delz = xj[j][2] - xi[i][2];
                  ireal r2 = delx * delx + dely * dely + delz * delz;
                  if (r2 < cutoff2_i[jtype] - TINY) {
                    total_int += 2;
                    ireal r2inv = 1.0 / r2;
                    ireal r6inv = r2inv * r2inv * r2inv;
                    ireal force = r6inv * (r6inv * c12_i[jtype] - c6_i[jtype]) * r2inv;

                    fi[i][0] -= delx * force;
                    fi[i][1] -= dely * force;
                    fi[i][2] -= delz * force;
                    
                    fj[j][0] += delx * force;
                    fj[j][1] += dely * force;
                    fj[j][2] += delz * force;
                    evdwl += r6inv * (r6inv * ec12_i[jtype] - ec6_i[jtype]) * 2;
		    virial += r2 * force;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  md->accu_local.epot += evdwl;
  md->accu_local.virial += virial;
#ifdef DEBUG_THIS_FILE
  areal evdwl_gbl;
  esmd_global_sum_scalar(md, &evdwl_gbl, evdwl);

  if (md->mpp.pid == 0){
    debug("%d %f %d %d\n", md->mpp.pid, evdwl_gbl, total_atoms, total_int);
  }
#endif
}
