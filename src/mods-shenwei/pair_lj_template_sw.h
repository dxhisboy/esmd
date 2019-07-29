#ifdef CPE

void CAT(pair_lj_force_cpe, VER_CODE)(esmd_t *gl_md) {
  lwpf_enter(PAIR_LJ);
  lwpf_start(ALL);
  dma_init();
  esmd_t md;
  pe_get(gl_md, &md, sizeof(esmd_t));
  dma_syn();
  box_t box;
  potential_conf_t pot_conf;
  swdata_t swdata;
  pe_get(md.platformdata, &swdata, sizeof(swdata_t));
  pe_get(md.box, &box, sizeof(box_t));
  pe_get(md.pot_conf, &pot_conf, sizeof(potential_conf_t));
  dma_syn();
  //box_t *box = &(md->box);
  lj_param_t *lj_param = &(pot_conf.param.lj);
  areal evdwl = 0;
  areal virial = 0;
  areal xi[CELL_SIZE][3], fi[CELL_SIZE][3];
  int ti[CELL_SIZE];
  areal xj[CELL_SIZE][3], fj[CELL_SIZE][3];
  int tj[CELL_SIZE];
  memset(fi, 0, sizeof(areal) * CELL_SIZE * 3);
  int *nlocal = box.nlocal;
  int *nall = box.nall;
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  int nallcells = nall[0] * nall[1] * nall[2];
  lwpf_start(FILL);
  //for (int icell = _MYID; icell < nallcells; icell += 64){
  /* cell_t cell_self; */
  /* pe_get(box.cells + self_off, &cell_self, sizeof(cell_t)); */
  /* dma_syn(); */
  for (int irep = _MYID; irep < swdata.nupdate; irep += 64){
    //celldata_t *data_self = box.celldata + icell;
    pe_put(swdata.frep + irep, fi, sizeof(areal) * 3 * CELL_SIZE);
    dma_syn();
  }
  //}
  lwpf_stop(FILL);
  athread_syn(ARRAY_SCOPE, 0xffff);
  //if (_MYID > 0) return;
  areal max_cut2 = pot_conf.cutoff * pot_conf.cutoff;
  unsigned long long lowerpe_mask = (1ULL << _MYID) - 1;
  lwpf_start(COMP);
  int *start = swdata.start[_MYID], *count = swdata.count[_MYID];
  /* for (int icell = _MYID; icell < ncells; icell += 64){ */
  /*   int kk = icell / (nlocal[0] * nlocal[1]); */
  /*   int jj = icell / nlocal[0] % nlocal[1]; */
  /*   int ii = icell % nlocal[0]; */
  for (int kk = start[2]; kk < start[2] + count[2]; kk ++){
    for (int jj = start[1]; jj < start[1] + count[1]; jj ++){
      for (int ii = start[0]; ii < start[0] + count[0]; ii ++){
        int self_off = get_cell_off((&box), ii, jj, kk);
        cell_t cell_self;
        lwpf_start(RI);
        pe_get(box.cells + self_off, &cell_self, sizeof(cell_t));
        dma_syn();
        pe_get(box.celldata[self_off].x, xi, cell_self.natoms * sizeof(areal) * 3);
        pe_get(box.celldata[self_off].type, ti, cell_self.natoms * sizeof(int));
        memset(fi, 0, cell_self.natoms * sizeof(areal) * 3);
        dma_syn();
        lwpf_stop(RI);

        for (int dx = -NCELL_CUT; dx <= 0; dx ++) {
          int dytop = (dx == 0) ? 0 : NCELL_CUT;
          for (int dy = -NCELL_CUT; dy <= dytop; dy ++) {
            int dztop = (dx == 0 && dy == 0) ? 0 : NCELL_CUT;
            for (int dz = -NCELL_CUT; dz <= dztop; dz ++) {
              int neigh_x = dx + ii;
              int neigh_y = dy + jj;
              int neigh_z = dz + kk;
              int neigh_off = get_cell_off((&box), ii + dx, jj + dy, kk + dz);
              int irep;

              cell_t cell_neigh;
              lwpf_start(RJ);
              pe_get(box.cells + neigh_off, &cell_neigh, sizeof(cell_t));
              dma_syn();
              asm ("ctpop %1, %0\n\t" : "=r"(irep): "r"(lowerpe_mask & cell_neigh.pemask));
              pe_get(box.celldata[neigh_off].x, xj, cell_neigh.natoms * sizeof(areal) * 3);
              pe_get(box.celldata[neigh_off].type, tj, cell_neigh.natoms * sizeof(int));
              if (md.mpp->pid == 0 && _MYID == 20){
                //printf("%d %d %d %p %d %d %p\n", ii + dx, jj + dy, kk + dz, cell_neigh.frep, irep, neigh_off, box.cells[5121].frep);
              }
              pe_get(cell_neigh.frep + irep, fj, cell_neigh.natoms * sizeof(areal) * 3);
              dma_syn();
              lwpf_stop(RJ);
              int self_interaction = 0;
              if (dx == 0 && dy == 0 && dz == 0) {
                self_interaction = 1;
              }
              areal (*bbox)[3] = cell_neigh.bbox_ideal;
              areal box_o[3], box_h[3];
              box_o[0] = 0.5 * (bbox[0][0] + bbox[1][0]);
              box_o[1] = 0.5 * (bbox[0][1] + bbox[1][1]);
              box_o[2] = 0.5 * (bbox[0][2] + bbox[1][2]);
              box_h[0] = 0.5 * (bbox[1][0] - bbox[0][0]);
              box_h[1] = 0.5 * (bbox[1][1] - bbox[0][1]);
              box_h[2] = 0.5 * (bbox[1][2] - bbox[0][2]);
              for (int i = 0; i < cell_self.natoms; i ++) {
                int jtop = self_interaction ? i : cell_neigh.natoms;
                //if (dsq_atom_box(xi[i], box_o, box_h) >= max_cut2) continue;
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
                  if (r2 < cutoff2_i[jtype]) {
                    ireal r2inv = 1.0 / r2;
                    ireal r6inv = r2inv * r2inv * r2inv;
                    ireal force = r6inv * (r6inv * c12_i[jtype] - c6_i[jtype]) * r2inv;

                    fi[i][0] -= delx * force;
                    fi[i][1] -= dely * force;
                    fi[i][2] -= delz * force;
                    
                    fj[j][0] += delx * force;
                    fj[j][1] += dely * force;
                    fj[j][2] += delz * force;
                    if (ENERGY){
                      evdwl += r6inv * (r6inv * ec12_i[jtype] - ec6_i[jtype]) * 2;
                    }
                    if (VIRIAL){
                      virial += r2 * force;
                    }
                  }
                }
              }
              if (self_interaction){
                for (int i = 0; i < cell_self.natoms; i ++){
                  fj[i][0] += fi[i][0];
                  fj[i][1] += fi[i][1];
                  fj[i][2] += fi[i][2];
                }
              }
              lwpf_start(WJ);
              pe_put(cell_neigh.frep + irep, fj, sizeof(areal) * 3 * cell_neigh.natoms);
              dma_syn();
              lwpf_stop(WJ);
              lwpf_start(SYN);
              //athread_syn(ARRAY_SCOPE, 0xffff);
              lwpf_stop(SYN);
            }
          }
        }
      }
    }
  }
  lwpf_stop(COMP);
  athread_syn(ARRAY_SCOPE, 0xffff);
  /* if (ncells % 64 != 0 && _MYID >= ncells % 64){ */
  /*   int scan_w = NCELL_CUT * 2 + 1; */
  /*   int ncell_neigh = (scan_w * scan_w * scan_w - 1) / 2 + 1; */
  /*   for (int i = 0; i < ncell_neigh; i ++){ */
  /*     athread_syn(ARRAY_SCOPE, 0xffff); */
  /*   } */
  /* } */
  lwpf_start(FINI);
  int ncellsall = box.nall[0] * box.nall[1] * box.nall[2];
  for (int icell = _MYID; icell < ncellsall; icell += 64){
    cell_t cell;
    pe_get(box.cells + icell, &cell, sizeof(cell_t));
    memset(fi, 0, CELL_SIZE * sizeof(areal) * 3);
    dma_syn();
    for (int irep = 0; irep < cell.nreplicas; irep ++){
      pe_get(cell.frep + irep, fj, cell.natoms * sizeof(areal) * 3);
      dma_syn();
      for (int i = 0; i < cell.natoms; i ++){
        fi[i][0] += fj[i][0];
        fi[i][1] += fj[i][1];
        fi[i][2] += fj[i][2];
      }
    }
    pe_put(box.celldata[icell].f, fi, cell.natoms * sizeof(areal) * 3);
    dma_syn();
  }

  /* athread_syn(ARRAY_SCOPE, 0xffff); */

  doublev4 ev = 0, ev4 = evdwl, vv4 = virial;
  ev = simd_vinsf0(ev4, ev);
  ev = simd_vinsf1(vv4, ev);
  ev = reg_reduce_doublev4(ev);
  if (_MYID == 0){
    gl_md->accu_local.epot += simd_vextf0(ev);
    gl_md->accu_local.virial += simd_vextf1(ev);
  }
  lwpf_stop(FINI);
  lwpf_stop(ALL);
  lwpf_exit(PAIR_LJ);
}
#endif
#ifdef MPE
extern void CAT(slave_pair_lj_force_cpe, VER_CODE)(esmd_t *gl_md);
#endif
