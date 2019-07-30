#ifdef CPE
static inline void CAT(pair_lj_force_b2b_self, VER_CODE)(cell_t *self, potential_conf_t *pot_conf, areal (*xi)[3], areal (*fi)[3], int *ti, areal (*fj)[3], areal *evdwl, areal *virial){
  int max_cut2 = pot_conf->cutoff * pot_conf->cutoff;
  lj_param_t *lj_param = &(pot_conf->param.lj);

  for (int i = 0; i < self->natoms; i ++) {
    ireal *cutoff2_i = lj_param->cutoff2[ti[i]];
    ireal *c6_i = lj_param->c6[ti[i]];
    ireal *c12_i = lj_param->c12[ti[i]];
    ireal *ec6_i = lj_param->ec6[ti[i]];
    ireal *ec12_i = lj_param->ec12[ti[i]];

    for (int j = 0; j < i; j ++) {
      int jtype = ti[j];
      ireal delx = xi[j][0] - xi[i][0];
      ireal dely = xi[j][1] - xi[i][1];
      ireal delz = xi[j][2] - xi[i][2];
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
	  *evdwl += r6inv * (r6inv * ec12_i[jtype] - ec6_i[jtype]) * 2;
	}
	if (VIRIAL){
	  *virial += r2 * force;
	}
      }
    }
  }
  for (int i = 0; i < self->natoms; i ++){
    fj[i][0] += fi[i][0];
    fj[i][1] += fi[i][1];
    fj[i][2] += fi[i][2];
  }
}
void CAT(pair_lj_force_b2b_neighbor, VER_CODE)(cell_t *self, cell_t *neigh, potential_conf_t *pot_conf, areal (*xi)[3], areal (*fi)[3], int *ti, areal (*xj)[3], areal (*fj)[3], int *tj, areal *evdwl, areal *virial){
  int max_cut2 = pot_conf->cutoff * pot_conf->cutoff;
  lj_param_t *lj_param = &(pot_conf->param.lj);
  areal (*bbox)[3] = self->bbox_ideal;
  areal box_o[3], box_h[3];
  box_o[0] = 0.5 * (bbox[0][0] + bbox[1][0]);
  box_o[1] = 0.5 * (bbox[0][1] + bbox[1][1]);
  box_o[2] = 0.5 * (bbox[0][2] + bbox[1][2]);
  box_h[0] = 0.5 * (bbox[1][0] - bbox[0][0]);
  box_h[1] = 0.5 * (bbox[1][1] - bbox[0][1]);
  box_h[2] = 0.5 * (bbox[1][2] - bbox[0][2]);
  doublev4 cutoff_prof[MAX_TYPES], c6_prof[MAX_TYPES], c12_prof[MAX_TYPES];
  doublev4 ec6_prof[MAX_TYPES], ec12_prof[MAX_TYPES];
  double fake_cut[MAX_TYPES];
  doublev4 fj_v4[CELL_SIZE][3], evdwl_v4 = 0, virial_v4 = 0;
  doublev4 c0_v4 = 0, c1_v4 = 1, c2_v4 = 2;
  lwpf_start(B2B);
  for (int j = 0; j < neigh->natoms; j ++){
    fj_v4[j][0] = c0_v4;
    fj_v4[j][1] = c0_v4;
    fj_v4[j][2] = c0_v4;
  }
  for (int itype = 0; itype < pot_conf->ntypes; itype ++){
    fake_cut[itype] = -1.0;
  }
  int jlist[CELL_SIZE];
  int jcnt = 0;
  for (int j = 0; j < neigh->natoms; j ++){
    if (dsq_atom_box(xj[j], box_o, box_h) < max_cut2)
      jlist[jcnt ++] = j;
  }
  int iptr = 0;
  //for (int i = 0; i < self.natoms; i += 4){
  while (iptr < self->natoms){
    lwpf_start(CONV);
    int ni = 0;
    int is[4];
    ireal *cut[4], *c6[4], *c12[4], *ec6[4], *ec12[4];
    lwpf_start(PROF);
    while (iptr < self->natoms && ni < 4){
      //if (dsq_atom_box(xi[iptr], box_o, box_h) < max_cut2){
      is[ni] = iptr;
      int itype = ti[iptr];
      cut [ni] = lj_param->cutoff2[itype];
      c6  [ni] = lj_param->c6  [itype];
      c12 [ni] = lj_param->c12 [itype];
      if (ENERGY){
	ec6 [ni] = lj_param->ec6 [itype];
	ec12[ni] = lj_param->ec12[itype];
      }
      ni ++;
      //}
      iptr ++;
    }
    lwpf_stop(PROF);
    for (int pad = ni; pad < 4; pad ++){
      int itype = 0;
      cut [pad] = fake_cut;
      c6  [pad] = lj_param->c6  [itype];
      c12 [pad] = lj_param->c12 [itype];
      if (ENERGY){
    	ec6 [pad] = lj_param->ec6 [itype];
    	ec12[pad] = lj_param->ec12[itype];
      }
      is[pad] = self->natoms;
    }

    doublev4 fixv4 = 0, fiyv4 = 0, fizv4 = 0;
    doublev4 xiv4, yiv4, ziv4, txiv4, tyiv4, tziv4;/* , tx0, tx1, tx2, tx3; */
    simd_load_4x3d(&xiv4, &yiv4, &ziv4, xi[is[0]]);
    
    for (int jtype = 0; jtype < pot_conf->ntypes; jtype ++){
      cutoff_prof[jtype] = simd_set_doublev4(cut[0][jtype], cut[1][jtype], cut[2][jtype], cut[3][jtype]);
      c6_prof[jtype] = simd_set_doublev4(c6[0][jtype], c6[1][jtype], c6[2][jtype], c6[3][jtype]);
      c12_prof[jtype] = simd_set_doublev4(c12[0][jtype], c12[1][jtype], c12[2][jtype], c12[3][jtype]);
      if (ENERGY){
    	ec6_prof[jtype] = simd_set_doublev4(ec6[0][jtype], ec6[1][jtype], ec6[2][jtype], ec6[3][jtype]);
    	ec12_prof[jtype] = simd_set_doublev4(ec12[0][jtype], ec12[1][jtype], ec12[2][jtype], ec12[3][jtype]);
      }
    }
    lwpf_stop(CONV);
    
    for (int jidx = 0; jidx < jcnt; jidx ++){
      int j = jlist[jidx];
      int jtype = tj[j];
      doublev4 xjv4 = xj[j][0];
      doublev4 yjv4 = xj[j][1];
      doublev4 zjv4 = xj[j][2];
      doublev4 delx_v4 = xjv4 - xiv4;
      doublev4 dely_v4 = yjv4 - yiv4;
      doublev4 delz_v4 = zjv4 - ziv4;
      doublev4 r2_v4 = delx_v4 * delx_v4 + dely_v4 * dely_v4 + delz_v4 * delz_v4;
      doublev4 cutoff_v4 = cutoff_prof[jtype];
      doublev4 ltcut = simd_vfcmplt(r2_v4, cutoff_v4);
      if (simd_vmatchd(ltcut, 0x40000000)){
	doublev4 r2inv_v4 = c1_v4 / r2_v4;
	doublev4 r6inv_v4 = r2inv_v4 * r2inv_v4 * r2inv_v4;
	doublev4 c6_v4   = c6_prof[jtype];
	doublev4 c12_v4  = c12_prof[jtype];
	doublev4 force_v4 = r6inv_v4 * (r6inv_v4 * c12_v4 - c6_v4) * r2inv_v4;
	force_v4 = simd_vseleq(ltcut, c0_v4, force_v4);
	fixv4 -= delx_v4 * force_v4;
	fiyv4 -= dely_v4 * force_v4;
	fizv4 -= delz_v4 * force_v4;
	fj_v4[j][0] += delx_v4 * force_v4;
	fj_v4[j][1] += dely_v4 * force_v4;
	fj_v4[j][2] += delz_v4 * force_v4;
	if (ENERGY) {
	  doublev4 ec6_v4  = ec6_prof[jtype];
	  doublev4 ec12_v4 = ec12_prof[jtype];
	  doublev4 vdw_tmp_v4 = r6inv_v4 * (r6inv_v4 * ec12_v4 - ec6_v4);
	  vdw_tmp_v4 = simd_vseleq(ltcut, c0_v4, vdw_tmp_v4);
	  evdwl_v4 += vdw_tmp_v4;
	}
	if (VIRIAL) {
	  virial_v4 += r2_v4 * force_v4;
	}
      }
    }
    lwpf_start(CONV);
    fi[is[0]][0] += simd_vextf0(fixv4);
    fi[is[0]][1] += simd_vextf0(fiyv4);
    fi[is[0]][2] += simd_vextf0(fizv4);
    fi[is[1]][0] += simd_vextf1(fixv4);
    fi[is[1]][1] += simd_vextf1(fiyv4);
    fi[is[1]][2] += simd_vextf1(fizv4);
    fi[is[2]][0] += simd_vextf2(fixv4);
    fi[is[2]][1] += simd_vextf2(fiyv4);
    fi[is[2]][2] += simd_vextf2(fizv4);
    fi[is[3]][0] += simd_vextf3(fixv4);
    fi[is[3]][1] += simd_vextf3(fiyv4);
    fi[is[3]][2] += simd_vextf3(fizv4);
    lwpf_stop(CONV);

  }
  lwpf_start(CONV);
  for (int j = 0; j < neigh->natoms; j ++){
    fj[j][0] += simd_vsumd(fj_v4[j][0]);
    fj[j][1] += simd_vsumd(fj_v4[j][1]);
    fj[j][2] += simd_vsumd(fj_v4[j][2]);
  }
  lwpf_stop(CONV);
  *evdwl += simd_vsumd(evdwl_v4 * 2);
  *virial += simd_vsumd(virial_v4);
  lwpf_stop(B2B);
}

void CAT(pair_lj_force_cpe, VER_CODE)(esmd_t *gl_md) {
  /* if (_MYID == 0 && gl_md->mpp->pid == 0){ */
  /*   test_4x3d(); */
  /* } */
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
  int ncellsall = nall[0] * nall[1] * nall[2];
  lwpf_start(FILL);

  for (int icell = _MYID; icell < ncellsall; icell += 64){
    cell_t cell;
    pe_get(box.cells + icell, &cell, sizeof(cell_t));
    dma_syn();
    int natoms_rep = cell.nreplicas * cell.natoms;
    areal (*gl_frep)[3] = cell.frep[0];
    for (int offset = 0; offset < natoms_rep; offset += CELL_SIZE) {
      pe_put(gl_frep + offset, fi, CELL_SIZE * sizeof(areal) * 3);
      dma_syn();
    }
  }

  lwpf_stop(FILL);
  athread_syn(ARRAY_SCOPE, 0xffff);
  //if (_MYID > 0) return;
  areal max_cut2 = pot_conf.cutoff * pot_conf.cutoff;
  unsigned long long lowerpe_mask = (1ULL << _MYID) - 1;
  lwpf_start(COMP);
  int *start = swdata.start[_MYID], *count = swdata.count[_MYID];

  for (int icellst = _MYID * BLK; icellst < ncells; icellst += NCPE * BLK){
    int icelled = icellst + BLK;
    if (icelled > ncells) icelled = ncells;
    for (int icell = icellst; icell < icelled; icell ++){
      int kk = icell / (nlocal[0] * nlocal[1]);
      int jj = icell / nlocal[0] % nlocal[1];
      int ii = icell % nlocal[0];
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
	    int self_interaction = (dx == 0 && dy == 0 && dz == 0);
            cell_t cell_neigh;
            lwpf_start(RJ);
            pe_get(box.cells + neigh_off, &cell_neigh, sizeof(cell_t));
            dma_syn();
            asm ("ctpop %1, %0\n\t" : "=r"(irep): "r"(lowerpe_mask & cell_neigh.pemask));
            pe_get(box.celldata[neigh_off].x, xj, cell_neigh.natoms * sizeof(areal) * 3);
            pe_get(box.celldata[neigh_off].type, tj, cell_neigh.natoms * sizeof(int));
	    areal (*gl_frep)[3] = cell_neigh.frep[0];
	    int roff = irep * cell_neigh.natoms;
            pe_get(gl_frep + roff, fj, cell_neigh.natoms * sizeof(areal) * 3);
            dma_syn();
            lwpf_stop(RJ);
	    if (self_interaction){
	      CAT(pair_lj_force_b2b_self, VER_CODE)(&cell_self, &pot_conf,
						    xi, fi, ti, fj,
						    &evdwl, &virial);
	    } else {
	      CAT(pair_lj_force_b2b_neighbor, VER_CODE)(&cell_self, &cell_neigh, &pot_conf,
							xi, fi, ti, xj, fj, tj,
							&evdwl, &virial);
	    }
            lwpf_start(WJ);
            pe_put(gl_frep + roff, fj, sizeof(areal) * 3 * cell_neigh.natoms);
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
  for (int icell = _MYID; icell < ncellsall; icell += 64){
    cell_t cell;
    pe_get(box.cells + icell, &cell, sizeof(cell_t));
    memset(fi, 0, CELL_SIZE * sizeof(areal) * 3);
    dma_syn();
    areal (*gl_frep)[3] = cell.frep[0];

    for (int irep = 0; irep < cell.nreplicas; irep ++){
      int roff = irep * cell.natoms;
      pe_get(gl_frep + roff, fj, cell.natoms * sizeof(areal) * 3);
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
