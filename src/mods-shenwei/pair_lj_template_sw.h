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
  double (*cut )[MAX_TYPES] = lj_param->cutoff2;
  double (*c6  )[MAX_TYPES] = lj_param->c6;
  double (*c12 )[MAX_TYPES] = lj_param->c12;
  double (*ec6 )[MAX_TYPES] = lj_param->ec6;
  double (*ec12)[MAX_TYPES] = lj_param->ec12;
  doublev4 ec6_prof[MAX_TYPES], ec12_prof[MAX_TYPES];
  double fake_cut[MAX_TYPES];
  doublev4 fj_v4[CELL_SIZE][3], evdwl_v4 = 0, virial_v4 = 0;
  doublev4 c0_v4 = 0, c1_v4 = 1, c2_v4 = 2;
  lwpf_start(B2B);
  lwpf_start(PROF);
  for (int j = 0; j < neigh->natoms; j ++){
    fj_v4[j][0] = c0_v4;
    fj_v4[j][1] = c0_v4;
    fj_v4[j][2] = c0_v4;
  }
  for (int itype = 0; itype < pot_conf->ntypes; itype ++){
    cut[pot_conf->ntypes][itype] = -1.0;
  }
  for (int i = self->natoms; i < self->natoms + 4; i ++) {
    ti[i] = pot_conf->ntypes;
  }
  int jlist[CELL_SIZE];
  int jcnt = 0;
  for (int j = 0; j < neigh->natoms; j ++){
    if (dsq_atom_box(xj[j], box_o, box_h) < max_cut2)
      jlist[jcnt ++] = j;
  }
  if (pot_conf->ntypes == 1){
    cutoff_prof[0] = cut[0][0];
    c6_prof[0] = c6[0][0];
    c12_prof[0] = c12[0][0];
    ec6_prof[0] = ec6[0][0];
    ec12_prof[0] = ec12[0][0];
  }
  lwpf_stop(PROF);
  for (int i = 0; i < self->natoms; i += 4){
    lwpf_start(CONV);
    int ni = 4;
    if (ni > self->natoms - i) ni = self->natoms - i;
    int it0 = ti[i + 0];
    int it1 = ti[i + 1];
    int it2 = ti[i + 2];
    int it3 = ti[i + 3];
    doublev4 fixv4 = 0, fiyv4 = 0, fizv4 = 0;
    doublev4 xiv4, yiv4, ziv4, txiv4, tyiv4, tziv4;/* , tx0, tx1, tx2, tx3; */
    simd_load_4x3d(&xiv4, &yiv4, &ziv4, xi[i]);
    if (pot_conf->ntypes > 1){
      for (int jt = 0; jt < pot_conf->ntypes; jt ++){
        cutoff_prof[jt] = simd_set_doublev4(cut[it0][jt], cut[it1][jt], cut[it2][jt], cut[it3][jt]);
        c6_prof[jt] = simd_set_doublev4(c6[it0][jt], c6[it1][jt], c6[it2][jt], c6[it3][jt]);
        c12_prof[jt] = simd_set_doublev4(c12[it0][jt], c12[it1][jt], c12[it2][jt], c12[it3][jt]);
        if (ENERGY){
          ec6_prof[jt] = simd_set_doublev4(ec6[it0][jt], ec6[it1][jt], ec6[it2][jt], ec6[it3][jt]);
          ec12_prof[jt] = simd_set_doublev4(ec12[it0][jt], ec12[it1][jt], ec12[it2][jt], ec12[it3][jt]);
        }
      }
    } else if (ni != 4){
      double *tp = (double*)cutoff_prof;
      for (int pad = ni; pad < 4; pad ++){
        tp[pad] = -1.0;
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
      int has_lt_cut;
      simd_vmatchd_m(has_lt_cut, ltcut, 0x40000000);
      if (has_lt_cut){
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
    /*
t0 = vshff( r2 , r1 , 0x44 )
t1 = vshff( r2 , r1 , 0xee )
t2 = vshff( r1 , r0 , 0xee )
tt0 = vshff( t0 , r0 , 0x84 )
tt2 = vshff( t2 , r2 , 0xde )
w0 = vshff( tt0 , tt0 , 0x78 )
w1 = vshff( t2 , t0 , 0x8d )
w2 = vshff( t1 , tt2 , 0xd8 )
['02', '12', '22', '32'] ['02', '12', '22', '32'] 0x44 => ['01', '11', '02', '12']
['02', '12', '22', '32'] ['02', '12', '22', '32'] 0xee => ['21', '31', '22', '32']
['01', '11', '21', '31'] ['01', '11', '21', '31'] 0xee => ['20', '30', '21', '31']
['01', '11', '02', '12'] ['01', '11', '02', '12'] 0x84 => ['00', '10', '01', '02']
['20', '30', '21', '31'] ['20', '30', '21', '31'] 0xde => ['22', '32', '30', '31']
['00', '10', '01', '02'] ['00', '10', '01', '02'] 0x78 => ['00', '01', '02', '10']
['20', '30', '21', '31'] ['20', '30', '21', '31'] 0x8d => ['11', '12', '20', '21']
['21', '31', '22', '32'] ['21', '31', '22', '32'] 0xd8 => ['22', '30', '31', '32']

     */
    fi[i + 0][0] += simd_vextf0(fixv4);
    fi[i + 0][1] += simd_vextf0(fiyv4);
    fi[i + 0][2] += simd_vextf0(fizv4);
    fi[i + 1][0] += simd_vextf1(fixv4);
    fi[i + 1][1] += simd_vextf1(fiyv4);
    fi[i + 1][2] += simd_vextf1(fizv4);
    fi[i + 2][0] += simd_vextf2(fixv4);
    fi[i + 2][1] += simd_vextf2(fiyv4);
    fi[i + 2][2] += simd_vextf2(fizv4);
    fi[i + 3][0] += simd_vextf3(fixv4);
    fi[i + 3][1] += simd_vextf3(fiyv4);
    fi[i + 3][2] += simd_vextf3(fizv4);
    lwpf_stop(CONV);

  }
  lwpf_start(CONV);
  for (int j = 0; j < neigh->natoms; j ++){
    simd_vsumd_m(fj_v4[j][0]);
    simd_vsumd_m(fj_v4[j][1]);
    simd_vsumd_m(fj_v4[j][2]);
    fj[j][0] += fj_v4[j][0];
    fj[j][1] += fj_v4[j][1];
    fj[j][2] += fj_v4[j][2];
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
  areal xj_cache[9][CELL_SIZE][3], fj_cache[9][CELL_SIZE][3];
  int tj_cache[9][CELL_SIZE][3];
  cell_t cj_cache[9];
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
  //int *start = swdata.start[_MYID], *count = swdata.count[_MYID];
  int *read_from = swdata.read_from, *store_to = swdata.store_to;
  for (int icellst = _MYID * BLK; icellst < ncells; icellst += NCPE * BLK){
    int icelled = icellst + BLK;
    if (icelled > ncells) icelled = ncells;
    for (int icell = icellst; icell < icelled; icell ++){
      int kk = icell / (nlocal[0] * nlocal[1]);
      int jj = icell / nlocal[0] % nlocal[1];
      int ii = icell % nlocal[0];
      int read_last = (icell != icellst && ii != 0);
      int stor_next = (ii != nlocal[0] - 1 && icell != icelled - 1);
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
      for (int dz = -NCELL_CUT; dz <= 0; dz ++) {
	int dytop = (dz == 0) ? 0 : NCELL_CUT;
	for (int dy = -NCELL_CUT; dy <= dytop; dy ++) {
	  int dxtop = (dz == 0 && dy == 0) ? 0 : NCELL_CUT;
	  for (int dx = -NCELL_CUT; dx <= dxtop; dx ++) {
	    int should_stor = stor_next && dx > -NCELL_CUT;
	    int should_read = read_last && dx < NCELL_CUT;
            int neigh_x = dx + ii;
            int neigh_y = dy + jj;
            int neigh_z = dz + kk;
            int neigh_off = get_cell_off((&box), ii + dx, jj + dy, kk + dz);
	    int dn = dx + dy * 3 + dz * 9 + 13;
            int irep;
	    int self_interaction = (dx == 0 && dy == 0 && dz == 0);
            cell_t cell_neigh;
            lwpf_start(RJ);
            pe_get(box.cells + neigh_off, &cell_neigh, sizeof(cell_t));
            dma_syn();
            asm ("ctpop %1, %0\n\t" : "=r"(irep): "r"(lowerpe_mask & cell_neigh.pemask));
	    areal (*gl_frep)[3] = cell_neigh.frep[0];
	    int roff = irep * cell_neigh.natoms;
            //pe_get(gl_frep + roff, fj, cell_neigh.natoms * sizeof(areal) * 3);
	    if (should_read && read_from[dn] != -1){
	      simd_cpyo(xj, xj_cache[read_from[dn]], cell_neigh.natoms * sizeof(areal) * 3);
	      simd_cpyo(fj, fj_cache[read_from[dn]], cell_neigh.natoms * sizeof(areal) * 3);
	      simd_cpyo(tj, tj_cache[read_from[dn]], cell_neigh.natoms * sizeof(int));
	    } else {
              pe_get(gl_frep + roff, fj, cell_neigh.natoms * sizeof(areal) * 3);
	      pe_get(box.celldata[neigh_off].x, xj, cell_neigh.natoms * sizeof(areal) * 3);
              pe_get(box.celldata[neigh_off].type, tj, cell_neigh.natoms * sizeof(int));
	    }
	    //
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
	    if (should_stor && store_to[dn] != -1){
	      simd_cpyo(xj_cache[store_to[dn]], xj, sizeof(areal) * 3 * cell_neigh.natoms);
              simd_cpyo(fj_cache[store_to[dn]], fj, sizeof(areal) * 3 * cell_neigh.natoms);
	      simd_cpyo(tj_cache[store_to[dn]], tj, sizeof(int) * cell_neigh.natoms);
	    } else {
              pe_put(gl_frep + roff, fj, sizeof(areal) * 3 * cell_neigh.natoms);
            }
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
