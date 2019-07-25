#include <data.h>
#ifdef CPE
#include <slave.h>
#include <dma_macros.h>
#include <assert.h>
#include <cal.h>
void initial_integrate_nve_cpe(esmd_t *gl_md){
  dma_init();
  esmd_t md;
  pe_get(gl_md, &md, sizeof(esmd_t));
  dma_syn();
  box_t box;
  pe_get(md.box, &box, sizeof(box_t));
  potential_conf_t pot_conf;
  pe_get(md.pot_conf, &pot_conf, sizeof(potential_conf_t));
  dma_syn();

  areal dt = md.dt;
  areal dtf = dt * 0.5;
  areal *rmass = pot_conf.rmass;

  areal v[CELL_SIZE][3], x[CELL_SIZE][3], f[CELL_SIZE][3];
  int type[CELL_SIZE];
  int *nlocal = box.nlocal;
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  for (int icell = _MYID; icell < ncells; icell += 64){
    int kk = icell / (nlocal[0] * nlocal[1]);
    int jj = icell / nlocal[0] % nlocal[1];
    int ii = icell % nlocal[0];

    int cell_off = get_cell_off((&box), ii, jj, kk);
    cell_t cell;
    pe_get(box.cells + cell_off, &cell, sizeof(cell_t));
    dma_syn();
    celldata_t *data = box.celldata + cell_off;
    pe_get(data->v, v, sizeof(areal) * 3 * cell.natoms);
    pe_get(data->x, x, sizeof(areal) * 3 * cell.natoms);
    pe_get(data->f, f, sizeof(areal) * 3 * cell.natoms);
    pe_get(data->type, type, sizeof(int) * cell.natoms);
    dma_syn();

    for (int i = 0; i < cell.natoms; i ++){
      /* if (type[i] != 0) { */
      /*   cal_locked_printf("%d %d %d %d %d\n", ii, jj, kk, i, type[i]); */
      /*   assert(0); */
      /* } */
      v[i][0] += dtf * f[i][0] * rmass[type[i]];
      v[i][1] += dtf * f[i][1] * rmass[type[i]];
      v[i][2] += dtf * f[i][2] * rmass[type[i]];
      x[i][0] += dt * v[i][0];
      x[i][1] += dt * v[i][1];
      x[i][2] += dt * v[i][2];
    }
    pe_put(data->v, v, sizeof(areal) * 3 * cell.natoms);
    pe_put(data->x, x, sizeof(areal) * 3 * cell.natoms);
    dma_syn();
  }
}

void final_integrate_nve_cpe(esmd_t *gl_md){
  dma_init();
  esmd_t md;
  pe_get(gl_md, &md, sizeof(esmd_t));
  dma_syn();
  box_t box;
  pe_get(md.box, &box, sizeof(box_t));
  potential_conf_t pot_conf;
  pe_get(md.pot_conf, &pot_conf, sizeof(potential_conf_t));
  dma_syn();

  areal dt = md.dt;
  areal dtf = dt * 0.5;
  areal *rmass = pot_conf.rmass;

  areal v[CELL_SIZE][3], x[CELL_SIZE][3], f[CELL_SIZE][3];
  int type[CELL_SIZE];
  int *nlocal = box.nlocal;
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  for (int icell = _MYID; icell < ncells; icell += 64){
    int kk = icell / (nlocal[0] * nlocal[1]);
    int jj = icell / nlocal[0] % nlocal[1];
    int ii = icell % nlocal[0];

    int cell_off = get_cell_off((&box), ii, jj, kk);
    cell_t cell;
    pe_get(box.cells + cell_off, &cell, sizeof(cell_t));
    dma_syn();
    celldata_t *data = box.celldata + cell_off;
    pe_get(data->v, v, sizeof(areal) * 3 * cell.natoms);
    pe_get(data->f, f, sizeof(areal) * 3 * cell.natoms);
    pe_get(data->type, type, sizeof(int) * cell.natoms);
    dma_syn();
    for (int i = 0; i < cell.natoms; i ++){
      v[i][0] += dtf * f[i][0] * rmass[type[i]];
      v[i][1] += dtf * f[i][1] * rmass[type[i]];
      v[i][2] += dtf * f[i][2] * rmass[type[i]];
    }
    pe_put(data->v, v, sizeof(areal) * 3 * cell.natoms);
    dma_syn();
  }
}

#endif
#ifdef MPE
#include <timer.h>
#include <athread.h>
extern void slave_initial_integrate_nve_cpe(esmd_t *md);
void initial_integrate_nve_sw(esmd_t *md) {
  timer_start("initial_integrate");
  athread_spawn(initial_integrate_nve_cpe, md);
  athread_join();
  timer_stop("initial_integrate");
}

extern void slave_final_integrate_nve_cpe(esmd_t *md);
void final_integrate_nve_sw(esmd_t *md) {
  timer_start("final_integrate");
  athread_spawn(final_integrate_nve_cpe, md);
  athread_join();
  timer_stop("final_integrate");
}

#endif
