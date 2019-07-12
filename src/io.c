#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include "io.h"
#include "data.h"
#define NATOMS_BUF 2048
void load_raw_x_atoms(esmd_t *md, const char *path){
  int fd = open(path, O_RDONLY);
  puts(path);
  size_t buffer_size = NATOMS_BUF * 3 * sizeof(double);
  double (*buffer)[3] = esmd_malloc(buffer_size, "load_raw_buffer");
  size_t count;

  box_t *box = &(md->box);
  double *rlcell = box->rlcell;
  while ((count = read(fd, buffer, buffer_size)) > 0) {
    int natoms = count / (3 * sizeof(double));
    for (int i = 0; i < natoms; i ++){
      if (buffer[i][0] < 0) buffer[i][0] += box->lglobal[0];
      if (buffer[i][1] < 0) buffer[i][1] += box->lglobal[1];
      if (buffer[i][2] < 0) buffer[i][2] += box->lglobal[2];
      if (buffer[i][0] > box->lglobal[0]) buffer[i][0] -= box->lglobal[0];
      if (buffer[i][1] > box->lglobal[1]) buffer[i][1] -= box->lglobal[1];
      if (buffer[i][2] > box->lglobal[2]) buffer[i][2] -= box->lglobal[2];
      int cellx = floor(buffer[i][0] * rlcell[0] + TINY);
      int celly = floor(buffer[i][1] * rlcell[1] + TINY);
      int cellz = floor(buffer[i][2] * rlcell[2] + TINY);
      int celloff = get_cell_off(box, cellx, celly, cellz);
      cell_t *cell = box->cells + celloff;
      celldata_t *celldata = box->celldata + celloff;
      int curatom = cell->natoms;
      celldata->x[curatom][0] = buffer[i][0];
      celldata->x[curatom][1] = buffer[i][1];
      celldata->x[curatom][2] = buffer[i][2];
      celldata->type[curatom] = 0;
      cell->natoms ++;
    }
  }
  esmd_free(buffer);
}

void load_raw_xv_atoms(esmd_t *md, const char *path){
  int fd = open(path, O_RDONLY);
  puts(path);
  size_t buffer_size = NATOMS_BUF * 6 * sizeof(double);
  double (*buffer)[6] = esmd_malloc(buffer_size, "load_raw_buffer");
  size_t count;

  box_t *box = &(md->box);
  double *rlcell = box->rlcell;
  while ((count = read(fd, buffer, buffer_size)) > 0) {
    int natoms = count / (6 * sizeof(double));
    for (int i = 0; i < natoms; i ++){
      if (buffer[i][0] < 0) buffer[i][0] += box->lglobal[0];
      if (buffer[i][1] < 0) buffer[i][1] += box->lglobal[1];
      if (buffer[i][2] < 0) buffer[i][2] += box->lglobal[2];
      if (buffer[i][0] > box->lglobal[0]) buffer[i][0] -= box->lglobal[0];
      if (buffer[i][1] > box->lglobal[1]) buffer[i][1] -= box->lglobal[1];
      if (buffer[i][2] > box->lglobal[2]) buffer[i][2] -= box->lglobal[2];
      int cellx = floor(buffer[i][0] * rlcell[0] + TINY);
      int celly = floor(buffer[i][1] * rlcell[1] + TINY);
      int cellz = floor(buffer[i][2] * rlcell[2] + TINY);
      int celloff = get_cell_off(box, cellx, celly, cellz);
      cell_t *cell = box->cells + celloff;
      celldata_t *celldata = box->celldata + celloff;
      int curatom = cell->natoms;
      celldata->x[curatom][0] = buffer[i][0];
      celldata->x[curatom][1] = buffer[i][1];
      celldata->x[curatom][2] = buffer[i][2];
      celldata->v[curatom][0] = buffer[i][3];
      celldata->v[curatom][1] = buffer[i][4];
      celldata->v[curatom][2] = buffer[i][5];

      celldata->type[curatom] = 0;
      cell->natoms ++;
    }
  }
  esmd_free(buffer);
}

/* void print_atoms_x(esmd_t *md){ */
/*   box_t *box = &(md->box); */
/*   for (int i = -NCELL_CUT; i < box->nlocal[0] + NCELL_CUT; i ++){ */
/*     for (int j = -NCELL_CUT; j < box->nlocal[1] + NCELL_CUT; j ++){ */
/*       for (int k = -NCELL_CUT; k < box->nlocal[2] + NCELL_CUT; k ++){ */
/*         cell_t *cell = box->cells + get_cell_off(box, i, j, k); */
/*         celldata_t *celldata = box->celldata + get_cell_off(box, i, j, k); */
/*         printf("%d %d %d\n", i, j, k); */
/*         for (int p = 0; p < cell->natoms; p ++){ */
/*           printf("%d %d %d %f %f %f\n", i, j, k, celldata->x[p][0], celldata->x[p][1], celldata->x[p][2]); */
/*         } */
/*       } */
/*     } */
/*   } */
/* } */

void print_atoms_f(esmd_t *md){
  box_t *box = &(md->box);
  for (int k = 0; k < box->nlocal[2]; k ++){
    for (int j = 0; j < box->nlocal[1]; j ++){
      for (int i = 0; i < box->nlocal[0]; i ++){
        cell_t *cell = box->cells + get_cell_off(box, i, j, k);
        celldata_t *celldata = box->celldata + get_cell_off(box, i, j, k);
        //printf("%d %d %d\n", i, j, k);
        for (int p = 0; p < cell->natoms; p ++){
          printf("%g %g %g %g %g %g\n", celldata->x[p][0], celldata->x[p][1], celldata->x[p][2], celldata->f[p][0], celldata->f[p][1], celldata->f[p][2]);
          //printf("%g %g %g\n", celldata->f[p][0], celldata->f[p][1], celldata->f[p][2]);
        }
      }
    }
  }
}
