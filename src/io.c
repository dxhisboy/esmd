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
  double rlcell = box->rlcell;
  while ((count = read(fd, buffer, buffer_size)) > 0) {
    int natoms = count / (3 * sizeof(double));
    for (int i = 0; i < natoms; i ++){
      if (buffer[i][0] < 0) buffer[i][0] += box->lglobal[0];
      if (buffer[i][1] < 0) buffer[i][1] += box->lglobal[1];
      if (buffer[i][2] < 0) buffer[i][2] += box->lglobal[2];
      if (buffer[i][0] > box->lglobal[0]) buffer[i][0] -= box->lglobal[0];
      if (buffer[i][1] > box->lglobal[1]) buffer[i][1] -= box->lglobal[1];
      if (buffer[i][2] > box->lglobal[2]) buffer[i][2] -= box->lglobal[2];
      int cellx = floor(buffer[i][0] * rlcell + TINY);
      int celly = floor(buffer[i][1] * rlcell + TINY);
      int cellz = floor(buffer[i][2] * rlcell + TINY);
      int celloff = get_cell_off(box, cellx, celly, cellz);
      cell_t *cell = box->cells + celloff;
      celldata_t *celldata = box->celldata + celloff;
      int curatom = cell->natoms;
      celldata->x[curatom][0] = buffer[i][0];
      celldata->x[curatom][1] = buffer[i][1];
      celldata->x[curatom][2] = buffer[i][2];
      cell->natoms ++;
    }
  }
  esmd_free(buffer);
}
