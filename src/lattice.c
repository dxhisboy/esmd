#include <math.h>
#include <lattice.h>
#include <loops.h>
#include <multiproc.h>
#include <thermo.h>
#include <box.h>

//#define DEBUG_THIS_FILE
#include <log.h>

#define IA   16807
#define IM   2147483647
#define AM   (1.0/IM)
#define IQ   127773
#define IR   2836
#define MASK 123459876

double minimd_random(int *seed){
  int k = *seed / IQ;
  *seed = IA * (*seed - k * IQ) - IR * k;
  if (*seed < 0){
    *seed += IM;
  }
  return AM * *seed;
}

double next_rand(int *seed){
  for (int i = 0; i < 5; i ++)
    minimd_random(seed);
  return minimd_random(seed);
}
areal offset_fcc[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.5, 0.0},
  {0.5, 0.0, 0.5},
  {0.0, 0.5, 0.5}
};

int weight_fcc = 4;

lattice_t lattices[] = {
  //fcc
  {offset_fcc, 1, 1, 1, 4}
};

void esmd_lattice_scale(esmd_t *md, lattice_conf_t *conf){
  if (md->utype == UNIT_LJ){
    enum lattice_type type = conf->type;
    int *atom_types = conf->atom_types;
    ireal *masses = md->pair_conf.mass;
    
    lattice_t *lat = lattices + type;
    
    int weight = 0;
    if (atom_types){
      for (int i = 0; i < lat->natoms; i ++){
	weight += masses[atom_types[i]];
      }
    } else {
      weight = masses[0] * lat->natoms;
    }
    areal volume = lat->lx * lat->ly * lat->lz;
    areal scale = pow(weight / (conf->dens * volume), 1 / 3.0);
    conf->scale = scale;    
  }
}
void esmd_set_box_size_by_lattice(esmd_t *md){
  lattice_conf_t *conf = &(md->lat_conf);
  enum lattice_type type = conf->type;
  lattice_t *lat = lattices + type;

  md->box->lglobal[0] = conf->nx * conf->scale * lat->lx;
  md->box->lglobal[1] = conf->ny * conf->scale * lat->ly;
  md->box->lglobal[2] = conf->nz * conf->scale * lat->lz;
}

void esmd_create_atoms_by_lattice(esmd_t *md) {
  /* int lat_lo_x, lat_hi_x, lat_lo_y, lat_hi_y, lat_lo_z, lat_hi_z; */
  lattice_conf_t *conf = &(md->lat_conf);
  box_t *box = md->box;
  areal scale = conf->scale;
  enum lattice_type type = conf->type;
  lattice_t *lat = lattices + type;
  areal (*offset)[3] = lat->offset;
  areal *rlcell = box->rlcell;
  
  int lat_lo_x = (int)floor(box->olocal[0] / (scale * lat->lx));
  int lat_lo_y = (int)floor(box->olocal[1] / (scale * lat->ly));
  int lat_lo_z = (int)floor(box->olocal[2] / (scale * lat->lz));
  int lat_hi_x = (int)ceil((box->olocal[0] + box->llocal[0]) / (scale * lat->lx));
  int lat_hi_y = (int)ceil((box->olocal[1] + box->llocal[1]) / (scale * lat->ly));
  int lat_hi_z = (int)ceil((box->olocal[2] + box->llocal[2]) / (scale * lat->lz));

  areal vtot[3];
  vtot[0] = 0;
  vtot[1] = 0;
  vtot[2] = 0;
  debug("%f %f %d %d\n", box->olocal[0], box->olocal[0] + box->llocal[0], lat_lo_x, lat_hi_x);
  for (int kk = lat_lo_z; kk < lat_hi_z; kk ++){
    for (int jj = lat_lo_y; jj < lat_hi_y; jj ++){
      for (int ii = lat_lo_x; ii < lat_hi_x; ii ++){
	for (int io = 0; io < lat->natoms; io ++){
	  areal x[3], v[3];
	  x[0] = (ii * lat->lx + offset[io][0]) * scale;
	  x[1] = (jj * lat->ly + offset[io][1]) * scale;
	  x[2] = (kk * lat->lz + offset[io][2]) * scale;
	  if (x[0] >= box->olocal[0] && x[0] < box->olocal[0] + box->llocal[0] &&
	      x[1] >= box->olocal[1] && x[1] < box->olocal[1] + box->llocal[1] &&
	      x[2] >= box->olocal[2] && x[2] < box->olocal[2] + box->llocal[2]){
	    //this part follows the seeding strategy of minimd, to be changed
	    int ri = (int)round((x[0] / scale) * 2);
	    int rj = (int)round((x[1] / scale) * 2);
	    int rk = (int)round((x[2] / scale) * 2);
	    int seed = rk * (2 * conf->ny) * (2 * conf->nx) + rj * (2 * conf->nx) + ri + 1;
	    //int seed = ((kk * conf->ny + jj) * conf->nx + ii) * lat->natoms + io + 1;
	    v[0] = next_rand(&seed);
	    v[1] = next_rand(&seed);
	    v[2] = next_rand(&seed);
	    vtot[0] += v[0];
	    vtot[1] += v[1];
	    vtot[2] += v[2];
	  }
	  box_add_atom(box, x, v, 0, 0);
	}
      }
    }
  }
  areal vavg[3];
  esmd_global_sum_vec(md, vavg, vtot);
  int natoms = lat->natoms * conf->nx * conf->ny * conf->nz;
  vavg[0] = vavg[0] / natoms;
  vavg[1] = vavg[1] / natoms;
  vavg[2] = vavg[2] / natoms;
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	celldata->v[i][0] -= vavg[0];
	celldata->v[i][1] -= vavg[1];
	celldata->v[i][2] -= vavg[2];
	if (celldata->x[i][0] == 0 && celldata->x[i][1] == 0 && celldata->x[i][2] == 0){
	  debug("%f %f %f\n", celldata->v[i][0], celldata->v[i][1], celldata->v[i][2]);
	}
      }
    });
  md->natoms = natoms;
  debug("%f\n", temperature(md));
}
