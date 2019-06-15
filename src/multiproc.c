#include <data.h>
inline void part1d(int n, int np, int pid, int *start, int *count){
  int pncell = n / np;
  int rncell = n % np;
  if (pid < rncell) {
    *start = (pncell + 1) * pid;
    *count = pncell + 1;
  } else {
    *start = pncell * pid + rncell;
    *count = pncell;
  }
}
void esmd_multiproc_part_cart(esmd_t *md, int npx, int npy, int npz, int pid){
  box_t *box = &(md->box);
  int pidx = pid / (npy * npz);
  int pidy = (pid - pidx * npy * npz) / npz;
  int pidz = pid - (pidx * npy + pidy) * npz;
  int pcellx = box->nglobal[0] / npx;
  int pcelly = box->nglobal[1] / npy;
  int pcellz = box->nglobal[2] / npz;
  md->mpp.pidx = pidx;
  md->mpp.pidy = pidy;
  md->mpp.pidz = pidz;
  md->mpp.npx  = npx ;
  md->mpp.npy  = npy ;
  md->mpp.npz  = npz ;
  int stx, cntx, sty, cnty, stz, cntz;
  part1d(box->nglobal[0], npx, pidx, &stx, &cntx);
  part1d(box->nglobal[1], npy, pidy, &sty, &cnty);
  part1d(box->nglobal[2], npz, pidz, &stz, &cntz);
  /* box->llocal[0] = cntx * box->lcell; */
  /* box->llocal[1] = cnty * box->lcell; */
  /* box->llocal[2] = cntz * box->lcell; */
  /* box->lall[0] = (cntx + halo * 2) * box->lcell; */
  /* box->lall[1] = (cnty + halo * 2) * box->lcell; */
  /* box->lall[2] = (cntz + halo * 2) * box->lcell; */
  
  box->nall[0] = cntx + NCELL_CUT * 2;
  box->nall[1] = cnty + NCELL_CUT * 2;
  box->nall[2] = cntz + NCELL_CUT * 2;
  
  box->nlocal[0] = cntx;
  box->nlocal[1] = cnty;
  box->nlocal[2] = cntz;

  box->offset[0] = stx;
  box->offset[1] = sty;
  box->offset[2] = stz;
}

