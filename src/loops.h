#define EXTRACT_CELLDATA(celldata, sf)		\
  double (*v ## sf)[3] = celldata->v;		\
  double (*x ## sf)[3] = celldata->x;		\
  double (*f ## sf)[3] = celldata->f;		\
  double *q  ## sf     = celldata->q;		\
  int    *export ## sf = celldata->export;	\
  int    *type   ## sf = celldata->type;	\

#define ESMD_CELL_ITER(box, callback, ...)			\
  for (int kk = 0; kk < box->nlocal[2]; kk ++)			\
    for (int jj = 0; jj < box->nlocal[1]; jj ++)		\
      for (int ii = 0; ii < box->nlocal[0]; ii ++){		\
	int celloff = get_cell_off(box, ii, jj, kk);		\
	cell_t *cell = box->cells + celloff;			\
	celldata_t *celldata = box->celldata + celloff;		\
	EXTRACT_CELLDATA(celldata, __VA_ARGS__);		\
	{							\
	  callback;						\
	}							\
      }
