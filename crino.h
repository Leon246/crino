#ifndef CRINO_H
#define CRINO_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "inttypes.h"

typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64

namespace Ui {
class crino;
}

class crino {
public:
    int me,nprocs;
    int gflag,rflag;

    int nlocal,nbond,nangle,ndihedral;
    int ilocal,ibond,iangle,idihedral;
    int ntype,nbondtype,nangletype,ndihedraltype;
    int typemetal1,typemetal2,typeS,typeC,typeH;

    double xlo,xhi,ylo,yhi,zlo,zhi;
    double xlh,ylh,zlh;
    double xhalf,yhalf,zhalf;

    double samylo,samyhi,samzlo,samzhi;

    int *type,*tag,*mol,*q;
    int **bond,**angle,**dihedral;
    double **x,**v;

    double *theta,*phi,*chi;
    double ave_theta,ave_phi,ave_chi;

    int sam_units,nsam;
    double **sam;

    double iphi,itheta,igamma; // rotation angles
    double r[3][3]; // rotation matrix

    double percent_hydracarbon;
    double percent_alloy;

    crino();
    ~crino();

    void init(double *,int);
    void seed(int);
    void output(std::string);
    void multipleSAM(double);
    void pair_bond(int,int);
    void pair_angle(int,int);
    void pair_dihedral(int,int);
    void Au100();
    void Au111();
    void shiftyz();
    void dumptextdump(std::string);

private:
    int single_mol(int,int);
    static const int nflagmax = 20;

    void rotate(double,double,double);

    /* ----------------------------------------------------------------------
     create a 1d array
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n;
      if (nbytes == 0) return NULL;
      array = (TYPE *) malloc(nbytes);
      if (array == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array " << name << std::endl;
        exit(1);
      }
      return array;
    }

    template <typename TYPE>
    TYPE **create(TYPE **&array, int n, const char *name)
    {
      std::cout << "Cannot create/grow a vector/array of pointers for " << name << std::endl;
      exit(1);
    }

    /* ----------------------------------------------------------------------
     grow or shrink 1d array
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    TYPE *grow(TYPE *&array, int n, const char *name)
    {
      if (array == NULL) return create(array,n,name);
      bigint nbytes = sizeof(TYPE) * n;
      if (nbytes == 0) {
        destroy(array);
      }
      array = (TYPE *) realloc(array,nbytes);
      if (array == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array " << name << std::endl;
        exit(1);
      }
      return array;
    }

    template <typename TYPE>
    TYPE **grow(TYPE **&array, int n, const char *name)
    {
      std::cout << "Cannot create/grow a vector/array of pointers for " << name << std::endl;
      exit(1);
    }

    /* ----------------------------------------------------------------------
     destroy a 1d array
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    void destroy(TYPE *array)
    {
      if (array == NULL) return;
      free(array);
    }

    /* ----------------------------------------------------------------------
     create a 2d array
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n1*n2;
      if (nbytes == 0) return NULL;
      TYPE *data = (TYPE *) malloc(nbytes);
      if (data == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array 1" << name << std::endl;
        exit(1);
      }
      nbytes = sizeof(TYPE *) * n1;
      if (nbytes == 0) return NULL;
      array = (TYPE **) malloc(nbytes);
      if (data == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array 2" << name << std::endl;
        exit(1);
      }
      int n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }

    template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, const char *name)
    {
      std::cout << "Cannot create/grow a vector/array of pointers for " << name << std::endl;
      exit(1);
    }

    /* ----------------------------------------------------------------------
     grow or shrink 1st dim of a 2d array
     last dim must stay the same
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    TYPE **grow(TYPE **&array, int n1, int n2, const char *name)
    {
      if (array == NULL) return create(array,n1,n2,name);

      bigint nbytes = sizeof(TYPE) * n1*n2;
      if (nbytes == 0) {
        destroy(array);
        return NULL;
      }
      TYPE *data = (TYPE *) realloc(array[0],nbytes);
      if (data == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array " << name << std::endl;
        exit(1);
      }
      nbytes = sizeof(TYPE *) * n1;
      if (nbytes == 0) {
        destroy(array);
        return NULL;
      }
      array = (TYPE **) realloc(array,nbytes);
      if (data == NULL) {
        std::cout << "Failed to allocate " << nbytes << "bytes for array " << name << std::endl;
        exit(1);
      }
      int n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }

    template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, const char *name)
    {
      std::cout << "Cannot create/grow a vector/array of pointers for " << name << std::endl;
      exit(1);
    }


    /* ----------------------------------------------------------------------
     destroy a 2d array
     ------------------------------------------------------------------------- */

    template <typename TYPE>
    void destroy(TYPE **array)
    {
      if (array == NULL) return;
      free(array[0]);
      free(array);
    }
};

#endif // CRINO_H
