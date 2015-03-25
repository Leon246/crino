#include "crino.h"

crino::crino()
{
    nbond = nangle = ndihedral = nlocal = 0;
    ibond = iangle = idihedral = ilocal = 0;
    nbondtype = nangletype = ndihedraltype = ntype = 0;
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0;
    type = tag = mol = NULL;
    bond = angle = dihedral = NULL;
    x = v = NULL;
    sam = NULL;
}

crino::~crino()
{
    destroy(x);
    destroy(v);
    destroy(tag);
    destroy(type);
}

void crino::init(double *box, int nflag)
{
    typemetal1 = 1;
    typemetal2 = 2;
    typeS = 3;
    typeC = 4;
    typeH = 5;

    sam_units = nflag;
    nsam = 0;
    ntype = 5;
    nbondtype = 3;
    nangletype = 5;
    ndihedraltype = 5;

    xlo = box[0];
    xhi = box[1];
    ylo = box[2];
    yhi = box[3];
    zlo = box[4];
    zhi = box[5];

    xlh = xhi - xlo; ylh = yhi - ylo; zlh = zhi - zlo;
    xhalf = xlh/2; yhalf = ylh/2; zhalf = zlh/2;

    // shift simulation box in y and z
    ylo = -yhalf;
    yhi = yhalf;
    zlo = -zhalf;
    zhi = zhalf;

    // specific restrictions to sam molecules
    double spacing = 5.0;
    int nysamlattice = 2*yhalf/(spacing*sqrt(3));
    samyhi = double(nysamlattice+1)*(spacing*sqrt(3)/2);
    samylo = -samyhi+spacing*sqrt(3)/2+0.5;
    int nzsamlattice = zhalf/(spacing);
    samzhi = double(nzsamlattice+1)*(spacing);
    samzlo = -samzhi+spacing/2;
}

/* ----------------------------------------------------------------------
 Create single molecule
 sflag == 0: hydrocarbon
 sflag == 1: alkan-thiol
 nunits: length of chain

 hydrocarbon:
 0..n-1 carbon atoms
 n,n+1:     hydrogen for carbon No. 1
 n+2,n+3:   hydrogen for carbon No. 2
 ...
 3n-3,3n-2: hydrogen for carbon No. n-1
 3n-1:      hydrogen for carbon No. n-1
 3n,3n+1,3n+2: hydrogen for carbon No. 0

 alkan-thiol
 0: sulfur atom
 1..n-1 carbon atoms
 n,n+1:     hydrogen for carbon No. 1
 n+2,n+3:   hydrogen for carbon No. 2
 ...
 3n-3,3n-2: hydrogen for carbon No. n-1
 3n-1:      hydrogen for carbon No. n-1

------------------------------------------------------------------------- */

int crino::single_mol(int sflag, int nflag)
{
//  double PI = atan(1.)*4.;
    double CC = 1.526;
    double SC = 1.810;
    double CH = 1.090;
    double CCC = 1.91113553;//1.9373;
    double CCH = 1.91113553;//1.9106;
    double HCH = 1.91113553;//1.8832;

    double cosCCC = cos(CCC);
    double sinCCC = sin(CCC);
    double cosCCH = cos(CCH);
    double sinCCH = sin(CCH);
    double sinHCHdiv2 = sin(HCH/2.);
    double cosHCHdiv2 = cos(HCH/2.);
    double costheta = -cosCCH/cos(HCH/2.);
    double sintheta = sqrt(1-costheta*costheta);
    double vec_x,vec_y;

    if (sflag != 0 && sflag != 1) {
//        std::cout << "Error setting sflag for single_mol" << std::endl;
        return 0;
    }

    int nunits = nflag;
    int nxs_sam = nunits*3+2;
    int nlocal_sam = 0;
    create(sam, nxs_sam, 4, "crino::sam");

    // STEP: 1
    // constructure base chain, xy plain
    // first atom
    sam[0][0] = sam[0][1] = sam[0][2] = 0.0;
    sam[0][3] = typeC;
    nlocal_sam++;
    // second atom
    sam[1][0] = CC;
    sam[1][1] = sam[1][2] = 0.0;
    sam[1][3] = typeC;
    nlocal_sam++;
    // the rest
    for (int i = 0; i < nunits-2; i++) {
      if (i%2 == 0) {
        vec_x = cosCCC*(sam[i][0]-sam[i+1][0])+sinCCC*(sam[i][1]-sam[i+1][1]);
        vec_y = cosCCC*(sam[i][1]-sam[i+1][1])-sinCCC*(sam[i][0]-sam[i+1][0]);
      } else {
        vec_x = cosCCC*(sam[i][0]-sam[i+1][0])-sinCCC*(sam[i][1]-sam[i+1][1]);
        vec_y = cosCCC*(sam[i][1]-sam[i+1][1])+sinCCC*(sam[i][0]-sam[i+1][0]);
      }
      sam[i+2][0] = sam[i+1][0]+vec_x;
      sam[i+2][1] = sam[i+1][1]+vec_y;
      sam[i+2][2] = 0.0;
      sam[i+2][3] = typeC;
      nlocal_sam++;
    }
    // STEP: 2
    // add hydrogen atoms for 2..n-1
    // CH2, from head to tail:
    for (int i = 0; i < nunits-1; i++) {
      if (i%2 == 0) {
        vec_x = costheta*(sam[i][0]-sam[i+1][0])+sintheta*(sam[i][1]-sam[i+1][1]);
        vec_y = costheta*(sam[i][1]-sam[i+1][1])-sintheta*(sam[i][0]-sam[i+1][0]);
      } else {
        vec_x = costheta*(sam[i][0]-sam[i+1][0])-sintheta*(sam[i][1]-sam[i+1][1]);
        vec_y = costheta*(sam[i][1]-sam[i+1][1])+sintheta*(sam[i][0]-sam[i+1][0]);
      }
      sam[nlocal_sam][0] = sam[i+1][0]-vec_x*CH*cosHCHdiv2/CC;
      sam[nlocal_sam][1] = sam[i+1][1]-vec_y*CH*cosHCHdiv2/CC;
      sam[nlocal_sam][2] = sam[i+1][2]+CH*sinHCHdiv2;
      sam[nlocal_sam][3] = typeH;
      nlocal_sam++;
      sam[nlocal_sam][0] = sam[nlocal_sam-1][0];
      sam[nlocal_sam][1] = sam[nlocal_sam-1][1];
      sam[nlocal_sam][2] = sam[i+1][2]-CH*sinHCHdiv2;
      sam[nlocal_sam][3] = typeH;
      nlocal_sam++;
      if (i == nunits-2) {
      if (i%2 == 0) {
        vec_x = cosCCH*(sam[i][0]-sam[i+1][0])+sinCCH*(sam[i][1]-sam[i+1][1]);
        vec_y = cosCCH*(sam[i][1]-sam[i+1][1])-sinCCH*(sam[i][0]-sam[i+1][0]);
      } else {
        vec_x = cosCCH*(sam[i][0]-sam[i+1][0])-sinCCH*(sam[i][1]-sam[i+1][1]);
        vec_y = cosCCH*(sam[i][1]-sam[i+1][1])+sinCCH*(sam[i][0]-sam[i+1][0]);
      }
      sam[nlocal_sam][0] = sam[i+1][0]+vec_x*CH/CC;
      sam[nlocal_sam][1] = sam[i+1][1]+vec_y*CH/CC;
      sam[nlocal_sam][2] = sam[i+1][2];
      sam[nlocal_sam][3] = typeH;
      nlocal_sam++;
    }
  }
  if (sflag == 0) {
    // STEP: 3
    // head hydrogen, bottom
    // first atom
    sam[nlocal_sam][0] = CH*cosCCH;
    sam[nlocal_sam][1] = -CH*sinCCH;
    sam[nlocal_sam][2] = 0.0;
    sam[nlocal_sam][3] = typeH;
    nlocal_sam++;
    // second atom
    vec_x = costheta*CC;
    sam[nlocal_sam][0] = -costheta*CH*cosHCHdiv2;
    sam[nlocal_sam][1] = sintheta*CH*cosHCHdiv2;
    sam[nlocal_sam][2] = CH*sinHCHdiv2;
    sam[nlocal_sam][3] = typeH;
    nlocal_sam++;
    // third atom
    sam[nlocal_sam][0] = -costheta*CH*cosHCHdiv2;
    sam[nlocal_sam][1] = sintheta*CH*cosHCHdiv2;
    sam[nlocal_sam][2] = -CH*sinHCHdiv2;
    sam[nlocal_sam][3] = typeH;
    nlocal_sam++;
  } else if (sflag == 1) {
    // extend sulfur
    sam[0][3] = typeS;
    double shift = SC-CC;
    for (int i = 1; i < nlocal_sam; i++) {
      sam[i][0] += shift;
    }
  }
  return nlocal_sam;
}

void crino::multipleSAM(double xlayer)
{
  double spacing = 5.;
  //double ratio;
  int m,n,maxm,maxn;
  double a[3],b[3],vec[3];
  double xinit,xinvert;

  xinit = abs(xlayer);
  xinvert = xinit == xlayer ? 1.0 : -1.0;

  a[0] = b[0] = 0.;
  a[2] = spacing;
  b[2] = spacing/2.;
  a[1] = 0.;
  b[1] = spacing*sqrt(3.)/2.;

  maxm = int(2.*(abs(ylo)+abs(yhi))/spacing);
  maxn = int(2.*(abs(zlo)+abs(zhi))/spacing);

  //rotate(0.,0.,-0.567239742);
  rotate(0.,0.,-0.602146327);

  nsam = 0;
  vec[0] = 0;

  for (m = -maxm; m <= maxn; m++) {
    for (n = -maxn; n <= maxn; n++) {
      //vec[0] = (double)m*a[0] + (double)n*b[0];
      vec[1] = (double)m*a[1] + (double)n*b[1];
      vec[2] = (double)m*a[2] + (double)n*b[2];
      //if (vec[0] > xhi || vec[0] < xlo) continue;
      if (vec[1] > samyhi || vec[1] < samylo) continue;
      if (vec[2] > samzhi || vec[2] < samzlo) continue;
      // check pbc

      nsam++;
      int sflag = 1, nflag = sam_units;
      int sam_atoms = single_mol(sflag,nflag);
      for (int i = 0; i < sam_atoms; i++) {
        if (ilocal == nlocal) {
          nlocal += 1000;
          grow(type, nlocal, "crino::type");
          grow(tag,  nlocal, "crino::tag");
          grow(mol,  nlocal, "crino::mol");
          grow(x, nlocal, 3, "crino::x");
          grow(v, nlocal, 3, "crino::v");
        }
        type[ilocal] = int(sam[i][3]);
        double xnew = r[0][0]*sam[i][0]+r[0][1]*sam[i][1]+r[0][2]*sam[i][2];
        double ynew = r[1][0]*sam[i][0]+r[1][1]*sam[i][1]+r[1][2]*sam[i][2];
        double znew = r[2][0]*sam[i][0]+r[2][1]*sam[i][1]+r[2][2]*sam[i][2];
        x[ilocal][0] = vec[0]+xnew*xinvert+xinit;
        x[ilocal][1] = vec[1]+ynew;
        x[ilocal][2] = vec[2]+znew;
        v[ilocal][0] = 0.0;
        v[ilocal][1] = 0.0;
        v[ilocal][2] = 0.0;
        mol[ilocal] = nsam;
        tag[ilocal] = ilocal+1;
        ilocal++;
      }
      pair_bond(sflag, nflag);
      pair_angle(sflag, nflag);
      pair_dihedral(sflag, nflag);
    }
  }
}

/* ----------------------------------------------------------------------
 // type 1: S-C
 // type 2: C-C
 // type 3: H-C
------------------------------------------------------------------------- */

void crino::pair_bond(int sflag, int nflag)
{
  int plocal = ilocal-nflag*3-2+sflag*3;
  int i,j;
  nbond += (nflag-sflag)*3+1;
  grow(bond, nbond, 3, "crino::bond");
  if (sflag == 1) {
    // S-C
    bond[ibond][0] = 1; // bond type
    bond[ibond][1] = tag[plocal];   // atom S
    bond[ibond][2] = tag[plocal+1]; // atom C
    ibond++;
  } else if (sflag == 0) {
    // C-C
    bond[ibond][0] = 2; // bond type
    bond[ibond][1] = tag[plocal];   // atom C
    bond[ibond][2] = tag[plocal+1]; // atom C
    ibond++;
    // H-C
    bond[ibond][0] = 3; // bond type
    bond[ibond][1] = tag[plocal];   // atom C
    bond[ibond][2] = tag[ilocal-1]; // atom H
    ibond++;
    // H-C
    bond[ibond][0] = 3; // bond type
    bond[ibond][1] = tag[plocal];   // atom C
    bond[ibond][2] = tag[ilocal-2]; // atom H
    ibond++;
    // H-C
    bond[ibond][0] = 3; // bond type
    bond[ibond][1] = tag[plocal];   // atom C
    bond[ibond][2] = tag[ilocal-3]; // atom H
    ibond++;
  }
  // C-C
  for (i = 1; i < nflag-1; i++) {
    bond[ibond][0] = 2; // bond type
    bond[ibond][1] = tag[plocal+i];   // atom C
    bond[ibond][2] = tag[plocal+i+1]; // atom C
    ibond++;
  }
  // H-C
  for (i = 1; i < nflag; i++) {
    bond[ibond][0] = 3; // bond type
    bond[ibond][1] = tag[plocal+i]; // atom C
    j = plocal+nflag+i*2-2;
    bond[ibond][2] = tag[j]; // atom H
    ibond++;
    bond[ibond][0] = 3; // bond type
    bond[ibond][1] = tag[plocal+i]; // atom C
    j = plocal+nflag+i*2-1;
    bond[ibond][2] = tag[j]; // atom H
    ibond++;
    if (i == nflag-1) {
      bond[ibond][0] = 3; // bond type
      bond[ibond][1] = tag[plocal+i]; // atom C
      j = plocal+nflag+i*2;
      bond[ibond][2] = tag[j]; // atom H
      ibond++;
    }
  }
  if (ibond != nbond) {
    //cout << "bond count error" << endl;
    exit(1);
  }
  return;
}

/* ----------------------------------------------------------------------
 // type 1: C-C-C
 // type 2: C-C-H, H-C-C
 // type 3: H-C-H
 // type 4: S-C-C, C-C-S
 // type 5: S-C-H, H-C-S
------------------------------------------------------------------------- */

void crino::pair_angle(int sflag, int nflag)
{
  int plocal = ilocal-nflag*3-2+sflag*3;
  int pbond = nbond-(nflag-sflag)*3-1;
  int pangle = iangle;
  int i,j,k,l,m,d;
  for (i = plocal; i < ilocal; i++) {
    // search for first pair
    for (j = plocal; j < ilocal; j++) {
      // avoid identical atoms
      if (j == i) continue;
      for (l = pbond; l < nbond; l++) {
        // found first pair
        if ((tag[i] == bond[l][1] && tag[j] == bond[l][2]) || \
            (tag[i] == bond[l][2] && tag[j] == bond[l][1])) {
          // search for second pair
          for (k = plocal; k < ilocal; k++) {
            // avoid identical atoms
            if (k == i || k == j) continue;
            for (m = pbond; m < nbond; m++) {
              // found second pair
              if ((tag[j] == bond[m][1] && tag[k] == bond[m][2]) || \
                  (tag[j] == bond[m][2] && tag[k] == bond[m][1])) {
                // grow angle array if necessary
                if (iangle == nangle) {
                  nangle += 10000;
                  grow(angle, nangle, 4, "crino::angle");
                }
                // determine angle type
                if (type[i] == typeC && type[j] == typeC && type[k] == typeC) {
                  // type 1: C-C-C
                  angle[iangle][0] = 1;
                } else if ((type[i] == typeC && type[j] == typeC && type[k] == typeH) || \
                           (type[i] == typeH && type[j] == typeC && type[k] == typeC)) {
                  // type 2: C-C-H, H-C-C
                  angle[iangle][0] = 2;
                } else if (type[i] == typeH && type[j] == typeC && type[k] == typeH) {
                  // type 3: H-C-H
                  angle[iangle][0] = 3;
                } else if ((type[i] == typeS && type[j] == typeC && type[k] == typeC) || \
                           (type[i] == typeC && type[j] == typeC && type[k] == typeS)) {
                  // type 4: S-C-C, C-C-S
                  angle[iangle][0] = 4;
                } else if ((type[i] == typeS && type[j] == typeC && type[k] == typeH) || \
                           (type[i] == typeH && type[j] == typeC && type[k] == typeS)) {
                  // type 5: S-C-H, H-C-S
                  angle[iangle][0] = 5;
                } else {
                  // debug
                  //cout << "BUG: cannot detect angle type" << endl;
                  exit(1);
                }
                angle[iangle][1] = tag[i];
                angle[iangle][2] = tag[j];
                angle[iangle][3] = tag[k];
                // avoid duplicate angles
                for (d = pangle; d < iangle; d++) {
                  if (angle[d][0] == angle[iangle][0] && \
                      angle[d][1] == angle[iangle][3] && \
                      angle[d][3] == angle[iangle][1] && \
                      angle[d][2] == angle[iangle][2]) {
                    iangle--;
                    break;
                  }
                }
                iangle++;
              }
            }
          }
        }
      }
    }
  }
  return;
}

/* ----------------------------------------------------------------------
 // type 1 C-C-C-H, H-C-C-C
 // type 2 H-C-C-H
 // type 3 C-C-C-C
 // type 4 S-C-C-C, C-C-C-S
 // type 5 S-C-C-H, H-C-C-S
------------------------------------------------------------------------- */

void crino::pair_dihedral(int sflag, int nflag)
{
  int plocal = ilocal-nflag*3-2+sflag*3;
  int pbond = ibond-(nflag-sflag)*3-1;
  int pdihedral = idihedral;
  int i,j,k,l,o,p,q,d;
  for (i = plocal; i < ilocal; i++) {
    // search for the first pair
    for (j = plocal; j < ilocal; j++) {
      // avoid identical atoms
      if (i == j) continue;
      for (o = pbond; o < ibond; o++) {
        // found first pair
        if ((tag[i] == bond[o][1] && tag[j] == bond[o][2]) || \
            (tag[i] == bond[o][2] && tag[j] == bond[o][1])) {
          // search for second pair
          for (k = plocal; k < ilocal; k++) {
            // avoid identical atoms
            if (k == i || k == j) continue;
            for (p = pbond; p < ibond; p++) {
              // found second pair
              if ((tag[j] == bond[p][1] && tag[k] == bond[p][2]) || \
                  (tag[j] == bond[p][2] && tag[k] == bond[p][1])) {
                // search for third pair
                for (l = plocal; l < ilocal; l++) {
                  // avoid identical atoms
                  if (l == i || l == j || l == k) continue;
                  for (q = pbond; q < ibond; q++) {
                    // found third pair
                    if ((tag[k] == bond[q][1] && tag[l] == bond[q][2]) || \
                        (tag[k] == bond[q][2] && tag[l] == bond[q][1])) {
                      // grow dihedral if necessary
                      if (idihedral == ndihedral) {
                        ndihedral += 10000;
                        grow(dihedral, ndihedral, 5, "crino::dihedral");
                      }
                      // determine dihedral type
                      if ((type[i] == typeC && type[j] == typeC && type[k] == typeC && type[l] == typeH) || \
                          (type[i] == typeH && type[j] == typeC && type[k] == typeC && type[l] == typeC)) {
                        // type 1: C-C-C-H, H-C-C-C
                        dihedral[idihedral][0] = 1;
                      } else if (type[i] == typeH && type[j] == typeC && type[k] == typeC && type[l] == typeH) {
                        // type 2: H-C-C-H
                        dihedral[idihedral][0] = 2;
                      } else if (type[i] == typeC && type[j] == typeC && type[k] == typeC && type[l] == typeC) {
                        // type 3: C-C-C-C
                        dihedral[idihedral][0] = 3;
                      } else if ((type[i] == typeS && type[j] == typeC && type[k] == typeC && type[l] == typeC) || \
                                 (type[i] == typeC && type[j] == typeC && type[k] == typeC && type[l] == typeS)) {
                        // type 4: S-C-C-C, C-C-C-S
                        dihedral[idihedral][0] = 4;
                      } else if ((type[i] == typeS && type[j] == typeC && type[k] == typeC && type[l] == typeH) || \
                                 (type[i] == typeH && type[j] == typeC && type[k] == typeC && type[l] == typeS)) {
                        // type 5: S-C-C-H, H-C-C-S
                        dihedral[idihedral][0] = 5;
                      } else {
                        // debug
                        //cout << "BUG: cannot detect dihedral type" << endl;
                        exit(1);
                      }
                      dihedral[idihedral][1] = tag[i];
                      dihedral[idihedral][2] = tag[j];
                      dihedral[idihedral][3] = tag[k];
                      dihedral[idihedral][4] = tag[l];
                      // avoid duplicate dihedral angle
                      for (d = pdihedral; d < idihedral; d++) {
                        if (dihedral[d][0] == dihedral[idihedral][0] && \
                            dihedral[d][1] == dihedral[idihedral][4] && \
                            dihedral[d][4] == dihedral[idihedral][1] && \
                            dihedral[d][2] == dihedral[idihedral][3] && \
                            dihedral[d][3] == dihedral[idihedral][2]) {
                          idihedral--;
                          break;
                        }
                      }
                      idihedral++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

/* ----------------------------------------------------------------------
 create Au surface with (100) plane, fcc lattice
 ------------------------------------------------------------------------- */

void crino::Au100()
{
  int i = 0;
  double alat = 4.08;
  double basis[3][4],vec[3];
  double shift = -1.0;
  int l,m,n,maxl,maxm,maxn;

  maxl = int(10.*(abs(xlo)+abs(xhi))/alat);
  maxm = int(10.*(abs(ylo)+abs(yhi))/alat);
  maxn = int(10.*(abs(zlo)+abs(zhi))/alat);

  // initialize lattice vector
  basis[0][0] = basis[1][0] = basis[2][0] = 0.;
  basis[1][1] = basis[2][1] = 0.5;
  basis[0][2] = basis[2][2] = 0.5;
  basis[0][3] = basis[1][3] = 0.5;
  basis[0][1] = basis[1][2] = basis[2][3] = 0.;

  for (l = -maxl; l <= maxl; l++) {
    for (m = -maxm; m <= maxm; m++) {
      for (n = -maxn; n <= maxn; n++) {
        for (i = 0; i < 4; i++) {

          vec[0] = (basis[0][i]+(double)l)*alat;
          vec[1] = (basis[1][i]+(double)m)*alat;
          vec[2] = (basis[2][i]+(double)n)*alat;

          if (vec[0] >= shift || vec[0] <= xlo) continue;
          if (vec[1] >= yhi || vec[1] <= ylo) continue;
          if (vec[2] >= zhi || vec[2] <= zlo) continue;

          if (nlocal == nlocal) {
            nlocal += 1000;
            grow(type, nlocal, "crino::type");
            grow(tag,  nlocal, "crino::tag");
            grow(mol,  nlocal, "crino::mol");
            grow(x, nlocal, 3, "crino::x");
            grow(v, nlocal, 3, "crino::v");
          }
          type[ilocal] = typemetal1; // Au
          x[ilocal][0] = vec[0];
          x[ilocal][1] = vec[1];
          x[ilocal][2] = vec[2];
          v[ilocal][0] = 0.0;
          v[ilocal][1] = 0.0;
          v[ilocal][2] = 0.0;
          mol[ilocal] = 0;
          tag[ilocal] = ilocal+1;
          ilocal++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
 create Au surface with (111) plane, fcc lattice
 ------------------------------------------------------------------------- */

void crino::Au111()
{
  int i = 0;
  double PI = atan(1.)*4.;
  double alat = 4.08;
  double basis[4][3],vec[3],vec0[3];
  double cut = -0.5, shiftx = 0.0;
  int l,m,n,maxl,maxm,maxn;

  maxl = int(5.*(abs(xlo)+abs(xhi))/alat);
  maxm = int(5.*(abs(ylo)+abs(yhi))/alat);
  maxn = int(5.*(abs(zlo)+abs(zhi))/alat);

  // initialize lattice vector
  basis[0][0] = basis[0][1] = basis[0][2] = 0.;
  basis[1][0] = basis[2][1] = basis[3][2] = 0.;
  basis[1][1] = basis[1][2] = 0.5;
  basis[2][0] = basis[2][2] = 0.5;
  basis[3][0] = basis[3][1] = 0.5;

  rotate(0.,PI/4.,atan(sqrt(2.)/2.));

  for (l = -maxl; l <= maxl; l++) {
    for (m = -maxm; m <= maxm; m++) {
      for (n = -maxn; n <= maxn; n++) {
        for (i = 0; i < 4; i++) {
          vec0[0] = (basis[i][0]+(double)l)*alat;
          vec0[1] = (basis[i][1]+(double)m)*alat;
          vec0[2] = (basis[i][2]+(double)n)*alat;
          vec[0] = r[0][0]*vec0[0]+r[0][1]*vec0[1]+r[0][2]*vec0[2];
          vec[1] = r[1][0]*vec0[0]+r[1][1]*vec0[1]+r[1][2]*vec0[2];
          vec[2] = r[2][0]*vec0[0]+r[2][1]*vec0[1]+r[2][2]*vec0[2];
          if (vec[0] >= cut || vec[0] <= xlo) continue;
          if (vec[2] >= yhi || vec[2] <= ylo) continue;
          if (vec[1] >= zhi || vec[1] <= zlo) continue;
          if (ilocal == nlocal) {
            nlocal += 1000;
            grow(type, nlocal, "crino::type");
            grow(tag,  nlocal, "crino::tag");
            grow(mol,  nlocal, "crino::mol");
            grow(x, nlocal, 3, "crino::x");
            grow(v, nlocal, 3, "crino::v");
          }
          type[ilocal] = typemetal1; // Au
          x[ilocal][0] = vec[0]+shiftx;
          x[ilocal][1] = vec[2];
          x[ilocal][2] = vec[1];
          v[ilocal][0] = 0.0;
          v[ilocal][1] = 0.0;
          v[ilocal][2] = 0.0;
          mol[ilocal] = 0;
          tag[ilocal] = ilocal+1;
          ilocal++;
        }
      }
    }
  }

  return;
}

void crino::shiftyz()
{
  for (int i = 0; i < ilocal; i++) {
    x[i][1] += yhalf;
    x[i][2] += zhalf;
  }
  ylo += yhalf;
  yhi += yhalf;
  zlo += zhalf;
  zhi += zhalf;
}


void crino::output(std::string fileDir)
{
  std::fstream file;
  std::string fileName = fileDir + "data.crino";

  file.open(fileName.c_str(),std::fstream::out);
  if (file.is_open()) { std::cout << "Creating output file..." << std::endl; }
  else { std::cout << "Failed to creat file " << fileName << std::endl; }

  file << "crino generated" << std::endl;
  file << std::endl;

  file.precision(8);
  file.setf(std::ios::right,std::ios::adjustfield);
  file.setf(std::ios::dec,std::ios::basefield);
  file.setf(std::ios::fixed,std::ios::floatfield);
  file.setf(std::ios_base::showbase | std::ios_base::uppercase);


  file.width(10); file << std::right << ilocal << '\t' << "atoms" << std::endl;
  file.width(10); file << std::right << ibond  << '\t' << "bonds" << std::endl;
  file.width(10); file << std::right << iangle << '\t' << "angles" << std::endl;
  file.width(10); file << std::right << idihedral << '\t' << "dihedrals" << std::endl;
  file << std::endl;

  file.width(10); file << std::right << ntype << '\t' << "atom types" << std::endl;

  file.width(10); file << std::right << nbondtype << '\t' << "bond types" << std::endl;
  file.width(10); file << std::right << nangletype << '\t' << "angle types" << std::endl;
  file.width(10); file << std::right << ndihedraltype << '\t' << "dihedral types" << std::endl;

  file << std::endl;

  file.width(17); file << xlo;
  file.width(17); file << xhi << " xlo xhi" << std::endl;
  file.width(17); file << ylo;
  file.width(17); file << yhi << " ylo yhi" << std::endl;
  file.width(17); file << zlo;
  file.width(17); file << zhi << " zlo zhi" << std::endl;
  file << std::endl;

  file << "Masses" << std::endl << std::endl;
  file << "1 196.97" << std::endl;
  file << "2 58.69" << std::endl;
  file << "3 32.065" << std::endl;
  file << "4 12.0107" << std::endl;
  file << "5 1.00794" << std::endl;

  file << std::endl;

  file << "Bond Coeffs" << std::endl << std::endl;
  file << "1 10.27632 1.810" << std::endl;  // type 1: S-C
  file << "2 13.44160 1.526" << std::endl;  // type 2: C-C
  file << "3 14.74240 1.090" << std::endl; // type 3: C-H
  file << std::endl;
  file << "Angle Coeffs" << std::endl << std::endl;
  file << "1 1.73440 109.50" << std::endl; // type 1: C-C-C
  file << "2 2.16800 109.50" << std::endl; // type 2: C-C-H
  file << "3 1.51760 109.50" << std::endl; // type 3: H-C-H
  file << "4 2.16800 108.60" << std::endl; // type 4: S-C-C
  file << "5 2.16800 109.50" << std::endl; // type 5: S-C-H
  file << std::endl;
  file << "Dihedral Coeffs" << std::endl << std::endl;
  file << "1 0.060704 1 3" << std::endl; // type 1 C-C-C-H
  file << "2 0.060704 1 3" << std::endl; // type 2 H-C-C-H
  file << "3 0.060704 1 3" << std::endl; // type 3 C-C-C-C
  file << "4 0.060704 1 3" << std::endl; // type 4 S-C-C-C
  file << "5 0.060704 1 3" << std::endl; // type 5 S-C-C-H
  file << std::endl;

  file << "Atoms" << std::endl << std::endl;

  for (int i = 0; i < ilocal; i++) {

    file.width(10); file << tag[i];
    file.width(6);
    file << mol[i];
    file.width(4);  file << type[i];
    file.width(17); file << x[i][0];
    file.width(17); file << x[i][1];
    file.width(17); file << x[i][2];
    file << std::endl;
  }

  file << std::endl;


    file << "Bonds" << std::endl << std::endl;

    for (int i = 0; i < ibond; i++) {
      file.width(10); file << i+1;         // tag
      file.width(4);  file << bond[i][0]; // bond type
      file.width(10); file << bond[i][1]; // x 1
      file.width(10); file << bond[i][2]; // x 2
      file << std::endl;
    }
    file << std::endl;

    file << "Angles" << std::endl << std::endl;

    for (int i = 0; i < iangle; i++) {
      file.width(10); file << i+1;          // tag
      file.width(4);  file << angle[i][0]; // angle type
      file.width(10); file << angle[i][1]; // x 1
      file.width(10); file << angle[i][2]; // x 2
      file.width(10); file << angle[i][3]; // x 3
      file << std::endl;
    }
    file << std::endl;

    file << "Dihedrals" << std::endl << std::endl;

    for (int i = 0; i < idihedral; i++) {
      file.width(10); file << i+1;             // tag
      file.width(4);  file << dihedral[i][0]; // dihedral type
      file.width(10); file << dihedral[i][1]; // x 1
      file.width(10); file << dihedral[i][2]; // x 2
      file.width(10); file << dihedral[i][3]; // x 3
      file.width(10); file << dihedral[i][4]; // x 4
      file << std::endl;
    }
    file << std::endl;

  /*
  if (gflag == 1) {
    file << "Velocities" << endl << endl;
    for (int i = 0; i < ilocal; i++) {
      file.width(10); file << tag[i];
      file.width(17); file << v[i][0];
      file.width(17); file << v[i][1];
      file.width(17); file << v[i][2];
      file << endl;
    }
  }*/
  file.close();
}

void crino::dumptextdump(std::string fileDir)
{
  std::fstream outfile;
  std::string fileName = fileDir + "dump.0";

  outfile.open(fileName.c_str(),std::fstream::out);
  if (outfile.is_open()) { std::cout << "Creating LAMMPS dump file for preview" << std::endl; }
  else { std::cout << "Failed to creat file " << fileName << std::endl; }

  outfile << "ITEM: TIMESTEP" << std::endl;
  outfile << "0" << std::endl;
  outfile << "ITEM: NUMBER OF ATOMS" << std::endl;
  outfile << ilocal << std::endl;
  outfile << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
  outfile << xlo << " " << xhi << std::endl;
  outfile << ylo << " " << yhi << std::endl;
  outfile << zlo << " " << zhi << std::endl;
    outfile << "ITEM: ATOMS id type x y z" << std::endl;
    for (int i = 0; i < ilocal; i++) {
      outfile << tag[i] << " " << type[i] << " ";
      outfile << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
    }

  outfile.close();
}

/* ----------------------------------------------------------------------
 rotation matrix
 Euler Angles: phi theta psi
             [ 1    0         0     ]
 Rx[phi]   = [ 0 cos(phi) -sin(phi) ]
             [ 0 sin(phi)  cos(phi) ]
             [ cos(theta) 0 sin(theta) ]
 Ry[theta] = [      0     1      0     ]
             [-sin(theta) 0 cos(theta) ]
             [ cos(psi) -sin(psi) 0 ]
 Rz[psi]   = [ sin(psi)  cos(psi) 0 ]
             [    0         0     1 ]
 R = Rz(psi)Ry(theta)Rx(phi)
------------------------------------------------------------------------- */

void crino::rotate(double phi, double theta, double psi)
{
  r[0][0] = cos(theta)*cos(psi);
  r[0][1] = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
  r[0][2] = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
  r[1][0] = cos(theta)*sin(psi);
  r[1][1] = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
  r[1][2] = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
  r[2][0] = -sin(theta);
  r[2][1] = sin(phi)*cos(theta);
  r[2][2] = cos(phi)*cos(theta);
}
