#include "samana.h"
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>
#include <omp.h>

#define PI 3.14159265359

samana::samana()
{
    alllist.reset();
    samlist.reset();
    nlocal = nsam = 0;
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0;
    fileid = "NULL.txt";
}

samana::~samana()
{
    alllist.release();
    samlist.release();
}

void samana::init(int nflag)
{
    nunits = nflag;
    if (nunits > 32) {
        std::cout << "ERROR: " << "sam units cannot exceed 32" << std::endl;
        exit(1);
    }
    if (nunits%2 == 0) {
        std::cout << "ERROR: " << "sam units must be odd numbered" << std::endl;
        exit(1);
    }
}

int samana::readfile(std::string filename, int vflag)
{
    int i, j, k;
    char line[256], *pch;
    std::fstream dumpfile;
    fileid = filename;
    dumpfile.open(filename.c_str(), std::fstream::in);
    if (dumpfile.is_open()) {
        std::cout << "Reading file " << filename << " ..." << std::endl;
        alllist.reset();
        samlist.reset();
    }
    else {
        std::cout << "File " << filename << " not found" << std::endl;
        return 0;
    }
    i = 0;
    j = 0;
    nsam = 0;
    while (dumpfile.good()) {
        dumpfile.getline(line, 256, '\n');
        i++;
        if (i == 2) {
            mdstep = atoi(line);
        }
        else if (i == 4) {
            nlocal = atoi(line);
            alllist = std::unique_ptr<atom []>(new atom[nlocal]);
        }
        else if (i == 6) {
            pch = strtok(line, " \t\n");
            xlo = atof(pch);
            pch = strtok(NULL, " \t\n");
            xhi = atof(pch);
        }
        else if (i == 7) {
            pch = strtok(line, " \t\n");
            ylo = atof(pch);
            pch = strtok(NULL, " \t\n");
            yhi = atof(pch);
        }
        else if (i == 8) {
            pch = strtok(line, " \t\n");
            zlo = atof(pch);
            pch = strtok(NULL, " \t\n");
            zhi = atof(pch);
        }
        else if (i > 9) {
            pch = strtok(line, " \t\n");
            alllist[j].id = atoi(pch);
            if (vflag > 0) {
                pch = strtok(NULL, " \t\n");
                k = atoi(pch);
                alllist[j].mid = k;
                nsam = nsam > k ? nsam : k;
                pch = strtok(NULL, " \t\n");
                alllist[j].type = atoi(pch);
            }
            else if (vflag == 0) {
                pch = strtok(NULL, " \t\n");
                alllist[j].type = atoi(pch);
                if (alllist[j].type == 3) nsam++;
                alllist[j].mid = alllist[j].type < 3 ? 0 : nsam;
            }
            // x y z
            pch = strtok(NULL, " \t\n");
            alllist[j].x = atof(pch);
            pch = strtok(NULL, " \t\n");
            alllist[j].y = atof(pch);
            pch = strtok(NULL, " \t\n");
            alllist[j].z = atof(pch);
            // vx vy vz
            pch = strtok(NULL, " \t\n");
            alllist[j].vx = atof(pch);
            pch = strtok(NULL, " \t\n");
            alllist[j].vy = atof(pch);
            pch = strtok(NULL, " \t\n");
            alllist[j].vz = atof(pch);
            if (vflag == 2) {
                // pe/atom
                pch = strtok(NULL, " \t\n");
                alllist[j].pea = atof(pch);
            }
            j++;
            if (j == nlocal){
                break;
            }
        }
    }
    dumpfile.close();
    return 1;
}

void samana::assemblesam()
{
    int i, j, k, itmp;
    double dtmp;
    std::unique_ptr<int[]> jj(new int[nsam]);
    samlist = std::unique_ptr<sam[]>(new sam[nsam]);
    std::cout << "Reassembling SAM molecules... " << std::endl;
    for (j = 0; j < nsam; j++) {
        jj[j] = 0;
    }
    for (i = 0; i < nlocal; i++) {
        k = alllist[i].mid;
        if (k > 0) {
            if (alllist[i].type != 5) {
                j = jj[k - 1];
                samlist[k - 1].mol[j] = alllist[i].mid;
                samlist[k - 1].type[j] = alllist[i].type;
                samlist[k - 1].tag[j] = alllist[i].id;
                samlist[k - 1].x[j] = alllist[i].x;
                samlist[k - 1].y[j] = alllist[i].y;
                samlist[k - 1].z[j] = alllist[i].z;
                samlist[k - 1].pe[j] = alllist[i].pea;
                jj[k - 1]++;
            }
        }
    }
    // bubble sort
    for (j = 0; j < nsam; j++) {
        for (k = 1; k < nunits; k++) {
            for (i = 0; i < k; i++) {
                if (samlist[j].tag[k] < samlist[j].tag[i]) {
                    itmp = samlist[j].tag[k];
                    samlist[j].tag[k] = samlist[j].tag[i];
                    samlist[j].tag[i] = itmp;
                    itmp = samlist[j].type[k];
                    samlist[j].type[k] = samlist[j].type[i];
                    samlist[j].type[i] = itmp;
                    itmp = samlist[j].mol[k];
                    samlist[j].mol[k] = samlist[j].mol[i];
                    samlist[j].mol[i] = itmp;
                    dtmp = samlist[j].x[k];
                    samlist[j].x[k] = samlist[j].x[i];
                    samlist[j].x[i] = dtmp;
                    dtmp = samlist[j].y[k];
                    samlist[j].y[k] = samlist[j].y[i];
                    samlist[j].y[i] = dtmp;
                    dtmp = samlist[j].z[k];
                    samlist[j].z[k] = samlist[j].z[i];
                    samlist[j].z[i] = dtmp;
                }
            }
        }
    }
    // sum pe/atom to pe/mol
    for (j = 0; j < nsam; j++) {
        samlist[j].totalpe = 0;
        for (k = 0; k < 17; k++) {
            samlist[j].totalpe += samlist[j].pe[k];
        }
    }
    // pbc check
    double ylh = yhi - ylo;
    double zlh = zhi - zlo;
    for (j = 0; j < nsam; j++) {
        for (k = 1; k < nunits; k++) {
            if (samlist[j].y[k] - samlist[j].y[0] > 0.5*ylh) {
                samlist[j].y[k] -= ylh;
            }
            else if (samlist[j].y[0] - samlist[j].y[k] > 0.5*ylh) {
                samlist[j].y[k] += ylh;
            }
            if (samlist[j].z[k] - samlist[j].z[0] > 0.5*zlh) {
                samlist[j].z[k] -= zlh;
            }
            else if (samlist[j].z[0] - samlist[j].z[k] > 0.5*zlh) {
                samlist[j].z[k] += zlh;
            }
        }
    }
}

void samana::compute()
{
    int i, j;
    double tmp;
    double a, b, c;
    double x0, y0, z0;
    double D1, D2, D3, twistangle;
    std::cout << "Computing angles..." << std::endl;
    //std::ofstream testfile;
    //testfile.open("testfile.txt", std::ofstream::out);
    for (i = 0; i < nsam; i++) {
        // tilt angle and tilt direction
        // determined by S and even number C atoms
        samlist[i].t[0] = 0;
        samlist[i].t[1] = 0;
        samlist[i].t[2] = 0;
        for (j = 1; j < floor(nunits/2); j++) {
            a = samlist[i].x[(j + 1) * 2] - samlist[i].x[2];
            b = samlist[i].y[(j + 1) * 2] - samlist[i].y[2];
            c = samlist[i].z[(j + 1) * 2] - samlist[i].z[2];
            tmp = a * a + b * b + c * c;
            tmp = sqrt(tmp);
            a /= tmp;
            b /= tmp;
            c /= tmp;
            samlist[i].t[0] += a;
            samlist[i].t[1] += b;
            samlist[i].t[2] += c;
        }
        a = samlist[i].t[0];
        b = samlist[i].t[1];
        c = samlist[i].t[2];
        tmp = a * a + b * b + c * c;
        tmp = sqrt(tmp);
        a /= tmp;
        b /= tmp;
        c /= tmp;
        samlist[i].t[0] = a;
        samlist[i].t[1] = b;
        samlist[i].t[2] = c;
        // tilt angle
        samlist[i].tilt = atan(sqrt(b*b + c*c) / a)*180/PI;
        // tilt direction
        if (b == 0) {
            samlist[i].tiltdir = c > 0 ? 90 : -90.0;
        } else {
            samlist[i].tiltdir = atan(c / b) * 180 / PI;
        }
        if (c > 0 && b < 0) {
            samlist[i].tiltdir = 180.0 + samlist[i].tiltdir;
        } else if (c < 0 && b < 0) {
            samlist[i].tiltdir = -180.0 + samlist[i].tiltdir;
        }

        // twist angle
        tmp = 0;
        for (j = 1; j < floor(nunits/2); j++) {
            x0 = samlist[i].x[2 * j + 1] - samlist[i].x[2*j];
            y0 = samlist[i].y[2 * j + 1] - samlist[i].y[2*j];
            z0 = samlist[i].z[2 * j + 1] - samlist[i].z[2*j];
            // point to plane
            D1 = (c*y0 - b*z0) / sqrt(c*c + b*b);
            // point to line
            D2 = sqrt(a*a*(z0*z0 + y0*y0) + b*b*(x0*x0 + z0*z0) + c*c*(x0*x0 + y0*y0)
                - 2 * a*b*x0*y0 - 2 * a*c*x0*z0 - 2 * b*c*y0*z0) / sqrt(a*a + b*b + c*c);
            // point to orthgnal plane ((bb+cc)x-aby-acz=0)
            D3 = ((b*b + c*c)*x0 - a*b*y0 - a*c*z0) / sqrt((b*b + c*c)*(b*b + c*c) + a*a*b*b + a*a*c*c);
            twistangle = asin(D1 / D2) * 180 / PI;
            if (D1 > 0 && D3 > 0) {
                twistangle = 180 - twistangle;
            }
            else if (D1 < 0 && D3 > 0) {
                twistangle = -180 - twistangle;
            }
            tmp += twistangle;
            //testfile << i << " " << j << " | ";
            //testfile << a << " " << b << " " << c << " | ";
            //testfile << x0 << " " << y0 << " " << z0 << " | ";
            //testfile << D1 << " " << D2 << " " << D3 << " " << twistangle << std::endl;
        }
        samlist[i].twist = tmp / (floor(nunits/2)-1);
    }
    //testfile.close();
}

void samana::output()
{
    int i;
    std::ofstream output;
    std::string filename;
    filename = fileid + ".output";
    output.open(filename.c_str(), std::ofstream::out);
    for (i = 0; i < nsam; i++) {
        output << samlist[i].x[0] << " ";
        output << samlist[i].y[0] << " ";
        output << samlist[i].z[0] << " ";
        output << samlist[i].t[0] << " ";
        output << samlist[i].t[1] << " ";
        output << samlist[i].t[2] << " ";
        output << samlist[i].tilt << " ";
        output << samlist[i].tiltdir << " ";
        output << samlist[i].twist << " ";
        output << samlist[i].totalpe << std::endl;
    }
    output.close();
}

void samana::dumpout()
{
    int i;
    unsigned int j;
    const unsigned int buffersize = 8192;
    char buffer[buffersize];
    char line[256];
    std::ofstream dumpout;
    std::string filename;
    filename = fileid + ".pem";
    dumpout.open(filename.c_str(), std::ofstream::out);
    // empty buffer
    for (j = 0; j < buffersize; j++) buffer[j] = '\0';
    strcpy(buffer,"ITEM: TIMESTEP\n");
    sprintf(line,"%d\n",mdstep);
    strcat(buffer,line);
    strcat(buffer,"ITEM: NUMBER OF ATOMS\n");
    sprintf(line,"%d\n",nlocal);
    strcat(buffer,line);
    strcat(buffer,"ITEM: BOX BOUNDS pp pp pp\n");
    sprintf(line,"%f %f\n",xlo,xhi);
    strcat(buffer,line);
    sprintf(line,"%f %f\n",ylo,yhi);
    strcat(buffer,line);
    sprintf(line,"%f %f\n",zlo,zhi);
    strcat(buffer,line);
    strcat(buffer,"ITEM: ATOMS id mol type x y z vx vy vz c_cpe c_cpem c_cpemnos tilt tiltdir twist\n");
    dumpout << buffer;
    // empty buffer
    for (j = 0; j < buffersize; j++) buffer[j] = '\0';
    j = 0;
    for (i = 0; i < nlocal; i++) {
        if (alllist[i].mid == 0) {
            sprintf(line,"%d %d %d %f %f %f %f %f %f %f %d %d %d %d %d\n",
                    i+1,alllist[i].mid,alllist[i].type,alllist[i].x,alllist[i].y,alllist[i].z,
                    alllist[i].vx,alllist[i].vy,alllist[i].vz,alllist[i].pea,0,0,0,0,0);
        }
        else {
            sprintf(line,"%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    i+1,alllist[i].mid,alllist[i].type,alllist[i].x,alllist[i].y,alllist[i].z,
                    alllist[i].vx,alllist[i].vy,alllist[i].vz,alllist[i].pea,samlist[alllist[i].mid - 1].totalpe,
                    samlist[alllist[i].mid - 1].totalpe - samlist[alllist[i].mid - 1].pe[0],
                    samlist[alllist[i].mid - 1].tilt,
                    samlist[alllist[i].mid - 1].tiltdir,
                    samlist[alllist[i].mid - 1].twist);
        }
        strcat(buffer,line);
        j++;
        if (j >= buffersize/256) {
            dumpout << buffer;
            // empty buffer
            for (j = 0; j < buffersize; j++) buffer[j] = '\0';
            j = 0;
        } else if (i == nlocal-1) {
            dumpout << buffer;
            // empty buffer
            for (j = 0; j < buffersize; j++) buffer[j] = '\0';
            j = 0;
        }
    }
    dumpout.close();
}

void samana::structurefactor(double rflag, int nflag)
{
    int i,j,l;
    const int k1d = nflag;
    const int k2d = nflag*nflag;
    double krange = rflag;
    double step = krange*2.0/double(k1d);
    double yhalf = (yhi-ylo)/2.0;
    double zhalf = (zhi-zlo)/2.0;
    typedef struct {
        double y;
        double z;
        double s;
    }sf;
    std::unique_ptr<sf []> k(new sf [k2d]);
    std::unique_ptr<sf []> rhoksulfur(new sf [k2d]);
    std::unique_ptr<sf []> rhokoddcarbon(new sf [k2d]);
    std::unique_ptr<sf []> rhokevencarbon(new sf [k2d]);
    std::unique_ptr<sf []> rhokall(new sf [k2d]);
    double alpha;
    std::cout << "Calculating Structure Factors..." << std::endl;
//    double factorsulfur = 1/sqrt(nsam);
//    double factoroddcarbon = 1/sqrt(nsam*double(nunits%2));
//    double factorevencarbon = 1/sqrt(nsam*double(nunits%2));
//    double factorall = 1/sqrt(nsam*double(nunits));
    double factorsulfur = 1/double(nsam);
    double factoroddcarbon = 1/(nsam*double(nunits%2));
    double factorevencarbon = 1/(nsam*double(nunits%2));
    double factorall = 1/(nsam*double(nunits));
    for (i = 0; i < k1d; i++) {
        for (j = 0; j < k1d; j++) {
            l = i*k1d+j;
            // ky
            k[l].y = (double(i)+0.5)*step-krange;
            // kz
            k[l].z = (double(j)+0.5)*step-krange;
        }
    }
#pragma omp parallel num_threads(4) private(i,alpha,j) \
    shared(rhokall,rhokevencarbon,rhokoddcarbon,rhoksulfur, \
           factorall,factorevencarbon,factoroddcarbon,factorsulfur, \
           yhalf,zhalf,k)
{
    #pragma omp for
    for (i = 0; i < k2d; i++) {
        rhoksulfur[i].y = rhoksulfur[i].z = 0.0;
        rhokoddcarbon[i].y = rhokoddcarbon[i].z = 0.0;
        rhokevencarbon[i].y = rhokevencarbon[i].z = 0.0;
        rhokall[i].y = rhokall[i].z = 0.0;
        for (j = 0; j < nsam; j++) {
            // rho(k)=sum_j=1^N(exp(-i(ky*ry+kz*rz)))
            alpha = (samlist[j].y[0]-yhalf)*k[i].y + (samlist[j].z[0]-zhalf)*k[i].z;
            // sulfur only
            rhoksulfur[i].y += cos(alpha);
            rhoksulfur[i].z += sin(alpha);
            // odd carbon
            for (l = 1; l < nunits; l+=2) {
                alpha = (samlist[j].y[l]-yhalf)*k[i].y + (samlist[j].z[l]-zhalf)*k[i].z;
                rhokoddcarbon[i].y += cos(alpha);
                rhokoddcarbon[i].z += sin(alpha);
            }/*
            // even carbon
            for (l = 2; l < nunits; l+=2) {
                alpha = samlist[j].y[l]*k[i].y + samlist[j].z[l]*k[i].z;
                rhokevencarbon[i].y += cos(alpha);
                rhokevencarbon[i].z += sin(alpha);
            }
            // all
            for (l = 0; l < nunits; l++) {
                alpha = samlist[j].y[l]*k[i].y + samlist[j].z[l]*k[i].z;
                rhokall[i].y += cos(alpha);
                rhokall[i].z += sin(alpha);
            }*/
        }
        rhoksulfur[i].s = (rhoksulfur[i].y*rhoksulfur[i].y+rhoksulfur[i].z*rhoksulfur[i].z)*factorsulfur;
        rhokoddcarbon[i].s = (rhokoddcarbon[i].y*rhokoddcarbon[i].y+rhokoddcarbon[i].z*rhokoddcarbon[i].z)*factoroddcarbon;
        //rhokevencarbon[i].s = (rhokevencarbon[i].y*rhokevencarbon[i].y+rhokevencarbon[i].z*rhokevencarbon[i].z)*factorevencarbon;
        //rhokall[i].s = (rhokall[i].y*rhokall[i].y+rhokall[i].z*rhokall[i].z)*factorall;
    }
}
    std::ofstream sfout;
    std::string filename;
    filename = fileid + ".sf";
    sfout.open(filename.c_str(), std::ofstream::out);
    double step2 = krange*0.01;
    double srhoksulfur,srhokoddcarbon,srhokevencarbon,srhokall;
    int ii,jj;
    int nsum = k1d/200;
    for (i = 0; i < 200; i++) {
        for (j = 0; j < 200; j++) {
            sfout << (double(i)+0.5)*step2-krange << " " << (double(j)+0.5)*step2-krange << " ";
            srhoksulfur = srhokoddcarbon = srhokevencarbon = srhokall = 0.0;
            for (ii = 0; ii < nsum; ii++) {
                for (jj = 0; jj < nsum; jj++) {
                    l = (i*nsum+ii)*k1d+j*nsum+jj;
                    srhoksulfur += rhoksulfur[l].s;
                    srhokoddcarbon += rhokoddcarbon[l].s;
                    srhokevencarbon += rhokevencarbon[l].s;
                    srhokall += rhokall[l].s;
                }
            }
            sfout << srhoksulfur << " " << srhokoddcarbon << " ";
            sfout << srhokevencarbon << " " << srhokall << std::endl;
        }
    }
    sfout.close();
}
