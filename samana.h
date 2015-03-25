#ifndef SAMANA_H
#define SAMANA_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <memory>

namespace Ui {
class samana;
}

typedef struct
{
    int id;
    int type;
    int mid;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double pea;
    double pem;
}atom;

typedef struct
{
    int tag[32];
    int mol[32];
    int type[32];
    double x[32];
    double y[32];
    double z[32];
    double pe[32];
    double t[3];
    double tilt;
    double tiltdir;
    double twist;
    double totalpe;
}sam;

class samana {
public:
    std::unique_ptr<atom []> alllist;
    std::unique_ptr<sam []> samlist;
    int nlocal;
    int nsam;
    int nunits;
    int mdstep;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    std::string fileid;
    samana();
    ~samana();
    void init(int);
    int readfile(std::string, int);
    void assemblesam();
    void compute();
    void output();
    void dumpout();
    void structurefactor(double,int);
protected:
private:
};

#endif // SAMANA_H
