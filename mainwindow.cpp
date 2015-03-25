#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QString>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_Search_start_clicked()
{
    int i = 0;

    //  read files
    QString fileDir = ui->search_dir->text();
    fileDir.append("dump.");
    QString fileName;
    QFile fin;
    std::unique_ptr<samana> obj(new samana());
    int vflag = 2;
    if (ui->isdumpfileoldsam->isChecked()) {
        vflag = 0;
    } else if (ui->isdumpfilewopea->isChecked()) {
        vflag = 1;
    } else if (ui->isdumpfilewpea->isChecked()) {
        vflag = 2;
    }
    obj->init(ui->samunits->value());
    for (i = ui->numstart->value(); i <= ui->numend->value(); i += 1000) {
        //step = i/1000;
        fileName = fileDir + QString::number(i);
        if (fin.exists(fileName)) {
            //fin.close();
            ui->label_2->setText(QString("Reading file ").append(fileName));
            obj->readfile(fileName.toStdString().c_str(),vflag);
            obj->assemblesam();
            obj->compute();
            obj->dumpout();
            obj->output();
        } else {
            //fin.close();
            std::cout << "File " << fileName.toStdString() << " not found" << std::endl;
        }
    }

}

// histogram
void MainWindow::on_Compute_Histogram_clicked()
{
    int i,j;
    int binsize = ui->Binsizevalue->value();
    int nstart = ui->numstart->value();
    int nend = ui->numend->value();
    int ntimesteps = (nend-nstart)/1000+1;
    int step = 0;
    // -90 ~ 90
    std::unique_ptr<int []> histotilt(new int [binsize*ntimesteps]);
    // -180 ~ 180
    std::unique_ptr<int []> histotiltdir(new int [2*binsize*ntimesteps]);
    std::unique_ptr<int []> histotwist(new int [2*binsize*ntimesteps]);
    for (i = 0; i < binsize*ntimesteps; i++) {
        histotilt[i] = 0;
        histotiltdir[i] = 0;
        histotwist[i] = 0;
    }
    for (i = binsize*ntimesteps; i < 2*binsize*ntimesteps; i++) {
        histotiltdir[i] = 0;
        histotwist[i] = 0;
    }
    double tilt,tiltdir,twist;
    int itilt,itiltdir,itwist;

    QString fileDir = ui->search_dir->text();
    fileDir.append("dump.");
    QString fileName;
    QFile fin;
    QString line;
    QStringList list;
    //int progressbarvalue = 0;
    for (i = nstart; i <= nend; i += 1000) {
        step = (i-nstart)/1000;
        fileName = fileDir + QString::number(i) + QString(".output");
        fin.setFileName(fileName);
        if (!fin.open(QIODevice::ReadOnly | QIODevice::Text)) continue;
        std::cout << "Reading file dump." << i << ".output ..." << std::endl;
        while (!fin.atEnd()) {
            line = fin.readLine();
            list = line.split(" ");
            tilt = list[6].toDouble();
            tiltdir = list[7].toDouble();
            twist = list[8].toDouble();
            itilt = int((tilt+90.0)/(180.0/double(binsize)));
            itiltdir = int((tiltdir+180.0)/(180.0/double(binsize)));
            itwist = int((twist+180.0)/(180.0/double(binsize)));
            histotilt[step*binsize+itilt]++;
            histotiltdir[step*binsize*2+itiltdir]++;
            histotwist[step*binsize*2+itwist]++;
        }
        fin.close();
        //progressbarvalue = int(double(step)/double(ntimesteps)*100.0);

    }

    // result
    std::ofstream output;
    std::string filename;
    filename = (ui->search_dir->text()).toStdString() + "histotilt.txt";
    output.open(filename.c_str(), std::ofstream::out);
    for (i = 0; i < ntimesteps; i++) {
        for (j = 0; j < binsize; j++) {
            output << i << " " << -90.0+(double(j)+0.5)*(180.0/double(binsize));
            output << " " << histotilt[i*binsize+j] << std::endl;
        }
    }
    output.close();
    filename = (ui->search_dir->text()).toStdString() + "histotiltdir.txt";
    output.open(filename.c_str(), std::ofstream::out);
    for (i = 0; i < ntimesteps; i++) {
        for (j = 0; j < binsize*2; j++) {
            output << i << " " << -180.0+(double(j)+0.5)*(180.0/double(binsize));
            output << " " << histotiltdir[i*binsize*2+j] << std::endl;
        }
    }
    output.close();
    filename = (ui->search_dir->text()).toStdString() + "histotwist.txt";
    output.open(filename.c_str(), std::ofstream::out);
    for (i = 0; i < ntimesteps; i++) {
        for (j = 0; j < binsize*2; j++) {
            output << i << " " << -180.0+(double(j)+0.5)*(180.0/double(binsize));
            output << " " << histotwist[i*binsize*2+j] << std::endl;
        }
    }
    output.close();
}

void MainWindow::on_Preview_clicked()
{
    double xlo, xhi, ylo, yhi, zlo, zhi;
    xlo = ui->xlo->text().toDouble();
    xhi = ui->xhi->text().toDouble();
    ylo = ui->ylo->text().toDouble();
    yhi = ui->yhi->text().toDouble();
    zlo = ui->zlo->text().toDouble();
    zhi = ui->zhi->text().toDouble();
    double xlh, ylh, zlh;
    xlh = xhi - xlo;
    ylh = yhi - ylo;
    zlh = zhi - zlo;
    ui->xlh->setText("xlh = " + QString::number(xlh));
    ui->ylh->setText("ylh = " + QString::number(ylh));
    ui->zlh->setText("zlh = " + QString::number(zlh));
}

void MainWindow::on_SAM_generate_start_clicked()
{
    std::cout << "Generating crino.data..." << std::endl;
    crino *ice = new crino();
    double zfactor = 4.08*sqrt(0.375)*2.0; //4.08*sqrt(0.375)/2.0;
    double yfactor = 4.08*sqrt(0.5)*3.0; //4.08*sqrt(0.5)/2.0;
    double box[6];

    box[0] = -60.;
    box[1] = 100.;
    box[2] = -17.6*yfactor;
    box[3] = 17.4*yfactor;
    box[4] = -30.3*zfactor;
    box[5] = 30.7*zfactor;

    box[0] = ui->xlo->text().toDouble();
    box[1] = ui->xhi->text().toDouble();
    box[2] = ui->ylo->text().toDouble();
    box[3] = ui->yhi->text().toDouble();
    box[4] = ui->zlo->text().toDouble();
    box[5] = ui->zhi->text().toDouble();
    int nunit = ui->samunits->text().toInt();

    ice->init(box, nunit);
    std::cout << "Running SAM..." << std::endl;
    ice->multipleSAM(300.0);
    //std::cout << "Running Au111..." << std::endl;
    //ice->Au111();
    ice->shiftyz();
    ice->output(ui->search_dir->text().toStdString());
    ice->dumptextdump(ui->search_dir->text().toStdString());
    delete ice;
    ui->SAM_generate_status->setText("Complete");
}

void MainWindow::on_Compute_SF_clicked()
{
    int i = 0;

    //  read files
    QString fileDir = ui->search_dir->text();
    fileDir.append("dump.");
    QString fileName;
    QFile fin;
    std::unique_ptr<samana> obj(new samana());
    int vflag = 2;
    if (ui->isdumpfileoldsam->isChecked()) {
        vflag = 0;
    } else if (ui->isdumpfilewopea->isChecked()) {
        vflag = 1;
    } else if (ui->isdumpfilewpea->isChecked()) {
        vflag = 2;
    }
    obj->init(ui->samunits->value());
    for (i = ui->numstart->value(); i <= ui->numend->value(); i += 1000) {
        //step = i/1000;
        fileName = fileDir + QString::number(i);
        if (fin.exists(fileName)) {
            //fin.close();
            ui->label_2->setText(QString("Reading file ").append(fileName));
            obj->readfile(fileName.toStdString().c_str(),vflag);
            obj->assemblesam();
            obj->compute();
            obj->structurefactor(ui->krange->value(),ui->k1d->value());
        } else {
            //fin.close();
            std::cout << "File " << fileName.toStdString() << " not found" << std::endl;
        }
    }
}

void MainWindow::on_Compute_Normal_clicked()
{
    int i,j;
    int binsize = ui->Binsizevalue_2->value();
    int nstart = ui->numstart->value();
    int nend = ui->numend->value();
    int ntimesteps = (nend-nstart)/1000+1;
    int step = 0;
    std::unique_ptr<int []> histotilt(new int [binsize]);
    std::unique_ptr<int []> histotiltdir(new int [binsize]);
    std::unique_ptr<int []> histotwist(new int [binsize]);
    /*
    QString fileDir = ui->search_dir->text();
    fileDir.append("dump.");
    QString fileName;
    QFile fin;
    QString line;
    QStringList list;
    //int progressbarvalue = 0;
    for (i = nstart; i <= nend; i += 1000) {
        step = (i-nstart)/1000;
        fileName = fileDir + QString::number(i) + QString(".output");
        fin.setFileName(fileName);
        if (!fin.open(QIODevice::ReadOnly | QIODevice::Text)) continue;
        std::cout << "Reading file dump." << i << ".output ..." << std::endl;
        while (!fin.atEnd()) {
            line = fin.readLine();
            list = line.split(" ");
            tilt = list[6].toDouble();
            tiltdir = list[7].toDouble();
            twist = list[8].toDouble();
            itilt = int((tilt+90.0)/(180.0/double(binsize)));
            itiltdir = int((tiltdir+90.0)/(180.0/double(binsize)));
            itwist = int((twist+180.0)/(180.0/double(binsize)));
            histotilt[step*binsize+itilt]++;
            histotiltdir[step*binsize+itiltdir]++;
            histotwist[step*binsize*2+itwist]++;
        }
        fin.close();
        //progressbarvalue = int(double(step)/double(ntimesteps)*100.0);

    }*/
}
