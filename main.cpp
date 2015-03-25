#include "mainwindow.h"
#include "crino.h"
#include "samana.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.setWindowTitle(QApplication::translate("crino","crino"));
    w.show();

    return a.exec();
}
