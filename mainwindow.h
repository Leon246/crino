#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "crino.h"
#include "samana.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_SAM_generate_start_clicked();

    void on_Search_start_clicked();

    void on_Compute_Histogram_clicked();

    void on_Preview_clicked();

    void on_Compute_SF_clicked();

    void on_Compute_Normal_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
