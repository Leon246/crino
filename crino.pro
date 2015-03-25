#-------------------------------------------------
#
# Project created by QtCreator 2014-05-08T16:39:53
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = crino
TEMPLATE = app

CONFIG += c++11
QMAKE_CXXFLAGS += -std=gnu++0x
QMAKE_CXXFLAGS += -fopenmp

LIBS += -fopenmp

SOURCES += main.cpp\
        mainwindow.cpp \
    crino.cpp \
    samana.cpp

HEADERS  += mainwindow.h \
    crino.h \
    samana.h

FORMS    += mainwindow.ui
