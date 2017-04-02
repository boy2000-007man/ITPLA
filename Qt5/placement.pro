#-------------------------------------------------
#
# Project created by QtCreator 2015-02-11T20:46:54
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = placement
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp \
        mainwidget.cpp

HEADERS  += mainwindow.h \
        mainwidget.h \
    ITPLA.h

FORMS    += mainwindow.ui

QMAKE_CXXFLAGS += -std=c++11 -O2 -frounding-math
LIBS += -lgmp -lCGAL -lCGAL_Core

win32 {
    INCLUDEPATH += D:/
    INCLUDEPATH += D:/Software/CGAL-4.5.2/include/
    LIBS += -L"D:/Box2D" -lBox2D
}

unix:!macx {
    INCLUDEPATH += /mnt/Zero_Data
    LIBS += -L"/mnt/Zero_Data/Box2D" -lBox2D
}
