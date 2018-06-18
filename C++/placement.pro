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

win32 {
    QMAKE_CXXFLAGS += -std=c++11 -O2 -frounding-math
    INCLUDEPATH += D:/Software/CGAL-4.5.2/include/
    LIBS += -lCGAL -lCGAL_Core -lgmp
    INCLUDEPATH += D:/
    LIBS += -L"D:/Box2D" -lBox2D
}

unix:!macx {
    QMAKE_CXXFLAGS += -std=c++11 -O2 -frounding-math
    LIBS += -lCGAL -lCGAL_Core -lgmp
    INCLUDEPATH += /mnt/Zero_Data
    LIBS += -L"/mnt/Zero_Data/Box2D" -lBox2D
}

macx {
    QMAKE_CXXFLAGS += -O2
    CONFIG += c++11
    INCLUDEPATH += /usr/local/include/
    LIBS += -lgmp
    LIBS += -L"/usr/local/lib" -lcgal -lcgal_Core
    INCLUDEPATH += /Users/zero/Projects/Box2D-master
    LIBS += -L"/Users/zero/Projects/Box2D-master/Box2D" -lBox2D
}
