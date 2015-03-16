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
    sdk.h \
    mainwidget.h

FORMS    += mainwindow.ui

QMAKE_CXXFLAGS += -frounding-math
LIBS += -lgmp -lCGAL -lCGAL_Core

win32 {
    INCLUDEPATH += D:/
    LIBS += -L"D:/Box2D" -lBox2D
}

unix:!macx {
    INCLUDEPATH += /mnt/Zero_Software
    LIBS += -L"/mnt/Zero_Software/Box2D" -lBox2D
}
