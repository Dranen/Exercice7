TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

QMAKE_CXXFLAGS_RELEASE += -O2
QMAKE_CXXFLAGS_RELEASE += -march=ivybridge


SOURCES += \
    Exercice7.cpp

HEADERS += \
    Exercice7_io.h \
    Exercice7_calcul.h

