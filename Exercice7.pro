TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

macx {
    QMAKE_CXX = g++-4.8
}

linux {
    QMAKE_CXX = g++-4.9
}

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp

linux {
    QMAKE_CXXFLAGS += -fopenmp-simd
    QMAKE_CXXFLAGS += -Wopenmp-simd
}

QMAKE_CXXFLAGS_DEBUG += -Wall -Wextra

QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -march=native

TARGET = ../Exercice7/Exercice7

SOURCES += \
    Exercice7.cpp

HEADERS += \
    Exercice7_io.h \
    Exercice7_calcul.h

