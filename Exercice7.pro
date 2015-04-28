TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

linux {
    QMAKE_CXX = g++-4.9
    QMAKE_CXXFLAGS_DEBUG += -Wall -Wextra

    QMAKE_CXXFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS_RELEASE += -s
    #QMAKE_CXXFLAGS_RELEASE += -march=nehalem
    #QMAKE_CXXFLAGS_RELEASE += -static -static-libgcc -static-libstdc++
    QMAKE_CXXFLAGS_RELEASE += -march=native

    QMAKE_CXXFLAGS += -std=c++11
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS +=  -fopenmp

    QMAKE_CXXFLAGS += -fopenmp-simd
    QMAKE_CXXFLAGS += -Wopenmp-simd
}

macx {
    QMAKE_CXXFLAGS += -std=c++11
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_CXXFLAGS_RELEASE += -O2
}

TARGET = ../Exercice7/Exercice7

SOURCES += \
    Exercice7.cpp

HEADERS += \
    Exercice7_io.h \
    Exercice7_calcul.h

