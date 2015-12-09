TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
#CONFIG -= qt

QMAKE_CXX += -g

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    math/vec3.cpp \
    math/random.cpp \
    io.cpp \
    potentials/potential.cpp \
    statisticssampler.cpp \
    integrators/eulercromer.cpp \
    unitconverter.cpp \
    potentials/lennardjones.cpp \
    thermostat.cpp \
    celllist.cpp

QMAKE_CXXFLAGS += -std=c++11 -fopenmp -O3
QMAKE_LFLAGS +=  -fopenmp

HEADERS += \
    atom.h \
    system.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    potentials/potential.h \
    potentials/lennardjones.h \
    statisticssampler.h \
    integrators/eulercromer.h \
    unitconverter.h \
    thermostat.h \
    celllist.h

