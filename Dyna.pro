QT += 3dcore 3drender 3dinput 3dextras

CONFIG += c++11

TARGET = Dyna

SOURCES += \
    Setup_Links.cpp \
    Dynamics.cpp \
    main.cpp \
    MakeMovieFrame.cpp

HEADERS += \
    Dynamics.h
