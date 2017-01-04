/****************************************************************************
**
** Copyright (C) 2016 William W. Armstrong
** Contact: William W. Armstrong
** 3624 - 108 Street NW, Edmonton, Alberta, Canada T6J 1B4
**
**
                   GNU LESSER GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.


  This version of the GNU Lesser General Public License incorporates
the terms and conditions of version 3 of the GNU General Public
License, supplemented by the additional permissions listed below.

  0. Additional Definitions.

  As used herein, "this License" refers to version 3 of the GNU Lesser
General Public License, and the "GNU GPL" refers to version 3 of the GNU
General Public License.

  "The Library" refers to a covered work governed by this License,
other than an Application or a Combined Work as defined below.

  An "Application" is any work that makes use of an interface provided
by the Library, but which is not otherwise based on the Library.
Defining a subclass of a class defined by the Library is deemed a mode
of using an interface provided by the Library.

  A "Combined Work" is a work produced by combining or linking an
Application with the Library.  The particular version of the Library
with which the Combined Work was made is also called the "Linked
Version".

  The "Minimal Corresponding Source" for a Combined Work means the
Corresponding Source for the Combined Work, excluding any source code
for portions of the Combined Work that, considered in isolation, are
based on the Application, and not on the Linked Version.

  The "Corresponding Application Code" for a Combined Work means the
object code and/or source code for the Application, including any data
and utility programs needed for reproducing the Combined Work from the
Application, but excluding the System Libraries of the Combined Work.

  1. Exception to Section 3 of the GNU GPL.

  You may convey a covered work under sections 3 and 4 of this License
without being bound by section 3 of the GNU GPL.

  2. Conveying Modified Versions.

  If you modify a copy of the Library, and, in your modifications, a
facility refers to a function or data to be supplied by an Application
that uses the facility (other than as an argument passed when the
facility is invoked), then you may convey a copy of the modified
version:

   a) under this License, provided that you make a good faith effort to
   ensure that, in the event an Application does not supply the
   function or data, the facility still operates, and performs
   whatever part of its purpose remains meaningful, or

   b) under the GNU GPL, with none of the additional permissions of
   this License applicable to that copy.

  3. Object Code Incorporating Material from Library Header Files.

  The object code form of an Application may incorporate material from
a header file that is part of the Library.  You may convey such object
code under terms of your choice, provided that, if the incorporated
material is not limited to numerical parameters, data structure
layouts and accessors, or small macros, inline functions and templates
(ten or fewer lines in length), you do both of the following:

   a) Give prominent notice with each copy of the object code that the
   Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the object code with a copy of the GNU GPL and this license
   document.

  4. Combined Works.

  You may convey a Combined Work under terms of your choice that,
taken together, effectively do not restrict modification of the
portions of the Library contained in the Combined Work and reverse
engineering for debugging such modifications, if you also do each of
the following:

   a) Give prominent notice with each copy of the Combined Work that
   the Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the Combined Work with a copy of the GNU GPL and this license
   document.

   c) For a Combined Work that displays copyright notices during
   execution, include the copyright notice for the Library among
   these notices, as well as a reference directing the user to the
   copies of the GNU GPL and this license document.

   d) Do one of the following:

       0) Convey the Minimal Corresponding Source under the terms of this
       License, and the Corresponding Application Code in a form
       suitable for, and under terms that permit, the user to
       recombine or relink the Application with a modified version of
       the Linked Version to produce a modified Combined Work, in the
       manner specified by section 6 of the GNU GPL for conveying
       Corresponding Source.

       1) Use a suitable shared library mechanism for linking with the
       Library.  A suitable mechanism is one that (a) uses at run time
       a copy of the Library already present on the user's computer
       system, and (b) will operate properly with a modified version
       of the Library that is interface-compatible with the Linked
       Version.

   e) Provide Installation Information, but only if you would otherwise
   be required to provide such information under section 6 of the
   GNU GPL, and only to the extent that such information is
   necessary to install and execute a modified version of the
   Combined Work produced by recombining or relinking the
   Application with a modified version of the Linked Version. (If
   you use option 4d0, the Installation Information must accompany
   the Minimal Corresponding Source and Corresponding Application
   Code. If you use option 4d1, you must provide the Installation
   Information in the manner specified by section 6 of the GNU GPL
   for conveying Corresponding Source.)

  5. Combined Libraries.

  You may place library facilities that are a work based on the
Library side by side in a single library together with other library
facilities that are not Applications and are not covered by this
License, and convey such a combined library under terms of your
choice, if you do both of the following:

   a) Accompany the combined library with a copy of the same work based
   on the Library, uncombined with any other library facilities,
   conveyed under the terms of this License.

   b) Give prominent notice with the combined library that part of it
   is a work based on the Library, and explaining where to find the
   accompanying uncombined form of the same work.

  6. Revised Versions of the GNU Lesser General Public License.

  The Free Software Foundation may publish revised and/or new versions
of the GNU Lesser General Public License from time to time. Such new
versions will be similar in spirit to the present version, but may
differ in detail to address new problems or concerns.

  Each version is given a distinguishing version number. If the
Library as you received it specifies that a certain numbered version
of the GNU Lesser General Public License "or any later version"
applies to it, you have the option of following the terms and
conditions either of that published version or of any later version
published by the Free Software Foundation. If the Library as you
received it does not specify a version number of the GNU Lesser
General Public License, you may choose any version of the GNU Lesser
General Public License ever published by the Free Software Foundation.

  If the Library as you received it specifies that a proxy can decide
whether future versions of the GNU Lesser General Public License shall
apply, that proxy's public statement of acceptance of any version is
permanent authorization for you to choose that version for the
Library.
****************************************************************************/
// file main.cpp

// A description of the recursive solution of the equations of motion of a tree-linkage is
// "The dynamics of articulated rigid bodies for purposes of animation"
// by William W. Armstrong and Mark W. Green, The Visual Computer (1985) 1: 231 - 240
// This new implementation uses Eigen, Boost, and Qt and views the motion using Qt3D
#include <QGuiApplication>
#include <Qt3DWindow>
#include <QPropertyAnimation>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DCore/QAspectEngine>
#include <Qt3DInput/QInputAspect>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QCameraLens>
#include <Qt3DRender/QRenderAspect>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QPlaneMesh>
#include <Qt3DExtras/QCylinderMesh>
#include "Dynamics.h"
//#include <iostream>
#include <fstream>
#define radianstodegrees 57.29578049
constexpr int movieSize = 1000; // Number of frames in the movie. If you want a longer movie, increase this.
double movie[movieSize][36];
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
// Get a specific frame of the movie (since the movie doesn't yet work)
int captureFrame = 999;
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void initialize();
void motion();
void storeFrame(int fn, double st);
QVector3D viewCentre;
QVector3D cameraPosition;
Qt3DCore::QEntity *createScene();
void run();
using namespace std;

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);
    initialize(); // Sets up the dynamics and graphics
    // We simulate the motion and store the data of a frame for later viewing in a movie
    for(int frame_number = 0; frame_number < movieSize ; frame_number++)
    {
        // Motion performs dynamic simulation for a short time, e.g. 20 timesteps, between frames of the movie
        motion();
        // We collect data about link positions and orientations for later viewing
        storeFrame(frame_number, SimTime);
    }
    // We view the object from 40 cm distance
    cameraPosition = viewCentre;
    cameraPosition.setZ(viewCentre.z()-40.0f);
    // The view
    Qt3DExtras::Qt3DWindow view;
    QSurfaceFormat format;
    format.setVersion(3, 2);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setDepthBufferSize(24);
    view.setFormat(format);
    // Create a camera to point at the scene in the view window
    Qt3DRender::QCamera *camera = view.camera();
    camera->lens()->setPerspectiveProjection(100.0f, 16.0f/9.0f, 0.1f, 1000.0f);
    // The camera position (now depends on the chosen frame, like the view centre)
    // Later it must follow the motion
    camera->setPosition(cameraPosition);
    // camera view centre
    camera->setViewCenter(viewCentre);
    view.setRootEntity(createScene());
    view.show();
    return app.exec();
}

void storeFrame(int fn, double st)
{
    // fn = frame_number, st = the simulation time of that frame
    Vector3d tempv;
    movie[fn][0] = st;
    // We store everything simply as doubles, which facilitates a conversion
    // of quaternions and vectors from Eigen format, used for dynamics,
    // to Qt format, used for graphics, via their components
    for(int n = 0; n < numlinks; n++)
    {
        movie[fn][1 + 4 * n] = qRI[n].w();
        movie[fn][2 + 4 * n] = qRI[n].x();
        movie[fn][3 + 4 * n] = qRI[n].y();
        movie[fn][4 + 4 * n] = qRI[n].z();
        tempv = p[n] + RI[n] * c[n];       // the c.g. is the reference point of cylinder meshes
        movie[fn][21 + 3 * n] = tempv(0);
        movie[fn][22 + 3 * n] = tempv(1);
        movie[fn][23 + 3 * n] = tempv(2);
    }
}

//Array of cylinder entities
Qt3DCore::QEntity *cylinderEntity[5];
// Array of cylinder meshes
Qt3DExtras::QCylinderMesh *cylinderMesh[5];
// Array of cylinder Transforms
Qt3DCore::QTransform *cylinderTransform[5];

Qt3DCore::QEntity *createScene()
{
    // The root entity contains components which make up the scene graph.
    Qt3DCore::QEntity *linkageFrame = new Qt3DCore::QEntity;
    // Material
    Qt3DExtras::QPhongMaterial *material = new Qt3DExtras::QPhongMaterial(linkageFrame);
    material->setDiffuse(QColor(QRgb(0x928327)));
    // Plane entity -- we put in a plane at floor level (y = 0)
    Qt3DCore::QEntity *planeEntity = new Qt3DCore::QEntity(linkageFrame);
    // Plane mesh
    // Note: plane size to be enlarged when Qt3D computes hidden surfaces correctly
    Qt3DExtras::QPlaneMesh *planeMesh = new Qt3DExtras::QPlaneMesh;
    planeMesh->setWidth(18.0f);
    planeMesh->setHeight(8.0f);
    // Plane transform
    Qt3DCore::QTransform *planeTransform = new Qt3DCore::QTransform;
    // Put the plane close to where hinge 0 of the linkage is (so we can see something!)
    planeTransform->setTranslation(QVector3D(0.0f, (float)p[0](1), 0.0f));
    //Add plane components
    planeEntity->addComponent(planeMesh);
    planeEntity->addComponent(planeTransform);
    planeEntity->addComponent(material);

    int fn = captureFrame;
    viewCentre = QVector3D(
                movie[fn][21],
                movie[fn][22],
                movie[fn][23]); // These are the coordinates of the c.g. of link 0 of the specified frame
    // Linkage image data
    for(int n = 0; n < numlinks; n++)
    {
        // Create the cylinder entity, mesh and transform for link n
        cylinderEntity[n] = new Qt3DCore::QEntity(linkageFrame);
        cylinderMesh[n] = new Qt3DExtras::QCylinderMesh;
        cylinderMesh[n]->setRadius(radius[n]);
        cylinderMesh[n]->setLength(length[n]);
        cylinderMesh[n]->setRings(20);
        cylinderMesh[n]->setSlices(20);
        cylinderTransform[n] = new Qt3DCore::QTransform;
        cylinderTransform[n]->setScale3D(QVector3D(1.0f, 1.0f, 1.0f));
        SimTime = movie[fn][0];
        cylinderTransform[n]->setRotation(QQuaternion(
            movie[fn][1 + 4 * n],
            movie[fn][2 + 4 * n],
            movie[fn][3 + 4 * n],
            movie[fn][4 + 4 * n]));
        cylinderTransform[n]->setTranslation(QVector3D(
                    movie[fn][21 + 3 * n],
                    movie[fn][22 + 3 * n],
                    movie[fn][23 + 3 * n]));
        cylinderEntity[n]->addComponent(cylinderMesh[n]);
        cylinderEntity[n]->addComponent(cylinderTransform[n]);
        cylinderEntity[n]->addComponent(material);
    }
    return linkageFrame;
}

void run() // not used =-- the movie doesn't work yet
{
    for(int fn = 0;fn < movieSize; fn++)
    {
        for(int n = 0; n < numlinks; n++)
        {
            SimTime = movie[fn][0];
            cylinderTransform[n]->setRotation(QQuaternion(
                movie[fn][1 + 4 * n],
                movie[fn][2 + 4 * n],
                movie[fn][3 + 4 * n],
                movie[fn][4 + 4 * n]));
            cylinderTransform[n]->setTranslation(QVector3D(
                movie[fn][21 + 3 * n],
                movie[fn][22 + 3 * n],
                movie[fn][23 + 3 * n]));
        }
    }
}



