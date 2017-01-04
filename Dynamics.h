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
#pragma once

// file DynaEigen.h

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Most of the following quantities are initialized in Setup_Links.cpp

// Integer quantities which define the structure of the tree-linkage.
extern const int numlinks;      // The links are numbered from 0 to numlinks - 1. set in #define LINKS
extern const int Parent[];      // The index of a link should always follow the index of its parent in this array.
extern const int length[];      // The length of the cylinders to display graphically
extern const int radius[];      // The radius of the cylinders
extern const int isParent[];    // This is 1 if the link is the parent of another link
extern const double elasticity; // floor elasticity
extern const double viscosity;  // floor viscosity
extern const double air_resistance; // air resistance to prevent excessive velocities in the x and y directions.
extern const double rotElastic; // control of links in elastic directions
extern const double rotDamp;    // rotational damping to diminish high speed rotations along long link axes
// The simulation time increment
extern const double deltat;
extern const double deltat_sq;
extern double SimTime; // time since the start of simulation

// Vector of scalars
extern const double m[];           // masses of the links, e.g.m[r] is the mass of link r
extern const double totalMass;
// The acceleration of gravity in the inertial frame,
// a right-handed coordinate system whose third axis (index 2) points out of the screen
extern Vector3d aG;                // The acceleration of gravity 980.7 cm/sec^2

// The following are vectors of 3-dimensional vectors, indexed by link number;
// e.g.fEdistal[r] is a 3-vector (represented as a column matrix) associated with link r
// It is in the inertial frame (of the viewer). Other quantities are represented in
// frames attached to the links because such a quantity is constant in that frame.
// The value in component i of fEdistal[r] is fEdistal[r](i), for i = 0,1,2.
// Using brackets in fEdistal[r](i) for the link number instead of fE(r)(i)
// means that bounds-checking is omitted to make execution faster.

// Representations in the inertial frame

extern vector<Vector3d> p;         // vector from the origin of the inertial frame to the proximal hinge of link r, which joins link r to its parent
extern vector<Vector3d> v;         // velocity of the proximal hinge of link r
extern vector<Vector3d> gE;        // external torque acting on link r ( for a rigid link, the position of application of the torque doesn't matter)
extern vector<Vector3d> pEproximal;// position at the proximal ends of links for applying external force.
extern vector<Vector3d> pEdistal;  // position at the distal ends of links for applying external force
extern vector<Vector3d> fEproximal; // external force acting on link r at the position given by vector pEproximal
extern vector<Vector3d> fEdistal;    // external force acting on link r at the position given by vector pEdistal
//extern vector<Vector3d> f1Eproximal;// timestep -1 of application of external force
//extern vector<Vector3d> f1Edistal;  // timestep -1 of application of external force
//extern vector<Vector3d> f2Eproximal;// timestep -2 of application of external force
//extern vector<Vector3d> f2Edistal;  // timestep -2 of application of external force
extern vector<Vector3d> vdistal;    // this is the velocity of the distal end of link r translated down one radius
extern vector<Vector3d> downRadius; // down one radius from the contact point of the link

// Representations in the frames of link r
extern vector<Vector3d> a;         // acceleration of the proximal hinge
extern vector<Vector3d> omega;     // angular velocity
extern vector<Vector3d> omegadot;  // rate of change of angular velocity(LINKS,0);
extern vector<Vector3d> deltau;    // the incremental angle turned
extern vector<Vector3d> c;         // vector from proximal hinge to centre of mass
extern vector<Vector3d> f;         // force on its parent at the proximal hinge
extern vector<Vector3d> g;         // torque of link on its parent
extern vector<Vector3d> pEproximal;// points of application of external force at the proximal ends of the links
extern vector<Vector3d> pEdistal;  // points of application of external force at the distal ends of the links
extern double           m_min;     // a mass less than half of the lightest link
extern double           m_max;     // a mass greater than the total mass of the linkage
extern vector<Matrix3d> J;         // moment of inertia matrices of link about its proximal hinge
extern vector<Vector3d> l;          //vector from proximal hinge of parent to proximal hinge of link r in the frame of the parent (or the inertial frame if r is 0)
// The following are vectors of 3x3 rotation matrices, indexed by links, e.g. R[2] is a 3x3 rotation matrix
extern vector<Matrix3d> R;         // converts from frame of link r to frame of parent, or to the inertial frame in the case r = 0
extern vector<Matrix3d> RT;        // inverse, i.e. transpose, of R matrix
extern vector<Matrix3d> RI;        // converts from frame of link r to the inertial frame
extern vector<Matrix3d> RIT;       // inverse, i.e. transpose, of RI matrix
extern vector<Quaterniond> qR;     // rotation vectors for the links
extern vector<Quaterniond> qRI;    // rotation vector of the root w.r.t. the inertial frame
extern vector<Quaterniond> qRdesired;      // quaternion to show where the link should be positioned
// The quantities below are of interest only for the internal workings of the dynamics computations.
// They are initialized when computed
extern vector<Matrix3d> Q;         // useful for slowband inbound part of solving the equations of motion
extern vector<Matrix3d> W;         //   "
extern vector<Matrix3d> T;         //   "
extern vector<Matrix3d> K;         //   "
extern vector<Vector3d> d;         //   "
extern vector<Matrix3d> M;         //   "
extern vector<Vector3d> fprime;    // the reactive force that results from accelerating a subtree of the linkage
extern vector<Vector3d> fdblprime; // the last above, but in parent coordinates
extern vector<Matrix3d> munitm;    // mass of link times the unit matrix
extern vector<Vector3d> mc;        // mass of link times vector to centre of mass
extern vector<Matrix3d> mctilde;   // mass times tilde matrix of vector to centre of mass (produces cross-product)
extern vector<Matrix3d> ltilde;    // tilde of l
extern vector<Vector3d> aC;        // centripetal acceleration in frame of the parent hinge
extern vector<Vector3d> maG;       // force of gravity on the link
extern vector<Vector3d> fsigma;    // a sum of forces useful in slowband
extern vector<Vector3d> gsigma;    // a sum of torques useful in slowband
extern vector<Vector3d> g1sigma;   // a sum of torques useful in fastband
extern vector<Vector3d> deltarot;  // an incremental rotation vector
extern Matrix3d         zerom;      // the 3x3 zero matrix
extern Vector3d         zerov;      // the 3d zero vector

