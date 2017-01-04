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

// file Setup_Links.cpp

// SEE POSSIBLE TEMPORARY VALUES FOR TESTING AT THE BOTTOM OF THIS FILE
#define PI_OVER_2 1.570796
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Dynamics.h"
using namespace std;
using namespace Eigen;
const double s7071  = sqrt(0.5);
void adjustLinks();
// Units used are CGS: centimeter, gram, second.
// Tree linkage structure:
// Constant: numlinks = the number of links in the linkage.
// The links are numbered 0,..., numlinks -1.
// The index of a parent link should always precede the index of the link in this order.
// The main link (the "body") has index 0.
// Each link s has exactly one parent r which has a lower index in the order of links,
// except for r = 0 whose parent has the out-of-bounds parent value numlinks.
// For a chain of 5 links, the initialization of Parent would be 5, 0, 1, 2, 3
// For an animal with four legs the Parent array would be 5, 0, 0, 0, 0
// You can change the tree-linkage below, e.g. to make a creature with head, neck, and other extremities
constexpr int numlinks = 5;
// Each of links has a parent except link 0 which has non-existing parent  = numlinks
const int Parent[numlinks] = {5,0,0,0,0};
const int length[numlinks] = {20, 20, 20, 20, 20};
const int radius[numlinks] = {5, 1, 1, 1, 1};
const int isParent[numlinks] = {1, 0, 0, 0, 0}; // 1 if the link is a parent of another link
// Simulation time step
constexpr double deltat = 0.001;
constexpr double deltat_sq = deltat * deltat;
double SimTime = 0.0; // elapsed time since the start of a simulation run.
// initialize the masses in grams
const double m[numlinks] = {10000.0, 1000.0, 1000.0, 1000.0, 1000.0};
const double totalMass = 14000.0; // not used

// initialize the viscosity, elasticity etc. of the floor
const double viscosity  = 4500.0;//  Consider a viscosity of 1000.0 dynes/(cm/sec)
// After entering the floor after falling a quarter of a second, say , the velocity could be
// say 1/32 * 980.7 cm/sec or 30.6469 cm/sec, giving rise to a force of 30646.9 dynes on each foot.
// If a mass of 14000.0 grams is pushed with this force, it will decelerate at 2.189 cm/ sec^2.
// So stopping in this case is too slow.  We try 4000 as viscosity.
// If the body has all 10 hinge ends below the floor, then the force is multiplied by 10 instead of 4.

const double air_resistance = 250.0;  // There is an unexpected velocity buildup in the x direction.
// This could be correct dynamics and nothing to slow down the motion of the linkage in x, or
// it could be due to the matrix M[0] being ill-conditioned, which could lead to a faulty acceleration a[0].

const double elasticity = 68650.0; // Note: Consider 34325.0 -- this elasticity, measured in dynes/cm, allows the four feet
// to penetrate around 10 cm before the force upward on four feet equals the weight of the linkage.
// With the weight compensated for, we don't need as much elasticity. 68650 supports 5 cm penetration

// Initialize quantities that will eventually control the linkage, set to 0.0 for floppy limbs
const double rotElastic = 0.0; // For a value 10,000,000.0: 1000 g at a distance of 10 cm under aG needs counter torque
// of over 9,807,000 dyne cm/radian.  This could hold a 20 cm limb sideways at an angle of about one radian down.
// The value 10 times greater, could hold the position at a deviation about 6 degrees.
const double rotDamp = 100.0;  // This is mainly to try to counteract the spinning of the legs in the long axis.
// There are also fast rotations in the other axes when the links hit the floor. Now this is "frigged" as
// they say in the aircraft simulator lingo in order to make the simulation more stable.

// The acceleration of gravity vector is in the inertial frame, a right-handed coordinate system
// whose x-axis [0] points to the right, whose y-axis [1] points up, and
// whose z-axis [2] points out of the screen towards the viewer.
// The acceleration of gravity is negative in the y direction.
// It is initialized here because of the quantities below which need it: mc, maG, J
Vector3d aG = Vector3d(0.0, -980.7, 0.0); // Acceleration of gravity (downward)

// The following are vectors of 3-dimensional vectors, indexed by link number;
// e.g.fEdistal[r] is a 3-vector (represented as a column matrix) associated with link r
// The value in component i of fEdistal[r] is fEdistal[r](i), for i = 0,1,2.

// Representations in the inertial frame
vector<Vector3d> p(numlinks);         // vector from the origin of the 3D model to the proximal hinge of link r, which joins link r to its parent
vector<Vector3d> v(numlinks);         // velocity of the proximal hinge of link r
vector<Vector3d> fEproximal(numlinks);        // external force acting on link r at the proximal end
vector<Vector3d> fEdistal(numlinks);        // external force acting on link r at the distal end
vector<Vector3d> gE(numlinks);        // external torque acting on link r ( for a rigid link, the position of application of the torque doesn't matter)
vector<Vector3d> lastpEproximal(numlinks); // something to track the contact point of the proximal ends of cylinders
vector<Vector3d> downRadius(numlinks); // in effect, adds balls to the end of all cylinders to which floor forces are applied
// Representations in the frames of link r
vector<Vector3d> a(numlinks);         // acceleration of the proximal hinge
vector<Vector3d> omega(numlinks);     // angular velocity
vector<Vector3d> omegadot(numlinks);  // rate of change of angular velocity
vector<Vector3d> deltau(numlinks);    // the incremental deltau turned
vector<Vector3d> c(numlinks);         // vector from proximal hinge to centre of mass (c.g.), also the centre of the meshes in Qt3D
vector<Vector3d> f(numlinks);         // force on its parent at the proximal hinge
vector<Vector3d> g(numlinks);         // torque of link on its parent
// The following are in the inertial frame, since is applied at a point straight down one radius from a link end
vector<Vector3d> pEproximal(numlinks);// position at the proximal ends of links for applying external force.
vector<Vector3d> pEdistal(numlinks);  // position of application of external force at the distal ends of links
vector<Vector3d> vdistal(numlinks);    // this is the velocity of the distal end of link r translated down one radius
double           m_min;                // a mass less than half of the lightest link
double           m_max;                // a mass greater than the total mass of the linkage
vector<Matrix3d> J(numlinks);          // moment of inertia matrices of link about its proximal hinge ("frigged" for added stability)

vector<Vector3d> l(numlinks); // Vector to the proximal hinge of link r from the proximal hinge of the parent of r
                              //  in the frame of the parent. (not used if r is 0, but p[0] is about the same thing)
// The following are vectors of 3x3 rotation matrices, indexed by numlinks, e.g. R[2] is a 3x3 rotation matrix
vector<Matrix3d> R(numlinks);         // converts from frame of link r to frame of parent, or to the inertial frame in the case r = 0
vector<Matrix3d> RT(numlinks);        // inverse, i.e. transpose, of R matrix
vector<Matrix3d> RI(numlinks);        // converts from frame of link r to the inertial frame
vector<Matrix3d> RIT(numlinks);       // inverse, i.e. transpose, of RI matrix
vector<Quaterniond> qR(numlinks);     // quaternion instead of R
vector<Quaterniond> qRI(numlinks);    // quaternion instead of RI
vector<Quaterniond> qRdesired(numlinks);  // quaternion to show where the link should be positioned for control

// The quantities below are of interest only for the internal workings of the dynamics computations.
// They are initialized when computed
vector<Matrix3d> Q(numlinks);         // useful for slowband inbound
vector<Matrix3d> W(numlinks);         // useful for slowband inbound
vector<Matrix3d> T(numlinks);         // useful matrix
vector<Matrix3d> K(numlinks);         // recursive coefficient
vector<Vector3d> d(numlinks);         // recursive coefficient
vector<Matrix3d> M(numlinks);         // recursive coefficient
vector<Vector3d> fprime(numlinks);    // recursive coefficient
vector<Vector3d> fdblprime(numlinks); // above in parent coordinates
vector<Matrix3d> munitm(numlinks);    // mass of link times the unit matrix
vector<Vector3d> mc(numlinks);        // mass of link times vector to centre of mass
vector<Matrix3d> mctilde(numlinks);   // mass times tilde matrix of vector to centre of mass
vector<Matrix3d> ltilde(numlinks);    // tilde of l
vector<Vector3d> aC(numlinks);        // centripetal acceleration in frame of the parent hinge
vector<Vector3d> maG(numlinks);       // force of gravity on the link
vector<Vector3d> fsigma(numlinks);    // a sum of forces useful in slowband
vector<Vector3d> gsigma(numlinks);    // a sum of torques useful in slowband
vector<Vector3d> g1sigma(numlinks);   // a sum of torques useful in fastband
vector<Vector3d> deltarot(numlinks);  // accumulator for the incremental rotation vector
Vector3d         zerov;               // the 3d zero vector
Matrix3d         zerom;               // zero matrix


void initialize()
{

    // When setting up the proximal hinge locations think of all the frames as the inertial one
    // where x is into the right, y up and z pointing out of the screen.
    // All cylinder axes of meshes initially point up and the proximal hinge is at the bottom end.
    // For Qt3D graphics, the centre of the cylinder is the position of the cylinder, which
    // we take as the centre of mass at c[r]. Links are rotated into position prior to starting simulation.
    p[0] = Vector3d(0.0, 22.0, 0.0); // this is  the initial position of link 0.
    // The falling linkage should be stopped by the floor by the application of appropriate forces.
    // Proximal hinge positions could be updated from v[0] and a[0] by integration,
    // but to prevent errors accumulating, it is done here by adjustLinks()
    // based on p[0] and rotations qR[r] of all links which, inverted, change the representation to the inertial frame.
    // adjustLinks() computes all other quaternions and all rotation matrices.

    l[0] = Vector3d(0,0,0); // This value is never used
    l[1] = Vector3d(0,0,-5.0); // The link vectors are in the frame of the PARENT(!)
    l[2] = Vector3d(0,0, 5.0);
    l[3] = Vector3d(0,20.0,-5.0);
    l[4] = Vector3d(0,20.0, 5.0);
    c[0] = Vector3d(0,length[0]/2.0,0);
    c[1] = Vector3d(0,length[1]/2.0,0);
    c[2] = Vector3d(0,length[2]/2.0,0);
    c[3] = Vector3d(0,length[3]/2.0,0);
    c[4] = Vector3d(0,length[4]/2.0,0);

    // Initialize quaternions
    // Rotate the legs ... degrees counter-clockwise in the frame of the PARENT(!)
    qR[1] = qR[2] = qR[3] = qR[4] = Quaterniond(s7071, 0.0, 0.0,s7071);
    // Rotate link 0, the body, ... degrees counter-clockwise in the inertial frame
    qR[0] = Quaternion<double>(s7071, 0.0f, 0.0f, s7071 );
    // This should rotate the body plus limbs together.
    // Changing the latter s7071 to -s7071 turns the "animal" upside down.

    // adjustLinks() computes all quaternions besides the qR's, then
    // derives from them all rotation matrices and the positions
    // (in the inertial frame) of the proximal hinges p[r]
    adjustLinks();

    // The following loop initializes everything with default values, which need to be changed for specific cases
    for( int r = 0; r < numlinks; r++)
    {
        // temporarily set the desired orientation of the child links to the original one
        qRdesired[r] = qR[r];
        // the distance vector from a cylinder end at the centre to a likely point of
        // contact with the floor
        downRadius[r] = zerov;
        downRadius[r](1) = - radius[r];
        int nRows, nr, nCols, nc;
        nRows = nCols = 3;
        //initialize the (column) vectors
        for(nr = 0; nr < nRows; ++nr)
        {
            // set masses times the vector to the c.o.g of the links
            mc[r](nr) = m[r] * c[r](nr);
            // motion variables
            v[r](nr) = f[r](nr) = fEdistal[r](nr) = fEproximal[r](nr)
                     = a[r](nr) = aC[r](nr) =  omega[r](nr) = omegadot[0](nr)
                     = g[r](nr) = gE[r](nr) =  gsigma[r](nr) = g1sigma[r](nr)  = 0;
            maG[r](nr) = m[r] * aG(nr);
            // Now the vector recursive coefficients for angular acceleration and force on link r
            d[r](nr) = fprime[r](nr) = fdblprime[r](nr) = zerov(nr) = 0;
        }
        // initialize some matrices for the links including the inertia matrix J
        // The moment of inertia of a cylindrical shell about one end is given by
        double sym = m[r]*(0.333333*length[r]*length[r] + 0.5 * (double)radius[r] * (double)radius[r]);
        // for (axes x and z in our case)
        // and
        double nonsym = m[r] * radius[r] * radius[r];
        // for the axis of symmetry (y = 1 in our case).  We have had to use hollow spheres to increase stability
        // but this is temporary. Using the correct values gets you lots of nans in the output file.
        for(nr = 0; nr < nRows; ++nr)
        {
            for(nc = 0; nc < nCols; ++nc)
            {
                // J[r](nr,nc)= (nr==nc)?((nr == 1)?sym:nonsym):0.0; we replace this by the moment of interia of a hollow sphere
                J[r](nr,nc)= (nr==nc)?0.6666 * m[r] * 100.0:0.0; // this should be replaced later by the above
                munitm[r](nr,nc) = (nr == nc)?m[r]:0;
                // Now the matrix recursive coefficients for angular acceleration and force on link r
                // This can be changed later
                M[r](nr,nc) = K[r](nr,nc) = Q[r](nr, nc) = W[r](nr,nc) = T[r](nr,nc) = zerom(nr,nc) = 0;
            }
        }
        // initialize more matrices
        mctilde[r] << 0,      -mc[r](2), mc[r](1),
                      mc[r](2), 0,       -mc[r](0),
                     -mc[r](1), mc[r](0), 0;
        // initialize ltilde
        ltilde[r] << 0,      -l[r](2),  l[r](1),
                     l[r](2), 0,       -l[r](0),
                    -l[r](1), l[r](0),  0;
        // the recursive coefficients of links r which have no children, i.e. the extremities,
        // are computed as follows. The quantities are constants.
        // If the links have children, the values will be recalculated in solving the equations of motion
        // so this initialization will be overridden.
        if(isParent[r] == 0) // we do the following for the extremities
        {
            T[r] = J[r].inverse();
            K[r] = -J[r].inverse() * mctilde[r];
            M[r] = - mctilde[r] * J[r].inverse() * mctilde[r] - munitm[r];
        }
        // Important quantities for computing the floor forces on the linkage
        pEproximal[r] = p[r] + downRadius[r];
        pEdistal[r] = pEproximal[r] + 2.0 * RI[r] * c[r];
    }
    m_min = 250.0; // these two values are bounds on the mass of the links and the linkage
    m_max = 14000.0; // they are used in computing floor force

}

/*  POSSIBLE VALUES FOR TESTING
 *  some random forces are input
 *
    What follows are possible values for testing with a call to adjust links at the end
    random_device rndDevice;
    mt19937 eng(rndDevice());
    uniform_int_distribution<int> dist(1,99);
    double scale = 0.01;
    for(int r = 0; r < numlinks; r++)
    {
       for(int nr = 0; nr < 3; nr++)
       {
            // we set all the external forces and control torques randomly
            // to see if the equations are being solved
            fEdistal[r](nr) = scale * (double) dist(eng);
            fEproximal[r](nr) = scale * (double) dist(eng);
            omega[r](nr) = scale * (double) dist(eng);
            g[r](nr) = scale * (double) dist(eng);
            gE[r](nr) = scale * (double) dist(eng);
       }
    }

    // Experiments for understanding quaternions
    Quaterniond quat1 = Quaterniond(4.2, -3.6, 9.2, 4.8);
    Quaterniond quat2 = Quaterniond(-7.4, 2.5, 119.0, -6);
    quat1.normalize();
    quat2.normalize();
    Matrix3d m1 = quat1.toRotationMatrix();
    Matrix3d m2 = quat2.toRotationMatrix();
    Matrix3d m12 = m1 * m2;
    Matrix3d m21 = m2 * m1;
    Quaterniond quat12 = quat1*quat2;
    Quaterniond quat21 = quat2*quat1;
    quat12.normalize();
    quat21.normalize();
    Matrix3d mq12 = quat12.toRotationMatrix();
    Matrix3d mq21 = quat21.toRotationMatrix();
    // Now use the debugger and stop below.
    // Conclusion: the quaternions are multiplied in the *same* order
    // as the rotation matrices.

    // aG = Vector3d(0.0, 0.0 , 0.0); This must be changed ABOVE because other things depend on it
    //omega[0] = Vector3d(0.0, 0.0, 3.1416e-2);
    //omega[1] = Vector3d(0.0, 0.0, 3.1416e-2);
    //omega[2] = Vector3d(0.0, 0.0, 3.1416e-2);
    //omega[3] = Vector3d(0.0, 0.0, 3.1416e-2);
    //omega[4] = Vector3d(0.0, 0.0, 3.1416e-2);
    //qR[0] = Quaternion<double>(.7071, 0.0, 0.0, .7071); // to test the floor only
    //qR[0].normalize();
    //p[0]  = Vector3d(0.0, 0.0, 0.0);
    // fEdistal[1](0) = 10.0; fEdistal[3](0) = 10.0; // inject some rotation to check momentum
    END OF LIST OF POSSIBLE VALUES FOR TESTING
*/



