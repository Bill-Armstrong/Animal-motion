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
// file Dynamics.cpp
#include <Eigen/Dense>
#include <vector>
#include "Dynamics.h"
//#include <iostream>
#include <fstream>
using namespace std;

Matrix3d accuWltilde;
Matrix3d accuW;
Matrix3d accuQ;
Matrix3d accuQltilde;
Vector3d accuRg;

void slowband_in();
void fastband_in();
void integration();
void floorSupport();
Vector3d floor_force(Vector3d p, Vector3d v, Vector3d p1, Vector3d p2,
                   Vector3d f1, Vector3d f2, double& m_eff);
void adjustLinks();
void linkControl();

std::ofstream ofs("link3.txt", std::ofstream::out);



void slowband_in()
{
    for( int r = numlinks -1; r >= 0; r--)
    {

        // initialize the accumulatorslowband_s
        accuW = accuQ = accuQltilde = zerom;
        accuRg = zerov;
        accuWltilde = J[r]; // start with this part of T[r]

        // Summing over the children of link r involves summing over all links and filtering for the r parent
        // Recall that all children s of link r have a higher index than r in the Parent
        for(int s = numlinks - 1; s > r ; s--)
        {
           if(isParent[r] == 0) continue; // we initialized the values for the non-parents
           if(Parent[s] == r)
           {
                aC[s] = omega[r].cross(omega[r].cross(l[s]));
                Q[s] = R[s] * M[s] * RT[s];
                W[s] = ltilde[s] * Q[s];
                accuQ += Q[s];
                accuQltilde += Q[s] * ltilde[s];
                accuW += W[s];
                accuWltilde += W[s]*ltilde[s];
                accuRg += R[s] * g[s];
           }
        }
        T[r] = accuWltilde.inverse();
        K[r] = T[r] * (accuW - mctilde[r]);
        M[r] = mctilde[r] * K[r] - munitm[r] + accuQ - accuQltilde * K[r];
        // Note: If link r has no children, then setting the sums over the non-existent children to zero we get
        // T[r] = J[r].inverse()
        // K[r] = -J[r].inverse() * mctilde[r]
        // M[r] = - mctilde[r] * J[r].inverse() * mctilde[r] - munitm[r]
        // This is what appears in initialize(), since all these values are constants for each link.

    }
}
    
void fastband_in()
{
    Vector3d accuRg;
    Vector3d acculfQa;
    Vector3d accufQaCld;
    for( int r = numlinks -1; r >= 0; --r)
    {
        if(isParent[r] == 0)
        {
            fdblprime[r] = zerov; // this needs to be initialized somewhere
            aC[r] = omega[Parent[r]].cross(omega[Parent[r]].cross(l[r]));
        }
        // We assume the external force is rapidly varying, e.g. when a link hits the floor
        // so the computations of g1sigma and fsigma have been movedinto the fastband
        // We have also chosen to use inertial coordinates for the positions of the external forces from the floor
        // since that is where it is easier to see what is going on
        g1sigma[r] = - omega[r].cross(J[r] * omega[r]) + RIT[r] * gE[r] +
                (2.0 * c[r] + RIT[r] * downRadius[r]).cross(RIT[r] * fEdistal[r])
                          + (RIT[r] * downRadius[r]).cross(RIT[r] * fEproximal[r])
                + mc[r].cross(RIT[r] * aG);
        fsigma[r] = - omega[r].cross(omega[r].cross(mc[r])) +
                    RIT[r] * (fEdistal[r] + fEproximal[r] + maG[r]);

        accuRg = acculfQa = zerov;
        for( int s = numlinks - 1; s > r ; s--)
        {
            if(Parent[s] == r)
            {
                accuRg += R[s] * g[s];
                acculfQa += l[s].cross(fdblprime[s] + Q[s] * aC[s]);
            }
        }
        gsigma[r] = g1sigma[r] - g[r] + accuRg;
        d[r] =T[r] * (gsigma[r] + acculfQa);
        accufQaCld = zerov;
        for( int s = numlinks - 1; s > r ; s--)
        {
            if(Parent[s] == r)
            {
                accufQaCld += fdblprime[s] + Q[s] * (aC[s] - l[s].cross(d[r]));
            }
        }
        fprime[r] = fsigma[r] + mc[r].cross(d[r]) + accufQaCld;
        fdblprime[r] = R[r] * fprime[r];
        // For links r which have no children, we just set all the accu... quantities to zero
        // and use T[r] from the slowband_in computation to get values for d[r] and fprime[r]
        // which, together with K[r] and M[r], complete the set of recursive coefficients
    }
}

void integration()
{

    // This takes care of updating the rotation quaternions and matrices of the links
    // as well as updating the velocities v[r] and omega[r] and position p[r]
    // adjustLinks takes care of the rest
    a[0] = - M[0].inverse() * fprime[0];
    omegadot[0] = K[0] * a[0] + d[0];
    // we assume accelerations are constant during the coming time step
    omega[0] += deltat * omegadot[0]; // omega, omegadot are all in the frame of link r
    deltau[0] = deltat * ( omega[0] - 0.5 * omegadot[0] * deltat); // analogous to p[0] done below
    v[0] += deltat * RI[0] * a[0]; // v is in the inertial frame
    p[0] += deltat * (v[0] - 0.5 * deltat * RI[0] * a[0]); // p is in the inertial frame; v[0] is AFTER update, so there is a minus sign
    f[0] = zerov; // link zero exerts 0 force on a non-existent parent

    for( int r = 1; r < numlinks; r++)
    {
        a[r] = RT[r] * (aC[r] + a[Parent[r]] - l[r].cross(omegadot[Parent[r]]));
        omegadot[r] = K[r] * a[r] + d[r];
        omega[r] += deltat * omegadot[r];
        deltau[r] = deltat * ( omega[r] - 0.5 * omegadot[r] * deltat); // analogous to p[0] done above
        // The following quantities may be needed (e.g.for checking or computing floor support)
        f[r] = M[r] * a[r] + fprime[r];
        v[r] += deltat * RI[r] * a[r]; // temporary: avoid accumulation of errors!
        vdistal[r] = v[r] + 2.0 *  RI[r]* (omega[r].cross(c[r]));
        // p[r] for r > 0 is computed using the l[r]'s and the rotation matrices,
        // thus avoiding accumulation of errors
     }
     // Updating the rotation quaternions of the links
    double magnitude, w, x, y, z;
    Quaterniond quaternionTemp1, quaternionTemp2;
    for(int r = 0; r < numlinks; r++)
    {
       // Angle of change of the link during the last time step
       magnitude = sqrt((deltau[r](0) * deltau[r](0)) + (deltau[r](1) * deltau[r](1)) + (deltau[r](2) * deltau[r](2)));
       if (magnitude <= 0.0)  // In this case, we don't update
       {
         w = 1.0; x = y = z = 0.0; // i.e. don't do any incremental rotation
       }
       else
       {
           w = cos(0.5 * magnitude);
           x = sin(0.5 * magnitude) * deltau[r](0)/magnitude;
           y = sin(0.5 * magnitude) * deltau[r](1)/magnitude;
           z = sin(0.5 * magnitude) * deltau[r](2)/magnitude;
       }
       quaternionTemp1 = Quaterniond(w, x, y ,z); // This is the incremental rotation on link r which is
       // now to be combined with the current quaternion qR[r]. Does the compiler optimize this?
       quaternionTemp2 = quaternionTemp1 * qR[r];
       qR[r] = quaternionTemp2;
    }
    // Now adjust the other quaternions, rotation matrices and positions
    adjustLinks();
}

void adjustLinks()
{
    // The inputs are the quaternions qR[r] for all links; all other rotation matrices are computed
    // The rotation matrices are needed for dynamics, and the qRI are used for the display transforms
    Vector3d velocityDiff;
    for(int r = 0; r < numlinks; r++)
    {
       qR[r].normalize();
       // Fix the rotation matrices, using the quaternions as basis
       if(r == 0)
       {
         qRI[0] = qR[0];
         R[0] = qR[0].toRotationMatrix();
         RT[0] = R[0].transpose();
         RI[0] = R[0];
         RIT[0] = RI[0].transpose();
       }
       else
       {
         R[r] = qR[r].toRotationMatrix();
         RT[r] = R[r].transpose();
         qRI[r] = qRI[Parent[r]]*qR[r];
         qRI[r].normalize(); // This may be unnecessary
         RI[r] = qRI[r].toRotationMatrix();
         RIT[r] = RI[r].transpose();
       }
       // Fix the proximal hinge positions, velocities and momenta in the inertial frame
       // The origin of the graphics is the reference point for positions.
       if(r > 0)
       {
           p[r] = p[Parent[r]] + RI[Parent[r]] * l[r]; // recall that l[r] is in the frame of the parent
           velocityDiff = -v[r];
           v[r] = v[Parent[r]] + RI[Parent[r]] * (omega[Parent[r]].cross(l[r]));
           velocityDiff += v[r];
       }
    }
}


void check()
{
    vector<Vector3d> accuRg(numlinks);
    vector<Vector3d> acculRf(numlinks);
    vector<Vector3d> accuRf(numlinks);
    vector<Vector3d> amDiff(numlinks); // error in angular momentum in solving the equations of motion
    vector<Vector3d> forceDiff(numlinks); // error in hinge force in solving the equations of motion
    vector<Vector3d> accelDiff(numlinks); // error hinge in acceleration for the proximal hinge of link s

    for(int r = 0; r < numlinks; r++)
    {
        accuRg[r]  = zerov;
        acculRf[r] = zerov;
        accuRf[r]  = zerov;
        for(int s = numlinks - 1; s > r; s--)
        {
            if(Parent[s] == r)
            {
                acculRf[r] += l[s].cross(R[s] * f[s]);
                accuRg[r]  += R[s] * g[s];
                accuRf[r]  += R[s] * f[s];
                accelDiff[s] = omega[r].cross(omega[r].cross(l[s])) + a[r] - l[s].cross(omegadot[r]);
                accelDiff[s] += - R[s] * a[s]; // Equation (5) in Visual Computer 1985 : 1
            }
        }
        amDiff[r] = -omega[r].cross(J[r] * omega[r])- g[r];
        amDiff[r] += accuRg[r] + RIT[r] * gE[r];
        amDiff[r] += mc[r].cross(RIT[r] * aG)
                     + (2.0* c[r] + RIT[r] * downRadius[r]).cross(RIT[r] * fEdistal[r])
                                  + (RIT[r] * downRadius[r]).cross(RIT[r] * fEproximal[r]);
        // up to this point, we have computed gsigma[r]
        amDiff[r] += -mc[r].cross(a[r]) + acculRf[r];
        amDiff[r] -= J[r] * omegadot[r] ; // Equation (1) in Visual Computer 1985 : 1
        forceDiff[r] = -m[r] * omega[r].cross(omega[r].cross(c[r]));
        forceDiff[r] += RIT[r]*(fEdistal[r] + fEproximal[r] + maG[r]);
        // up to here we have the fsigma[r] part
        forceDiff[r] += - m[r] * a[r] + mc[r].cross(omegadot[r]) + accuRf[r];
        forceDiff[r] -= f[r]; // Equatio0n (3) in Visual Computer 1985 : 1

        // Examine the amDiff, forceDiff and accelDiff in the debugger
        // while pushing the f5 key.
        // All of them should keep very close to zero.
    }
}

char dummy; // for diagnosis
void floorSupport()
{
    double masspenetration; // this divides the floor support proportionally to the mass of the link times its penetration
    double forceperunit; // the total upward force (not counting elasticity and viscosity) is a multiple of the weicht
    masspenetration = 0.0;
    forceperunit =  0.28 * -14000.0*aG[1]; // the total force applied in the floor is a fraction of what
    // is required to support the weight of the whole object.  This allows us to attenuate elasticity.
    // This force is applied to the links in the floor, and depends on their masses and penetrations
    // The rest of support comes from elasticity.

    // aG[1] is negative, so forceperunit is positive
    for(int r = 0; r < numlinks; r++)
    {
        // The ends of the cylinders are rounded with spheres, and the bottom of the sphere is the point of application of the external force.
        pEproximal[r] = p[r] + downRadius[r]; // a force is to be located near the proximal hinge on the links (pE's are in the INERTIAL FRAME!)
        pEdistal[r] =   pEproximal[r] + 2.0 * RI[r] * c[r]; // a force is to be located near the most distal point on the links.
        // Properties of the floor generating external forces can be added below
        // e.g. elasticity, viscosity, etc.  The floor must be "bouncy" since the
        // recursive solution method works only for tree-linkages, i.e. acyclic graphs.
        // We weight the forces applied to the links in proportion to the product of their mass and penetration
        if(pEproximal[r](1) < 0.0) masspenetration += m[r] * pEproximal[r](1);
        if(pEdistal[r](1) < 0.0) masspenetration += m[r] * pEdistal[r](1);
    }

    if(masspenetration > -0.01)// masspenetration is always a negative or zero quantity
    {
        // almost nothing has penetrated the floor: apply zero force
        // this avoids division by penetration == zero
        for(int r = 0; r < numlinks; r++)
        {
            fEproximal[r] = -air_resistance * v[r];
            fEdistal[r] =  -air_resistance * v[r];
        }
    }
    else
    {
        // We have (significant) penetration of the floor
        for(int r = 0; r < numlinks; r++)
        {
            if(pEproximal[r](1) < 0.0)
            {
                fEproximal[r] = - v[r] * viscosity;
                // The following adds enough upward force to push the link back to the surface
                fEproximal[r](1) -= pEproximal[r](1) * elasticity;
                // The following adds enough upward force to support the weight
                fEproximal[r](1) += forceperunit * pEproximal[r](1) * m[r]/masspenetration;
            }
            else
            {
                fEproximal[r] = -air_resistance * v[r];
            }

            if(pEdistal[r](1) < 0.0)
            {
                fEdistal[r] = -(v[r] + 2.0 * RI[r] * (omega[r].cross(c[r]))) * viscosity; // blows up if omega gets huge!
                // The following adds enough upward force to push the link back to the surface
                fEdistal[r](1) -= pEdistal[r](1) * elasticity;
                // The following adds enough upward force to support the weight
                fEdistal[r](1) += forceperunit * pEdistal[r](1) * m[r]/masspenetration;
            }
            else
            {
                fEdistal[r] =  -air_resistance * (v[r] + 2.0 * RI[r] * (omega[r].cross(c[r])));
            }
        } // end for r
    } // end if penetration > -0.01 ... else
    if(int r = 3)
    {
        //#include <fstream> for output file vs cout
        //std::ofstream ofs("bounce.output", std::ofstream::out);
        // For writing to a file: Change cout to ofs
        // ofs << std::setfill(' ') <<
        ofs << SimTime
        << " fEproximal[" << r << "] " << fEproximal[r](0) << " " << fEproximal[r](1) << " "<< fEproximal[r](2) << " "
        << " fEdistal[" << r << "] " << fEdistal[r](0) << " " << fEdistal[r](1) << " "<< fEdistal[r](2) << " "
        << " p[" << r << "] " << p[r](0) << " " <<  p[r](1) << " " << p[r](2)
        << " a[" << r << "] " << (RI[r] * a[r])(0) << " " << (RI[r] * a[r])(1) << " " << (RI[r] * a[r])(2)
        << " v[" << r << "] " << (RI[r] * v[r])(0) << " " << (RI[r] * v[r])(1) << " " << (RI[r] * v[r])(2)
        << " o[" << r << "] " << (RI[r] * omega[r])(0) << " " << (RI[r] * omega[r])(1) << " " << (RI[r] * omega[r])(2) << std::endl;
    }
}
/*
Vector3d floor_force(Vector3d p, Vector3d v, Vector3d p1, Vector3d p2,
                      Vector3d f1, Vector3d f2, double& m_eff)
// Computes the force for the current timestep of duration tau
// p is the position at the start of the timestep.
// p1 is the position one timestep ago and
// p2 is the position two timesteps ago.
// f1 was applied during the last timestep and
// f2 during the one before that.
{
    Vector3d f_eff;               // effective force in the middle interval of duration tau in the past two timesteps
    Vector3d accel;               // (usually) upward acceleration caused by f_eff which has been added to gravity
    Vector3d f;                   // force value to be returned

    f = zerov; // We often just change the y component of this f[1] below before returning f.
    m_eff =  0.0;
    f_eff = zerov;

    // This routine has access to p, p1, p2, f1, and f2 and may need to store m_eff for smoothing purposes later.
    if(p[1] <= 0.0)
    {
        // the mass is below floor level 0, has an adequate test force been applied?
        if(0.5 * (f1[1]+f2[1]) < -m_min * aG[1])
        {
            // an insufficient test force was applied during the last two timesteps, so
            // we add an upward test force to zerov for the current timestep (aG[1] is negative)
            f[1] = -m_min * aG[1];
        }
        else
        {
            // The mass has entered the floor and an adequate upward force
            // has been applied to try to stop it.
            // The velocity has gone from (p1-p2)/tau to (p-p1)/tau in time tau
            // in the middle of the last two timesteps,
            // The part of the acceleration *not* caused by gravity was
            accel = (p + p2 - 2 * p1)/deltat_sq - aG; // acceleration due to the forces f1, f2
            f_eff = 0.5 * (f1 + f2);               // effective force causing that acceleration
            if(accel[1] < 50.0) // must be a sizeable fraction of the acceleration of gravity in magnitude
            {
                // We need the applied force to add a measurable upward acceleration to gravity,
                // but not enough to cancel gravity completely (then we get a division by zero), so
                // we apply an upward test force for the next timestep. Thus, we will
                // have less negative acceleration for this link than the acceleration of gravity
                //(which could happen when a link is pushed from above by another link).
                // We continue adding the minimal force to f.
                f[1] = -m_min * aG[1];
            }
            else
            {
                // We compute an effective mass (which is correct for one body, but maybe not for a linkage).
                // For a linkage, we may have to smooth this over several timesteps while p[1] < 0.0.
                m_eff =  0.25* m_eff + 0.75 * f_eff[1] / accel[1];
                // For the next time step we apply a sum of forces to
                // 1. counteract gravity,
                f = - m_eff * aG;
                // 2. push the mass upward by an elastic force
                //f[1] -= p[1] * elasticity;
                // and 3. slow the mass enough to counteract the current velocity
                // by applying a force so momentum is removed (in all directions)
                // The 0.5 takes the momentum away over 2 timesteps.
                f = -m_eff * v/deltat;
                // We also add viscosity to hold it back
                //f -= viscosity * v;
            }
        }
        f= zerov;
        f[1] = 3432450.0;  // this works !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    else
    {
        // the link end is not under floor level
        f = zerov;
        m_eff = 0.0;
    }
    //double atest, btest, ctest;
    //atest =f[0]; btest = f[1]; ctest = f[2];
    return f;
}


// Update the remembered inputs to floor force p2E <- p1E, p1E <- pE, f2E <- f1E, f1E <- fE
void update_floor_force_inputs()
{
    for(int r = 0; r < numlinks; r++)
    {
       p2Eproximal[r] = p1Eproximal[r];
       p2Edistal[r]   = p1Edistal[r];
       p1Eproximal[r] = pEproximal[r];
       p1Edistal[r]   = pEdistal[r];
       f2Eproximal[r] = f1Eproximal[r];
       f2Edistal[r]   = f1Edistal[r];
       f1Eproximal[r] = fEproximal[r];
       f1Edistal[r]   = fEdistal[r];
    }
}
*/

Quaterniond rotDev;
Vector3d axis;
void linkControl()
{
    // Question: can we just set the joint torque using g[]???
    g[0] = zerov;
    for(int r = 1; r < numlinks; r++) // only links which have a parent generate a torque
    {
        rotDev = qR[r]*qRdesired[r].inverse(); // the rotDev (rotation deviation) is the difference from desired
        // acos returns an angle between 0 and pi radians.
        rotDev.normalize();
        if(rotDev.w() > 0.9999)
        {
            g[r] = zerov; // There is very little deviation from the desired orientation
        }
        else if(rotDev.w() < - 0.9999)
        {
            g[r] = zerov; // Maximal 180 degrees deviation, and we choose a torque in the y direction
            g[r](1) = -rotElastic * 3.1415926;
        }
        else //
        {
           axis = Vector3d(rotDev.x(), rotDev.y(), rotDev.z());
           g[r] = -rotElastic * acos(rotDev.w()) * RI[r] * axis.normalized();
        }
        // We introduce rotational damping as an external torque to hinder
        // a link spinning. Mostly the problem is exponential speedup along the cylinder axis.
        // We don't know what causes the speedup, so the rotDamp has to be set experimentally.
        g[r] -= rotDamp * RI[r] * omega[r];
    }
}

// Comments on floor support
//The total mass of the linkage is 14000 grams, so the weight on each foot is 3500 grams force, and
// that times aG gives 3432450 dynes on each foot needed to support if all is balanced.
// So if the foot has penetrated 1 cm, the force upward should be of that magnitude
// Hence the elasticity should be 343245 dynes/cm to stop the links from going deeper than 10 cm.
// There are two external forces on each link, one at the proximal end, fEproximal[r],
// and one at the distal end fEdistal[r]. These are applied at positions pEproximal[r] and pEdistal[r],
// where the fE and pE are in inertial coordinates.
// If one wishes to have other external forces, these could be added into the program
// at places where fEdistal and pEdistal etc. are mentioned
// The floor must be "bouncy" since the
// recursive solution method works only for tree-linkages, i.e. acyclic graphs.
// Our next idea is to apply an average force between the old elastic force and the new elastic force
// over successive time intervals.

