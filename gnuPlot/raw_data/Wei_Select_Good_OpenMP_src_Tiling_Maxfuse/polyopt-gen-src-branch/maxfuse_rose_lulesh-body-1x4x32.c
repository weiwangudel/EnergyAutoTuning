/*

                 Copyright (c) 2010.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 1.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lulesh.h"
#include "lulesh-others.c"
#define LULESH_SHOW_PROGRESS 0 
/* Stuff needed for boundary conditions */
/* 2 BCs on each of 6 hexahedral faces (12 bits) */

void VoluDer(const Real_t x0,const Real_t x1,const Real_t x2,const Real_t x3,const Real_t x4,const Real_t x5,const Real_t y0,const Real_t y1,const Real_t y2,const Real_t y3,const Real_t y4,const Real_t y5,const Real_t z0,const Real_t z1,const Real_t z2,const Real_t z3,const Real_t z4,const Real_t z5,Real_t *dvdx,Real_t *dvdy,Real_t *dvdz)
{
  const Real_t twelfth = (1.0 / 12.0);
   *dvdx = (((((((y1 + y2) * (z0 + z1)) - ((y0 + y1) * (z1 + z2))) + ((y0 + y4) * (z3 + z4))) - ((y3 + y4) * (z0 + z4))) - ((y2 + y5) * (z3 + z5))) + ((y3 + y5) * (z2 + z5)));
   *dvdy = ((((((-(x1 + x2) * (z0 + z1)) + ((x0 + x1) * (z1 + z2))) - ((x0 + x4) * (z3 + z4))) + ((x3 + x4) * (z0 + z4))) + ((x2 + x5) * (z3 + z5))) - ((x3 + x5) * (z2 + z5)));
   *dvdz = ((((((-(y1 + y2) * (x0 + x1)) + ((y0 + y1) * (x1 + x2))) - ((y0 + y4) * (x3 + x4))) + ((y3 + y4) * (x0 + x4))) + ((y2 + y5) * (x3 + x5))) - ((y3 + y5) * (x2 + x5)));
   *dvdx *= twelfth;
   *dvdy *= twelfth;
   *dvdz *= twelfth;
}

void CalcElemVolumeDerivative(Real_t dvdx[8UL],Real_t dvdy[8UL],Real_t dvdz[8UL],const Real_t x[8UL],const Real_t y[8UL],const Real_t z[8UL])
{
  VoluDer(x[1],x[2],x[3],x[4],x[5],x[7],y[1],y[2],y[3],y[4],y[5],y[7],z[1],z[2],z[3],z[4],z[5],z[7],(dvdx + 0),(dvdy + 0),(dvdz + 0));
  VoluDer(x[0],x[1],x[2],x[7],x[4],x[6],y[0],y[1],y[2],y[7],y[4],y[6],z[0],z[1],z[2],z[7],z[4],z[6],(dvdx + 3),(dvdy + 3),(dvdz + 3));
  VoluDer(x[3],x[0],x[1],x[6],x[7],x[5],y[3],y[0],y[1],y[6],y[7],y[5],z[3],z[0],z[1],z[6],z[7],z[5],(dvdx + 2),(dvdy + 2),(dvdz + 2));
  VoluDer(x[2],x[3],x[0],x[5],x[6],x[4],y[2],y[3],y[0],y[5],y[6],y[4],z[2],z[3],z[0],z[5],z[6],z[4],(dvdx + 1),(dvdy + 1),(dvdz + 1));
  VoluDer(x[7],x[6],x[5],x[0],x[3],x[1],y[7],y[6],y[5],y[0],y[3],y[1],z[7],z[6],z[5],z[0],z[3],z[1],(dvdx + 4),(dvdy + 4),(dvdz + 4));
  VoluDer(x[4],x[7],x[6],x[1],x[0],x[2],y[4],y[7],y[6],y[1],y[0],y[2],z[4],z[7],z[6],z[1],z[0],z[2],(dvdx + 5),(dvdy + 5),(dvdz + 5));
  VoluDer(x[5],x[4],x[7],x[2],x[1],x[3],y[5],y[4],y[7],y[2],y[1],y[3],z[5],z[4],z[7],z[2],z[1],z[3],(dvdx + 6),(dvdy + 6),(dvdz + 6));
  VoluDer(x[6],x[5],x[4],x[3],x[2],x[0],y[6],y[5],y[4],y[3],y[2],y[0],z[6],z[5],z[4],z[3],z[2],z[0],(dvdx + 7),(dvdy + 7),(dvdz + 7));
}

void CalcElemFBHourglassForce(Real_t *xd,Real_t *yd,Real_t *zd,Real_t *hourgam0,Real_t *hourgam1,Real_t *hourgam2,Real_t *hourgam3,Real_t *hourgam4,Real_t *hourgam5,Real_t *hourgam6,Real_t *hourgam7,Real_t coefficient,Real_t *hgfx,Real_t *hgfy,Real_t *hgfz)
{
  Index_t i00 = 0;
  Index_t i01 = 1;
  Index_t i02 = 2;
  Index_t i03 = 3;
  Real_t h00 = ((((((((hourgam0[i00] * xd[0]) + (hourgam1[i00] * xd[1])) + (hourgam2[i00] * xd[2])) + (hourgam3[i00] * xd[3])) + (hourgam4[i00] * xd[4])) + (hourgam5[i00] * xd[5])) + (hourgam6[i00] * xd[6])) + (hourgam7[i00] * xd[7]));
  Real_t h01 = ((((((((hourgam0[i01] * xd[0]) + (hourgam1[i01] * xd[1])) + (hourgam2[i01] * xd[2])) + (hourgam3[i01] * xd[3])) + (hourgam4[i01] * xd[4])) + (hourgam5[i01] * xd[5])) + (hourgam6[i01] * xd[6])) + (hourgam7[i01] * xd[7]));
  Real_t h02 = ((((((((hourgam0[i02] * xd[0]) + (hourgam1[i02] * xd[1])) + (hourgam2[i02] * xd[2])) + (hourgam3[i02] * xd[3])) + (hourgam4[i02] * xd[4])) + (hourgam5[i02] * xd[5])) + (hourgam6[i02] * xd[6])) + (hourgam7[i02] * xd[7]));
  Real_t h03 = ((((((((hourgam0[i03] * xd[0]) + (hourgam1[i03] * xd[1])) + (hourgam2[i03] * xd[2])) + (hourgam3[i03] * xd[3])) + (hourgam4[i03] * xd[4])) + (hourgam5[i03] * xd[5])) + (hourgam6[i03] * xd[6])) + (hourgam7[i03] * xd[7]));
  hgfx[0] = (coefficient * ((((hourgam0[i00] * h00) + (hourgam0[i01] * h01)) + (hourgam0[i02] * h02)) + (hourgam0[i03] * h03)));
  hgfx[1] = (coefficient * ((((hourgam1[i00] * h00) + (hourgam1[i01] * h01)) + (hourgam1[i02] * h02)) + (hourgam1[i03] * h03)));
  hgfx[2] = (coefficient * ((((hourgam2[i00] * h00) + (hourgam2[i01] * h01)) + (hourgam2[i02] * h02)) + (hourgam2[i03] * h03)));
  hgfx[3] = (coefficient * ((((hourgam3[i00] * h00) + (hourgam3[i01] * h01)) + (hourgam3[i02] * h02)) + (hourgam3[i03] * h03)));
  hgfx[4] = (coefficient * ((((hourgam4[i00] * h00) + (hourgam4[i01] * h01)) + (hourgam4[i02] * h02)) + (hourgam4[i03] * h03)));
  hgfx[5] = (coefficient * ((((hourgam5[i00] * h00) + (hourgam5[i01] * h01)) + (hourgam5[i02] * h02)) + (hourgam5[i03] * h03)));
  hgfx[6] = (coefficient * ((((hourgam6[i00] * h00) + (hourgam6[i01] * h01)) + (hourgam6[i02] * h02)) + (hourgam6[i03] * h03)));
  hgfx[7] = (coefficient * ((((hourgam7[i00] * h00) + (hourgam7[i01] * h01)) + (hourgam7[i02] * h02)) + (hourgam7[i03] * h03)));
  h00 = ((((((((hourgam0[i00] * yd[0]) + (hourgam1[i00] * yd[1])) + (hourgam2[i00] * yd[2])) + (hourgam3[i00] * yd[3])) + (hourgam4[i00] * yd[4])) + (hourgam5[i00] * yd[5])) + (hourgam6[i00] * yd[6])) + (hourgam7[i00] * yd[7]));
  h01 = ((((((((hourgam0[i01] * yd[0]) + (hourgam1[i01] * yd[1])) + (hourgam2[i01] * yd[2])) + (hourgam3[i01] * yd[3])) + (hourgam4[i01] * yd[4])) + (hourgam5[i01] * yd[5])) + (hourgam6[i01] * yd[6])) + (hourgam7[i01] * yd[7]));
  h02 = ((((((((hourgam0[i02] * yd[0]) + (hourgam1[i02] * yd[1])) + (hourgam2[i02] * yd[2])) + (hourgam3[i02] * yd[3])) + (hourgam4[i02] * yd[4])) + (hourgam5[i02] * yd[5])) + (hourgam6[i02] * yd[6])) + (hourgam7[i02] * yd[7]));
  h03 = ((((((((hourgam0[i03] * yd[0]) + (hourgam1[i03] * yd[1])) + (hourgam2[i03] * yd[2])) + (hourgam3[i03] * yd[3])) + (hourgam4[i03] * yd[4])) + (hourgam5[i03] * yd[5])) + (hourgam6[i03] * yd[6])) + (hourgam7[i03] * yd[7]));
  hgfy[0] = (coefficient * ((((hourgam0[i00] * h00) + (hourgam0[i01] * h01)) + (hourgam0[i02] * h02)) + (hourgam0[i03] * h03)));
  hgfy[1] = (coefficient * ((((hourgam1[i00] * h00) + (hourgam1[i01] * h01)) + (hourgam1[i02] * h02)) + (hourgam1[i03] * h03)));
  hgfy[2] = (coefficient * ((((hourgam2[i00] * h00) + (hourgam2[i01] * h01)) + (hourgam2[i02] * h02)) + (hourgam2[i03] * h03)));
  hgfy[3] = (coefficient * ((((hourgam3[i00] * h00) + (hourgam3[i01] * h01)) + (hourgam3[i02] * h02)) + (hourgam3[i03] * h03)));
  hgfy[4] = (coefficient * ((((hourgam4[i00] * h00) + (hourgam4[i01] * h01)) + (hourgam4[i02] * h02)) + (hourgam4[i03] * h03)));
  hgfy[5] = (coefficient * ((((hourgam5[i00] * h00) + (hourgam5[i01] * h01)) + (hourgam5[i02] * h02)) + (hourgam5[i03] * h03)));
  hgfy[6] = (coefficient * ((((hourgam6[i00] * h00) + (hourgam6[i01] * h01)) + (hourgam6[i02] * h02)) + (hourgam6[i03] * h03)));
  hgfy[7] = (coefficient * ((((hourgam7[i00] * h00) + (hourgam7[i01] * h01)) + (hourgam7[i02] * h02)) + (hourgam7[i03] * h03)));
  h00 = ((((((((hourgam0[i00] * zd[0]) + (hourgam1[i00] * zd[1])) + (hourgam2[i00] * zd[2])) + (hourgam3[i00] * zd[3])) + (hourgam4[i00] * zd[4])) + (hourgam5[i00] * zd[5])) + (hourgam6[i00] * zd[6])) + (hourgam7[i00] * zd[7]));
  h01 = ((((((((hourgam0[i01] * zd[0]) + (hourgam1[i01] * zd[1])) + (hourgam2[i01] * zd[2])) + (hourgam3[i01] * zd[3])) + (hourgam4[i01] * zd[4])) + (hourgam5[i01] * zd[5])) + (hourgam6[i01] * zd[6])) + (hourgam7[i01] * zd[7]));
  h02 = ((((((((hourgam0[i02] * zd[0]) + (hourgam1[i02] * zd[1])) + (hourgam2[i02] * zd[2])) + (hourgam3[i02] * zd[3])) + (hourgam4[i02] * zd[4])) + (hourgam5[i02] * zd[5])) + (hourgam6[i02] * zd[6])) + (hourgam7[i02] * zd[7]));
  h03 = ((((((((hourgam0[i03] * zd[0]) + (hourgam1[i03] * zd[1])) + (hourgam2[i03] * zd[2])) + (hourgam3[i03] * zd[3])) + (hourgam4[i03] * zd[4])) + (hourgam5[i03] * zd[5])) + (hourgam6[i03] * zd[6])) + (hourgam7[i03] * zd[7]));
  hgfz[0] = (coefficient * ((((hourgam0[i00] * h00) + (hourgam0[i01] * h01)) + (hourgam0[i02] * h02)) + (hourgam0[i03] * h03)));
  hgfz[1] = (coefficient * ((((hourgam1[i00] * h00) + (hourgam1[i01] * h01)) + (hourgam1[i02] * h02)) + (hourgam1[i03] * h03)));
  hgfz[2] = (coefficient * ((((hourgam2[i00] * h00) + (hourgam2[i01] * h01)) + (hourgam2[i02] * h02)) + (hourgam2[i03] * h03)));
  hgfz[3] = (coefficient * ((((hourgam3[i00] * h00) + (hourgam3[i01] * h01)) + (hourgam3[i02] * h02)) + (hourgam3[i03] * h03)));
  hgfz[4] = (coefficient * ((((hourgam4[i00] * h00) + (hourgam4[i01] * h01)) + (hourgam4[i02] * h02)) + (hourgam4[i03] * h03)));
  hgfz[5] = (coefficient * ((((hourgam5[i00] * h00) + (hourgam5[i01] * h01)) + (hourgam5[i02] * h02)) + (hourgam5[i03] * h03)));
  hgfz[6] = (coefficient * ((((hourgam6[i00] * h00) + (hourgam6[i01] * h01)) + (hourgam6[i02] * h02)) + (hourgam6[i03] * h03)));
  hgfz[7] = (coefficient * ((((hourgam7[i00] * h00) + (hourgam7[i01] * h01)) + (hourgam7[i02] * h02)) + (hourgam7[i03] * h03)));
}

void CalcFBHourglassForceForElems(Real_t *determ,Real_t *x8n,Real_t *y8n,Real_t *z8n,Real_t *dvdx,Real_t *dvdy,Real_t *dvdz,Real_t hourg)
{
/*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/
  Index_t numElem = m_numElem;
  Index_t numElem8 = (numElem * 8);
  Real_t *fx_elem = Allocate(numElem8);
  Real_t *fy_elem = Allocate(numElem8);
  Real_t *fz_elem = Allocate(numElem8);
  Real_t gamma[4UL][8UL];
{
    gamma[0][0] = 1.0;
    gamma[0][1] = 1.0;
    gamma[0][2] = (-1.0);
    gamma[0][3] = (-1.0);
    gamma[0][4] = (-1.0);
    gamma[0][5] = (-1.0);
    gamma[0][6] = 1.0;
    gamma[0][7] = 1.0;
    gamma[1][0] = 1.0;
    gamma[1][1] = (-1.0);
    gamma[1][2] = (-1.0);
    gamma[1][3] = 1.0;
    gamma[1][4] = (-1.0);
    gamma[1][5] = 1.0;
    gamma[1][6] = 1.0;
    gamma[1][7] = (-1.0);
    gamma[2][0] = 1.0;
    gamma[2][1] = (-1.0);
    gamma[2][2] = 1.0;
    gamma[2][3] = (-1.0);
    gamma[2][4] = 1.0;
    gamma[2][5] = (-1.0);
    gamma[2][6] = 1.0;
    gamma[2][7] = (-1.0);
    gamma[3][0] = (-1.0);
    gamma[3][1] = 1.0;
    gamma[3][2] = (-1.0);
    gamma[3][3] = 1.0;
    gamma[3][4] = 1.0;
    gamma[3][5] = (-1.0);
    gamma[3][6] = 1.0;
    gamma[3][7] = (-1.0);
  }
/*************************************************/
/*    compute the hourglass modes */
  Index_t i2;
  
#pragma omp parallel for firstprivate ( numElem, hourg )
  for (i2 = 0; i2 < numElem; ++i2) {
    Real_t *fx_local;
    Real_t *fy_local;
    Real_t *fz_local;
    Real_t hgfx[8UL];
    Real_t hgfy[8UL];
    Real_t hgfz[8UL];
    Real_t coefficient;
    Real_t hourgam0[4UL];
    Real_t hourgam1[4UL];
    Real_t hourgam2[4UL];
    Real_t hourgam3[4UL];
    Real_t hourgam4[4UL];
    Real_t hourgam5[4UL];
    Real_t hourgam6[4UL];
    Real_t hourgam7[4UL];
    Real_t xd1[8UL];
    Real_t yd1[8UL];
    Real_t zd1[8UL];
    const Index_t *elemToNode = (m_nodelist + (8 * i2));
    Index_t i3 = (8 * i2);
    Real_t volinv = (1.0 / determ[i2]);
    Real_t ss1;
    Real_t mass1;
    Real_t volume13;
    for (Index_t i1 = 0; i1 < 4; ++i1) {
      Real_t hourmodx = ((((((((x8n[i3] * gamma[i1][0]) + (x8n[i3 + 1] * gamma[i1][1])) + (x8n[i3 + 2] * gamma[i1][2])) + (x8n[i3 + 3] * gamma[i1][3])) + (x8n[i3 + 4] * gamma[i1][4])) + (x8n[i3 + 5] * gamma[i1][5])) + (x8n[i3 + 6] * gamma[i1][6])) + (x8n[i3 + 7] * gamma[i1][7]));
      Real_t hourmody = ((((((((y8n[i3] * gamma[i1][0]) + (y8n[i3 + 1] * gamma[i1][1])) + (y8n[i3 + 2] * gamma[i1][2])) + (y8n[i3 + 3] * gamma[i1][3])) + (y8n[i3 + 4] * gamma[i1][4])) + (y8n[i3 + 5] * gamma[i1][5])) + (y8n[i3 + 6] * gamma[i1][6])) + (y8n[i3 + 7] * gamma[i1][7]));
      Real_t hourmodz = ((((((((z8n[i3] * gamma[i1][0]) + (z8n[i3 + 1] * gamma[i1][1])) + (z8n[i3 + 2] * gamma[i1][2])) + (z8n[i3 + 3] * gamma[i1][3])) + (z8n[i3 + 4] * gamma[i1][4])) + (z8n[i3 + 5] * gamma[i1][5])) + (z8n[i3 + 6] * gamma[i1][6])) + (z8n[i3 + 7] * gamma[i1][7]));
      hourgam0[i1] = (gamma[i1][0] - (volinv * (((dvdx[i3] * hourmodx) + (dvdy[i3] * hourmody)) + (dvdz[i3] * hourmodz))));
      hourgam1[i1] = (gamma[i1][1] - (volinv * (((dvdx[i3 + 1] * hourmodx) + (dvdy[i3 + 1] * hourmody)) + (dvdz[i3 + 1] * hourmodz))));
      hourgam2[i1] = (gamma[i1][2] - (volinv * (((dvdx[i3 + 2] * hourmodx) + (dvdy[i3 + 2] * hourmody)) + (dvdz[i3 + 2] * hourmodz))));
      hourgam3[i1] = (gamma[i1][3] - (volinv * (((dvdx[i3 + 3] * hourmodx) + (dvdy[i3 + 3] * hourmody)) + (dvdz[i3 + 3] * hourmodz))));
      hourgam4[i1] = (gamma[i1][4] - (volinv * (((dvdx[i3 + 4] * hourmodx) + (dvdy[i3 + 4] * hourmody)) + (dvdz[i3 + 4] * hourmodz))));
      hourgam5[i1] = (gamma[i1][5] - (volinv * (((dvdx[i3 + 5] * hourmodx) + (dvdy[i3 + 5] * hourmody)) + (dvdz[i3 + 5] * hourmodz))));
      hourgam6[i1] = (gamma[i1][6] - (volinv * (((dvdx[i3 + 6] * hourmodx) + (dvdy[i3 + 6] * hourmody)) + (dvdz[i3 + 6] * hourmodz))));
      hourgam7[i1] = (gamma[i1][7] - (volinv * (((dvdx[i3 + 7] * hourmodx) + (dvdy[i3 + 7] * hourmody)) + (dvdz[i3 + 7] * hourmodz))));
    }
{
/* compute forces */
/* store forces into h arrays (force arrays) */
      ss1 = m_ss[i2];
      mass1 = m_elemMass[i2];
    }
    volume13 = cbrt(determ[i2]);
    Index_t n0si2 = elemToNode[0];
    Index_t n1si2 = elemToNode[1];
    Index_t n2si2 = elemToNode[2];
    Index_t n3si2 = elemToNode[3];
    Index_t n4si2 = elemToNode[4];
    Index_t n5si2 = elemToNode[5];
    Index_t n6si2 = elemToNode[6];
    Index_t n7si2 = elemToNode[7];
{
      xd1[0] = m_xd[n0si2];
      xd1[1] = m_xd[n1si2];
      xd1[2] = m_xd[n2si2];
      xd1[3] = m_xd[n3si2];
      xd1[4] = m_xd[n4si2];
      xd1[5] = m_xd[n5si2];
      xd1[6] = m_xd[n6si2];
      xd1[7] = m_xd[n7si2];
      yd1[0] = m_yd[n0si2];
      yd1[1] = m_yd[n1si2];
      yd1[2] = m_yd[n2si2];
      yd1[3] = m_yd[n3si2];
      yd1[4] = m_yd[n4si2];
      yd1[5] = m_yd[n5si2];
      yd1[6] = m_yd[n6si2];
      yd1[7] = m_yd[n7si2];
      zd1[0] = m_zd[n0si2];
      zd1[1] = m_zd[n1si2];
      zd1[2] = m_zd[n2si2];
      zd1[3] = m_zd[n3si2];
      zd1[4] = m_zd[n4si2];
      zd1[5] = m_zd[n5si2];
      zd1[6] = m_zd[n6si2];
      zd1[7] = m_zd[n7si2];
      coefficient = ((((-hourg * 0.01) * ss1) * mass1) / volume13);
    }
    CalcElemFBHourglassForce(xd1,yd1,zd1,hourgam0,hourgam1,hourgam2,hourgam3,hourgam4,hourgam5,hourgam6,hourgam7,coefficient,hgfx,hgfy,hgfz);
{
      fx_local = (fx_elem + i3);
      fx_local[0] = hgfx[0];
      fx_local[1] = hgfx[1];
      fx_local[2] = hgfx[2];
      fx_local[3] = hgfx[3];
      fx_local[4] = hgfx[4];
      fx_local[5] = hgfx[5];
      fx_local[6] = hgfx[6];
      fx_local[7] = hgfx[7];
      fy_local = (fy_elem + i3);
      fy_local[0] = hgfy[0];
      fy_local[1] = hgfy[1];
      fy_local[2] = hgfy[2];
      fy_local[3] = hgfy[3];
      fy_local[4] = hgfy[4];
      fy_local[5] = hgfy[5];
      fy_local[6] = hgfy[6];
      fy_local[7] = hgfy[7];
      fz_local = (fz_elem + i3);
      fz_local[0] = hgfz[0];
      fz_local[1] = hgfz[1];
      fz_local[2] = hgfz[2];
      fz_local[3] = hgfz[3];
      fz_local[4] = hgfz[4];
      fz_local[5] = hgfz[5];
      fz_local[6] = hgfz[6];
      fz_local[7] = hgfz[7];
    }
  }
{
    Index_t numNode = m_numNode;
    Index_t gnode;
    
#pragma omp parallel for firstprivate ( numNode )
    for (Index_t gnode = 0; gnode < numNode; ++gnode) {
      Index_t count = m_nodeElemCount[gnode];
      Index_t start = m_nodeElemStart[gnode];
      Real_t fx = 0.0;
      Real_t fy = 0.0;
      Real_t fz = 0.0;
      for (Index_t i = 0; i < count; ++i) {
        Index_t elem = m_nodeElemCornerList[start + i];
        fx += fx_elem[elem];
        fy += fy_elem[elem];
        fz += fz_elem[elem];
      }
{
        m_fx[gnode] += fx;
        m_fy[gnode] += fy;
        m_fz[gnode] += fz;
      }
    }
  }
  Release(&fz_elem);
  Release(&fy_elem);
  Release(&fx_elem);
}

void CalcHourglassControlForElems(Real_t determ[],Real_t hgcoef)
{
  Index_t numElem = m_numElem;
  Index_t numElem8 = (numElem * 8);
  Real_t *dvdx = Allocate(numElem8);
  Real_t *dvdy = Allocate(numElem8);
  Real_t *dvdz = Allocate(numElem8);
  Real_t *x8n = Allocate(numElem8);
  Real_t *y8n = Allocate(numElem8);
  Real_t *z8n = Allocate(numElem8);
  Index_t i;
  Index_t j;
  Index_t k;
// #pragma omp parallel for firstprivate(numElem)
  Real_t pfx[8UL];
  Real_t pfy[8UL];
  Real_t pfz[8UL];
  Index_t ii;
  Index_t jj;
{
    int c2;
    int c4;
    int c5;
    int c1;
    int c0;
#pragma omp parallel for private(c1, c5, c4, c2)
    for (c0 = 0; c0 <= 127; c0++) {
      for (c1 = 0; c1 <= 31; c1++) {
        for (c2 = 0; c2 <= 3; c2++) {
          for (c4 = 4 * c1; c4 <= 4 * c1 + 3; c4++) {
#pragma ivdep
#pragma vector always
#pragma simd
            for (c5 = 32 * c2; c5 <= 32 * c2 + 31; c5++) {
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = (1.0 / 12.0 * (((((((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1])) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = (1.0 / 12.0 * ((((((-(m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1])) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = (1.0 / 12.0 * ((((((-(m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1])) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = (1.0 / 12.0 * (((((((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1])) - ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = (1.0 / 12.0 * ((((((-(m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1])) + ((m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = (1.0 / 12.0 * ((((((-(m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1])) + ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = (1.0 / 12.0 * (((((((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4])) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = (1.0 / 12.0 * ((((((-(m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4])) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = (1.0 / 12.0 * ((((((-(m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4])) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = (1.0 / 12.0 * (((((((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129])) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = (1.0 / 12.0 * ((((((-(m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129])) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = (1.0 / 12.0 * ((((((-(m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129])) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = (1.0 / 12.0 * (((((((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1])) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = (1.0 / 12.0 * ((((((-(m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1])) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = (1.0 / 12.0 * ((((((-(m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1])) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = (1.0 / 12.0 * (((((((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129])) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = (1.0 / 12.0 * ((((((-(m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129])) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = (1.0 / 12.0 * ((((((-(m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129])) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = (1.0 / 12.0 * (((((((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129])) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = (1.0 / 12.0 * ((((((-(m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129])) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = (1.0 / 12.0 * ((((((-(m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129])) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1]))) + ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))) - ((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129]))));
              dvdx[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = (1.0 / 12.0 * (((((((m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1])) - ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) + ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))));
              dvdy[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = (1.0 / 12.0 * ((((((-(m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1])) + ((m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))));
              dvdz[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = (1.0 / 12.0 * ((((((-(m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1])) + ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1]) * (m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1] + m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129]))) - ((m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]) * (m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1] + m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1]))) + ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))) - ((m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129] + m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]) * (m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129] + m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4]))));
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = m_x[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = m_x[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = m_x[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1];
              x8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = m_x[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = m_y[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = m_y[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = m_y[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1];
              y8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = m_y[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 0] = m_z[(((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 1] = m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 1];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 2] = m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129) + 1];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 3] = m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 4] = m_z[((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 5] = m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 1];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 6] = m_z[((((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129) + 1];
              z8n[(8 * ((((c0 * 128) * 128) + (c4 * 128)) + c5)) + 7] = m_z[(((((((c0 * 128) * 128) + (c4 * 128)) + c5) + (c0 * (2 * 128 + 1))) + c4) + 129 * 129) + 129];
              determ[(((c0 * 128) * 128) + (c4 * 128)) + c5] = (m_volo[(((c0 * 128) * 128) + (c4 * 128)) + c5] * m_v[(((c0 * 128) * 128) + (c4 * 128)) + c5]);
            }
          }
        }
      }
    }
  }
  if (hgcoef > 0.0) {
    CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef);
  }
  Release(&z8n);
  Release(&y8n);
  Release(&x8n);
  Release(&dvdz);
  Release(&dvdy);
  Release(&dvdx);
}

int main(int argc,char *argv[])
{
//   Index_t edgeElems = 45 ;
//edgeNodes = edgeElems+1;
  Real_t tx;
  Real_t ty;
  Real_t tz;
  Index_t nidx;
  Index_t zidx;
  Index_t domElems;
{
/* get run options to measure various metrics */
/* ... */
/****************************/
/*   Initialize Sedov Mesh  */
/****************************/
/* construct a uniform box for this processor */
    m_sizeX = 128;
    m_sizeY = 128;
    m_sizeZ = 128;
    m_numElem = (128 * 128 * 128);
    m_numNode = (129 * 129 * 129);
    domElems = m_numElem;
  }
/* allocate field memory */
  AllocateElemPersistent(m_numElem);
  AllocateElemTemporary(m_numElem);
  AllocateNodalPersistent(m_numNode);
  AllocateNodesets((129 * 129));
{
/* initialize nodal coordinates */
    nidx = 0;
    tz = 0.0;
  }
  for (Index_t plane = 0; plane < 129; ++plane) {{
      ty = 0.0;
    }
    for (Index_t row = 0; row < 129; ++row) {{
        tx = 0.0;
      }
      for (Index_t col = 0; col < 129; ++col) {
        m_x[nidx] = tx;
        m_y[nidx] = ty;
        m_z[nidx] = tz;
        ++nidx;
        tx = ((1.125 * (col + 1)) / 128);
      }
{
        ty = ((1.125 * (row + 1)) / 128);
      }
    }
{
      tz = ((1.125 * (plane + 1)) / 128);
    }
  }
{
/* embed hexehedral elements in nodal point lattice */
    nidx = 0;
    zidx = 0;
  }
  for (Index_t plane = 0; plane < 128; ++plane) {
    for (Index_t row = 0; row < 128; ++row) {
      for (Index_t col = 0; col < 128; ++col) {
//        printf("<zidx,nidx>:=<%d,%d>\n", zidx,nidx); 
//        printf("<zidx,nidx>:=<%d,%d>\n", plane*edgeElems*edgeElems+row*edgeElems+col, plane*edgeElems*edgeElems+row*edgeElems+col+plane*(2*edgeElems+1)+row);
        Index_t *localNode = (m_nodelist + (8 * zidx));
        localNode[0] = nidx;
        localNode[1] = (nidx + 1);
        localNode[2] = ((nidx + 129) + 1);
        localNode[3] = (nidx + 129);
        localNode[4] = (nidx + 129 * 129);
        localNode[5] = ((nidx + 129 * 129) + 1);
        localNode[6] = (((nidx + 129 * 129) + 129) + 1);
        localNode[7] = ((nidx + 129 * 129) + 129);
        ++zidx;
        ++nidx;
      }
{
        ++nidx;
      }
    }
{
      nidx += 129;
    }
  }
  AllocateNodeElemIndexes();
/* Create a material IndexSet (entire same material for now) */
  for (Index_t i = 0; i < domElems; ++i) {
    m_matElemlist[i] = i;
  }
{
/* initialize material parameters */
    m_dtfixed = (-1.0e-7);
    m_deltatime = 1.0e-7;
    m_deltatimemultlb = 1.1;
    m_deltatimemultub = 1.2;
    m_stoptime = 1.0e-6;
    m_dtcourant = 1.0e+20;
    m_dthydro = 1.0e+20;
    m_dtmax = 0.01;
    m_time = 0.0;
    m_cycle = 0;
    m_e_cut = 1.0e-7;
    m_p_cut = 1.0e-7;
    m_q_cut = 1.0e-7;
    m_u_cut = 1.0e-7;
    m_v_cut = 1.0e-10;
    m_hgcoef = 3.0;
    m_ss4o3 = (4.0 / 3.0);
    m_qstop = 1.0e+12;
    m_monoq_max_slope = 1.0;
    m_monoq_limiter_mult = 2.0;
    m_qlc_monoq = 0.5;
    m_qqc_monoq = (2.0 / 3.0);
    m_qqc = 2.0;
    m_pmin = 0.0;
    m_emin = (-1.0e+15);
    m_dvovmax = 0.1;
    m_eosvmax = 1.0e+9;
    m_eosvmin = 1.0e-9;
    m_refdens = 1.0;
  }
/* initialize field data */
  for (Index_t i = 0; i < domElems; ++i) {
    Real_t x_local[8UL];
    Real_t y_local[8UL];
    Real_t z_local[8UL];
    Index_t *elemToNode = (m_nodelist + (8 * i));
    for (Index_t lnode = 0; lnode < 8; ++lnode) {
      Index_t gnode = elemToNode[lnode];
      x_local[lnode] = m_x[gnode];
      y_local[lnode] = m_y[gnode];
      z_local[lnode] = m_z[gnode];
    }
// volume calculations
    Real_t volume = CalcElemVolume_3(x_local,y_local,z_local);
{
      m_volo[i] = volume;
      m_elemMass[i] = volume;
    }
    for (Index_t j = 0; j < 8; ++j) {
      Index_t idx = elemToNode[j];
      m_nodalMass[idx] += (volume / 8.);
    }
  }
{
/* deposit energy */
    m_e[0] = 3.948746e+7;
/* set up symmetry nodesets */
    nidx = 0;
  }
  for (Index_t i = 0; i < 129; ++i) {
    Index_t planeInc = ((i * 129) * 129);
    Index_t rowInc = (i * 129);
    for (Index_t j = 0; j < 129; ++j) {
      m_symmX[nidx] = (planeInc + (j * 129));
      m_symmY[nidx] = (planeInc + j);
      m_symmZ[nidx] = (rowInc + j);
      ++nidx;
    }
  }
{
/* set up elemement connectivity information */
    m_lxim[0] = 0;
  }
  for (Index_t i = 1; i < domElems; ++i) {
    m_lxim[i] = (i - 1);
    m_lxip[i - 1] = i;
  }
{
    m_lxip[domElems - 1] = (domElems - 1);
  }
  for (Index_t i = 0; i < 128; ++i) {
    m_letam[i] = i;
    m_letap[(domElems - 128) + i] = ((domElems - 128) + i);
  }
  for (Index_t i = 128; i < domElems; ++i) {
    m_letam[i] = (i - 128);
    m_letap[i - 128] = i;
  }
  for (Index_t i = 0; i < 128 * 128; ++i) {
    m_lzetam[i] = i;
    m_lzetap[(domElems - 128 * 128) + i] = ((domElems - 128 * 128) + i);
  }
  for (Index_t i = (128 * 128); i < domElems; ++i) {
    m_lzetam[i] = (i - 128 * 128);
    m_lzetap[i - 128 * 128] = i;
  }
/* set up boundary condition information */
  for (Index_t i = 0; i < domElems; ++i) {
/* clear BCs by default */
    m_elemBC[i] = 0;
  }
/* faces on "external" boundaries will be */
/* symmetry plane or free surface BCs */
  for (Index_t i = 0; i < 128; ++i) {
    Index_t planeInc = ((i * 128) * 128);
    Index_t rowInc = (i * 128);
    for (Index_t j = 0; j < 128; ++j) {
      m_elemBC[planeInc + (j * 128)] |= 1;
      m_elemBC[((planeInc + (j * 128)) + 128) - 1] |= 8;
      m_elemBC[planeInc + j] |= 16;
      m_elemBC[((planeInc + j) + 128 * 128) - 128] |= 128;
      m_elemBC[rowInc + j] |= 0x100;
      m_elemBC[((rowInc + j) + domElems) - 128 * 128] |= 0x800;
    }
  }
/* timestep to solution */
  while(m_time < m_stoptime){
    TimeIncrement();
    LagrangeLeapFrog();
#if LULESH_SHOW_PROGRESS
#endif
  }
  return 0;
}
