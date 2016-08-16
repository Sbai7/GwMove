/*****************************************************************************
*
* This file is part of GwMove hydrogeological software developed by
* Dr. M. A. Sbai
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*
*****************************************************************************/

#define _USE_MATH_DEFINES

// Project headers
#include <ptracker.h>
#include <myMath.h>

// C++ headers
#include <algorithm>
#include <cmath>



ParticleTracker::ParticleTracker(StrucMesh *nw_mesh,   	// Mesh pointer
	std::vector<Soil> soils,        // soil types
	std::vector<Tpoint3D> coords,		// coordinates of initial points
	double *Head,					// nodal groundwater heads          
	double *Flow)					// nodal flow rates
{
	mesh = nw_mesh; 
	head = Head; 
	Flow = flow;
	initialCoords = coords; 
	Soils = soils; 

	time = 0.0;
	unsigned int ne = mesh->GetNelements();
	dt = new double [ne];

	// calculate coordinates bounds of the computational domain and locally 
	ComputeElemCoordBounds();

	// compute element-wise local time steps 
	ComputeTimeSteps(); 

}

ParticleTracker::~ParticleTracker()
{
	// free memory
	if (dt)  delete[] dt;
	initialCoords.clear();
	flowline.clear(); 
	xemin.clear(); yemin.clear(); zemin.clear();
	xemax.clear(); yemax.clear(); zemax.clear();
	Soils.clear();
}

void ParticleTracker::ComputeElemCoordBounds()
{
	unsigned int i, ie;
	for (i = 0; i < 3; i++) {
		xmin[i] = 0.;
		xmax[i] = 0.; 
	}

	unsigned int ne = mesh->GetNelements();

	// initialise vectors 
	xemin.reserve(ne); yemin.reserve(ne); zemin.reserve(ne);
	xemax.reserve(ne); yemax.reserve(ne); zemax.reserve(ne);

	for (ie = 0; ie < ne; ie++) {
		double vmin, vmax;
		vmin = 1.E+99;
		vmax = 1.E-99;
		xemin.push_back(vmin); yemin.push_back(vmin); zemin.push_back(vmin);
		xemax.push_back(vmax); yemax.push_back(vmin); zemax.push_back(vmin);
	}

	// calculate min & max values
	for (ie = 0; ie < ne; ie++) {
		double coord[3][8];
		mesh->GetNeighboorNodesCoord(ie, coord);
		for (i = 1; i < 8; i++) {
			// search for minimal values
			xemin[ie] = std::min(xemin[ie], coord[0][i]);
			yemin[ie] = std::min(yemin[ie], coord[1][i]);
			zemin[ie] = std::min(zemin[ie], coord[2][i]);

			// search for maximal values
			xemax[ie] = std::max(xemax[ie], coord[0][i]);
			yemax[ie] = std::max(yemax[ie], coord[1][i]);
			zemax[ie] = std::max(zemax[ie], coord[2][i]);

			// update xmin & xmax
			xmin[0] = std::min(xmin[0], xemin[ie]);
			xmin[1] = std::min(xmin[1], xemin[ie]);
			xmin[2] = std::min(xmin[2], xemin[ie]);
			xmax[0] = std::max(xmax[0], xemax[ie]);
			xmax[1] = std::max(xmax[1], xemax[ie]);
			xmax[2] = std::max(xmax[2], xemax[ie]);
		}
		
	}
}

void ParticleTracker::ComputeTimeSteps()
{
	unsigned int i,ie,j;

	unsigned int ne = mesh->GetNelements();

	for (ie = 0; ie < ne; ie++) {
		// construct hexahedral cell i 

		// calculate local conductance matrix
		unsigned int connec[8];
		double coord[3][8];
		double K[3]; // hydraulic conductivities of element ie
		mesh->GetConnectivity(ie, connec);
		mesh->GetNeighboorNodesCoord(ie, coord);
		Soil soil = Soils[mesh->soils[ie] - 1];
		soil.GetK(K);
		Hex8 hex(K, coord, connec);
		hex.dimension = 3; // to setup space dimension	

		// get calculated groundwater heads at nodes of that element
		double h[8]; 
		for (i = 0; i < 8; i++) {
			h[i] = head[connec[i]];
		}

		// calculate velocity at cell center of cell ie 
		Tpoint3D v = hex.CellCentroidVelocity(h);

		double v_total = 0.;
		double step    = 0.;
		double time_min = 1.E+99; 
		for (j = 0; j < 3; j++) {
			if (v[j] == 0.) continue;
			if (j == 0)			step = abs((xemax[ie] - xemin[ie]) / v[j]);
			else if (j == 1)	step = abs((yemax[ie] - yemin[ie]) / v[j]);
			else if (j == 2)	step = abs((zemax[ie] - zemin[ie]) / v[j]);
			time_min = std::min(step, time_min);
		}
		dt[ie] = 0.2 * time_min; // 0.2 multiplier to customize
	}
}

bool ParticleTracker::NaiveFindElement(const Tpoint3D pt, unsigned int &elem, Tpoint3D & pt_local)
{
	// search sequentially among all hex elements in the mesh 
	unsigned int ie;
	unsigned int ne = mesh->GetNelements();
	for (ie = 0; ie < ne; ie++) {
		if (IsInside(pt, ie, pt_local)) {
			elem = ie;
			return true;
		}
	}
	return false;
}

bool ParticleTracker::FindElement(const Tpoint3D pt, unsigned int &elem, Tpoint3D & pt_local)
{
	unsigned int quadElem = 0;

	/*
	 * We decompose the costly full 3D search into 2 components: 
	 * 1) A 2D randomized search of the quadrilateral element (in a 2D slice of the mesh) 
	 *    containing the 2D point having the same X,Y coordinates of the 3D point pt. 
	 * 2) Once the quad element is found we loop over the restricted set of hex elements 
	 *    (of all laters) whose top/bottom quad face contain the projected 3D point.
	 * The algorithm search complexity is redced from n^3 to n^2 or even n^a (1<a<2). 
	 */

	// Do a 2D search ... 
	if (FindQuadElement(pt, quadElem)) {
		// now do a 1D search in the vertical Z direction 
		if (FindLayer(pt, quadElem, elem, pt_local)) return true; 
	}

	return false;
}

bool ParticleTracker::FindQuadElement(const Tpoint3D pt, unsigned int &quadElem)
{
	// ... 
	Tpoint3D quadPolygon[5]; 
	unsigned int ie, iquad;
	unsigned int ne = (mesh->GetNi() - 1) * (mesh->GetNj() - 1); // # of quads in one slice

	for (ie = 0; ie < ne; ie++) {
		unsigned int connec[8];
		double coord[3][8]; 
		mesh->GetConnectivity(ie, connec);
		mesh->GetNeighboorNodesCoord(ie, coord);
		// extract top face quadrilateral
		for (iquad = 0; iquad < 4; iquad++) {
			quadPolygon[iquad].x = coord[0][iquad+4];
			quadPolygon[iquad].y = coord[1][iquad+4];
		}
		quadPolygon[4] = quadPolygon[0]; 
		if (IsInsidePolygon(pt, quadPolygon, 5)) {
			quadElem = ie; 
			return true; 
		}
	}
	return false;
}

bool ParticleTracker::FindLayer(const Tpoint3D pt, const unsigned int quadElem, unsigned int &elem, Tpoint3D &pt_local)
{
	unsigned int l, ie = quadElem;
	unsigned int ni = mesh->GetNi();
	unsigned int nj = mesh->GetNj();
	unsigned int nl = mesh->GetNslices() - 1;

	for (l = 0; l < nl; l++) {
		ie += (ni - 1) * (nj - 1) * l;
		if (IsInside(pt, ie, pt_local)) {
			elem = ie;
			return true;
		}
	}

	return false;
}

bool ParticleTracker::IsInside(const Tpoint3D pt, const unsigned int ie, Tpoint3D & pt_local)
{
	unsigned int i, j, k, l;
	double xr[3], x[3];
	double mat[8][8];


	unsigned int connec[8];
	mesh->GetConnectivity(ie, connec);
	double coord[3][8];
	mesh->GetNeighboorNodesCoord(ie, coord);

	// a 'fake' hex element
	double K[3];
	for (i = 0; i < 3; i++) K[i] = 0.;
	Hex8 hex(K, coord, connec);
	hex.dimension = 3;	
	int lc[3][8];
	hex.GetLocalCoord(lc); 

	// calculate local interpolation coefficients 
	double denom = xemax[ie] - xemin[ie];
	for (i = 0; i < 8; i++) {
		xr[0] = (coord[0][i] - xemin[ie]) / denom;
		xr[1] = (coord[1][i] - xemin[ie]) / denom;
		xr[2] = (coord[2][i] - xemin[ie]) / denom;

		mat[0][i] = 1.;
		mat[1][i] = xr[0];
		mat[2][i] = xr[1];
		mat[3][i] = xr[2];
		mat[4][i] = xr[0] * xr[1];
		mat[5][i] = xr[1] * xr[2];
		mat[6][i] = xr[2] * xr[0];
		mat[7][i] = xr[0] * xr[1] * xr[2];

		// invert the local interpolation matrix 
		if (!InvertMatrix(mat)) return false;

		// now, we're ready to calculate all interpolation coefficients 
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 8; k++) {
				weight[k][j] = 0.;
				for (l = 0; l < 8; l++) {
					weight[k][j] += mat[l][k] * lc[j][l];
				}
			}
		}
		
	}

	/* Section 2 - Now check if point 'pt' is inside element ie */

	// initialize local coordinates point : 'pt_local' 
	pt_local = Tpoint3D(0., 0., 0.); 

	// do quick check 
	if (pt.x < xemin[ie]) return false;
	if (pt.y < yemin[ie]) return false;
	if (pt.z < zemin[ie]) return false;
	if (pt.x > xemax[ie]) return false;
	if (pt.y > yemax[ie]) return false;
	if (pt.z > zemax[ie]) return false;

	// do detailed check 
	xr[0] = (pt.x - xemin[ie]) / denom;
	xr[1] = (pt.y - xemin[ie]) / denom;
	xr[2] = (pt.z - xemin[ie]) / denom;

	x[0] = 1.;
	x[1] = xr[0];
	x[2] = xr[1];
	x[3] = xr[2];
	x[4] = xr[0] * xr[1];
	x[5] = xr[1] * xr[2];
	x[6] = xr[2] * xr[0];
	x[7] = xr[0] * xr[1] * xr[2];
	for (j = 0; j < 8; j++) {
		pt_local.x += x[j] * weight[j][0];
		pt_local.y += x[j] * weight[j][1];
		pt_local.z += x[j] * weight[j][2];
	}

	if (std::abs(pt_local.x) > (1. + pointPrecision)) return false;
	if (std::abs(pt_local.y) > (1. + pointPrecision)) return false;
	if (std::abs(pt_local.z) > (1. + pointPrecision)) return false;

	// we arrive here => the point is inside ! 
	return true;
}

bool ParticleTracker::IsSourceOrSink(const unsigned int elem, Tpoint3D &pt)
{
	unsigned int k;
	unsigned int nr;
	double Qw;

	for (k = 0; k < 8; k++) {
		nr = mesh->GetNode(elem, k);
		if (isForward)	Qw =  flow[nr];
		else			Qw = -flow[nr];

		/*** source or sink detection algorithm ***/
		if (Qw < 0.0) {
			/* In case of unsaturated flow we should calculate the
			 * water content (wc) at this element. Otherwise, in fully
			 * saturated conditions it equals one. 
			 *
			 * pressure = head[nr]-pt_nr.z; 
			 * unsat(ie,pressure,wc,uc);
			 */

			// calculate radius of influence 
			double base   = 3 * (-Qw) * dt[elem] / (4 * M_PI);  // 3 * (-Qw) * dt[nr] / (4 * M_PI * wc )
			double radius = pow(base, 1. / 3.);

			// calculate distance between current particle position
			// and the kth node of this element
			Tpoint3D pt_nr; 
			mesh->GetCoord(nr, pt_nr);
			double distance = Distance(pt, pt_nr);

			// Attraction ?
			if (distance < radius) {
				// calculate final travel time
				time += (4 * M_PI * pow(distance,3.)) / (3 * (-Qw)); // (4 * M_PI * wc * pow(distance,3.)) / (3 * (-Qw))

				// update particle position
				pt = pt_nr;

				return true;
			}
		}
	}
	return false;
}

bool ParticleTracker::Solve()
{
	// empty flowline vector from previous calculations
	flowline.clear(); 

	// reserve memory 
	unsigned int nlines = initialCoords.size(); 
	flowline.reserve(nlines);

	unsigned int fl, ie, ie_new, i;
	Tpoint3D pt_local; 

	// loop over all flowlines ... 
	for (fl = 0; fl < nlines; fl++) {
		Tpoint3D pt = initialCoords[fl]; 
		FlowLine Line(fl + 1, isForward, pt);

		// initialize travel time
		time = 0.0; 

		// add initial position/time to pathline 
		Line.point.push_back(pt); Line.time.push_back(time);

		do {
			// trace the next flowline if false returned 
			if (!FindElement(pt, ie, pt_local)) continue;

			// test wether the particle has arrived to a source/sink position 
			if (IsSourceOrSink(ie, pt)) continue;

			/* Predictor step of the algorithm */

			// calculate local conductance matrix
			unsigned int connec[8];
			double coord[3][8];
			double K[3]; // hydraulic conductivities of element ie
			mesh->GetConnectivity(ie, connec);
			mesh->GetNeighboorNodesCoord(ie, coord);
			Soil soil = Soils[mesh->soils[ie] - 1];
			soil.GetK(K);
			Hex8 hex(K, coord, connec);
			hex.dimension = 3; // to setup space dimension	

			// get calculated groundwater heads at nodes of that element
			double h[8];
			for (i = 0; i < 8; i++) h[i] = head[connec[i]];

			// calculate velocity at cell center of cell ie 
			Tpoint3D v = hex.CellCentroidVelocity(h);
			Tpoint3D pt_new = v*dt[ie] / 2.0;
			if (!isForward) pt_new = pt_new*(-1.0);
			pt_new = pt + pt_new;

			/* Correction step */
			if (!FindElement(pt_new, ie_new, pt_local)) {
				pt_new = v*dt[ie];
				if (!isForward) pt_new = pt_new*(-1.0);
				pt = pt + pt_new;
				time += dt[ie];

				// add position/time to pathline object
				Line.point.push_back(pt); Line.time.push_back(time);
				// add id of current soil type id
				Line.soilId.push_back(soil.GetId()); 
			}
			else {
				mesh->GetConnectivity(ie_new, connec);
				mesh->GetNeighboorNodesCoord(ie_new, coord);
				Soil soil = Soils[mesh->soils[ie_new] - 1]; soil.GetK(K);
				Hex8 hex(K, coord, connec); hex.dimension = 3;
				for (i = 0; i < 8; i++) h[i] = head[connec[i]];
				v = hex.CellCentroidVelocity(h);
				pt_new = v*dt[ie_new];
				if (!isForward) pt_new = pt_new*(-1.0);
				pt = pt + pt_new;
				time += dt[ie_new];

				// add position/time to pathline object
				Line.point.push_back(pt); Line.time.push_back(time);
				// add id of current soil type 
				Line.soilId.push_back(soil.GetId());
			}

		} while (time < maxTime);

		// Store current flow path
		flowline.push_back(Line);
	}
	 
	return true; 
}

