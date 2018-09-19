/*
MIT License

Copyright (c) 2018 Blur Studio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once
#include <vector>
#include <array>
#include <memory>
#include <utility>
#include <cmath>
#include <numeric>
#include "mpkdtree.h"


// Template overrides for possible implementations different from std
template <typename PointArray>
inline void resize(PointArray &a, unsigned size) {
	a.resize(size);
}

template <typename PointArray>
inline unsigned size(const PointArray &a) {
	return a.size();
}

template <typename Vector, typename Float = double>
inline Float dot(const Vector &a, const Vector &b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename Vector, typename Float = double>
inline Float length(const Vector &a) {
	return sqrt(dot<Vector, Float>(a, a));
}

template <typename Vector, typename Float = double>
inline Vector normalized(const Vector &a) {
	return a / length<Vector, Float>(a);
}

template <typename Vector>
inline Vector cross(const Vector &a, const Vector &b) {
	return Vector((a[1] * b[2] - a[2] * b[1]), (a[0] * b[2] - a[2] * b[0]), (a[0] * b[1] - a[1] * b[0]));
}


// Functions for solving automatic tangents
template <typename Vector, typename Float = double>
Vector smoothTangent(const Vector &preTfm, const Vector &curTfm, const Vector &nxtTfm, Float weight){
	// This code calculates only the smooth tangent.
	// It doesn't do anything about the linear, or user tans
	Vector smo;
	Vector preLeg = preTfm - curTfm;
	Vector nextLeg = nxtTfm - curTfm;
	Float preLegLen = length(preLeg);
	Float nextLegLen = length(nextLeg);
	Vector preNorm = preLeg / preLegLen;
	Vector postNorm = nextLeg / nextLegLen;

	Float dprod = dot<Vector, Float>(preNorm, postNorm);
	if (abs(dprod) >= 0.999999999 || preLegLen == 0.0) { // Linear case
		smo = nextLeg / 3.0;
	}
	else { // Nonlinear
		Vector bin = normalized<Vector, Float>(cross(preNorm, postNorm));
		smo = normlaized<Vector, Float>(cross(bin, preNorm) + cross(bin, postNorm));
		smo *= nextLegLen / 3.0;
	}
	smo *= weight;
	return smo;
}

template <typename Vector, typename Float = double>
Vector userTangent(
		const Vector &userTfm, const Vector &curTfm,
		const Vector &smoothTfm, const Vector &nlt,
		Float smooth, Float autoTan
	){
	Vector lin = normalized(nlt - curTfm) * length(smoothTfm);
	Vector result = (1.0 - smooth) * lin + smooth * smoothTfm;
	Vector freeLeg = userTfm - curTfm;
	result = (1.0 - autoTan) * freeLeg + autoTan * result;
	result += curTfm;
	return result;
}

// A helper function for quickly finding an index and segment percentage for an input tValue
template <typename Float=double>
void linearIndex(Float t, const std::vector<Float> &samp, Float &segT, size_t &segIdx) {
	// Remap a length T-Value to a param TValue
	auto ubIt = std::upper_bound(samp.begin(), samp.end(), t);
	// get the new T value and the segment index
	if (ubIt == samp.begin()) {
		segIdx = 0;
		if (samp.size() > 1) {
			segT = t / samp[1];
		}
		else {
			segT = 0.0;
		}
	}
	else if (ubIt == samp.end()) {
		// subtract 1 because (size-1) is the end index
		// subtract another 1 because we're indexing into
		// the segments between the samples
		segIdx = samp.size() - 2;
		if (samp.size() > 1) {
			segT = (t - samp[segIdx]) / (samp[segIdx + 1] - samp[segIdx]);
		}
		else {
			segT = 1.0;
		}
	}
	else {
		auto lbIt = ubIt - 1;
		segT = (t - *lbIt) / (*ubIt - *lbIt);
		segIdx = lbIt - samp.begin(); 
		if (segIdx < samp.size() - 2 && segT > 0.99999999999) {
			// Snap values close to the end to the next segment
			// This removes a lot of jitter
			segT = 0.0;
			++segIdx;
		}
	}
}

template <typename Float = double>
void multiLinearIndexes(const std::vector<Float> &params, const std::vector<Float> &samp, std::vector<Float> &segTs, std::vector<size_t> &segIdxs) {
	// The params aren't necessarily sorted, so to do the "merge", I have to sort them first.
	// but, so I can 'unsort' them later, I'm doing this the annoying way, by index
	std::vector<size_t> pOrder(params.size());
	std::iota(pOrder.begin(), pOrder.end(), 0); // like python's range()
	sort(pOrder.begin(), pOrder.end(), [&params](size_t i1, size_t i2) {return params[i1] < params[i2]; });

	size_t segIdx = 0;
	segTs.resize(params.size());
	segIdxs.resize(params.size());
	bool capped = false; // are we past the end of the remap param values?
	for (size_t i = 0; i < params.size(); ++i) {
		Float t = params[pOrder[i]];
		if (!capped && t >= remap[segIdx + 1]) {
			// If I'm past the end of this segment, go on to the next
			// unless I'm at the very end, in which case, cap me
			if (segIdx + 1 >= params.size()) capped = true;
			++segIdx;
		}
		segIdxs[pOrder[i]] = segIdx;
		segTs[pOrder[i]] = (t - remap[segIdx - 1]) / (remap[segIdx] - remap[segIdx - 1]);
	}
}


// Prototype
template <typename PointArray, typename Point, typename VectorArray, typename Vector, typename QuatArray, typename Quat, typename Float = double>
class TwistSpline;

// This class doesn't own the Point objects, they just look at them
// This class should not do any memory management as it will be
// completely managed by it's parent TwistSpline instance
template <typename PointArray, typename Point, typename VectorArray, typename Vector, typename QuatArray, typename Quat, typename Float = double>
class TwistSplineSegment {
private:
	std::array<Point*, 4> verts; // Convenience pointers to the transforms
	std::array<Quat*, 4> quats; // Pointer to the position/orientation of the verts

	std::array<Vector, 3> d1verts; // first derivative segments
	std::array<Vector, 2> d2verts; // second derivative segments

	PointArray points; // The point lookup table along the curve
	VectorArray tangents;
	VectorArray rnormals; // Raw normals frame-walked starting from the inorm
	VectorArray rbinormals;
	VectorArray tnormals; // Normals twisted around the spline by the user
	VectorArray tbinormals;
	std::vector<Float> twistVals; // the angle, sin, and cos for each twist along the spline

	VectorArray units; // unit vectors along the curve. Makes closest point lookup faster
	Vector iNorm;
	std::vector<Float> sampleLengths; // The sampleLength param for each point in the points
	size_t lutSteps;

public:
	TwistSplineSegment(Point *p0, Point *p1, Point *p2, Point *p3, Quat *q0, Quat *q1, Quat *q2, Quat *q3, const Vector &iNorm, size_t lutSteps) {
		verts[0] = p0; verts[1] = p1; verts[2] = p2; verts[3] = p3;
		quats[0] = q0; quats[1] = q1; quats[2] = q2; quats[3] = q3;
		this->iNorm = iNorm;
		this->lutSteps = lutSteps;

		// resize the point-wise arrays
		size_t pointSteps = lutSteps + 1;
		resize(points, pointSteps);
		resize(tangents, pointSteps);
		resize(rnormals, pointSteps);
		resize(rbinormals, pointSteps);
		resize(sampleLengths, pointSteps);

		// resize the line-wise arrays
		resize(units, lutSteps);

		buildDverts();
		buildLut();
	}

	TwistSplineSegment(TwistSplineSegment const &old, std::array<Point*, 4> &verts, std::array<Quat*, 4> &quats){
		this->d1verts = old.d1verts;
		this->d2verts = old.d2verts;
		this->points = old.points;
		this->tangents = old.tangents;
		this->rnormals = old.rnormals;
		this->rbinormals = old.rbinormals;
		this->tnormals = old.tnormals;
		this->tbinormals = old.tbinormals;
		this->twistVals = old.twistVals;
		this->units = old.units;
		this->iNorm = old.iNorm;
		this->sampleLengths = old.sampleLengths;
		this->lutSteps = old.lutSteps;
		this->verts = verts;
		this->quats = quats;
	}

	~TwistSplineSegment() {}

	size_t getLutSteps() { return lutSteps; }
	PointArray& getPoints() { return points; }
	VectorArray& getTangents() { return tangents; }
	VectorArray& getRawNormals() { return rnormals; }
	VectorArray& getRawBinormals() { return rbinormals; }
	VectorArray& getTwistNormals() { return tnormals; }
	VectorArray& getTwistBinormals() { return tbinormals; }
	std::vector<Float>& getTwistVals() { return twistVals; }

	VectorArray& getUnits() { return units; }
	std::vector<Float>& getSampleLengths() { return sampleLengths; }
	Float getLength() const {
		size_t ss = sampleLengths.size();
		if (ss) {
			return sampleLengths[ss - 1];
		}
		return 0.0;
	}

	Point compute(Float t) const {
		if (t <= 0.0) return *(verts[0]);
		if (t >= 1.0) return *(verts[3]);
		// When calculating a single t-value, just use the bernstein because the pre-computation is heavier than one step
		Float mt = 1.0 - t;
		Float a = mt * mt * mt;
		Float b = mt * mt * t * 3.0;
		Float c = mt * t * t * 3.0;
		Float d = t * t * t;
		return a * (*verts[0]) + b * (*verts[1]) + c * (*verts[2]) + d * (*verts[3]);
	}

	Vector tangent(Float t) const {
		if (t <= 0.0) return d1verts[0];
		if (t >= 1.0) return d1verts[2];
		// When calculating a single t-value, just use the bernstein because the pre-computation is heavier than one step
		Float mt = 1.0 - t;
		Float a = mt * mt;
		Float b = mt * t * 2.0;
		Float c = t * t;

		return a * d1verts[0] + b * d1verts[1] + c * d1verts[2];
	}


	Vector rejectNormal(const Vector &onto, const Vector &n) const {
		// This is the only normal we need to calculate 99% of the time
		// Get the -y axis, then reject it onto the first derivative
		return normalized<Vector, Float>(n - ((dot<Vector, Float>(n, onto) / dot<Vector, Float>(onto, onto)) * onto));
	}

	Vector rejectInitialNormal(const Vector &n) const {
		// This is the only normal we need to calculate 99% of the time
		// Get the -y axis, then reject it onto the first derivative
		return rejectNormal(tangent(0.0), n);
	}

	Vector initialNormal() const {
		// This is the only normal we need to calculate 99% of the time
		// Get the -y axis, then reject it onto the first derivative
		Vector y = Vector(0.0, 1.0, 0.0);
		Vector n = y.rotateBy(*(quats[0]));
		return rejectInitialNormal(n);
	}

	Float postAngle() const {
		// Get the angle between the last normal and the y-axis of the last CV
		Vector y = Vector(0.0, 1.0, 0.0);
		Vector n = y.rotateBy(*(quats[3]));
		const Vector &tLast = tangents[lutSteps];
		const Vector &bLast = rbinormals[lutSteps];
		const Vector &fLast = rnormals[lutSteps];
		Vector cvLast = rejectNormal(tLast, n);

		Float dd = dot<Vector, Float>(cvLast, fLast);
		dd = std::max(std::min(dd, 1.0), -1.0);
		Float angle = acos(dd);
		Float disc = dot<Vector, Float>(bLast, cvLast);
		// One of these may be better? do some testing
		//if (disc >= 0) return -angle;
		if (disc >= 0) return -angle;
		return angle;
	}

	void buildDverts() {
		for (size_t i = 0; i < 3; ++i)
			d1verts[i] = 3 * (*(verts[i + 1]) - (*(verts[i])));

		for (size_t i = 0; i < 2; ++i)
			d2verts[i] = 2 * (d1verts[i + 1] - d1verts[i]);
	}

	void buildLut() {
		/////////////////////////////////////////////////////////////
		// Fancy math transforms the evaluation of a bernstein polynomial into evaluation of a finite Taylor Series
		// Which means (in this case) that we can pre-compute a bunch of stuff and make our loop nothing but some additions
		//
		// The intuition of this algorithm:
		// The der3 value is some constant, and to step between der2 values, you add (der3 * stepLength)
		// So it follows that stepping to der1 values is just some linear function of the der2 values
		// Same with the output. There's complicated proof, but I think of it like forces pulling each other around

		// Do the pre-calculations
		Vector a = (*(verts[3])) - 3 * (*(verts[2])) + 3 * (*(verts[1])) - (*(verts[0]));
		Vector b = 3 * (*(verts[2])) - 6 * (*(verts[1])) + 3 * (*(verts[0]));
		Vector c = 3 * (*(verts[1])) - 3 * (*(verts[0]));
		Vector d = (*(verts[0]));

		points[0] = *verts[0];
		Vector tan = tangent(0.0); // only place I ever need the non-normalized tangent
		tangents[0] = normalized<Vector, Float>(tan);

		Float h = 1.0 / (Float)lutSteps;
		Float h2 = h * h;
		Float h3 = h2 * h;

		Vector fd = a * h3 + b * h2 + c * h;
		Vector fd2 = 6 * a*h3 + 2 * b*h2;
		Vector fd3 = 6 * a*h3;

		Vector td = 3 * a*h2 + 2 * b*h;
		Vector td2 = 6 * a*h2;

		// Wonderfully simple loop
		size_t pointSteps = lutSteps + 1;
		for (size_t i = 1; i < pointSteps; i++) {
			d += fd;
			fd += fd2;
			fd2 += fd3;
			points[i] = d;

			tan += td;
			td += td2;
			tangents[i] = normalized<Vector, Float>(tan);
		}

		/////////////////////////////////////////////////////////////
		// Now compute the rnormals, rbinormals, and arc-length
		// Algorithm taken from "Computation of Rotation Minimizing Frames"
		// by W. Wang et. al.
		// This is called the "Double Refleciton Method". It takes one frame
		// Reflects it to the next sampled position, then reflects it again
		// so that the reflected tangent matches the computed one. The resulting
		// up-vector is a surprisingly good approximation of the true RMF

		rnormals[0] = rejectInitialNormal(iNorm);
		rbinormals[0] = cross(tangents[0], rnormals[0]);
		sampleLengths[0] = 0.0;

		// A chunk of this loop could be parallelized
		// but that might not play well with maya's parallel DAG evaluation
		for (size_t i = 0; i < lutSteps; i++) {
			Vector v1 = points[i + 1] - points[i];
			Float v12 = dot<Vector, Float>(v1, v1);
			Float c1 = 2.0 / v12;
			Vector tLi = tangents[i] - c1 * dot<Vector, Float>(v1, tangents[i]) * v1;
			Vector v2 = tangents[i + 1] - tLi;
			Float dv2 = dot<Vector, Float>(v2, v2);
			if (abs(dv2) < 1.0e-15) {
				rnormals[i + 1] = rnormals[i];
				rbinormals[i + 1] = rbinormals[i];
				continue;
			}
			Float c2 = 2.0 / dv2;

			// While I'm in here, I can calculate the sampleLengths and unit parameters
			Float len = sqrt(v12);
			sampleLengths[i + 1] = sampleLengths[i] + len;
			units[i] = v1 / len;

			// non-parallelizable Section
			Vector rLi = rnormals[i] - c1 * dot<Vector, Float>(v1, rnormals[i]) * v1;
			rnormals[i + 1] = rLi - c2 * dot<Vector, Float>(v2, rLi) * v2;
			rbinormals[i + 1] = cross(tangents[i + 1], rnormals[i + 1]);
		}
	}


	void twistedMatrixAtParam(Float rawT, Vector &tan, Vector &norm, Vector &binorm, Point &tran, Point &scl, Float &twist) const {
		matrixAtParam(rawT, tan, norm, binorm, tran, scl, twist, tnormals, tbinormals);
	}

	void rawMatrixAtParam(Float rawT, Vector &tan, Vector &norm, Vector &binorm, Point &tran, Point &scl, Float &twist) const {
		matrixAtParam(rawT, tan, norm, binorm, tran, scl, twist, rnormals, rbinormals);
	}

	void matrixAtParam(Float rawT, Vector &tan, Vector &norm, Vector &binorm, Point &tran, Point &scl, Float &twist, const VectorArray &normals, const VectorArray &binormals) const {
		// rawT is 0-1 normalized length
		// lenT is unnormalized length
		// segT is 0-1 param value
		// t is 0-1 param of a lut segment
		Float t, segT, lenT = rawT * getLength();
		size_t i;
		linearIndex(lenT, sampleLengths, t, i);
		segT = (t + i) / lutSteps;
		if (abs(t) < 1.0e-11) { // t == 0
			tan = tangents[i];
			norm = normals[i];
			binorm = binormals[i];
			tran = points[i];
			twist = twistVals[i];
			scl = Point(1.0, 1.0, 1.0); // TODO
			return;
		}
		if (t < 0.0) {
			tran = (normalized(tangents[0]) * lenT) + points[0];
			tan = tangents[0];
			norm = normals[0];
			binorm = binormals[0];
			twist = twistVals[0];
			scl = Point(1.0, 1.0, 1.0); // TODO
			return;
		}
		if (t > 0.99999999999) {
			tran = (normalized(tangents[size(tangents) - 1]) * (lenT - getLength())) + points[size(points) - 1];
			tan = tangents[size(tangents) - 1];
			norm = normals[size(normals) - 1];
			binorm = tbinormals[size(tbinormals) - 1];
			twist = twistVals[twistVals.size() - 1];
			scl = Point(1.0, 1.0, 1.0); // TODO
			return;
		}

		tran = compute(segT);
		tan = normalized(tangent(segT));

		// Interpolating the orientation could be done:
		// by lerp+normalization: Bad (can end up over-extrapolating in extreme cases)
		// by the double reflection method: Bad, too complicated/ slow
		// by Simple interpolation and rebuild the basis: Good enough
		auto n1 = normals[i];
		auto n2 = normals[i + 1];
		auto nn = n1 * (1-t) + n2 * t;
		binorm = normalized(cross(tan, nn));
		norm = cross(binorm, tan);

		auto tw1 = twistVals[i];
		auto tw2 = twistVals[i + 1];
		twist = tw1 * (1 - t) + tw2 * t;
		scl = Point(1.0, 1.0, 1.0); // TODO
	}

	// By here I've already figured out if the params are extrapolating off either end
	// AND I've taken care of the floating point weirdness
	//void matricesAtParams(Float rawT, Vector &tan, Vector &norm, Vector &binorm, Point &tran, Point &scl, Float &twist, const VectorArray &normals, const VectorArray &binormals) const {
	void matricesAtParams(
			std::vector<Float> segTs, VectorArray &norms, VectorArray &binorms, PointArray &trans,
			PointArray &scls, std::vector<Float> &twists, const VectorArray &normals, const VectorArray &binormals
		) const {

		// It's a little confusing, but I'm repurposing segTs for
		// all the different param types from the single version
		// so I don't have to reallocate memory
		std::vector<Float> ts;
		std::vector<size_t> segIdxs;
		Float len = getLength();
		for (auto &segT : segTs){
			segT *= len;
		}

		linearIndex(segTs, sampleLengths, ts, segIdxs);

		for (size_t i=0; i<ts.size(); ++i){
			Float segT = (ts[i] + segIdxs[i]) / lutSteps;

			tran = compute(segT);
			tan = normalized(tangent(segT));

			auto n1 = normals[segIdxs[i]];
			auto n2 = normals[segIdxs[i] + 1];
			auto nn = n1 * (1-t) + n2 * t;
			binorm = normalized(cross(tan, nn));
			norm = cross(binorm, tan);

			auto tw1 = twistVals[segIdxs[i]];
			auto tw2 = twistVals[segIdxs[i] + 1];
			twist = tw1 * (1 - t) + tw2 * t;
			scl = Point(1.0, 1.0, 1.0);

			norms.push_back(norm);
			binorms.push_back(binorm);
			trans.push_back(tran);
			scls.push_back(scl);
			twists.push_back(twist);
		}
	}

	void matricesAtPreParams() const {
		/*
		tran = (normalized(tangents[0]) * lenT) + points[0];
		tan = tangents[0];
		norm = normals[0];
		binorm = binormals[0];
		twist = twistVals[0];
		scl = Point(1.0, 1.0, 1.0);
		return;
		*/
	}

	void matricesAtPostParams() const {
		/*
		tran = (normalized(tangents[size(tangents) - 1]) * (lenT - getLength())) + points[size(points) - 1];
		tan = tangents[size(tangents) - 1];
		norm = normals[size(normals) - 1];
		binorm = tbinormals[size(tbinormals) - 1];
		twist = twistVals[twistVals.size() - 1];
		scl = Point(1.0, 1.0, 1.0);
		return;
		*/
	}




	void applyTwist(Float startAngle, Float endAngle){ // inRadians
		resize(tnormals, size(rnormals));
		resize(tbinormals, size(rbinormals));
		resize(twistVals, size(rnormals));
		Float len = getLength();
		
		for (size_t i=0; i <= lutSteps; ++i){
			Float perc = sampleLengths[i] / len; // parameterize by length
			Float angle = ((endAngle - startAngle) * perc) + startAngle;
			twistVals[i] = angle;
			const Vector &x = rnormals[i];
			const Vector &y = rbinormals[i];
			const Vector &n = tangents[i];

			// Everything should already be normalized
			//Float ca = cos(angle);
			//Float sa = sin(angle);
			//tnormals[i] = x * ca + n * dot<Vector, Float>(n, x) * (1 - ca) + cross(x, n) * sa;
			// dot(n, x) is 0 b/c they're perpendicular, and I've already got cross(n, x) as the binormal
			tnormals[i] = x * cos(angle) + y * sin(angle);
			tbinormals[i] = cross(n, tnormals[i]);
		}
	}
};



// For later:
// I'll bet I can check the derivatives of the spline
// to find where it's even possible to be closer, and turn a closest
// point lookup into something closer to a binary search.

template <typename PointArray, typename Point, typename VectorArray, typename Vector, typename QuatArray, typename Quat, typename Float = double>
class TwistSpline {
private:
	std::vector<std::unique_ptr<TwistSplineSegment<PointArray, Point, VectorArray, Vector, QuatArray, Quat, Float>>> segments;
	PointArray verts; // all verts, ordered, including tangents
	QuatArray quats; // The orientations of the verts

	std::vector<Float> lockPositions; // The locking t-values
	std::vector<Float> lockValues; // The [0-1] lock property for t-values
	std::vector<Float> userTwists; // The user defined twist offsets
	std::vector<Float> twistLocks; // The [0-1] lock property for user twist offsets
	std::vector<Float> orientLocks; // The [0-1] lock property for cv quaternion twist
	std::vector<Float> remap; // The remapped segment endpoint u-values

	size_t projSteps; // The number of steps for a closest point lookup
	size_t lutSteps; // The number of sub-segments per spline segment
	Float totalLength;
	KDNode<Point, Vector, Float> *kdTree;

public:
	TwistSpline() { lutSteps = 20; }
	~TwistSpline() {}

	PointArray getVerts() const { return verts; }
	QuatArray getQuats() const { return quats; }

	std::vector<Float> getLockPositions() const { return lockPositions; }
	std::vector<Float> getLockValues() const { return lockValues; }

	std::vector<Float> getTwistLocks() const { return twistLocks; }
	std::vector<Float> getOrientLocks() const { return orientLocks; }
	std::vector<Float> getUserTwists() const { return userTwists; }
	std::vector<Float> getRemap() const { return remap; }
	Float getTotalLength() const { return totalLength; }
	


	// copy-ish constructor
	TwistSpline(TwistSpline const &old){
		this->verts = old.verts;
		this->quats = old.quats;
		this->lockPositions = old.lockPositions;
		this->lockValues = old.lockValues;
		this->userTwists = old.userTwists;
		this->twistLocks = old.twistLocks;
		this->orientLocks = old.orientLocks;
		this->remap = old.remap;
		this->projSteps = old.projSteps;
		this->lutSteps = old.lutSteps;
		this->totalLength = old.totalLength;
		// specifically skipping the kdtree for now. Maybe later

		size_t numSegs = ((size(verts) - 1) / 3);
		segments.resize(numSegs);
		for (size_t i=0; i<numSegs; ++i){
			std::array<Point*, 4> vv;
			std::array<Quat*, 4> qq;
			vv = {&(verts[3*i]), &(verts[3*i + 1]), &(verts[3*i + 2]), &(verts[3*i + 3])};
			qq = {&(quats[3*i]), &(quats[3*i + 1]), &(quats[3*i + 2]), &(quats[3*i + 3])};
			segments[i] = std::unique_ptr<TwistSplineSegment<PointArray, Point, VectorArray, Vector, QuatArray, Quat, Float>>(
				new TwistSplineSegment<PointArray, Point, VectorArray, Vector, QuatArray, Quat, Float>(*(old.segments[i]), vv, qq)
			); 
		}
	}

	PointArray getPoints()const {
		PointArray pLut;
		resize(pLut, segments.size() * (lutSteps + 1));
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getPoints();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}
	VectorArray getTangents()const {
		VectorArray pLut;
		resize(pLut, segments.size() * (lutSteps + 1));
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getTangents();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}
	VectorArray getNormals()const {
		VectorArray pLut;
		resize(pLut, segments.size() * (lutSteps + 1));
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getTwistNormals();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}
	VectorArray getBinormals()const {
		VectorArray pLut;
		resize(pLut, segments.size() * (lutSteps + 1));
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getTwistBinormals();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}
	VectorArray getUnits()const {
		VectorArray pLut;
		resize(pLut, segments.size() * lutSteps);
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getUnits();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}
	std::vector<Float> getSampleLengths()const {
		std::vector<Float> pLut;
		resize(pLut, segments.size() * lutSteps);
		size_t c = 0;
		for (auto &seg : segments) {
			auto &pa = seg->getSampleLengths();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut[c++] = pa[j];
			}
		}
		return pLut;
	}

	void clearKDTree() {
		if (kdTree) delete kdTree;
	}

	void buildSegments() {
		size_t numSegs = ((size(verts) - 1) / 3);
		segments.resize(numSegs);
		for (size_t i = 0; i < numSegs; ++i) {
			Vector iNorm;
			if (i == 0) {
				Vector y = Vector(0.0, 1.0, 0.0);
				iNorm = y.rotateBy(quats[3 * i]);
			}
			else {
				const auto &pre = segments[i-1];
				size_t last = pre->getLutSteps() - 1;
				Vector d = pre->getRawNormals()[last];
				Vector a = pre->getTangents()[last];
				Vector b = normalized(verts[3 * i + 1] - verts[3 * i]);
				Vector ab = a + b;
				Float ll = length(ab);
				Vector n;
				if (ll == 0.0) {
					n = b;
				}
				else {
					n = normalized(a + b);
				}
				iNorm = d - 2 * dot<Vector, Float>(d, n) * n;
			}
			segments[i] = std::unique_ptr<TwistSplineSegment<PointArray, Point, VectorArray, Vector, QuatArray, Quat, Float>>(
				new TwistSplineSegment<PointArray, Point, VectorArray, Vector, QuatArray, Quat, Float>(
					&(verts[3*i]), &(verts[3*i + 1]), &(verts[3*i + 2]), &(verts[3*i + 3]),
					&(quats[3*i]), &(quats[3*i + 1]), &(quats[3*i + 2]), &(quats[3*i + 3]),
					iNorm, lutSteps
				)
			);
		}
	}

	void setVerts(
			const PointArray &verts,
			const QuatArray &quats,
			const std::vector<Float> &lockPositions,
			const std::vector<Float> &lockValues,
			const std::vector<Float> &userTwists,
			const std::vector<Float> &twistLocks,
			const std::vector<Float> &orientLocks) {
		this->verts = verts;
		this->quats = quats;
		this->lockPositions = lockPositions;
		this->lockValues = lockValues;
		this->userTwists = userTwists;
		this->twistLocks = twistLocks;
		this->orientLocks = orientLocks;
		buildSegments();

		std::vector<Float> ulens(1);
		Float runner = 0;
		for (auto &seg : segments) {
			runner += seg->getLength();
			ulens.push_back(runner);
		}
		totalLength = runner;

		solveParamMatrix(lockPositions, ulens, lockValues, remap);
		solveTwist();
	}

	void solveParamMatrix(
			const std::vector<Float> &rv /*restVals*/,
			const std::vector<Float> &cv /*currentVals*/,
			const std::vector<Float> &lv /*lockVals*/,
			std::vector<Float> &res) {
		// first, build the param matrix
		std::vector<std::array<Float, 3>> mat;
		size_t len = rv.size();
		if (len == 1) return;
		mat.resize(len);
		res.resize(len);

		// Build the matrix
		// Start case
		mat[0][0] = 0.0;
		mat[0][1] = -1.0;
		mat[0][2] = 1.0 - lv[0];
		res[0] = -lv[0] * rv[0] + (1.0 - lv[0]) * (cv[1] - cv[0]);

		// Mid Cases
		Float rvn = rv[rv.size() - 1];
		for (size_t i = 1; i < len - 1; ++i) {
			Float A = (cv[i] - cv[i - 1]) / (cv[i + 1] - cv[i - 1]);
			mat[i][0] = (1.0 - lv[i]) * (1.0 - A);
			mat[i][1] = -1.0;
			mat[i][2] = (1.0 - lv[i]) * A;
			res[i] = -lv[i] * rv[i];
		}

		// End case
		size_t e = len - 1;
		mat[e][0] = 1.0 - lv[e];
		mat[e][1] = -1.0;
		mat[e][2] = 0.0;
		res[e] = -lv[e] * rv[e] - (1.0 - lv[e]) * (cv[e] - cv[e - 1]);

		// Then pass it to the solver	
		solveTridiagonalMatrix(mat, res);
	}

	void solveTwistParamMatrix(
			const std::vector<Float> &rv /*restVals*/,
			const std::vector<Float> &cv /*currentVals*/,
			const std::vector<Float> &lv /*lockVals*/,
			std::vector<Float> &res) {
		// first, build the param matrix
		std::vector<std::array<Float, 3>> mat;
		size_t len = rv.size();
		if (len == 1) return;
		mat.resize(len);
		res.resize(len);

		// Build the matrix
		// Start case
		mat[0][0] = 0.0;
		mat[0][1] = -1.0;
		mat[0][2] = 1.0 - lv[0];
		res[0] = -lv[0] * rv[0];

		// Mid Cases
		for (size_t i = 1; i < len - 1; ++i) {
			Float A = (cv[i] - cv[i - 1]) / (cv[i + 1] - cv[i - 1]);
			mat[i][0] = (1.0 - lv[i]) * (1.0 - A);
			mat[i][1] = -1.0;
			mat[i][2] = (1.0 - lv[i]) * A;
			res[i] = -lv[i] * rv[i];
		}

		// End case
		size_t e = len - 1;
		mat[e][0] = 1.0 - lv[e];
		mat[e][1] = -1.0;
		mat[e][2] = 0.0;
		res[e] = -lv[e] * rv[e];

		// Then pass it to the solver	
		solveTridiagonalMatrix(mat, res);
	}

	static void solveTridiagonalMatrix(std::vector<std::array<Float, 3>> &mat, std::vector<Float> &res) {
		// Solve a tridiagonal matrix in linear time.
		// If you set up parameter values in a specific way, they can be thought of as a matrix
		// That only has numbers along the 3 most central diagonals, which can be solved efficiently

		// Run the first special case substitution
		mat[0][2] = mat[0][2] / mat[0][1];
		res[0] = res[0] / mat[0][1];

		// Do the middle forward substitutions
		for (size_t i = 1; i < res.size(); ++i) {
			mat[i][2] = mat[i][2] / (mat[i][1] - mat[i][0] * mat[i - 1][2]);
			res[i] = (res[i] - mat[i][0] * res[i - 1]) / (mat[i][1] - mat[i][0] * mat[i - 1][2]);
		}

		// Do the back substitution
		for (size_t i = res.size() - 1; i-- > 0; ) {
			res[i] = res[i] - mat[i][2] * res[i + 1];
		}
	}

	void buildKDTree() {
		// Only need the units for building the kdtree. No reason to build them
		// outside of this method

		// group the lines together and build the kdtree
		std::vector<Point *> pLut;
		pLut.reserve(segments.size() * (lutSteps + 1));

		std::vector<std::tuple<Point *, Point *, Vector *, size_t, size_t>> lines;
		lines.reserve(segments.size() * lutSteps);

		for (size_t i = 0; i < segments.size(); ++i) {
			auto &seg = segments[i];
			auto &pa = seg->getPoints();
			auto &un = seg->getUnits();
			for (size_t j = 0; j < size(pa); ++j) {
				pLut.push_back(&pa[j]);
				if (j > 0) {
					lines.push_back(std::make_tuple(&pa[j - 1], &pa[j], &un[j - 1], i, j));
				}
			}
		}

		size_t sLut = pLut.size() - 1;
		Point maxBounds(*(pLut[sLut])), minBounds(*(pLut[sLut]));
		for (auto &p : pLut) {
			auto &pp = *p;
			for (size_t d = 0; d < 3; ++d) {
				Float bb = pp[d];
				if (bb < minBounds[d]) minBounds[d] = bb;
				if (bb > maxBounds[d]) maxBounds[d] = bb;
			}
		}

		kdTree = KDNode<Point, Vector, Float>::split(pLut, lines, maxBounds, minBounds);
	}

	Point getClosestPoint(const Point &q) const {
		Point closest;
		Float curD2 = std::numeric_limits<Float>::max();
		std::array<Float, 3> off;
		kdTree->closestPoint(q, closest, curD2, 0.0, off);
		return closest;
	}

	Float getClosestParam(const Point &q) const {
		return 0.0; // TODO
	}

	void matrixAtParam(Float t, Vector &tan, Vector &norm, Vector &binorm, Point &tran, Point &scl, Float &twist, bool twisted) const {
		// Return the matrix for the given t-value
		// map the arc-length t-value through the remap vector
		Float segT;
		size_t segIdx;

		linearIndex(t, remap, segT, segIdx);

		// Then ask the segment about this new t-value
		if (twisted)
			segments[segIdx]->twistedMatrixAtParam(segT, tan, norm, binorm, tran, scl, twist);
		else
			segments[segIdx]->rawMatrixAtParam(segT, tan, norm, binorm, tran, scl, twist);
	}

	void matricesAtParams(const std::vector<Float> &params,
			VectorArray &tans, VectorArray &norms, VectorArray &binorms,
			PointArray &trans, PointArray &scls, std::vector<Float> &twists, bool twisted) {

		std::vector<Float> segTs;
		std::vector<size_t> segIdxs;
		multiLinearIndexes(params, remap, segTs, segIdxs);

		resize(tans, params.size());
		resize(norms, params.size());
		resize(binorms, params.size());
		resize(trans, params.size());
		resize(scls, params.size());
		twists.resize(params.size());

		for (size_t i = 0; i < params.size(); ++i) {
			if (twisted) {
				segments[segIdxs[i]]->twistedMatrixAtParam(segTs[i], tans[i], norms[i], binorms[i], trans[i], scls[i], twists[i]);
			}
			else {
				segments[segIdxs[i]]->rawMatrixAtParam(segTs[i], tans[i], norms[i], binorms[i], trans[i], scls[i], twists[i]);
			}
		}
	}

	void solveTwist() {
		// Get a running value of the lengths
		
		std::vector<Float> segLens(segments.size() + 1);
		std::vector<Float> orientVals(segments.size() + 1);
		for (size_t i=0; i<segments.size(); ++i){
			segLens[i+1] = segments[i]->getLength() + segLens[i];
			//restVals[i+1] = segments[i]->postAngle() + restVals[i];
			orientVals[i + 1] = segments[i]->postAngle();
		}
	
		// Get the angle difference between the last
		// frame of a spline and the final CV of the spline
		// Those are the orientVals
		std::vector<Float> oriMap;
		solveTwistParamMatrix(orientVals, segLens, orientLocks, oriMap);

		std::vector<Float> twistMap;
		solveTwistParamMatrix(userTwists, segLens, twistLocks, twistMap);

		for (size_t i=0; i<segments.size(); ++i){
			segments[i]->applyTwist(oriMap[i]+ twistMap[i], oriMap[i+1]+ twistMap[i+1]);
		}
	}
};
