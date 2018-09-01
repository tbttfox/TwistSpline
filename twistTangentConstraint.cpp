
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

#include <vector>

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MGlobal.h>
#include <maya/MTypes.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MMatrix.h>
#include <maya/MFnMatrixData.h>
#include <maya/MVector.h>

#include "twistTangentConstraint.h"

#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

MTypeId	TwistTangentConstraint::id(0x001226FA);

MObject TwistTangentConstraint::aOutLinearTarget;
MObject TwistTangentConstraint::aOutLinearTargetX;
MObject TwistTangentConstraint::aOutLinearTargetY;
MObject TwistTangentConstraint::aOutLinearTargetZ;

MObject TwistTangentConstraint::aSmoothTan;
MObject TwistTangentConstraint::aSmoothTanX;
MObject TwistTangentConstraint::aSmoothTanY;
MObject TwistTangentConstraint::aSmoothTanZ;

MObject TwistTangentConstraint::aOut;
MObject TwistTangentConstraint::aOutX;
MObject TwistTangentConstraint::aOutY;
MObject TwistTangentConstraint::aOutZ;

MObject TwistTangentConstraint::aInTangent;
MObject TwistTangentConstraint::aPrevVertex;
MObject TwistTangentConstraint::aCurrentVertex;
MObject TwistTangentConstraint::aNextVertex;
MObject TwistTangentConstraint::aNextLinearTarget;
MObject TwistTangentConstraint::aNextLinearTargetX;
MObject TwistTangentConstraint::aNextLinearTargetY;
MObject TwistTangentConstraint::aNextLinearTargetZ;
MObject TwistTangentConstraint::aAuto;
MObject TwistTangentConstraint::aSmooth;
MObject TwistTangentConstraint::aWeight;
MObject TwistTangentConstraint::aEndpoint;

MStatus TwistTangentConstraint::initialize() {
	MFnMatrixAttribute mAttr;
	MFnNumericAttribute nAttr;
	MFnEnumAttribute eAttr;
	MFnUnitAttribute uAttr;
	MStatus status;

	//----------------- Outputs -----------------

	aOutX = nAttr.create("outX", "ox", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aOutX")
	aOutY = nAttr.create("outY", "oy", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aOutY")
	aOutZ = nAttr.create("outZ", "oz", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aOutZ")
	aOut = nAttr.create("out", "v", aOutX, aOutY, aOutZ, &status);
	CHECKSTAT("aOut")
	addAttribute(aOut);

	aSmoothTanX = nAttr.create("smoothTanX", "stx", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aSmoothTanX")
	aSmoothTanY = nAttr.create("smoothTanY", "sty", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aSmoothTanY")
	aSmoothTanZ = nAttr.create("smoothTanZ", "stz", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aSmoothTanZ")
	aSmoothTan = nAttr.create("smoothTan", "st", aSmoothTanX, aSmoothTanY, aSmoothTanZ, &status);
	CHECKSTAT("aSmoothTan")
	addAttribute(aSmoothTan);

	aOutLinearTargetX = nAttr.create("outLinearTargetX", "ltx", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aLinearTargetX")
	aOutLinearTargetY = nAttr.create("outLinearTargetY", "lty", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aLinearTargetY")
	aOutLinearTargetZ = nAttr.create("outLinearTargetZ", "ltz", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aLinearTargetZ")
	aOutLinearTarget = nAttr.create("outLinearTarget", "lt", aOutLinearTargetX, aOutLinearTargetY, aOutLinearTargetZ, &status);
	CHECKSTAT("aLinearTarget")
	addAttribute(aOutLinearTarget);

	//----------------- Matrices -----------------
	aInTangent = mAttr.create("inTangent", "it");
	mAttr.setHidden(true);
	mAttr.setDefault(MMatrix::identity);
	addAttribute(aInTangent);

	aPrevVertex = mAttr.create("previousVertex", "pv");
	mAttr.setHidden(true);
	mAttr.setDefault(MMatrix::identity);
	addAttribute(aPrevVertex);

	aCurrentVertex = mAttr.create("currentVertex", "cv");
	mAttr.setHidden(true);
	mAttr.setDefault(MMatrix::identity);
	addAttribute(aCurrentVertex);

	aNextVertex = mAttr.create("nextVertex", "nv");
	mAttr.setHidden(true);
	mAttr.setDefault(MMatrix::identity);
	addAttribute(aNextVertex);

	aNextLinearTargetX = nAttr.create("inLinearTargetX", "nltx", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aNextLinearTargetX")
	aNextLinearTargetY = nAttr.create("inLinearTargetY", "nlty", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aNextLinearTargetY")
	aNextLinearTargetZ = nAttr.create("inLinearTargetZ", "nltz", MFnNumericData::kDouble, 0.0, &status);
	CHECKSTAT("aNextLinearTargetZ")
	aNextLinearTarget = nAttr.create("inLinearTarget", "nlt", aNextLinearTargetX, aNextLinearTargetY, aNextLinearTargetZ, &status);
	CHECKSTAT("aNextLinearTarget")
	addAttribute(aNextLinearTarget);

	//----------------- Weights -----------------
	aAuto = nAttr.create("auto", "a", MFnNumericData::kDouble, 0.0);
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);
	nAttr.setKeyable(true);
	addAttribute(aAuto);

	aSmooth = nAttr.create("smooth", "s", MFnNumericData::kDouble, 1.0);
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);
	nAttr.setKeyable(true);
	addAttribute(aSmooth);

	aWeight = nAttr.create("weight", "w", MFnNumericData::kDouble, 1.0);
	nAttr.setMin(0.0);
	nAttr.setMax(3.0);
	nAttr.setKeyable(true);
	addAttribute(aWeight);
	
	attributeAffects(aPrevVertex, aSmoothTan);
	attributeAffects(aPrevVertex, aSmoothTanX);
	attributeAffects(aPrevVertex, aSmoothTanY);
	attributeAffects(aPrevVertex, aSmoothTanZ);
	attributeAffects(aCurrentVertex, aSmoothTan);
	attributeAffects(aCurrentVertex, aSmoothTanX);
	attributeAffects(aCurrentVertex, aSmoothTanY);
	attributeAffects(aCurrentVertex, aSmoothTanZ);
	attributeAffects(aNextVertex, aSmoothTan);
	attributeAffects(aNextVertex, aSmoothTanX);
	attributeAffects(aNextVertex, aSmoothTanY);
	attributeAffects(aNextVertex, aSmoothTanZ);
	attributeAffects(aWeight, aSmoothTan);
	attributeAffects(aWeight, aSmoothTanX);
	attributeAffects(aWeight, aSmoothTanY);
	attributeAffects(aWeight, aSmoothTanZ);

	attributeAffects(aSmooth, aOutLinearTarget);
	attributeAffects(aSmooth, aOutLinearTargetX);
	attributeAffects(aSmooth, aOutLinearTargetY);
	attributeAffects(aSmooth, aOutLinearTargetZ);
	attributeAffects(aSmoothTan, aOutLinearTarget);
	attributeAffects(aSmoothTan, aOutLinearTargetX);
	attributeAffects(aSmoothTan, aOutLinearTargetY);
	attributeAffects(aSmoothTan, aOutLinearTargetZ);
	attributeAffects(aCurrentVertex, aOutLinearTarget);
	attributeAffects(aCurrentVertex, aOutLinearTargetX);
	attributeAffects(aCurrentVertex, aOutLinearTargetY);
	attributeAffects(aCurrentVertex, aOutLinearTargetZ);
	attributeAffects(aWeight, aOutLinearTarget);
	attributeAffects(aWeight, aOutLinearTargetX);
	attributeAffects(aWeight, aOutLinearTargetY);
	attributeAffects(aWeight, aOutLinearTargetZ);
	attributeAffects(aPrevVertex, aOutLinearTarget);
	attributeAffects(aPrevVertex, aOutLinearTargetX);
	attributeAffects(aPrevVertex, aOutLinearTargetY);
	attributeAffects(aPrevVertex, aOutLinearTargetZ);
	attributeAffects(aCurrentVertex, aOutLinearTarget);
	attributeAffects(aCurrentVertex, aOutLinearTargetX);
	attributeAffects(aCurrentVertex, aOutLinearTargetY);
	attributeAffects(aCurrentVertex, aOutLinearTargetZ);
	attributeAffects(aNextVertex, aOutLinearTarget);
	attributeAffects(aNextVertex, aOutLinearTargetX);
	attributeAffects(aNextVertex, aOutLinearTargetY);
	attributeAffects(aNextVertex, aOutLinearTargetZ);

	attributeAffects(aInTangent, aOut);
	attributeAffects(aInTangent, aOutX);
	attributeAffects(aInTangent, aOutY);
	attributeAffects(aInTangent, aOutZ);
	attributeAffects(aCurrentVertex, aOut);
	attributeAffects(aCurrentVertex, aOutX);
	attributeAffects(aCurrentVertex, aOutY);
	attributeAffects(aCurrentVertex, aOutZ);
	attributeAffects(aSmoothTan, aOut);
	attributeAffects(aSmoothTan, aOutX);
	attributeAffects(aSmoothTan, aOutY);
	attributeAffects(aSmoothTan, aOutZ);
	attributeAffects(aNextLinearTarget, aOut);
	attributeAffects(aNextLinearTarget, aOutX);
	attributeAffects(aNextLinearTarget, aOutY);
	attributeAffects(aNextLinearTarget, aOutZ);
	attributeAffects(aSmooth, aOut);
	attributeAffects(aSmooth, aOutX);
	attributeAffects(aSmooth, aOutY);
	attributeAffects(aSmooth, aOutZ);
	attributeAffects(aAuto, aOut);
	attributeAffects(aAuto, aOutX);
	attributeAffects(aAuto, aOutY);
	attributeAffects(aAuto, aOutZ);
	attributeAffects(aWeight, aOut);
	attributeAffects(aWeight, aOutX);
	attributeAffects(aWeight, aOutY);
	attributeAffects(aWeight, aOutZ);
	attributeAffects(aPrevVertex, aOut);
	attributeAffects(aPrevVertex, aOutX);
	attributeAffects(aPrevVertex, aOutY);
	attributeAffects(aPrevVertex, aOutZ);
	attributeAffects(aCurrentVertex, aOut);
	attributeAffects(aCurrentVertex, aOutX);
	attributeAffects(aCurrentVertex, aOutY);
	attributeAffects(aCurrentVertex, aOutZ);
	attributeAffects(aNextVertex, aOut);
	attributeAffects(aNextVertex, aOutX);
	attributeAffects(aNextVertex, aOutY);
	attributeAffects(aNextVertex, aOutZ);
	return MS::kSuccess;
}

MStatus	TwistTangentConstraint::compute(const MPlug& plug, MDataBlock& data) {
	// Don't care what plug it is, just compute everything and set the outputs clean
	MStatus status;
	// TODO: do this stuff with mVectors or mPoints instead of mMatrixes
	if (plug == aSmoothTan || plug == aSmoothTanX || plug == aSmoothTanY || plug == aSmoothTanZ) {
		// Calculate the weighted smooth tangents explicitly
		// Get the first matrix
		MDataHandle preH = data.inputValue(aPrevVertex);
		MDataHandle curH = data.inputValue(aCurrentVertex);
		MDataHandle nextH = data.inputValue(aNextVertex);
		MDataHandle weightH = data.inputValue(aWeight);

		MTransformationMatrix preTMat(preH.asMatrix());
		MTransformationMatrix curTMat(curH.asMatrix());
		MTransformationMatrix nextTMat(nextH.asMatrix());

		MVector preTfm = preTMat.getTranslation(MSpace::kWorld);
		MVector curTfm = curTMat.getTranslation(MSpace::kWorld);
		MVector nextTfm = nextTMat.getTranslation(MSpace::kWorld);
		double weight = weightH.asDouble();

		MVector preLeg = preTfm - curTfm;
		MVector nextLeg = nextTfm - curTfm;
		double preLegLen = preLeg.length();
		double nextLegLen = nextLeg.length();

		// We're pointing our tangent from pre->post 
		// This is an auto-tan point, so get the half-angle between
		// the perpendiculars of the legs
		MVector smo;
		MVector preNorm = preLeg / preLegLen;
		MVector postNorm = nextLeg / nextLegLen;
		double dot = preNorm * postNorm;
		if (abs(dot) == 1.0 || preLegLen == 0.0) { // Linear case
			smo = nextLeg / 3.0;
		}
		else { // Nonlinear
			MVector bin = preNorm ^ postNorm;
			bin.normalize();
			smo = ((bin ^ preNorm) + (bin ^ postNorm)).normal();
			smo *= nextLegLen / 3.0;
		}
		smo *= weight;
		MDataHandle matH = data.outputValue(aSmoothTan);
		matH.setMVector(smo);
		matH.setClean();
	}
	else if (plug == aOutLinearTarget || plug == aOutLinearTargetX || plug == aOutLinearTargetY || plug == aOutLinearTargetZ) {
		// Calculate 
		MDataHandle inH = data.inputValue(aSmoothTan);
		MVector smo = inH.asVector();
		MDataHandle smoothH = data.inputValue(aSmooth);
		double smooth = smoothH.asDouble();
		smo *= smooth;

		MDataHandle curH = data.inputValue(aCurrentVertex);
		MTransformationMatrix curTMat(curH.asMatrix());
		MVector curTfm = curTMat.getTranslation(MSpace::kWorld);

		MDataHandle matH = data.outputValue(aOutLinearTarget);
		matH.setMVector(smo + curTfm);
		matH.setClean();
	}
	else if (plug == aOut || plug == aOutX || plug == aOutY || plug == aOutZ) {
		// Get the first matrix
		MDataHandle inH = data.inputValue(aInTangent);
		MTransformationMatrix inTMat(inH.asMatrix());
		MVector inTfm = inTMat.getTranslation(MSpace::kWorld);

		MDataHandle curH = data.inputValue(aCurrentVertex);
		MTransformationMatrix curTMat(curH.asMatrix());
		MVector curTfm = curTMat.getTranslation(MSpace::kWorld);

		MDataHandle smoH = data.inputValue(aSmoothTan);
		MVector smo = smoH.asVector();

		MDataHandle nltH = data.inputValue(aNextLinearTarget);
		MVector nlt = nltH.asVector();

		MDataHandle smoothH = data.inputValue(aSmooth);
		double smooth = smoothH.asDouble();

		MDataHandle weightH = data.inputValue(aWeight);
		double weight = weightH.asDouble();

		MDataHandle autoH = data.inputValue(aAuto);
		double autoTan = autoH.asDouble();

		MVector result;
		MVector lin = (nlt - curTfm).normal() * smo.length();

		if (smooth == 0.0){
			result = lin;
		}
		else if (smooth == 1.0){
			result = smo;
		}
		else {
			result = lin + smooth*(smo - lin);
		}

		MVector freeLeg = inTfm - curTfm;
		result = freeLeg + autoTan*(result - freeLeg); // LERP with the free tangent
		result += curTfm;

		MDataHandle matH = data.outputValue(aOut);
		matH.setMVector(result);
		matH.setClean();
	}
	else {
		return MS::kUnknownParameter;
	}
	return MS::kSuccess;
}




/*
class TwistTangentConstraintCommand : public MPxConstraintCommand {

	virtual const MObject& constraintInstancedAttribute() const;
	virtual const MObject& constraintOutputAttribute() const;
	virtual const MObject& constraintTargetInstancedAttribute() const;
	virtual const MObject& constraintTargetAttribute() const;
	virtual const MObject& constraintTargetWeightAttribute() const;
	virtual const MObject& objectAttribute() const;

	virtual MStatus connectObjectAndConstraint( MDGModifier& modifier );

	virtual MStatus parseArgs(const MArgList &argList);
};
*/


MTypeId TwistTangentConstraintCommand::constraintTypeId() const {
	return TwistTangentConstraint::id;
}

void TwistTangentConstraintCommand::createdConstraint(MPxConstraint *constraint) {
	if (constraint) {
		TwistTangentConstraint *c = (TwistTangentConstraint*) constraint;
	}
	else {
		MGlobal::displayError("Failed to get created constraint.");
	}
}

MStatus TwistTangentConstraintCommand::doIt(const MArgList &argList) {
	MStatus stat = parseArgs(argList);
	if (stat = MS::kFailure) return stat; // Fail on failure
	return MS::kUnknownParameter; // Unknown on success so Maya handles the rest
}

// TODO: Should this use the return of some method? Like constraintTargetInstancedAttribute() ?
MStatus TwistTangentConstraintCommand::connectTarget(void *opaqueTarget, int index) {
	//MStatus status = connectTargetAttribute(opaqueTarget, index, TwistTangentConstraint::targetGeometry);
	if (!status) { status.perror("connectTarget"); return status;}
	return MS::kSuccess;
}

MStatus geometrySurfaceConstraintCommand::appendSyntax() {
	MStatus stat;
	MSyntax theSyntax = syntax(&stat);

	if (MS::kSuccess != stat) {
		MGlobal::displayError("Could not get the parent's syntax");
		return stat;
	}

	// Add our command flags
	theSyntax.addFlag(kPrevVertexFlag, kPrevVertexFlagLong);
	theSyntax.addFlag(kCurVertexFlag, kCurVertexFlagLong);
	theSyntax.addFlag(kNextVertexFlag, kNextVertexFlagLong);
	theSyntax.addFlag(kAutoFlag, kAutoFlagLong);
	theSyntax.addFlag(kSmoothFlag, kSmoothFlagLong);
	return stat;
}


MStatus geometrySurfaceConstraintCommand::parseArgs(const MArgList &argList) {
	MStatus stat;
	MArgDatabase argData(syntax(), argList, &stat);
	if (stat.error()) return MS::kFailure;

	MString userTangentName, preVertexName, curVertexName, nextVertexName; 
	double autoVal, smoothVal;

	if (argData.isFlagSet(kPrevVertexFlag)){
		stat = MString argData.getFlagArgument(kPrevVertexFlag, 0, preVertexName);
	}
	if (argData.isFlagSet(kCurVertexFlag)){
		stat = MString argData.getFlagArgument(kCurVertexFlag, 0, curVertexName);
	}
	if (argData.isFlagSet(kNextVertexFlag)){
		stat = MString argData.getFlagArgument(kNextVertexFlag, 0, nextVertexName);
	}

	if (argData.isFlagSet(kAutoFlag)){
		stat = MString argData.getFlagArgument(kAutoFlag, 0, autoVal);
	}
	if (argData.isFlagSet(kSmoothFlag)){
		stat = MString argData.getFlagArgument(kSmoothFlag, 0, smoothVal);
	}

	stat = getCommandArgument(0, userTangentName);

	// Settings only work at creation time. Would need an
	// attribute on the node in order to push this state
	// into the node at any time.
	ConstraintType typ;
	if (argData.isFlagSet(kConstrainToLargestWeightFlag))
		typ = geometrySurfaceConstraintCommand::kLargestWeight;
	else if (argData.isFlagSet(kConstrainToSmallestWeightFlag))
		typ = geometrySurfaceConstraintCommand::kSmallestWeight;
	else
		typ = geometrySurfaceConstraintCommand::kLargestWeight;
	weightType = typ;

	//Need maya to process
	return MS::kUnknownParameter;
}


