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

#include <maya/MStatus.h>
#include <maya/MObject.h>
#include <maya/MTypeId.h> 
 
#include <maya/MPxConstraint.h>
#include <maya/MPxConstraintCommand.h>

#define kPrevVertexFlagLong "-prevVertex"
#define kPrevVertexFlag "-pv"
#define kCurVertexFlagLong "-curVertex"
#define kCurVertexFlag "-cv"
#define kNextVertexFlagLong "-nextVertex"
#define kNextVertexFlag "-nv"
#define kAutoFlagLong "-auto"
#define kAutoFlag "-a"
#define kSmoothFlagLong "-smooth"
#define kSmoothFlag "-s"

class TwistTangentConstraintCommand : public MPxConstraintCommand {
public:
	TwistTangentConstraintCommand() {}
	~TwistTangentConstraintCommand() {}
	static void* creator() { return new TwistTangentConstraintCommand(); }
	
	virtual MStatus doIt(const MArgList &argList);
	virtual MStatus appendSyntax();

	virtual MTypeId constraintTypeId() const;
	virtual MPxConstraintCommand::TargetType targetType() const;

	virtual const MObject& constraintInstancedAttribute() const;
	virtual const MObject& constraintOutputAttribute() const;
	virtual const MObject& constraintTargetInstancedAttribute() const;
	virtual const MObject& constraintTargetAttribute() const;
	virtual const MObject& constraintTargetWeightAttribute() const;
	virtual const MObject& objectAttribute() const;

	virtual MStatus connectTarget(void *opaqueTarget, int index);
	virtual MStatus connectObjectAndConstraint( MDGModifier& modifier );

	virtual void createdConstraint(MPxConstraint *constraint);

protected:
	virtual MStatus parseArgs(const MArgList &argList);
};






class TwistTangentConstraint : public MPxConstraint {
public:
	TwistTangentConstraint() {}
	virtual	~TwistTangentConstraint() {}
	virtual void postConstructor() {}
	static	void*	creator() { return new TwistTangentConstraint(); }

	virtual	MStatus	compute(const MPlug& plug, MDataBlock& data);
	static	MStatus	initialize();

	virtual const MObject targetAttribute() const { return aCurrentVertex; }
	virtual const MObject constraintRotateOrderAttribute() const { return aRotateOrder; }
	virtual void getOutputAttributes(MObjectArray& attributeArray) {
		attibuteArray.clear();
		attributeArray.append(aOut);
	};

public:
	// outputs
	static MObject aOutLinearTarget;
	static MObject aOutLinearTargetX;
	static MObject aOutLinearTargetY;
	static MObject aOutLinearTargetZ;

	static MObject aSmoothTan;
	static MObject aSmoothTanX;
	static MObject aSmoothTanY;
	static MObject aSmoothTanZ;

	static MObject aOut;
	static MObject aOutX;
	static MObject aOutY;
	static MObject aOutZ;

	// inputs
	static MObject aRotateOrder;
	static MObject aInTangent;
	static MObject aPrevVertex;
	static MObject aCurrentVertex;
	static MObject aNextVertex;
	static MObject aNextLinearTarget;
	static MObject aNextLinearTargetX;
	static MObject aNextLinearTargetY;
	static MObject aNextLinearTargetZ;
	static MObject aAuto;
	static MObject aSmooth;

	static MTypeId	id;
};



