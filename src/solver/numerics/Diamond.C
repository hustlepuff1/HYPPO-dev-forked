/*=========================================================================================================*/
/*------------------------------------------Functions definitions------------------------------------------*/
/*=========================================================================================================*/

#include "Diamond.H"

// --- Computes the gradient with a diamond scheme
void Diamond::diamondFace(fvMesh& mesh, const volVectorField& C, volVectorField& U, volScalarField& T, const faceList& faces, const pointField& points, pointVectorField& Up, pointScalarField& Tp, const unallocLabelList& owner, const unallocLabelList& neighbour, const surfaceScalarField& magSf, surfaceVectorField& GradU, surfaceScalarField& GradT)
{
	forAll(faces,face)
	{
		if (mesh.isInternalFace(face)) // --- Internal face is found
		{
			const label own = owner[face];
			const label nei = neighbour[face];

			VecFace1 = points[faces[face][1]] - points[faces[face][0]];			// --- Defines two non colinear vectors on a face
			VecFace2 = points[faces[face][2]] - points[faces[face][0]];
			VecCenters = C[nei] - C[own];  										// --- Defines vector linking two cell centers
	
			U1 = Up[faces[face][1]] - Up[faces[face][0]];  						// --- Difference of value on VecFace 1 and 2
			U2 = Up[faces[face][2]] - Up[faces[face][0]];
			UC = U[nei] - U[own];  												// --- Difference of value on VecCenters

			T1 = Tp[faces[face][1]] - Tp[faces[face][0]];  						// --- Same for the temperature
			T2 = Tp[faces[face][2]] - Tp[faces[face][0]];
			TC = T[nei] - T[own];
;
			A.xx()=VecCenters[0]; A.xy()=VecCenters[1]; A.xz()=VecCenters[2]; 	// --- Tensor with vectors inside
			A.yx()=VecFace1[0]; A.yy()=VecFace1[1]; A.yz()=VecFace1[2];
			A.zx()=VecFace2[0]; A.zy()=VecFace2[1]; A.zz()=VecFace2[2]; 

			b.xx()=UC[0]; b.xy()=UC[1]; b.xz()=UC[2];							// --- Tensor with value of U (because U is a vector)
			b.yx()=U1[0]; b.yy()=U1[1]; b.yz()=U1[2];
			b.zx()=U2[0]; b.zy()=U2[1]; b.zz()=U2[2];		
 						
			bT.x()=TC; bT.y()=T1; bT.z()=T2; 									// --- Vector with value of T (because T is a scalar)

			GU_ = inv(A) & b;  													// --- We compute the velocity gradient
			GT_ = inv(A) & bT;  												// --- We compute the temperature gradient

			normalToFace =  mesh.Sf()[face]/magSf[face];  						// --- Compute the normal to the considered face
			GradU[face] = GU_.T() & normalToFace;  								// --- Project GU and GT normal to the face
			GradT[face] = GT_ & normalToFace;
		}
	}
};


// --- Computes the gradient on boundaries
void Diamond::diamondBoundary(fvMesh& mesh, const volVectorField& C, volVectorField& U, volScalarField& T, const faceList& faces, const pointField& points, pointVectorField& Up, pointScalarField& Tp, const unallocLabelList& owner, const unallocLabelList& neighbour, const surfaceScalarField& magSf, surfaceVectorField& GradU, surfaceScalarField& GradT,const polyBoundaryMesh& boundaryMesh, const surfaceVectorField& Cf)
{
	forAll(mesh.boundary(), patch)
	{
		forAll(mesh.boundary()[patch], facei)
		{
			const label& face = boundaryMesh[patch].start() + facei;        // --- Face index
			const label& own = mesh.owner()[face];

			const vectorField& Upb = U.boundaryField()[patch];  // --- Value at the face center at boundary
			const scalarField& Tpb = T.boundaryField()[patch];

			// --- The next is the same as before except that we use a face center for boundary
			VecFace1 = points[faces[face][1]] - points[faces[face][0]];
			VecFace2 = points[faces[face][2]] - points[faces[face][0]];
			VecCenters = Cf[face] - C[own];

			U1 = Up[faces[face][1]] - Up[faces[face][0]];
			U2 = Up[faces[face][2]] - Up[faces[face][0]];
			UC = Upb[facei] - U[own];

			T1 = Tp[faces[face][1]] - Tp[faces[face][0]];
			T2 = Tp[faces[face][2]] - Tp[faces[face][0]];
			TC = Tpb[facei] - T[own];

			A.xx()=VecCenters[0]; A.xy()=VecCenters[1]; A.xz()=VecCenters[2]; 	// --- Tensor with vectors inside
			A.yx()=VecFace1[0]; A.yy()=VecFace1[1]; A.yz()=VecFace1[2];
			A.zx()=VecFace2[0]; A.zy()=VecFace2[1]; A.zz()=VecFace2[2];

			b.xx()=UC[0]; b.xy()=UC[1]; b.xz()=UC[2];							// --- Tensor with value of U (because U is a vector)
			b.yx()=U1[0]; b.yy()=U1[1]; b.yz()=U1[2];
			b.zx()=U2[0]; b.zy()=U2[1]; b.zz()=U2[2];

			bT.x()=TC; bT.y()=T1; bT.z()=T2; 

			normalToFace = mesh.Sf()[face]/mag(mesh.Sf()[face]);

			GT_ = inv(A) & bT;
			GradT.boundaryField()[patch][facei] = GT_ & normalToFace;

			GU_ = inv(A) & b;
			GradU.boundaryField()[patch][facei] = GU_.T() & normalToFace;
		}
	}
};
