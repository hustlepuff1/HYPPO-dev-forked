/*=========================================================================================================*/
/*-------------------------------------Diamond computation method class------------------------------------*/
/*=========================================================================================================*/

class Diamond
{
	private:

		// --- Two non colinear vectors on a face
		vector VecFace1;
		vector VecFace2;

		// --- Links two cells centers
		vector VecCenters;

		// --- Difference of value on VecFace 1 
		vector U1;
		scalar T1;
	
		// --- Difference of value on VecFace 2
		vector U2;
		scalar T2;

		// --- Difference of value on VecCenters
		vector UC;
		scalar TC;
		
		// --- Vector normal to the face
		vector normalToFace;

		// --- Tensor containing vectors Vec...
		tensor A;
	
		// --- Tensor with value of U
		tensor b;

		// --- Vector with value of T
		vector bT;

		// --- Velocity gradient
		tensor GU_;
	
		// --- Temperature gradient
		vector GT_;

	public:
		// --- Computes the gradient with a diamond scheme
		void diamondFace(fvMesh& mesh, const volVectorField& C, volVectorField& U, volScalarField& T, const faceList& faces, const pointField& points, pointVectorField& Up, pointScalarField& Tp, const unallocLabelList& owner, const unallocLabelList& neighbour, const surfaceScalarField& magSf, surfaceVectorField& GradU, surfaceScalarField& GradT);

		// --- Computes the gradient on boundaries
		void diamondBoundary(fvMesh& mesh, const volVectorField& C, volVectorField& U, volScalarField& T, const faceList& faces, const pointField& points, pointVectorField& Up, pointScalarField& Tp, const unallocLabelList& owner, const unallocLabelList& neighbour, const surfaceScalarField& magSf, surfaceVectorField& GradU, surfaceScalarField& GradT, const polyBoundaryMesh& boundaryMesh, const surfaceVectorField& Cf );
};
