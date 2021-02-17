/*
 * RMSD.cpp
 *
 *  Created on: Feb 6, 2019
 *      Author: s2982206
 */

#include <geometry/RMSD.h>

namespace NSPgeometry {

TransForm::~TransForm(){
//	cout << "delete transForm" << endl;
}

void jRotate(double** vecA, int i, int j, int k, int l, double* gh, double s, double tau) {
    gh[0] = vecA[i][j];
    gh[1] = vecA[k][l];

    vecA[i][j] = gh[0] - s*(gh[1] + gh[0]*tau);
    vecA[k][l] = gh[1] + s*(gh[0] - gh[1]*tau);
}

void computeJacobi(double** vecA, double* vecD, double** vecV, int n){


    //  Computes all eigenvalues and eigenvectors of vecA real symmetric matrix
    //  vecA[0..n-1][0..n-1]. On output, elements of vecA above the diagonal are
    //  destroyed. vecD[0..n-1] returns the eigenvalues of vecA. vecV[0..n-1][0..n-1]
    //  is vecA matrix whose columns contain, on output, the normalized
    //  eigenvectors of vecA. iNumRot returns the number of Jacobi rotations
    //  that were required.  Note that vecD and vecV must be allocated prior
    //  to calling this function

    int max_sweeps = 50;

    // Declare iterators
    int j=0,         // iterator in inner loop for applying rotations
    iQ=0,             // column/q iterator
    iP=0,             // row/p iterator
    i=0;              // iterator over sweeps

    double threshold=0,    // Threshold for performing vecA pq rotation
    theta=0,          //rotation angle
    tau=0,            // sin(theta)/(1 + cos(theta))
    t=0,              // 1/(2*theta)
    sum=0,            // sum of off-diagonal elements
    s=0,              // sin(theta)
    h=0,              // a_pp - a_qq, which is used to determine theta
    g=0,              // 100*current off-diagonal element (used to check threshold)
    c=0;              // cos(theta)

    double gh[2];

    double b[n];
    double z[n]; // This array will accumulate terms
                   // of the form tapq as in equation (11.1.14)

    // Initialize eigenvectors to identity matrix
    for (iP = 0; iP < n; ++iP)
    {
        for (iQ = 0; iQ < n; ++iQ)
            vecV[iP][iQ]=0.0;
        vecV[iP][iP]=1.0;
    }

    for (iP = 0; iP < n; ++iP)
    {
        b[iP]=vecD[iP]=vecA[iP][iP];  // Initialize b and vecD to the diagonal of vecA
        z[iP] = 0.0;                  // initialize z to 0
    }

    // Start counting rotations
    int iNumRot = 0;
    for (i = 0; i < max_sweeps; ++i) // loop over sweeps
    {
        // Sum off-diagonal elements (of one triangle)
        sum = 0.0;
        for (iP=0;iP<n-1;++iP)
        {
            for (iQ = iP+1; iQ< n; ++iQ)
                sum += abs(vecA[iP][iQ]);
        }

        // Normal return: when off-diagonal elements sum to zero,
        //                within machine precision, we consider
        //                the matrix diagonalized.  For non-degenerate
        //                eigenvalues, we should have quadratic convergence
        //                to this state.

        if (sum == 0.0)
            return ; // success
        if (i < 3)
            threshold = 0.2*sum/(n*n);   //...on the 1st three sweeps.
        else
            threshold = 0.0;                //...thereafter.

        // Begin sweep over all possible pq rotations
        for (iP = 0; iP < n-1; ++iP)
        {
            for (iQ = iP+1; iQ < n; ++iQ)
            {
                g = 100.0 * abs(vecA[iP][iQ]);

                // After four sweeps, skip the rotation if
                // the off-diagonal element is small.
                if((i > 3) &&
                    (abs(vecD[iP])+g == abs(vecD[iP])) &&
                    (abs(vecD[iQ])+g == abs(vecD[iQ])))
                {
                    vecA[iP][iQ] = 0.0;
                }
                // perform pq rotation
                else if (abs(vecA[iP][iQ]) > threshold)
                {
                    // check how big the rotation needs to be
                    h = vecD[iQ] - vecD[iP];

                    // Perform pre-calculations and apply rotation so as to
                    // minimize roundoff error.
                    if (abs(h)+g == abs(h))
                        t = (vecA[iP][iQ]) / h;
                    else
                    {
                        theta = 0.5*h / (vecA[iP][iQ]); // Equation (11.1.10).
                        t = 1.0 / (abs(theta) + sqrt(1.0 + theta*theta) );
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1.0 + t*t);  // cos(theta)
                    s = t * c;              // sin(theta)
                    tau = s / (1.0 + c);
                    h = t * vecA[iP][iQ];
                    z[iP] -= h;
                    z[iQ] += h;
                    vecD[iP] -= h;
                    vecD[iQ] += h;
                    vecA[iP][iQ] = 0.0;
                    for (j = 0; j < iP; ++j) // Case of rotations 0 <= j <p.
                        jRotate(vecA,j,iP,j,iQ,gh,s,tau);
                    for (j = iP+1; j < iQ; ++j) // Case of rotations p <j<q.
                    	jRotate(vecA,iP,j,j,iQ,gh,s,tau);
                    for (j = iQ+1; j < n; ++j) // Case of rotations q <j < n.
                    	jRotate(vecA,iP,j,iQ,j,gh,s,tau);
                    for (j = 0; j < n; ++j) // update eigenvectors
                    	jRotate(vecV,j,iP,j,iQ,gh,s,tau);
                    ++iNumRot; // count this rotation
                }
            }
        }

        // Tidy up after current sweep
        for (iP = 0; iP < n; ++iP)
        {
            b[iP] += z[iP];
            vecD[iP] = b[iP];     // Update vecD with the sum of tapq,
            z[iP]= 0.0;     // and reinitialize z.
        }
    }
}

TransForm buildRotation(vector<XYZ>& points1, vector<XYZ>& points2){
    XYZ m,p;
    double** Q = new double*[4]; //all zero
    double** R = new double*[4]; //..
    double* ev = new double[4];
    for(int i=0;i<4;i++){
    	ev[i] = 0;
    	Q[i] = new double[4];
    	R[i] = new double[4];
    	for(int j=0;j<4;j++){
    		Q[i][j] = 0.0;
    		R[i][j] = 0.0;
    	}
    }

    int len = points1.size();

    // Calculate the P matrix from each pair of points
    for(int i=0;i<len && i<len; i++)
    {

        m = points1[i] - (points2[i]);    // vector between the points
        p = points1[i] + (points2[i]);    // sum of position vectors

        Q[0][0] += m.x_*m.x_+m.y_*m.y_+m.z_*m.z_;
        Q[1][1] += p.y_*p.y_+p.z_*p.z_+m.x_*m.x_;
        Q[2][2] += p.x_*p.x_+p.z_*p.z_+m.y_*m.y_;
        Q[3][3] += p.x_*p.x_+p.y_*p.y_+m.z_*m.z_;
        Q[0][1] += p.y_*m.z_-m.y_*p.z_; Q[1][0] = Q[0][1];
        Q[0][2] += m.x_*p.z_-p.x_*m.z_; Q[2][0] = Q[0][2];
        Q[0][3] += p.x_*m.y_-m.x_*p.y_; Q[3][0] = Q[0][3];
        Q[1][2] += m.x_*m.y_-p.x_*p.y_; Q[2][1] = Q[1][2];
        Q[1][3] += m.x_*m.z_-p.x_*p.z_; Q[3][1] = Q[1][3];
        Q[2][3] += m.y_*m.z_-p.y_*p.z_; Q[3][2] = Q[2][3];
    }

    computeJacobi(Q, ev, R, 4);

    // Each eigenvalue corresponds to the deviation of one structure
    // from the other, after the rotation. The smallest eigenvalue,
    // therefore, must be the best alignment.

    int j = 0;    // generic iterator for eigenvectors/values
    for (int i = 0; i < 4; i++)
        if (ev[i] < ev[j])
            j = i;    // i is a smaller eigenvalue

    // Step 3: Build the rotation matrix from the eigenvector

    XYZ q(R[1][j],R[2][j],R[3][j]);

    TransForm tf(q, R[0][j]);

    delete[] ev;
    for(int i=0;i<4;i++){
    	delete[] Q[i];
    	delete[] R[i];
    }
    delete[] Q;
    delete[] R;


    return tf;
}

double simpleRMSD(const vector<XYZ>& points1, const vector<XYZ>& points2){
	if(points1.size() != points2.size()) return -1;
	int len= points1.size();
	if(len == 0) return 0;
	double dd = 0;
	for(int i=0;i<len;i++){
		dd += squareDistance(points1[i], points2[i]);
	}
	return sqrt(dd/len);
}

XYZ getCOG(const vector<XYZ>& points){
	double x=0, y=0, z=0;
	int len = points.size();
	for(int i=0;i<len;i++){
		x += points[i].x_;
		y += points[i].y_;
		z += points[i].z_;
	}
	double w = 1.0/len;
	return XYZ(x*w, y*w, z*w);
}


double rmsd(const vector<XYZ>& points1, const vector<XYZ>& points2){
	XYZ Acog = getCOG(points1);
	XYZ Bcog = getCOG(points2);
	if(points1.size() != points2.size()) return -1;
	int len= points1.size();
	if(len == 0) return 0;

	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<len;i++){
		XYZ a = points1[i] - Acog;
		XYZ b = points2[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}

	TransForm tf = buildRotation(listA, listB);
	vector<XYZ> newListB;
	for(int i=0;i<len;i++) {
		XYZ c = tf.transform(listB[i]);
		newListB.push_back(c);
	}
	return simpleRMSD(listA, newListB);
}

double subRMSD(const vector<XYZ>& points1, const vector<XYZ>& points2, const vector<XYZ>& subList1, const vector<XYZ>& subList2){
	XYZ Acog = getCOG(points1);
	XYZ Bcog = getCOG(points2);
	if(points1.size() != points2.size()) return -1;
	int len= points1.size();
	if(len == 0) return 0;

	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<len;i++){
		XYZ a = points1[i] - Acog;
		XYZ b = points2[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}

	vector<XYZ> subA;
	vector<XYZ> subB;
	for(int i=0;i<subList1.size();i++){
		XYZ a = subList1[i] - Acog;
		XYZ b = subList2[i] - Bcog;
		subA.push_back(a);
		subB.push_back(b);
	}

	TransForm tf = buildRotation(listA, listB);

	vector<XYZ> newSubB;
	for(int i=0;i<subList1.size();i++){
		XYZ c = tf.transform(subB[i]);
		newSubB.push_back(c);
	}

	return simpleRMSD(subA, newSubB);

}

double rmsd2(const vector<XYZ>& points1, const vector<XYZ>& points2){
	cout << "rmsd2 " << endl;
	XYZ Acog = getCOG(points1);
	XYZ Bcog = getCOG(points2);
	if(points1.size() != points2.size()) return -1;
	int len= points1.size();
	if(len == 0) return 0;

	cout << "len: " << len << endl;

	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<len;i++){
		XYZ a = points1[i] - Acog;
		XYZ b = points2[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}

	for(int i=0;i<points1.size();i++) {
		cout << points1[i].toString() << endl;
	}
	cout << "cogA " << Acog.toString() << endl;

	for(int i=0;i<points2.size();i++) {
		cout << points2[i].toString() << endl;
	}
	cout << "cogB " << Bcog.toString() << endl;


	cout << listA[0].toString() << endl;
	cout << listB[0].toString() << endl;

	cout << "tf" << endl;
	TransForm tf = buildRotation(listA, listB);
	cout << "tf2" << endl;
	tf.print();
	cout << listB[0].toString() << endl;


	vector<XYZ> newListB;
	for(int i=0;i<len;i++) {
		cout << listB[i].toString() << endl;
		XYZ t = tf.transform(listB[i]);
		cout << i << " " << t.toString() << endl;
		newListB.push_back(tf.transform(listB[i]));
	}

	for(int i=0;i<listA.size();i++) {
		printf("%s %s %s\n",listA[i].toString().c_str(), listB[i].toString().c_str(), newListB[i].toString().c_str());
	}
	return simpleRMSD(listA, newListB);
}


} /* namespace NSPgeometry */
