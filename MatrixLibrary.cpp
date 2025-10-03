#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
using namespace std;
template<typename T>
class Matrix {
    private : 
    size_t rows,cols;
    T* vector ;
    T zero = 1e-12; // tolerance
    /////////////////////////////////////////////////////////////////////////Helper functions  /////////////////////////////////
    template<size_t M,size_t N>
    void copyData(const T (&src)[M][N]){                 //data copy
        for(size_t i = 0; i < rows;i++){
            for(size_t j = 0; j < cols ; j++){
                vector[i*cols + j] = src[i][j];
            }  
        }
    }
    
    ////////////////////////// helpers of detrminant function
    friend  Matrix minorMatrix(const Matrix& A , size_t r , size_t c){   //minor matrix
        Matrix minA;
        minA.rows = A.rows - 1;  minA.cols = A.cols - 1;
        minA.vector = new T[minA.rows*minA.cols];
        
        for(size_t i = 0; i < minA.rows ; i++ ){
            for(size_t j = 0 ; j < minA.cols ; j++){
                // minA[i,j] = A[x,y]// x = f(i,r)  // y = f(j,c) , makes sure we fill the minor with the correct elemets of A
                size_t  x = (i >= r) ? (i + 1) : i ;
                size_t  y = (j >= c) ? (j + 1) : j ;
                minA.vector[i*minA.cols + j] = A.vector[x * A.cols + y] ; 
            }
        }
        return minA;
    }
    
    friend Matrix cofI(const Matrix& matrix){           //cofactor just for computing inverse
        Matrix COF;
        COF.rows  = matrix.rows; COF.cols = matrix.cols;
        COF.vector = new T[COF.rows*COF.cols];
        for(size_t i = 0; i < COF.rows ; i++ ){
            for(size_t j = 0; j < COF.cols ; j++){
                COF.vector[i*COF.cols + j] = pow(-1,i+j)*DET(minorMatrix(matrix,i,j)); 
            }
        }
        return COF;
    }
    
    friend Matrix adjI(const Matrix& matrix){              //adjoint function just for computing inverse
        // assuming this is an n*n matrix that is not null
        Matrix ADJ = TR(cofI(matrix));
        return ADJ;
        
    }
    //////////////////// helpers of inverse funtion
    friend Matrix augmentedH(const Matrix& M){         // rrefAUG = [ A | I ]
        //creating augmented matrix for inverse computation only
        
        //create identity matrix compatible with the sqaure matrix M
        Matrix Identity(M.rows,M.cols);   for( size_t i = 0; i < Identity.rows ; i++ ) { Identity.vector[i*Identity.cols + i] = 1;}
        Matrix AUG; AUG.rows = M.rows; AUG.cols = 2*M.cols; AUG.vector = new T[AUG.rows*AUG.cols];
        
        for(size_t i = 0; i < AUG.rows ; i++ ){
            for(size_t j = 0; j < AUG.cols ; j++){
                if(j < M.cols){   //means we are in the matrix M
                    AUG.vector[i*AUG.cols+j] =  M.vector[i*M.cols+j] ;
                } else if(j > M.cols - 1){ // means we are in the identity matrix
                    AUG.vector[i*AUG.cols+j] =  Identity.vector[i*Identity.cols+j - M.cols] ; // becuse the column index of the identity elements in the AUG matrix are the indices of the indentity matrix shifted to the right by the number of M's columns
                }
            }
        }
        
        return AUG;
    }
    friend Matrix extractInvH(const Matrix& rrefAUG){    // rrefAUG = [I | A^-1]
        //creating augmented matrix for inverse computation only
        Matrix inverse; inverse.rows = rrefAUG.rows ; inverse.cols = (rrefAUG.cols)/2;
        inverse.vector = new T[inverse.rows*inverse.cols];
        
        for(size_t i = 0; i < inverse.rows; i++){
            for(size_t j = 0 ; j < inverse.cols ; j++){    // start form the columns of AUG that contain the inverse 
                inverse.vector[ i*inverse.cols  + j] = rrefAUG.vector[ i*rrefAUG.cols +(j + (rrefAUG.cols)/2)]; // the elements of th inverse in the AUG marix are shifted by the number of cols of the actual inverse we are filling by 1/2 the cols of AUG
            }
        }
        return inverse;
    }
    
    /////////////// ELEMENTARY ROW OPERATIONS HELPER FUNCTIONS + thier respective elementary matrices
    friend Matrix rSwapH(const Matrix& matrix, size_t r1 , size_t r2  ){ // no need to check if i,j are inbound or if matrix is null
        Matrix result(matrix);
        T temp{};
        for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
            temp =  result.vector[ r1*result.cols + k];
            result.vector[r1*result.cols + k] = result.vector[r2*result.cols + k];
            result.vector[r2*result.cols + k] = temp;
        }
        return result;
    }
    
    friend Matrix swapE(size_t r1, size_t r2,const Matrix& A){
        // define an identity matrix  m*m matrix compatible with A
        Matrix Identity(A.rows,A.rows);   for( size_t i = 0; i < Identity.rows ; i++ ) { Identity.vector[i*Identity.cols + i] = 1;}
        Matrix Ematrix = rSwapH(Identity,r1,r2);
        return Ematrix;
    }
    // ------------------------------------------------------------------------   
    
    friend Matrix rAddH(const Matrix& matrix, size_t r1 , size_t r2 , T scalar ){ // Rr1 -----> Rr1 + scalar * Rr2
        Matrix result(matrix);
        for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
            result.vector[r1*result.cols + k] = result.vector[r1*result.cols + k] + scalar*result.vector[r2*result.cols + k];
        }
        return result;
    }
    
    friend Matrix addE(size_t r1, size_t r2,T scalar,const Matrix& A){
        // define an identity matrix  m*m matrix compatible with A
        Matrix Identity(A.rows,A.rows);   for( size_t i = 0; i < Identity.rows ; i++ ) { Identity.vector[i*Identity.cols + i] = 1;}
        
        Matrix Ematrix = rAddH(Identity,r1,r2,scalar);
        return Ematrix;
    }
    
    //------------------------------------------------------------------
    friend Matrix rScaleH(const Matrix& matrix, size_t r , T scalar ){ // Rr -----> scalar *Rr
        Matrix result(matrix);
        for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
            result.vector[r*result.cols + k] = scalar*result.vector[r*result.cols + k];
        }
        return result;
    }
    
    friend Matrix scaleE(size_t r , T scalar,const Matrix& A){
        // define an identity matrix  m*m matrix compatible with A
        Matrix Identity(A.rows,A.rows);   for( size_t i = 0; i < Identity.rows ; i++ ) { Identity.vector[i*Identity.cols + i] = 1;}
        
        Matrix Ematrix = rScaleH(Identity,r,scalar);
        return Ematrix;
    }
    
    //------------------------------------------------------------------------
    
    public :
    //////////////////////////////////////////////////constructors //////////////////////////////////
    
    Matrix() : rows(0),cols(0),vector(nullptr){}     //default constructor
    
    Matrix(size_t rows ,size_t cols){                       //parametrized constructor
        if((rows > 0)&&(cols > 0)){
            vector = new(std::nothrow) T[rows*cols](); 
            if(vector == nullptr){
                cout << "Matrix declaration has failed!, creating a NULL Matrix..." << endl;
                this->rows = 0; this->cols = 0;
            } 
            this->rows = rows;  this->cols = cols;  
        } else {
            vector = nullptr;  this->rows = 0; this->cols = 0;
            cout << "The number of rows and columns must be strictly greater than zero,creating NULL Matrix... " << endl;
        }
    }
    
    Matrix(const Matrix& other_matrix){                  //copy constructor
        if(other_matrix.vector == nullptr){
            rows = 0;  cols = 0;  //if the others vector is null, then i the size of this must be zero
            vector = nullptr;
        } else {  
            rows = other_matrix.rows ; cols = other_matrix.cols ;  
            vector = new T[rows*cols];    //we dont need to check the allocation here, it succeded in last matrix so it should now
            for(size_t i = 0; i < rows*cols ; i++){
                vector[i]  = other_matrix.vector[i];
            } 
        }
    }
    
    Matrix& operator=(const Matrix& other_matrix){       // copy assignment (overrides the left side)
        if(this != &other_matrix){
            if(other_matrix.vector == nullptr){
                delete[] vector;
                vector = nullptr; cols = 0; rows = 0;
                
            } else {  
                delete[] vector;    //we could be deleting a null ptr here
                cols = other_matrix.cols; rows = other_matrix.rows; 
                vector = new T[rows*cols];
                for(size_t i = 0; i < rows*cols ; i++){
                    vector[i]  = other_matrix.vector[i];
                }
            }  
        } // else is the left and right operand are the same , then (A = B) returns A
        return *this;
    }
    
    ~Matrix(){                       //destructor
        delete[] vector;
    }
    
    ///////////////////////////// operations on matrices //////////// in all operations, the operands and arguments are never modified, a temporary object is always returned
    
    Matrix operator+(const Matrix& other){                                       //matrix addition
        if((this->vector == nullptr)||(other.vector == nullptr)){
            cout << "One of the operands is a Null matrix "<< endl;
            Matrix nullmatrix;
            return nullmatrix;  
        } else if ((this->rows != other.rows)||(this->cols != other.cols)){
            cout << "The dimentions of the operand matrices must match ! "<< endl;
            Matrix nullmatrix;
            return nullmatrix; 
        } else {
            size_t resultrows = this->rows ; size_t resultcols  =  this->cols ;
            Matrix sum; 
            sum.rows = resultrows ;sum.cols = resultcols ; //because its the same as the others
            sum.vector = new T[resultrows*resultcols];  //no need to check here, if it worked for two arrays it should work for the sum
            for(size_t i = 0; i < resultrows*resultcols; i++){
                sum.vector[i] = this->vector[i] + other.vector[i];
            }
            return sum ;  
        }
        
    }
    
    Matrix operator-(const Matrix& other){                           //matrix subtraction
        if((this->vector == nullptr)||(other.vector == nullptr)){
            cout << "One of the operands is a Null matrix "<< endl;
            Matrix nullmatrix;
            return nullmatrix;  
        } else if ((this->rows != other.rows)||(this->cols != other.cols)){
            cout << "The dimentions of the operand matrices must match ! "<< endl;
            Matrix nullmatrix;
            return nullmatrix; 
        } else {
            size_t resultrows = this->rows ; size_t resultcols  =  this->cols ;
            Matrix sum; // this will be destroyed once the operation ends ??????  (yes but a temporary object is created that the copy constructor copies info into)
            sum.rows = resultrows ;sum.cols = resultcols ; //because its the same as the others
            sum.vector = new T[resultrows*resultcols];  //no need to check here, if it worked for two arrays it should work for the sum
            for(size_t i = 0; i < resultrows*resultcols; i++){
                sum.vector[i] = this->vector[i] - other.vector[i];
            }
            return sum ;  
        }
        
    }
    
    Matrix operator*(const Matrix& other){                          // matrix multiplication
        if((this->vector == nullptr)||(other.vector == nullptr)){
            cout << "One of the operands is a Null matrix "<< endl;
            Matrix nullmatrix;
            return nullmatrix;  
        } else if ((this->cols != other.rows)){
            cout << "The number of columns of the first operand must match the second operand's number of rows ! "<< endl;
            Matrix nullmatrix;
            return nullmatrix; 
        } else {
            size_t resultrows = this->rows ; size_t resultcols  =  other.cols ;
            Matrix product; 
            product.rows = resultrows ;product.cols = resultcols ; 
            product.vector = new T[resultrows*resultcols];  //no need to check here, if it worked for two arrays it should work for the product
            T temp{}; //zero element
            for(size_t i = 0; i < this->rows  ; i++){
                for(size_t j = 0 ; j < other.cols ; j++){
                    
                    for(size_t k = 0; k < this->cols; k++){
                        temp = temp + (this->vector[i*(this->cols) + k])*(other.vector[k*(other.cols) + j]); 
                    }
                    product.vector[ i*(product.cols) + j] = temp;
                    temp = T{};   //clear for another element
                }
            }
            return product ;  
        }
    }
    
    // scalar multiplication
    
    friend Matrix operator*(const Matrix& matrix,T scalar) {    // M*scalar
        if(matrix.vector ==  nullptr){
            cout << "the operand matrix is null " << endl;
            return matrix;
        } else {
            Matrix scaledMatrix;
            scaledMatrix.rows = matrix.rows;  scaledMatrix.cols = matrix.cols;
            scaledMatrix.vector = new T[scaledMatrix.rows*scaledMatrix.cols];
            for(size_t i = 0; i < (matrix.rows)*(matrix.cols) ; i++){
                scaledMatrix.vector[i] = matrix.vector[i]*scalar;
            }
            return scaledMatrix;
        }
    }
    
    friend Matrix operator*(T scalar , const Matrix& matrix) {    // scalar*M
        if(matrix.vector ==  nullptr){
            cout << "the operand matrix is null " << endl;
            return matrix;
        } else {
            Matrix scaledMatrix;
            scaledMatrix.rows = matrix.rows;  scaledMatrix.cols = matrix.cols;
            scaledMatrix.vector = new T[scaledMatrix.rows*scaledMatrix.cols];
            for(size_t i = 0; i < (matrix.rows)*(matrix.cols) ; i++){
                scaledMatrix.vector[i] = matrix.vector[i]*scalar;
            }
            return scaledMatrix;
        }
    }
    
    friend Matrix TR(const Matrix& matrix){                                 //matrix transpose
        if(matrix.vector == nullptr){
            cout << "Cannot transpose a null matrix" << endl;
            return matrix;
        } else {
            Matrix transposed;
            transposed.rows = matrix.cols; transposed.cols = matrix.rows;
            transposed.vector = new T[ transposed.rows*transposed.cols];
            for(size_t i = 0; i < transposed.rows ; i++ ){
                for(size_t j = 0; j < transposed.cols ; j++){
                    transposed.vector[i*transposed.cols + j]  = matrix.vector[j*matrix.cols + i];
                }
            }
            return transposed;
        }
    } 
    
    friend T DET(const Matrix& matrix){       // determinant using laplace expansion (this shit is TOO slow O(n!) complexity) you must use gaussian elemination instead O(n^3)( idk how yet)
        if(matrix.vector == nullptr){
            cout << "Cannot compute determinant of null matrix " << endl;
            T zero{};
            return  zero;
        } else if(matrix.rows != matrix.cols) { 
            cout <<" Cannot compute determinant of non-square matrix " << endl;
            T zero{};
            return  zero;
        } else if(matrix.rows == 2){
            T det = matrix.vector[0]*matrix.vector[matrix.cols + 1] - matrix.vector[1]*matrix.vector[matrix.cols];
            return det;
        } else {
            T det{};
            for(size_t j = 0; j < matrix.cols ; j++){   // expanding along the 0th row i = 0
                det =  det + (matrix.vector[j])*pow(-1,j)*DET(minorMatrix(matrix,0,j));
            }
            return det; 
        }
        
    }  
    
    
    /*friend Matrix inv(const Matrix& matrix){          inverse using determinant adjoint method , and det = 0 invertibility test
    T DET {};
    if(matrix.vector == nullptr){
    cout << "Cannot invert a null matrix"<<endl;
    Matrix nullmatrix;
    return nullmatrix;
    } else if (matrix.rows != matrix.cols){
    cout << "Cannot invert a non-square matrix"<<endl;
    Matrix nullmatrix;
    return nullmatrix;
    } else if (abs(DET = det(matrix)) < 1e-12){       // this may never work with complecx types, you must use rank to determine invertibiblity
    cout << "det(matrix) = zero, the matrix is not invertible " << endl;
    Matrix nullmatrix;
    return nullmatrix;
    } else {
    Matrix INV = (adjI(matrix)*(1/DET));
    return INV;
    }
    }*/
    
    
    friend Matrix COF(const Matrix& matrix){         //// cofactor matrix
        if(matrix.vector == nullptr){ 
            cout << "Cannot compute the cofactor matrix of a null matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot compute the cofactor matrix of a non-square matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix cof = cofI(matrix);
            return cof;
        }
    }
    
    friend Matrix ADJ(const Matrix& matrix){                ////// adjoint matrix 
        if(matrix.vector == nullptr){ 
            cout << "Cannot compute the adjoint matrix of a null matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot compute the adjoint matrix of a non-square matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix adj = adjI(matrix);
            return adj;
        }
    }
    
    /////////////////////////// elementary column operations //////////////////  must be maid avaible to the user
    
    /*  ///1 - elementary column operations, ( you must remove H and make just one version)
    friend Matrix cSwapH(const Matrix& matrix, size_t c1 , size_t c2  ){ // no need to check if i,j are inbound or if matrix is null
    Matrix result(matrix);   
    T temp{};
    for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
    temp =  result.vector[ k*result.cols + c1];
    result.vector[ k*result.cols + c1] = result.vector[ k*result.cols + c2];
    result.vector[ k*result.cols + c2] = temp;
    }
    return result;
    }
    
    friend Matrix cAddH(const Matrix& matrix, size_t c1 , size_t c2 , T scalar ){ // Cc1 -----> Cc1 + scalar * Cc2
    Matrix result(matrix);
    for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
    result.vector[k*result.cols + c1] = result.vector[k*result.cols + c1] + scalar*result.vector[ k*result.cols + c2];
    }
    return result;
    }
    
    friend Matrix cScaleH(const Matrix& matrix, size_t c , T scalar ){ // Cc -----> scalar * Cc
    Matrix result(matrix);
    for(size_t k = 0; k < result.cols ; k++){   // [i*result.cols + j] = [i][j]
    result.vector[k*result.cols + c] = scalar * result.vector[k*result.cols + c];
    }
    return result;
    }
    */
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////// ROW ECHELON FORM
    
    friend Matrix REF(const Matrix& M){        //non normalized pivots,   // row echelon form  //testert
        if(M.vector == nullptr) {
            cout << "Cannot find row echelon form of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            // we have 3 function matrices that represent all row operations on any row with any scalar 
            // addE(r1,r2,s,M)   scaleE(r,s,M)   swapE(r1,r2,M)  , multiplying by xxE*M, applies the elementary row operation on M
            /* here we always take the pivot to be the first non zero entry in the column just for simplification purposes*/ 
            
            Matrix RE(M);      size_t m = RE.rows ; size_t n = RE.cols;   
            T zero = 1e-12; // the tolerance
            size_t cR = 0 , cC = 0;  //start at the top left corner
            size_t pivot = 0;  // row index of the first nonzer entry in a column
            size_t f = 0;   // increments when the FIRST non zero entry in a column is found
            /////////////vector[i*n + j] = RE[i][j]
            while((cR < m)&&(cC < n)){
                pivot = 0; f = 0;
                for(size_t i = cR; i < m ; i++){ // search for the first nonzero entry in the current column , starting from first row
                    if((abs(RE.vector[i*n + cC]) > zero)&&(f < 1)){    //(RE.vector[i*n + cC] != zero)
                        pivot = i ; f = 1;
                    }
                }
                if(f == 1){   //means we found a pivot
                    if(pivot != cR){ RE = rSwapH(RE,cR,pivot);} // raise the pivot to the current row
                    if((abs(RE.vector[cR*n + cC]) > zero)&&(RE.vector[cR*n + cC] != 1)){    //RE.vector[cR*n + cC] != zero
                        RE = rScaleH(RE,cR,1/(RE.vector[cR*n + cC])); // if its not 0 or 1 , normalize the pivot (that we have raised to the current row)
                    }
                    for(size_t i = cR + 1; i < m ; i++){   //elemenating below the pivot (only executes if cR >= m -2)
                        RE = rAddH(RE,i,cR,-RE.vector[i*n + cC]);
                        // RE.vector[i*n + cC] = T{}; makes sure below the pivot is zero
                    }
                    cR++ ; cC++; //go to create the next pivot
                    // RE.printMatrix(3,3); for debuggign purposes
                    cout << endl;
                } else { /// skip to the next column since no pivots where found
                    if(cC < n){ 
                        cC++ ; // so that we dont exceed the matrix bounds if the last column is a zero column
                    }
                }
            }
            return RE;
        }
        
    } 
    friend Matrix RREF(const Matrix& M){        /// reduced row echelon form
        if(M.vector == nullptr) {
            cout << "Cannot find row echelon form of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            T zero = 1e-12; // tolerance
            Matrix RE(REF(M)); //start from ref of the matrix
            size_t m = RE.rows ; size_t n = RE.cols;
            //size_t cR = m - 1; //here we dont need cC, it is replaced with i, and we start the operation from the bottom left corner of the matrix
            size_t f = 0 , pivot = 0; //this is the column index of the pivot found in the current row
            
            //N = a - b = (m-1) - (-1) = m times to scan all of the m rows of the matrix
            
            for(int cR = m - 1 ; cR > -1 ; cR--  ){ // each time the loop repeates , the algo researches for a pivot in the row above the last one /// this loop will run exactly m times
                pivot = 0; f = 0; //used to capture the first pivot and say if we found pivot
                for(size_t j = 0; j < n ; j++){ //search the current row for pivots ( the first 1 in this row) travel form left to right
                    if((abs( 1 - RE.vector[cR*n+j] ) < zero)&&(f < 1)){  
                        f = 1;
                        pivot = j; // located at column j
                    }  
                }     // analyse search results 
                if(f == 1){ // found a pivot
                    // eliminating above the pivot 
                    
                    // N = a - b  + 1 = (cR - 1) -(0) + 1 = cR times to clear all entries above the pivot
                    //we clear entries starting form the one above the current row, and end at the oth row
                    for(int  i = cR -1 ; i >= 0 ; i-- ) {   // when we reach i = 0 we dont execute, because there is nothing to eleminate, we start above the pivot , i is the index of the row we are eleminating in
                        RE = rAddH(RE,i,cR,-RE.vector[ i*n+ pivot]);  // eleminate even the last row  ( i> =0) but cant use size_t here
                    } 
                } else {   // if no pivots found 
                    // search in the upper row
                }
            }
            return RE;
        }
    }
    
    friend size_t RANK(const Matrix& M){ 
        if(M.vector == nullptr) {
            cout << "Cannot find row echelon form of a null matrix " << endl;
            return 0;
        } else {
            T zero = 1e-12; // tolerance
            Matrix RE(REF(M)); //start from ref of the matrix
            size_t m = RE.rows ; size_t n = RE.cols;
            //size_t cR = m - 1; //here we dont need cC, it is replaced with i, and we start the operation from the bottom left corner of the matrix
            size_t f = 0 , pivot = 0; //this is the column index of the pivot found in the current row
            
            size_t pivotRows = 0;  //help calculate rank for the next function
            // cR > 0 because we must eleminate above all the pivot rows exept the topmost one
            
            for(int cR = m - 1 ; cR > -1 ; cR--  ){ // each time the loop repeates , the algo researches for a pivot in the row above the last one /// this loop will run exactly m times
                pivot = 0; f = 0; //used to capture the first pivot and say if we found pivot
                for(size_t j = 0; j < n ; j++){ //search the current row for pivots ( the first 1 in this row) travel form left to right
                    if((abs( 1 - RE.vector[cR*n+j] ) < zero)&&(f < 1)){  
                        f = 1;
                        pivot = j; // located at column j
                        
                    }  
                }     // analyse search results 
                if(f == 1){ // found a pivot
                    // eliminating above the pivot 
                    for(int  i = cR -1 ; i >= 0 ; i-- ) {   // when we reach i = 0 we dont execute, because there is nothing to eleminate, we start above the pivot , i is the index of the row we are eleminating in
                        RE = rAddH(RE,i,cR,-RE.vector[ i*n+ pivot]);  // eleminate even the last row  ( i> =0) but cant use size_t here
                    } 
                    pivotRows++;
                    
                } else {   // if no pivots found 
                    // search in the upper row
                }
            }
            return pivotRows;  // return number of pivot rows found which is the number of nonzero rows , hence the rank
        }
    }
    
    friend size_t NULLITY(const Matrix& M){ 
        if(M.vector == nullptr) {
            cout << "Cannot find nullity of a null matrix " << endl;
            return 0;
        } else {
            size_t nullity = M.cols - RANK(M);
            return nullity;
            //rank nullity theorem  n = rank + nullity
        }
    }
    
    friend Matrix INV(const Matrix& matrix){           /// matrix  inverse ( too slow too, you have to use guassian elemination to check invertibility and return the inverse, do not use determinants alltogether IF RANK IS FULL THEN IT IS INVERTIBLE, dont have to use if (RREF(A) == Identity matrix)
        if(matrix.vector == nullptr){
            cout << "Cannot invert a null matrix"<<endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot invert a non-square matrix"<<endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (RANK(matrix) != matrix.cols){       // this may never work with complex types, you must use rank to determine invertibiblity
            cout << "The matrix is not full rank , hence not invertible " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            size_t m = matrix.rows ; size_t n = matrix.cols;
            Matrix inverse(n,n);
            //create an augmented matrix
            Matrix AUG = augmentedH(matrix);
            //get rref of aug
            AUG = RREF(AUG);
            //extract the inverse
            inverse = extractInvH(AUG);
            return inverse;
        }
    }
    
    //the problem was a pivot ,"1" sometimes isnt exactry 1 even though it was printed as 1, because 62.02565 - 61.02565 doesnt always yield the exact 1 but 0.999999999--->19 zeroes
    
    
    /////////////////////////////////////////////// setters and getters //////////////////////////////
    template<size_t M,size_t N>
    void fillMatrixA(const T (&arr)[M][N]){               //setter using array 
        if(vector == nullptr){  //create memory to the matrix
            rows = M; cols = N;
            vector = new T[rows*cols];
            copyData(arr);
        } else if((rows == M) && (cols == N)) {
            copyData(arr);
        }  else {
            cout << "The size of the array must match the size of the Matrix" << endl;
        }
    }
    
    void fillMatrixUI(){                                             //fill matrix with user input                             //// fill matrix with user input, if null, create a new size //// if not null fill the existing dimentions, you cant make it  larger or smaller
        if(vector == nullptr){
            cout << "Enter the number of rows" << endl;
            cin >> rows;
            cout << "Enter the number of columns" <<endl;
            cin >> cols;
            if((rows > 0)&&(cols > 0)){
                vector = new(std::nothrow) T[rows*cols]; 
                if(vector == nullptr){
                    cout << "Matrix initialization failed! setting matrix to NULL..." << endl;
                    rows = 0; cols = 0;
                } 
                for(size_t i = 0; i < rows;i++){
                    for(size_t j = 0; j < cols ; j++){
                        cout << "Enter element a(" << i << "," << j << "): ";
                        cin >>  vector[i*cols + j] ;
                    }  
                }
            } else {
                vector = nullptr;  rows = 0; cols = 0;
                cout << "The number of rows and columns must be strictly greater than zero, setting to NULL... " << endl;
            }
        } else { 
            for(size_t i = 0; i < rows;i++){
                for(size_t j = 0; j < cols ; j++){
                    cout << "Enter element a(" << i << "," << j << "): ";
                    cin >>  vector[i*cols + j] ;
                }  
            }
        }
    }
    template<typename Functype>   //this is so cool, becuase (int)func_a(int,int) and (void) func_b(double) and any other combination is a certain datatype the compiler decides automatically with auto or in this case the typename 
    void fillMatrixF(Functype customFiller){ //any callable works , but when declaing the function you must write Matrix<type>* in the argument list
        if(vector == nullptr){
            cout << "cannot fill a null matrix" << endl;
        } else {
            customFiller(this); //to fill the matrix , the custom handler must use setElement ,and getElement, and rowNbr , and colNbr
        }
        
        
        
    }
    T getElement (size_t r,size_t c) const{                //////// getter function ( get individual elements)
        if(vector == nullptr){
            cout << "This is a NULL Matrix" << endl ;
            T zero{}; //creating a zero element
            return zero;
        } else if( ( r >= 0 && r <= rows-1) && ( c >= 0 && c <= cols-1)){
            return vector[r*cols + c];
        } else {
            cout << "Choose a valid index for the elements of this matrix" <<endl;
            T zero{};
            return zero;
        }
    }                
    void setElement(size_t r,size_t c, T value){                     //// set individual elements 
        if(vector == nullptr){
            cout << "Cannot set the elements of a null matrix" << endl ;
        } else if( ( r >= 0 && r <= rows-1) && ( c >= 0 && c <= cols-1)){
            vector[r*cols + c] = value;
        } else {
            cout << "Choose a valid index for the elements of this matrix" <<endl;
        }
    }
    
    size_t rowNbr(){        //get number of rows
        return rows;    
    }
    size_t colNbr(){        //get number of columns
        return cols;
    }
    /*T* getVector(){         //this shit is dangarous, only if you know what you are doing 
        return this->vector;
    }*/
    /// //////////////////////////////////////////////diplayers ////////////////////////////////////
    
    void  printMatrix(size_t fieldlen = 3,size_t precision = 3) const {    //number of  total digits of each number
        if(precision > fieldlen){ cout << " precision should be smaller than the field lenght " << endl;}
        
        if(vector == nullptr){
            cout << "Cannot print an empty matrix" << endl;
        } else {
            size_t cursor = 0;
            for(size_t i = 0; i < rows*cols ; i++ ){
                if(abs(vector[i]) < 1e-12){cout << setw(fieldlen) << setprecision(precision) << T{} <<"";} else {
                    cout << setw(fieldlen) << setprecision(precision) << vector[i] <<"";
                }
                cursor++;
                if(cursor == cols){
                    cout << endl ;
                    cursor = 0;
                }
            }
        }     
    }
    
};


int main(){
    
    /*
    template<typename Functype>   
    void fillMatrixF(Functype customFiller){ 
    if(vector == nullptr){
    cout << "cannot fill a null matrix" << endl;
    } else {
    customHandler(this); 
    }
    */
    Matrix<int> A(3,3);
    A.fillMatrixF([](Matrix<int>* M ){
        size_t m = M->rowNbr(); size_t n = M->colNbr(); // try to get rid of it by passing them in the fill matrixF definition, this is highly unstable, you best bet is friend functions, (because you will have to write too many arguments for matrix*, vector, row...)
        for(size_t i = 0; i < m;i++){
            for(size_t j = 0; j < n ; j++){
                M->setElement(i,j,5);
            }  
        }
    });
    
    A.printMatrix();
    
    
    
    return 0;   // and then systems of linear equations, polynomial class 
}   









/*//MATRIX-VECTOR MULTIPLICATION
Matrix<double> A(2,3);
Matrix<double> x(3,1);
Matrix<double> b;

A.fillMatrixUI();
x.fillMatrixUI();

cout << "The matrix A is :" <<endl ;
A.printMatrix();
cout << " the vector x is :" << endl;
x.printMatrix();
b = A*x;
cout << "the result Ax = b is :" << endl;
b.printMatrix();*/


//A.fillMatrixA(data_array);
/*
cout << "the A is :"<< endl;
A.printMatrix(2,1);

cout << "the determinant of A is : " << DET(A) << endl;
cout <<"the rank of A is : " << RANK(A) << endl;
cout <<"the nullity of A is : " << NULLITY(A) << endl;

cout << "inverse of A is :" << endl;
INV(A).printMatrix(9,4);

cout << "transpose of A is :" << endl;
TR(A).printMatrix(9,4);

cout << "the cofactor matrix of A is :"<< endl;
COF(A).printMatrix(9,4);

cout << "the adjoint matrix of A is :" << endl;
ADJ(A).printMatrix(9,4);

cout << "row echelon form of A is :" << endl;
REF(A).printMatrix(9,4);

cout << "reduced row echelon form of A is :" << endl;
RREF(A).printMatrix(9,4);

(A*INV(A)).printMatrix(2,1);
(INV(A)*A).printMatrix(2,1);

*/  //test program
/*double data_array[6][6] = {
{ 1,  2,  3,  4,  5,  6 },
{ 0,  1,  4,  7, 10, 13 },
{ 2,  0,  1,  3,  5,  7 },
{ 1,  1,  0,  2,  4,  6 },
{ 3,  5,  2,  1,  0,  8 },
{ 4,  2,  1,  0,  3,  9 }
};*/
