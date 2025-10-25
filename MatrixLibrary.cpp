//declaring dependencies
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
    T tol = 1e-12; // tolerance
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
    friend Matrix augmentedH(const Matrix& M){         //AUG = [ A | I ]
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
            cout << "The number of rows and columns must be strictly greater than tol,creating NULL Matrix... " << endl;
        }
    }
    
    Matrix(const Matrix& other_matrix){                  //copy constructor
        if(other_matrix.vector == nullptr){
            rows = 0;  cols = 0;  //if the others vector is null, then i the size of this must be tol
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
            T temp{}; //tol element
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
    
    Matrix operator^(int k){  // matrix power A^k
        
        Matrix nullmatrix;
        if(vector == nullptr){
            cout << "the operand is a null matrix" << endl;
            return nullmatrix;
        } else if ( rows != cols){
            cout << "matrix powers are not defined for non-square matrices"<<endl;
            return nullmatrix;
        } else if(k < 0){
            cout << " the power must be a positive integer" << endl;
            return nullmatrix;
        } else if( k == 0){ 
            Matrix Identity(rows,cols);   for( size_t i = 0; i < Identity.rows ; i++ ) { Identity.vector[i*Identity.cols + i] = 1;}
            return Identity;
        } else if (k == 1){
            return *this;
        } else {
            Matrix result = (*this);
            for(size_t i = 0 ; i < k - 1; i++){
                result = (*this)*result; 
            }
            return result;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////// ROW ECHELON FORM
    
    friend Matrix REF(const Matrix& M){        //normalized pivots,   // row echelon form  //testert
        if(M.vector == nullptr) {
            cout << "Cannot find row echelon form of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            // we have 3 function matrices that represent all row operations on any row with any scalar 
            // addE(r1,r2,s,M)   scaleE(r,s,M)   swapE(r1,r2,M)  , multiplying by xxE*M, applies the elementary row operation on M
            /* here we always take the pivot to be the first non tol entry in the column just for simplification purposes*/ 
            
            Matrix RE(M);      size_t m = RE.rows ; size_t n = RE.cols;   
            T tol = 1e-12; // the tolerance
            size_t cR = 0 , cC = 0;  //start at the top left corner
            size_t pivot = 0;  // row index of the first nonzer entry in a column
            size_t f = 0;   // increments when the FIRST non tol entry in a column is found
            /////////////vector[i*n + j] = RE[i][j]
            while((cR < m)&&(cC < n)){
                pivot = 0; f = 0;
                for(size_t i = cR; i < m ; i++){ // search for the first nonzero entry in the current column , starting from current row
                    if((abs(RE.vector[i*n + cC]) > tol)&&(f < 1)){    //(RE.vector[i*n + cC] != 0)
                        pivot = i ; f = 1;
                    }
                }
                if(f == 1){   //means we found a pivot
                    if(pivot != cR){ RE = rSwapH(RE,cR,pivot);} // raise the pivot to the current row   Rpiovt <-----> Rcr
                    if((abs(RE.vector[cR*n + cC]) > tol)&&(RE.vector[cR*n + cC] != 1)){    //RE.vector[cR*n + cC] != tol
                        RE = rScaleH(RE,cR,1/(RE.vector[cR*n + cC])); // if its not 0 or 1 , normalize the pivot (that we have raised to the current row)
                    }  // Rcr ----> 1/(RE[cR][cC]) * Rcr
                    for(size_t i = cR + 1; i < m ; i++){   //elemenating below the pivot (only executes if cR >= m -2)
                        RE = rAddH(RE,i,cR,-RE.vector[i*n + cC]);  // Ri ---> Ri - RE[i][cC]* Rcr
                        // RE.vector[i*n + cC] = T{}; makes sure below the pivot is tol
                    }
                    cR++ ; cC++; //go to create the next pivot
                    
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
            T tol = 1e-12; // tolerance
            Matrix RE(REF(M)); //start from ref of the matrix
            size_t m = RE.rows ; size_t n = RE.cols;
            //size_t cR = m - 1; //here we dont need cC, it is replaced with i, and we start the operation from the bottom left corner of the matrix
            size_t f = 0 , pivot = 0; //this is the column index of the pivot found in the current row
            
            //N = a - b = (m-1) - (-1) = m times to scan all of the m rows of the matrix
            
            for(int cR = m - 1 ; cR > -1 ; cR--  ){ // each time the loop repeates , the algo researches for a pivot in the row above the last one /// this loop will run exactly m times
                pivot = 0; f = 0; //used to capture the first pivot and say if we found pivot
                for(size_t j = 0; j < n ; j++){ //search the current row for pivots ( the first 1 in this row) travel form left to right
                    if((abs( 1 - RE.vector[cR*n+j] ) < tol)&&(f < 1)){  
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
            T tol = 1e-12; 
            Matrix RE(REF(M)); 
            size_t m = RE.rows ; size_t n = RE.cols;
            size_t f = 0 , pivot = 0; 
            size_t pivotRows = 0;  //help calculate rank 
            for(int cR = m - 1 ; cR > -1 ; cR--  ){ 
                pivot = 0; f = 0; 
                for(size_t j = 0; j < n ; j++){
                    if((abs( 1 - RE.vector[cR*n+j] ) < tol)&&(f < 1)){  
                        f = 1;
                        pivot = j; 
                    }  
                }     
                if(f == 1){ // found a pivot
                    for(int  i = cR -1 ; i >= 0 ; i-- ) {   
                        RE = rAddH(RE,i,cR,-RE.vector[ i*n+ pivot]);
                    } 
                    pivotRows++;
                    
                } else {   // if no pivots found 
                    
                }
            }
            return pivotRows;  // return number of pivot rows found which is the number of nontol rows , hence the rank
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
    
    friend T DETl(const Matrix& matrix){       // determinant using laplace expansion (this shit is TOO slow O(n!) complexity) you must use gaussian elemination instead O(n^3)( idk how yet)
        if(matrix.vector == nullptr){
            cout << "Cannot compute determinant of null matrix " << endl;
            T tol{};
            return  tol;
        } else if(matrix.rows != matrix.cols) { 
            cout <<" Cannot compute determinant of non-square matrix " << endl;
            T tol{};
            return  tol;
        } else if(matrix.rows == 2){
            T det = matrix.vector[0]*matrix.vector[matrix.cols + 1] - matrix.vector[1]*matrix.vector[matrix.cols];
            return det;
        } else {
            T det{};
            for(size_t j = 0; j < matrix.cols ; j++){   // expanding along the 0th row i = 0
                det =  det + (matrix.vector[j])*pow(-1,j)*DETl(minorMatrix(matrix,0,j));
            }
            return det; 
        }
        
    }
    
    friend T DET(const Matrix& matrix){       // determinant using ref guassian elemination    O(n^3) complexity
        if(matrix.vector == nullptr){
            cout << "Cannot compute determinant of null matrix " << endl;
            T tol{};
            return  tol;
        } else if(matrix.rows != matrix.cols) { 
            cout <<" Cannot compute determinant of non-square matrix " << endl;
            T tol{};
            return  tol;
        } else if(RANK(matrix) < matrix.cols){  // singular matrix
            T det{};
            return det;
        } else {     //if( RANK(matrix) == matrix.cols) 
            T det{};
            T x = 1; // product of all scaling factors k used in row scaling operations to get ref(matrix)
            size_t p = 0; //number of row swap operations
            //////////////////////
            
            Matrix RE(matrix);      size_t m = RE.rows ; size_t n = RE.cols;   
            T tol = 1e-12; // the tolerance
            size_t cR = 0 , cC = 0; 
            size_t pivot = 0; 
            size_t f = 0;  
            
            while((cR < m)&&(cC < n)){  
                pivot = 0; f = 0;
                for(size_t i = cR; i < m ; i++){ 
                    if((abs(RE.vector[i*n + cC]) > tol)&&(f < 1)){ 
                        pivot = i ; f = 1;
                    }
                }
                if(f == 1){  
                    if(pivot != cR){ RE = rSwapH(RE,cR,pivot); p++;}  //swap opp
                    if((abs(RE.vector[cR*n + cC]) > tol)&&(RE.vector[cR*n + cC] != 1)){   
                        x = x*(RE.vector[cR*n + cC]);    // we will mutiply the row by k = 1/pivot1, then by 1/pivot2 ... so x = 1/p1*1/p2*..., and the det will be euqal to det(A) = det(REF(A))*1/x*(-1)^p, and det(REF(A)) = 1 becuase this is a nonsingluar matrix, so just mutiply by x = p1*p2*..
                        RE = rScaleH(RE,cR,1/(RE.vector[cR*n + cC])); 
                        for(size_t i = cR + 1; i < m ; i++){  
                            RE = rAddH(RE,i,cR,-RE.vector[i*n + cC]);
                        }
                        cR++ ; cC++;  
                        
                    } else { 
                        if(cC < n){ 
                            cC++ ; 
                        }
                    }
                }  
            }
            det = x*pow(-1,p);
            return det; 
        }  
        
    }
    
    friend Matrix INVd(const Matrix& matrix){        //  inverse using determinant adjoint method , and det = 0 invertibility test
        T DET {};
        if(matrix.vector == nullptr){
            cout << "Cannot invert a null matrix"<<endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot invert a non-square matrix"<<endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (abs(DET = DETl(matrix)) < 1e-12){       // this may never work with complecx types, you must use rank to determine invertibiblity
            cout << "det(matrix) = tol, the matrix is not invertible " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix INV = (adjI(matrix)*(1/DET));
            return INV;
        }
    }
    
    
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
    
    
    /////////////////////////// elementary column operations //////////////////  must be maid avaible to the user (add bounds checking)
    
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
    
    //the problem was a pivot ,"1" sometimes isnt exactlzzzzzfzy 1 even though it was printed as 1, because 62.02565 - 61.02565 doesnt always yield the exact 1 but 0.999999999--->19 toles
    
    
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
                cout << "The number of rows and columns must be strictly greater than tol, setting to NULL... " << endl;
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
            T tol{}; //creating a tol element
            return tol;
        } else if( ( r >= 0 && r <= rows-1) && ( c >= 0 && c <= cols-1)){
            return vector[r*cols + c];
        } else {
            cout << "Choose a valid index for the elements of this matrix" <<endl;
            T tol{};
            return tol;
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
    
    Matrix getRow(size_t r_indx){      /// row getters (dont forget to transpose)
        if(vector == nullptr){
            cout << "cannot get rows of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (!((r_indx >= 0) && (r_indx <= rows)) ) {
            cout << "choose a valid row index " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix rowVector ; rowVector.rows = 1; rowVector.cols = cols;    //you get a horizontal vector 1*n , you gotta transpose it
            rowVector.vector = new T[rowVector.rows*rowVector.cols];
            for(size_t j = 0; j < rowVector.cols; j++){
                rowVector.vector[0 + j] = vector[r_indx*cols + j];
            }
            return rowVector;
        }
    }
    
    Matrix getCol(size_t c_indx){   //column getters
        if(vector == nullptr){
            cout << "cannot get columns of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (!((c_indx >= 0) && (c_indx <= cols)) ) {
            cout << "choose a valid column index " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix colVector ; colVector.rows = rows; colVector.cols = 1;    //you get a vectrical vector m*1 , 
            colVector.vector = new T[colVector.rows*colVector.cols];
            for(size_t i = 0; i < colVector.rows; i++){
                colVector.vector[i*colVector.cols + 0] = vector[i*cols + c_indx];
            }
            return colVector;
        }
    }
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
    
    /////////////////////// aplications of linear algbra /////////////////
    //1- solving systems of linear equations Ax = b
    // the solve fucntion accepts A, the (m*n) coefficient matrix and b is the result vector (n*1)  and creates the augmented matrix AUG = [A | b] that is m*(n+1)
    // it computed rank(A | b)
    // if rank(A|b) != rank(A) => there exists no solutions
    //if rank(A|b) = rank(A) and nullity = 0 => there exists a unique solution
    //if  rank(A|b) = rank(A) and nullity != 0 => there exists an infinity of solutions  
    
    //2- the solve function returns a structure of type Solution_Struct
    //a 1*n matrix for the solution vector is there exists a single solution (A is n*n and invertible, use x = A^-1b)
    //a null matrix if there exists no solutions
    // if b = 0, return a n*k matrix with the basis vectors of a nullspace
    // if b != 0  a particular solution (basic = 1, free = 0) that is a n*1 matrix, and a n*k matrix with  vectors of a certain basis of the nullspace
    
    // create a struct that contains the solution
    struct Solution_Struct {
        bool hasSolutions;
        bool isUnique;
        Matrix X; //the unique solution
        Matrix Xp; //a n*1 matrix,  particluar soluion
        Matrix N_A; //a n*k matrix , certain nullspace basis
        
        // default constructor 
        Solution_Struct() : hasSolutions(0) ,isUnique(0) ,X(Matrix()) ,Xp(Matrix()) ,N_A(Matrix()) {}
        // or just Solution_Struct Solution{}; (aggregate initialization)
    };
    
    
    
    friend Solution_Struct Solve(const Matrix& A,const Matrix& b){    //Ax = b
        //declare the solution variable
        Solution_Struct Solution;
        
        if((A.vector == nullptr)||(b.vector == nullptr)){
            cout << "can solve a null system of equations " << endl;
            return Solution;
        } else if((b.cols != 1)||(A.rows != b.rows)) {   //the number of b coreffs must be the same as the number of equations, and b must be a vector
            cout << "The dimentions of the coefficeint matrix or vector are invalid  " << endl;
            return Solution;
        } else {  // now we can try to solve the system
            
            //construct the augmented matrix
            Matrix Aug; Aug.rows = A.rows; Aug.cols = A.cols + 1;   //Aug = [A | b]
            Aug.vector = new T[Aug.rows*Aug.cols];
            
            for(size_t i = 0; i < Aug.rows ; i++ ){
                for(size_t j = 0; j < Aug.cols ; j++){
                    
                    if(j < A.cols){   //means we are in the matrix A
                        Aug.vector[i*Aug.cols+j] =  A.vector[i*A.cols+j] ;
                    } else if(j > A.cols - 1){ // means we are in the vector b
                        Aug.vector[i*Aug.cols+j] =  b.vector[i*b.cols +j - A.cols] ;
                    }
                }
            }
            //declare solvabilty cases
            enum Solvability {NoSOLUTIONS = 0, UNIQUE = 1 , INFINITE = 2} ;
            Solvability solutionSet ;  //we dont have to write struct brfore the type in c++
            // determine solvability
            if(RANK(A) != RANK(Aug)) {
                solutionSet = NoSOLUTIONS;
                Solution.hasSolutions = 0;
            } else {
                Solution.hasSolutions = 1;
                if(NULLITY(A) == 0) {
                    solutionSet = UNIQUE;
                    Solution.isUnique = 1;
                    
                } else {
                    solutionSet =  INFINITE;
                    Solution.isUnique = 0;
                }
            }
            switch(solutionSet){
                
                case NoSOLUTIONS :{
                    cout << "the system of equations has no solutions" << endl;
                    return Solution;
                    // no need to break , when you return , the code exists the solve function
                }
                
                case UNIQUE : {
                    cout << "the system of equations has a unique solution" << endl;
                    Solution.X = INV(A)*b;//we are guarenteed that b will be n*1 becuase if the system has a unique solution then A is n*n and invertible
                    return Solution;
                } 
                
                case INFINITE :{
                    cout << "the system has infinitly many solutions " << endl;
                    
                    
                    Matrix rrefAug = RREF(Aug);//extract particular solution
                    Matrix R = RREF(A); //form a basis for the nullspace of A
                    Matrix RE(REF(A)); //search for pivots and free variables ( for pivot search algorithm) // but you could replace this with RREF
                    //1-FIND PARTICULAR SOLUTION
                    
                    // find the pivots , each time you do so , set Xp(c_index) <--- b(r_index)  where c_index and r_index the indices of the found pivot entries
                    Solution.Xp.rows = A.cols; Solution.Xp.cols = 1; Solution.Xp.vector = new T[Solution.Xp.rows*1](); //n*1 particula solution vector filled with tols
                    ///////////////////////////////
                    
                    T tol = 1e-12; 
                    
                    size_t m = RE.rows ; size_t n = RE.cols;
                    size_t f = 0 , c_index = 0 ,  r_index  =0; // column and row indecies of the found pivot
                    
                    for(int cR = m - 1 ; cR > -1 ; cR--  ){ 
                        c_index = 0; f = 0; 
                        for(size_t j = 0; j < n ; j++){
                            if((abs( 1 - RE.vector[cR*n+j] ) < tol)&&(f < 1)){  
                                f = 1;
                                c_index = j;  // here we saved the row and column indecies of the pivot
                                r_index = cR;
                            }  
                        }     
                        if(f == 1){ // found a pivot 
                            // set Xp(c_indexOfpivot) <---- b(r_indexOfpivot)
                            Solution.Xp.vector[(c_index)*Solution.Xp.cols + (0)] = rrefAug.vector[(r_index)*rrefAug.cols + (rrefAug.cols - 1)] ;//to get reduced b fix i to last column index of augmatrix (its number of cols - 1)
                            
                        } 
                    }
                    
                    //2-FORM A BASIS FOR THE NULLSPACE OF A
                    
                    //construct a n*k matrix   k = nullity , each column is a basis vector
                    size_t k  = NULLITY(A);
                    Solution.N_A.rows = A.cols ;    Solution.N_A.cols = k;
                    Solution.N_A.vector =  new T[Solution.N_A.rows*Solution.N_A.cols]() ; // fill it with zeros
                    
                    
                    ///----------------------> /*algo to find free variable indices (the indices of linearly dependent cols in the RREF(A))*/
                    
                    bool* isFree = new bool[A.cols]; //create a n elemtent array to store the state of each var index (false for basic , true for free)
                    //initialize it to true
                    for(int s = 0; s < A.cols ; s++) isFree[s] = true;
                    size_t* freeVar_index = new size_t[k](); //create an k element array to store the indices of the free vars (nullity is the number of basic vars)
                    
                    //run the pivot search algo on the rref(A) and set free indices to false in isFree[]
                    c_index = 0 ,  r_index  =0; //reset the vars form the last part just for safety
                    for(int cR = m - 1 ; cR > -1 ; cR--  ){  
                        c_index = 0; f = 0; 
                        for(size_t j = 0; j < n ; j++){
                            if((abs( 1 - RE.vector[cR*n+j] ) < tol)&&(f < 1)){   //the first entry = 1 in that current row is a pivot
                                f = 1;
                                c_index = j; 
                                r_index = cR;
                            } 
                        }     
                        if(f == 1){ // found a pivot
                            isFree[c_index] = false; //if found a pivot, then its colmun index is not that of a free variable
                        } 
                    }
                    
                    //fill the index array with the indices of the free vars
                    size_t z  = 0; //variable to help fill the correct indices of freeVar_index ( cuz n > k, we could access an invalid position in the array)
                    for(size_t s = 0; s < A.cols ; s++){
                        if(isFree[s]){freeVar_index[z] = s;  z++;} 
                        //to set the next element to the next free index
                    }
                    
                    //----------------------------------->
                    
                    // i is the index of the currect basic vector we are constructing in N_A
                    for(size_t i = 0 ; i <  Solution.N_A.cols ; i++ ){  //each iteration makes a basic vector
                        
                        Solution.N_A.vector[(freeVar_index[i])*Solution.N_A.cols + (i)] = 1;   // set a free var to 1 , other free vars are already initialized to 0
                        
                        /// set one basic variables of the current basic vector each time a pivot is found
                        c_index = 0 ,  r_index  =0; //reset the vars form the last part just for safety
                        for(int cR = m - 1 ; cR > -1 ; cR--  ){  
                            c_index = 0; f = 0; 
                            for(size_t j = 0; j < n ; j++){
                                if((abs( 1 - RE.vector[cR*n+j] ) < tol)&&(f < 1)){   //the first entry = 1 in that current row is a pivot
                                    f = 1;
                                    c_index = j;  // here we saved the row and column indecies of the pivot
                                    r_index = cR;
                                } 
                            }     
                            if(f == 1){ // found a pivot
                                // set N_A[c_index][index of the basis vector we are constructing now] <----- -R[r_index][freeVar_index (the col in R we are extracting values form)]
                                Solution.N_A.vector[(c_index)*Solution.N_A.cols + (i)] = -R.vector[(r_index)*R.cols + (freeVar_index[i])];
                            } 
                        }
                        // in the next iteration we fill the next basis vector in N_A
                    }
                    delete[] freeVar_index;
                    delete[] isFree;
                    
                    return Solution;
                    
                } 
                default : {
                    return Solution;
                }  //just to make the compiler happy (switch statements are not exhaustive)
            }
        } 
        
    }
};
//just an example to use the custom filler funtion
void filler(Matrix<double>* A){
    
    size_t m = A->rowNbr(); size_t n = A->colNbr(); 
    for(size_t i = 0; i < m;i++){
        for(size_t j = 0; j < n ; j++){
            if(i == j){
                A->setElement(i,j,1);
            } else {
                A->setElement(i,j,0);
            }
        }  
    }
    
}


int main(){  
    
    return 0; 
}   


// and then systems of linear equations,eigen vectors and eigen values , polynomial class


/* 
//MATRIX-VECTOR MULTIPLICATION DEMO
Matrix<double> A(4,3);
Matrix<double> x(3,1);
Matrix<double> b;
cout << "enter the matrix A" << endl;
A.fillMatrixUI();
cout << "enter the vector b" << endl;
x.fillMatrixUI();
cout << "The matrix A is :" <<endl ;
A.printMatrix();
cout << " the vector x is :" << endl;
x.printMatrix();
b = A*x;
cout << "the result Ax = b is :" << endl;
b.printMatrix();
*/


/*
//OPERATIONS ON MATRICES DEMO

Matrix<long double> A;
long double m[15][15] = {
{  3, -7,  12,  0, -5,  9,  8,  2, -11, 14,  1, -6,  4, 10, -3},
{ -2,  5, -9,  13,  7, -8,  0, 11,  6, -4, 15,  3, -10, 12,  1},
{  4,  6, -2,  8, -12,  3, -9,  7,  5, -1, 10, -14, 11, -7, 13},
{  9,  0, -6,  2,  1, -11,  13,  5, -8,  7, -3,  4, 14, -2, 12},
{ -8,  12,  0, -7,  6,  2, -5, 10, -4,  1,  9, -13, 11, -6,  3},
{  7, -3,  11,  5, -9,  14,  0, -2,  8, -12,  6,  1, -10,  13, -4},
{  1,  9, -5,  4,  0, -7,  12,  3, -2, 11, -8, 10, -6, 14, -13},
{ -4,  13,  7, -11,  3,  6, -1,  8,  2, -9,  5, -10,  15,  0, -12},
{  6, -10,  2,  14, -3,  0,  9, -5,  11, -7,  13, -1,  8, -4, 12},
{  5,  1, -8,  12,  4, -6,  0,  14, -9,  3, -11,  7, -13,  2, 10},
{ -11,  7,  0, -4,  10,  2, -12,  6,  1,  9, -3, 13, -7,  5, -8},
{  8, -2,  14, -9,  6, -13,  3,  0,  7, -5, 11, -1, 12, -10,  4},
{  2,  10, -7,  0,  13, -5,  9, -3,  8, 14, -6,  4, -12,  11, -1},
{ -6,  4,  1, -8,  11,  7, -2, 13,  0, -10,  12, -9,  5, -14,  3},
{  0,  15, -13,  6, -2,  8, -11,  4,  9, -7,  3, -12,  10,  1, -5}
};

A.fillMatrixA(m);
cout << "the matrix A is :"<< endl;
A.printMatrix(9,4);

cout << "the determinant of A is : " << DET(A) << endl;
cout <<"the rank of A is : " << RANK(A) << endl;
cout <<"the nullity of A is : " << NULLITY(A) << endl;

cout << "inverse of A is :" << endl;
INV(A).printMatrix(9,4);

cout << "transpose of A is :" << endl;
TR(A).printMatrix(9,4);

cout << "row echelon form of A is :" << endl;
REF(A).printMatrix(9,4);

cout << "reduced row echelon form of A is :" << endl;
RREF(A).printMatrix(9,4);

cout << "the cofactor matrix of A is :"<< endl;
COF(A).printMatrix(10,3);

cout << "the adjoint matrix of A is :" << endl;
ADJ(A).printMatrix(10,3);

// (A*INV(A)).printMatrix(2,1);
// (INV(A)*A).printMatrix(2,1);
*/ 


/* 
//DIAGNOLIZATION DEMO
Matrix<double> A(3,3),P(3,3);
A.fillMatrixUI();
P.fillMatrixUI();
(INV(P)*A*P).printMatrix(8,4);
*/


/*
// SOLVING SYSTMES OF LINEAR EQUATIONS DEMO
long double O[8][8] = {
{1,  2,  0,  1,  0,  3,  2,  0},
{2,  4,  0,  2,  0,  6,  4,  0},  
{0,  0,  1,  2,  1,  0,  1,  3},
{0,  0,  2,  4,  2,  0,  2,  6},   
{1,  1,  1,  0,  2,  1,  3,  1},
{2,  2,  2,  0,  4,  2,  6,  2},   
{0,  0,  0,  0,  1,  1,  0,  1},
{0,  0,  0,  0,  2,  2,  0,  2}   
};
long double p[8][1] = {
{3},
{6},  
{5},
{10},  
{4},
{8},   
{2},
{4}  
};
Matrix<long double> A,b;
A.fillMatrixA(O);
b.fillMatrixA(p);

Matrix<long double>::Solution_Struct SystemSol = Solve(A,b);  //the user defined type "Solution_Struct" must be qualified, means we must tell the compiler using the scope resolution specifier that the definition of this type is found in the matrix<long double> class, as each time you use the template class with a different datatype , the compiler creates a completely new version but with T replaced with long double
if(SystemSol.hasSolutions){
if(SystemSol.isUnique){   //when the solution exists and is unique
(SystemSol.X).printMatrix(20,19); 
} else {    //when there exists an infinity of solutions

cout << "a particular solution is :" << endl;
(SystemSol.Xp).printMatrix(20,19); 
cout << "a basis for the nullspace is :" << endl;
(SystemSol.N_A).printMatrix(8,4); // if we have k vectors that belong to N(A), and they are linearly indep, then if k = nullity then they are a valid basis
long double c1 = 4.0 , c2 = 100 , c3 = 40.7 , c4 = -125.0;
Matrix<long double> Xh = c1*SystemSol.N_A.getCol(0) + c2*SystemSol.N_A.getCol(1) +  c3*SystemSol.N_A.getCol(2)  ;
Matrix<long double> X = SystemSol.Xp + Xh;
(A*X).printMatrix(); //test if the new solution is correct
}
}else {  //when there are no solutions
(SystemSol.Xp).printMatrix(20,19); 
}

*/


/* 
// SOLVING SYSTMES OF LINEAR EQUATIONS (USER INPUT)
Matrix<double> A(3,4),b(3,1);
cout << "enter the coefficient matrix" << endl;
A.fillMatrixUI();
cout << "enter values of b" << endl;
b.fillMatrixUI();

Matrix<double>::Solution_Struct SystemSol = Solve(A,b);  //the user defined type "Solution_Struct" must be qualified, means we must tell the compiler using the scope resolution specifier that the definition of this type is found in the matrix<long double> class, as each time you use the template class with a different datatype , the compiler creates a completely new version but with T replaced with long double
if(SystemSol.hasSolutions){
if(SystemSol.isUnique){   //when the solution exists and is unique
(SystemSol.X).printMatrix(8,4); 
} else {    //when there exists an infinity of solutions

cout << "a particular solution is :" << endl;
(SystemSol.Xp).printMatrix(8,4); 
cout << "a basis for the nullspace is :" << endl;
(SystemSol.N_A).printMatrix(8,4); 
double c1 = 4.0 , c2 = 100 , c3 = 40.7 , c4 = -125.0;
Matrix<double> Xh = c1*SystemSol.N_A.getCol(0) + c2*SystemSol.N_A.getCol(1) ;
Matrix<double> X = SystemSol.Xp + Xh;
cout << "an other solution is " << endl;
X.printMatrix(8,4);
cout << " A*X "<<endl;
(A*SystemSol.X).printMatrix(8,4); //test if the new solution is correct


}
}else {  //when there are no solutions
(SystemSol.Xp).printMatrix(20,19); 
}
*/
