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
    /////////////////////////////////////////////////////////////////////////Helper functions 
    template<size_t M,size_t N>
    void copyData(const T (&src)[M][N]){                 //data copy
        for(size_t i = 0; i < rows;i++){
            for(size_t j = 0; j < cols ; j++){
                vector[i*cols + j] = src[i][j];
            }  
        }
    }
    
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
                COF.vector[i*COF.cols + j] = pow(-1,i+j)*det(minorMatrix(matrix,i,j)); 
            }
        }
        return COF;
    }
    
    friend Matrix adjI(const Matrix& matrix){              //adjoint function just for computing inverse
        // assuming this is an n*n matrix that is not null
        Matrix ADJ = tr(cofI(matrix));
        return ADJ;
        
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
        } 
        return *this;
    }
    
    ~Matrix(){                       //destructor
        delete[] vector;
    }
    
    ///////////////////////////// operations on matrices ////////////
    
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
    
    friend Matrix tr(const Matrix& matrix){                                 //matrix transpose
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
    
    friend T det(const Matrix& matrix){                 // determinant using laplace expansion
        if(matrix.vector == nullptr){
            cout << "Cannot compute determinant of null matrix " << endl;
            T zero{};
            return  zero;
        } else if(matrix.rows != matrix.cols) { 
            cout <<" Cannot compute determinant of non-square matrix " << endl;
            T zero{};
            return  zero;
        } else if(matrix.rows == 2){
            T DET = matrix.vector[0]*matrix.vector[matrix.cols + 1] - matrix.vector[1]*matrix.vector[matrix.cols];
            return DET;
        } else {
            T DET{};
            for(size_t j = 0; j < matrix.cols ; j++){   // expanding along the 0th row i = 0
                DET =  DET + (matrix.vector[j])*pow(-1,j)*det(minorMatrix(matrix,0,j));
            }
            return DET; 
        }
        
    }  
    
    friend Matrix inv(const Matrix& matrix){           /// matrix  inverse
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
    }
    
    friend Matrix cof(const Matrix& matrix){         //// cofactor matrix
        if(matrix.vector == nullptr){ 
            cout << "Cannot compute the cofactor matrix of a null matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot compute the cofactor matrix of a non-square matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix COF = cofI(matrix);
            return COF;
        }
    }
    
    friend Matrix adj(const Matrix& matrix){                ////// adjoint matrix 
        if(matrix.vector == nullptr){ 
            cout << "Cannot compute the adjoint matrix of a null matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else if (matrix.rows != matrix.cols){
            cout << "Cannot compute the adjoint matrix of a non-square matrix" << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            Matrix ADJ = adjI(matrix);
            return ADJ;
        }
    }
    
    
    
    
    /////////////////////////// elementary column operations //////////////////
    
    
    ///1 - elementary column operations, ( you must remove H and make just one version)
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
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////// ROW ECHELON FORM
    
    
    friend Matrix REF(const Matrix& M){        //non normalized pivots,   // row echelon form
        if(M.vector == nullptr) {
            cout << "Cannot find row echelon form of a null matrix " << endl;
            Matrix nullmatrix;
            return nullmatrix;
        } else {
            // we have 3 function matrices that represent all row operations on any row with any scalar 
            // addE(r1,r2,s,M)   scaleE(r,s,M)   swapE(r1,r2,M)  , multiplying by xxE*M, applies the elementary row operation on M
            /* here we always take the pivot to be the first non zero entry in the column just for simplification purposes*/ 
            
            Matrix RE(M);      size_t m = RE.rows ; size_t n = RE.cols;   
            /*  ----- -> */  T zero = 1e-12;
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
        

        friend Matrix()
    }
    /////////////////////////////////////////////// setters and getters //////////////////////////////
    template<size_t M,size_t N>
    void fillMatrixA(const T (&arr)[M][N]){               //setter using array 
        if(vector == nullptr){
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
    long double CrazyMatrix[20][11] = {
        {  1.223, -3.445,  2.118, -0.9912,  4.331, -2.771,  0.6732,  1.889, -2.114,  3.221, -0.5567 },
        { -2.781,  4.112, -1.008,  0.6621, -3.771,  2.119, -0.9876,  1.332, -2.445,  0.8732,  3.221 },
        {  3.998, -1.223,  0.4432, -2.118,  1.887, -0.6611,  4.555, -3.009,  1.672, -0.2299,  2.445 },
        { -0.7782,  1.332, -3.889,  4.111, -2.441,  0.5561, -1.223,  3.887, -0.6644,  2.998, -4.222 },
        {  2.771, -0.9932,  1.221, -4.332,  3.876, -2.118,  0.7781, -1.887,  4.119, -3.442,  1.003 },
        { -3.009,  2.665, -0.4412,  1.765, -2.333,  4.221, -0.7789,  3.001, -1.442,  0.6711, -2.871 },
        {  1.111, -2.876,  3.998, -1.442,  2.889, -0.5567,  4.332, -2.761,  0.8733, -3.221,  1.229 },
        { -4.112,  3.777, -2.991,  0.6611, -1.887,  4.009, -3.221,  1.776, -0.5543,  2.998, -0.9987 },
        {  0.6621, -1.223,  2.887, -3.111,  4.665, -2.888,  0.4422, -1.998,  3.333, -0.7781,  1.119 },
        { -2.665,  0.9981, -1.887,  4.441, -3.118,  0.5569, -2.771,  1.223, -4.112,  3.876, -0.2298 },
        {  3.765, -2.441,  0.8897, -1.223,  2.661, -3.777,  4.119, -0.6623,  1.445, -2.111,  0.7781 },
        { -0.8899,  1.119, -2.761,  3.221, -0.4456,  4.321, -3.222,  0.9981, -1.777,  2.887, -0.6621 },
        {  2.445, -3.221,  1.003, -0.7783,  4.119, -2.118,  0.6621, -1.998,  3.887, -0.2291,  1.772 },
        { -1.776,  0.5543, -3.887,  2.221, -0.6622,  1.442, -4.443,  3.111, -0.7789,  2.998, -1.009 },
        {  4.221, -2.111,  0.5567, -1.221,  3.887, -0.9932,  2.331, -4.001,  1.772, -0.6642,  3.118 },
        { -3.009,  1.998, -0.7781,  2.665, -4.441,  3.223, -1.118,  0.6623, -2.771,  1.889, -0.2294 },
        {  0.6622, -1.332,  2.887, -3.998,  4.119, -2.765,  1.223, -0.8891,  3.221, -1.776,  0.5543 },
        { -2.118,  3.765, -1.223,  0.7781, -3.887,  2.998, -0.6624,  1.119, -4.111,  3.887, -0.2291 },
        {  1.333, -0.7782,  4.119, -2.221,  0.6621, -1.118,  3.887, -0.4453,  2.665, -3.009,  1.776 },
        { -4.119,  2.887, -0.6611,  3.223, -1.009,  0.7784, -2.998,  4.221, -1.333,  0.6621, -3.887 }
    };
    
    long double StressTestMatrix[15][15] = {
        { 1.0, 0.5, 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667 },
        { 0.5, 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625 },
        { 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705 },
        { 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555 },
        { 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842 },
        { 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05 },
        { 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616 },
        { 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456 },
        { 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216 },
        { 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664 },
        { 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04 },
        { 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04, 0.038461538461538464 },
        { 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04, 0.038461538461538464, 0.037037037037037035 },
        { 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04, 0.038461538461538464, 0.037037037037037035, 0.03571428571428571 },
        { 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04, 0.038461538461538464, 0.037037037037037035, 0.03571428571428571, 0.034482758620689655 }
    };
    
    
    // Matrix<long double> A(4,3) ;A.fillMatrixUI();
    
    Matrix<long double> A;//A.fillMatrixA(CrazyMatrix); 
    A.fillMatrixA(StressTestMatrix);
    
    A.printMatrix(9,4);
    cout << endl;
    // (REF(A)).printMatrix(10,4);
    //cout << det(REF(A)) << endl;
    REF(A).printMatrix(3,1);
    cout << endl;
    //  cout << det(A);
    
    
    
    
    
    
    
    return 0;   // row echelon form and rref(A) , rank, nullity , and then systems of linear equations
}   








//cout << "Enter the matrix A" << endl;
//A.fillMatrixUI();


/*Matrix<double> A(2,3);
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
