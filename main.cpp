#include <iostream>
#include <vector>
#include <stdexcept>

using namespace std;


// HELPERS

class Tens {
  
  
  public:
    vector<int> data;
    vector<size_t> shape;


  Tens(vector<int> d, vector<size_t> s) : data(d), shape(s) {
        
        size_t _cumprod = 1;
        for (int i = 0; i  < s.size(); i++){
            
            _cumprod *= s[i];  
            
            if (s[i] < 1) {
              throw invalid_argument("Shape members must be greater than 0.\n");
            }

            if ((d.size() % s[i]) != 0) {
              throw invalid_argument("Shape members must divide d.size()");
            }
            
          }

        if (d.size() != _cumprod && d.size() > 0) {
            throw invalid_argument("Product of shape members must equal d.size()");
          }
    }

  Tens operator+(const Tens &b) {
    /*

    Component-wise Addition
    
    */


    if (data.size() != b.data.size()) {
      const string what_arg = "\'" + string(__func__) + "\' args a and b to must be the same size."; 
      throw runtime_error(what_arg);
    }
    
    vector<int> v(b.data.size());
    
    for (int i = 0; i < b.data.size(); i ++) {
      v[i] = data[i] + b.data[i];
    }
    
    return Tens(v, b.shape);
  }

  Tens operator*(const Tens &b) {
  /*

  If args a and b are respectively nxm and mxp matrices, then this returns the nxp matrix using the naive algorithm.
  
  */
  
  size_t m = shape[1];
  if (m != b.shape[0]) {
    const string what_arg = "\'" + string(__func__) + "\': Tensors must have compatible shapes for multiplication.";
    throw runtime_error(what_arg);
    }

  size_t n = shape[0];
  size_t p = b.shape[1];
  
  vector<int> ans(n*p);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {  
      int tmp = data[i*m + j];        
      for (int k = 0; k < p; k++) {
        ans[i*p + k] += tmp*b.data[j*p + k];
      } 
    }
  }
  vector<size_t> ans_shp = {n, p};
  return Tens(ans, ans_shp);
}

  friend ostream& operator<<(ostream& os, Tens t) 
  {   
      // TODO: 
      //  -Generalize to more than 2D
      //  -Line up each column
      string s = "[[ ";

      for (int i = 0; i < (t.data.size() - 1); i++) {

        if (t.data[i] > 0) {
          
          s += " ";     
        }
        s += to_string(t.data[i]) + " ";
        
        if ((i+1) % t.shape[1] == 0) {
          
          s += "],\n [ ";
        }
      }


      int last =  t.data[t.data.size() - 1];
      
      if (last > 0) {
        s += " ";
        }
      s += to_string(last) + " ]]\n\n";
      os << s;

      return os;
  }


};




// LINEAR ALGEBRA

Tens hprod(const Tens &a, const Tens &b) {
    /*

    Haddamard product - Component-wise mult
    
    */

    if (a.shape != b.shape) {
      const string what_arg = "\'" + string(__func__) + "\' args a and b to must be the same shape.";
      
      throw runtime_error(what_arg);
    }
    
    vector<int> v(a.data.size());
    
    for (int i = 0; i < a.data.size(); i ++) {
      v[i] = a.data[i]*b.data[i];
    }
    
    return Tens(v, b.shape);
}

int vDot(const vector<int> &a, const vector<int> &b) {

  /*

  Eventually want to overload this with types other than vectors

  */

    if (a.size() != b.size()) {
      const string what_arg = "\'" + string(__func__) + "\' args a and b to must be the same size.";
      
      throw runtime_error(what_arg);
    }
    
    int ans = 0;
    
    for (int i = 0; i < a.size(); i ++) {
      ans += a[i]*b[i];
    }
    
    return ans;
}

Tens mMul(const Tens &a, const Tens &b) {

  /*

  If args a and b are respectively nxm and mxp matrices, then this returns the nxp matrix using the naive algorithm.
  
  */
  
  size_t m = a.shape[1];
  if (m != b.shape[0]) {
    const string what_arg = "\'" + string(__func__) + "\': Tensors must have compatible shapes for multiplication.";
    throw runtime_error(what_arg);
    }

  size_t n = a.shape[0];
  size_t p = b.shape[1];
  
  vector<int> ans(n*p);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {  
      int tmp = a.data[i*m + j];        
      for (int k = 0; k < p; k++) {
        ans[i*p + k] += tmp*b.data[j*p + k];
      } 
    }
  }
  vector<size_t> ans_shp = {n, p};
  return Tens(ans, ans_shp);
}

Tens transpose(const Tens &a) {

  /*

  If arg a is an nxm matrix, then this returns the mxn transpose of a

  TODO: Make this optionally inplace
  
  */

  size_t n = a.shape[0];
  size_t m = a.data.size()/a.shape[0];
  
  vector<int> ans(a.data.size());
  
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {  
      ans[j*n + i] = a.data[i*m + j];        
    }
  }
  vector<size_t> ans_shp = {m, n};
  
  return Tens(ans, ans_shp);
}


int trace(const Tens &a) {

  /*

  Calculates trace given a matrix a
  
  */

  int ans = 0;
  
  int minn = min(a.shape[0], a.shape[1]);
  
  int tmp = a.shape[1] + 1;

  for (int i = 0; i < minn; i ++) {
    ans += a.data[tmp*i];
  }
  
  return ans;
}



/*
TODO: 

--General Development:

Set up a logger


--Linear Algebra Algorithms:

Mult:
Strassen

Inverse Algos:
Cramer's Rule

Decomps:
SVD
Cholsky
QR

Eq Test:
Freivald's


--Analysis Algorithms:

Transforms:
FT
FFT

Auto-grad

-- Probability Algorithms:

Gamma Function


Later:

- IMgui
- Logger

*/



int main() {

  vector<int> inp1 = {-1, 2, 3, -4, 10, 10,-4, 10, 10};
  vector<int> inp2 = {-1, -2, -3, 4, 5, 6,-4, 10, 10};
  
  // Num Cols in inp1 or num rows in inp2
  size_t m = 3;
  size_t n = inp1.size()/m;
  vector<size_t> s = {n,m};  
  Tens mat1 = Tens(inp1, s);
  
  size_t p = inp2.size()/m;
  vector<size_t> s2 = {m,p};  
  Tens mat2 = Tens(inp2, s2);


  cout << mat1;
  cout << mat2;
  cout << mat2 + mat1;

  //cout << to_string(inp1);
  //Tens sm = mat1 + mat2;
  //vector<int> ans = mMul(inp1, inp2, m);
  
  //int m2 = inp2.size()/m;

  //cout <<  "Transpose: " << "\n";
  //vector<int> ansT = tpose(ans, m2);

  //printMat(ansT, m);
  //printMat(ansT, m, true);
  //cout << trace(ans, m2) << "\n";


  return 0;
}