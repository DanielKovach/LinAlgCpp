#include <iostream>
#include <vector>
#include <stdexcept>


// HELPERS


class Tens {
  
  
  public:
    std::vector<int> data;
    std::vector<size_t> shape;


  Tens(std::vector<int> d, std::vector<size_t> s) : data(d), shape(s) {
        
        size_t _cumprod = 1;
        for (int i = 0; i  < s.size(); i++){
            
            _cumprod *= s[i];  
            
            if (s[i] < 1) {
              throw std::invalid_argument("Shape members must be greater than 0.\n");
            }

            if ((d.size() % s[i]) != 0) {
              throw std::invalid_argument("Shape members must divide d.size()");
            }
            
        }

        if (d.size() != _cumprod && d.size() > 0) {
            throw std::invalid_argument("Product of shape members must equal d.size()");
          }
    }

  Tens operator+(const Tens &b) {


  if (data.size() != b.data.size()) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' args a and b to must be the same size."; 
    throw std::runtime_error(what_arg);
  }
  
  std::vector<int> v(b.data.size());
  
  for (int i = 0; i < b.data.size(); i ++) {
    v[i] = data[i] + b.data[i];
  }
  
  return Tens(v, b.shape);
}



  friend std::ostream& operator<<(std::ostream& os, Tens t) 
  {   
      
      std::string s = "[[ ";

      for (int i = 0; i < (t.data.size() - 1); i++) {

        if (t.data[i] > 0) {
          
          s += " ";     
        }
        s += std::to_string(t.data[i]) + " ";
        
        if ((i+1) % t.shape[1] == 0) {
          
          s += "],\n [ ";
        }
      }


      int last =  t.data[t.data.size() - 1];
      
      if (last > 0) {
        s += " ";
        }
      s += std::to_string(last) + " ]]\n\n";
      os << s;

      return os;
  }


};




// LINEAR ALGEBRA



std::vector<int> hProd(const std::vector<int> &a, const std::vector<int> &b) {

    /*

    Haddamard product - Component-wise mult
    
    */

  if (a.size() != b.size()) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' args a and b to must be the same size.";
    
    throw std::runtime_error(what_arg);
  }
  
  std::vector<int> ans(a.size());
  
  for (int i = 0; i < a.size(); i ++) {
    ans[i] = a[i]*b[i];
  }
  
  return ans;
}


int vDot(const std::vector<int> &a, const std::vector<int> &b) {

  /*

  Eventually want to overload this with types other than vectors

  */

    if (a.size() != b.size()) {
      const std::string what_arg = "\'" + std::string(__func__) + "\' args a and b to must be the same size.";
      
      throw std::runtime_error(what_arg);
    }
    
    int ans = 0;
    
    for (int i = 0; i < a.size(); i ++) {
      ans += a[i]*b[i];
    }
    
    return ans;
}




std::vector<int> mMul(const std::vector<int> &a, const std::vector<int> &b, const int &m) {

  /*

  If args a and b are respectively nxm and mxp matrices, then this returns the nxp matrix using the naive algorithm.
  
  */

  if (m < 1) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' int arg m must be the greater than 0.";
    throw std::runtime_error(what_arg);
    }
  if ((a.size() % m) != 0) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' int arg m must divide the length of \'a\'.";
    throw std::runtime_error(what_arg);
    }
  if ((b.size() % m) != 0) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' int arg m must divide the length of \'b\'.";
    throw std::runtime_error(what_arg);
    }

  int n = a.size()/m;
  int p = b.size()/m;
  
  std::vector<int> ans(n*p);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {  
      int tmp = a[i*m + j];        
      for (int k = 0; k < p; k++) {
        ans[i*p + k] += tmp*b[j*p + k];
      } 
    }
  }
    
    return ans;
}

std::vector<int> tpose(const std::vector<int> &a, const int &m) {

  /*

  If arg a is an nxm matrix, then this returns the mxn transpose of a, optionally inplace
  
  */

  if (m < 1) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' int arg m must be the greater than 0.";
    throw std::runtime_error(what_arg);
    }
  if ((a.size() % m) != 0) {
    const std::string what_arg = "\'" + std::string(__func__) + "\' int arg m must divide the length of \'a\'.";
    throw std::runtime_error(what_arg);
    }


  int n = a.size()/m;
  
  
  std::vector<int> ans(a.size());
  
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {  
      ans[j*n + i] = a[i*m + j];        
    }
  }
    
    return ans;
}


int trace(const std::vector<int> &a, const int &numCols) {

  /*

  Calculates trace given a matrix a
  
  */

  int ans = 0;
  int n = a.size()/numCols;
  int min = std::min(n, numCols);
  
  int tmp = numCols + 1;

  for (int i = 0; i < min; i ++) {
    ans += a[tmp*i];
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

  std::vector<int> inp1 = {-1, 2, 3, -4, 10, 10,-4, 10, 10};
  std::vector<int> inp2 = {-1, -2, -3, 4, 5, 6,-4, 10, 10};
  
  // Num Cols in inp1 or num rows in inp2
  size_t m = 3;
  size_t n = inp1.size()/m;
  std::vector<size_t> s = {n,m};  
  Tens mat1 = Tens(inp1, s);
  
  size_t p = inp2.size()/m;
  std::vector<size_t> s2 = {m,p};  
  Tens mat2 = Tens(inp2, s2);


  std::cout << mat1;
  std::cout << mat2;
  std::cout << mat2 + mat1;

  //std::cout << std::to_string(inp1);
  //Tens sm = mat1 + mat2;
  //std::vector<int> ans = mMul(inp1, inp2, m);
  
  //int m2 = inp2.size()/m;

  //std::cout <<  "Transpose: " << "\n";
  //std::vector<int> ansT = tpose(ans, m2);

  //printMat(ansT, m);
  //printMat(ansT, m, true);
  //std::cout << trace(ans, m2) << "\n";


  return 0;
}