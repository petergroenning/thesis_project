#include <TMB.hpp>
using namespace density;

//////////// helper fun: find NA locations in vector ///////////
template <class Type>
vector<Type> is_not_na(vector<Type> x){
	vector<Type> y(x.size());
	y.fill(Type(1.0));
    for(int i=0; i<x.size(); i++){
		if( R_IsNA(asDouble(x(i))) ){
			y(i) = Type(0.0);
		}
	}
	return y;
}

//////////// helper fun: extract non-NAs from vector ///////////
template<class Type>
vector<Type> remove_nas__(vector<Type> data_vector, int number_of_datas, vector<Type> na_bool){
  int ii = 0;
	vector<Type> y_reduced(number_of_datas);
	for(int i=0; i < data_vector.size(); i++){
		if(na_bool(i) == Type(1.0)){
			y_reduced(ii) = data_vector(i);
			ii++;
		}
	}
  return y_reduced;
}

//////////// helper fun: construct permutation matrix ///////////
template <class Type>
matrix<Type> construct_permutation_matrix(int number_of_datas, int m, vector<Type> na_bool){
	matrix<Type> E(number_of_datas,m);
	E.setZero();
	/**/
	int j=0;
	for(int i=0; i < m; i++){
		if(na_bool(i) == Type(1.0)){ /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
			E(j,i) = Type(1.0);
			j += 1;
		}
	}
	return E;
}

//////////// loss function ///////////
template<class Type>
Type lossfunction__(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
  Type loss;
  if(lossFunc==1){
    Type a = tukeypars(0);
    Type b = tukeypars(1);
    Type c = tukeypars(2);
    Type d = tukeypars(3);
    loss = d * ( (Type(1.0)/(Type(1.0)+exp(-a*(x-b)))) + c );
  } else if (lossFunc==2){
    Type c_squared = pow(huber_c,2);
    loss = c_squared * (sqrt(1 + (x / c_squared)) - 1);
  } else {
    loss = x;
  }
  return(loss);
}

//////////// MAP estimation helper ///////////
template<class Type>
vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	vector<Type> ans(sum_mapints);
	int j=0;
	for(int i=0;i<mapints.size();i++){
		if(mapints(i)==1){
			ans(j) = parvec(i);
			j += 1;
		}
	}
	return(ans);
}

//////////// UKF sigma points ///////////
template<class Type>
matrix<Type> construct_Xsp__(vector<Type> x0, matrix<Type> s0){
  int n = 16;
  int nn = 33;
  matrix<Type> Xsp(n,nn);
  vector<Type> Si;
  Xsp.col(0) = x0;
  for(int i=1; i<n+1; i++){
    Si = s0.col(i-1);
    Xsp.col(i) = x0 + sqrt(1+n) * Si;
    Xsp.col(i+n) = x0 - sqrt(1+n) * Si;
  }
  return Xsp;
}

//////////// UKF function ///////////
template<class Type>
matrix<Type> Phi__(matrix<Type> M){
  matrix<Type> K(M.col(0).size(),M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
}

//////////// drift function ///////////
template<class Type>
vector<Type> f__(Type tau_x0, Type tau_x1, Type tau_x10, Type tau_x11, Type tau_x12, Type tau_x13, Type tau_x14, Type tau_x2, Type tau_x3, Type tau_x4, Type tau_x5, Type tau_x6, Type tau_x7, Type tau_x8, Type tau_x9, Type x0, Type x1, Type x10, Type x11, Type x12, Type x13, Type x14, Type x15, Type x2, Type x3, Type x4, Type x5, Type x6, Type x7, Type x8, Type x9){
	vector<Type> ans(16);
	ans(0) = (x1 - x0)/tau_x0;
	ans(1) = (x2 - x1)/tau_x1 - (x1 - x0)/tau_x0;
	ans(2) = (x3 - x2)/tau_x2 - (x2 - x1)/tau_x1;
	ans(3) = (x4 - x3)/tau_x3 - (x3 - x2)/tau_x2;
	ans(4) = (x5 - x4)/tau_x4 - (x4 - x3)/tau_x3;
	ans(5) = (x6 - x5)/tau_x5 - (x5 - x4)/tau_x4;
	ans(6) = (x7 - x6)/tau_x6 - (x6 - x5)/tau_x5;
	ans(7) = (x8 - x7)/tau_x7 - (x7 - x6)/tau_x6;
	ans(8) = (x9 - x8)/tau_x8 - (x8 - x7)/tau_x7;
	ans(9) = (x10 - x9)/tau_x9 - (x9 - x8)/tau_x8;
	ans(10) = (x11 - x10)/tau_x10 - (x10 - x9)/tau_x9;
	ans(11) = (x12 - x11)/tau_x11 - (x11 - x10)/tau_x10;
	ans(12) = (x13 - x12)/tau_x12 - (x12 - x11)/tau_x11;
	ans(13) = (x14 - x13)/tau_x13 - (x13 - x12)/tau_x12;
	ans(14) = (x15 - x14)/tau_x14 - (x14 - x13)/tau_x13;
	ans(15) = -((x15 - x14)/tau_x14);
	return ans;
}

//////////// UKF sigma points drift function ///////////
template<class Type>
matrix<Type> construct_F__(matrix<Type> Xsp, Type tau_x0, Type tau_x1, Type tau_x10, Type tau_x11, Type tau_x12, Type tau_x13, Type tau_x14, Type tau_x2, Type tau_x3, Type tau_x4, Type tau_x5, Type tau_x6, Type tau_x7, Type tau_x8, Type tau_x9){
  int n = 16;
  int nn = 33;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0(0), x0(1), x0(10), x0(11), x0(12), x0(13), x0(14), x0(15), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8), x0(9));
  }
  return F;
}

//////////// jacobian of drift function ///////////
template<class Type>
matrix<Type> dfdx__(Type tau_x0, Type tau_x1, Type tau_x10, Type tau_x11, Type tau_x12, Type tau_x13, Type tau_x14, Type tau_x2, Type tau_x3, Type tau_x4, Type tau_x5, Type tau_x6, Type tau_x7, Type tau_x8, Type tau_x9){
	matrix<Type> ans(16,16);
	ans(0,0) = -(1/tau_x0);
	ans(0,1) = 1/tau_x0;
	ans(0,2) = 0;
	ans(0,3) = 0;
	ans(0,4) = 0;
	ans(0,5) = 0;
	ans(0,6) = 0;
	ans(0,7) = 0;
	ans(0,8) = 0;
	ans(0,9) = 0;
	ans(0,10) = 0;
	ans(0,11) = 0;
	ans(0,12) = 0;
	ans(0,13) = 0;
	ans(0,14) = 0;
	ans(0,15) = 0;
	ans(1,0) = 1/tau_x0;
	ans(1,1) = -(1/tau_x0 + 1/tau_x1);
	ans(1,2) = 1/tau_x1;
	ans(1,3) = 0;
	ans(1,4) = 0;
	ans(1,5) = 0;
	ans(1,6) = 0;
	ans(1,7) = 0;
	ans(1,8) = 0;
	ans(1,9) = 0;
	ans(1,10) = 0;
	ans(1,11) = 0;
	ans(1,12) = 0;
	ans(1,13) = 0;
	ans(1,14) = 0;
	ans(1,15) = 0;
	ans(2,0) = 0;
	ans(2,1) = 1/tau_x1;
	ans(2,2) = -(1/tau_x1 + 1/tau_x2);
	ans(2,3) = 1/tau_x2;
	ans(2,4) = 0;
	ans(2,5) = 0;
	ans(2,6) = 0;
	ans(2,7) = 0;
	ans(2,8) = 0;
	ans(2,9) = 0;
	ans(2,10) = 0;
	ans(2,11) = 0;
	ans(2,12) = 0;
	ans(2,13) = 0;
	ans(2,14) = 0;
	ans(2,15) = 0;
	ans(3,0) = 0;
	ans(3,1) = 0;
	ans(3,2) = 1/tau_x2;
	ans(3,3) = -(1/tau_x2 + 1/tau_x3);
	ans(3,4) = 1/tau_x3;
	ans(3,5) = 0;
	ans(3,6) = 0;
	ans(3,7) = 0;
	ans(3,8) = 0;
	ans(3,9) = 0;
	ans(3,10) = 0;
	ans(3,11) = 0;
	ans(3,12) = 0;
	ans(3,13) = 0;
	ans(3,14) = 0;
	ans(3,15) = 0;
	ans(4,0) = 0;
	ans(4,1) = 0;
	ans(4,2) = 0;
	ans(4,3) = 1/tau_x3;
	ans(4,4) = -(1/tau_x3 + 1/tau_x4);
	ans(4,5) = 1/tau_x4;
	ans(4,6) = 0;
	ans(4,7) = 0;
	ans(4,8) = 0;
	ans(4,9) = 0;
	ans(4,10) = 0;
	ans(4,11) = 0;
	ans(4,12) = 0;
	ans(4,13) = 0;
	ans(4,14) = 0;
	ans(4,15) = 0;
	ans(5,0) = 0;
	ans(5,1) = 0;
	ans(5,2) = 0;
	ans(5,3) = 0;
	ans(5,4) = 1/tau_x4;
	ans(5,5) = -(1/tau_x4 + 1/tau_x5);
	ans(5,6) = 1/tau_x5;
	ans(5,7) = 0;
	ans(5,8) = 0;
	ans(5,9) = 0;
	ans(5,10) = 0;
	ans(5,11) = 0;
	ans(5,12) = 0;
	ans(5,13) = 0;
	ans(5,14) = 0;
	ans(5,15) = 0;
	ans(6,0) = 0;
	ans(6,1) = 0;
	ans(6,2) = 0;
	ans(6,3) = 0;
	ans(6,4) = 0;
	ans(6,5) = 1/tau_x5;
	ans(6,6) = -(1/tau_x5 + 1/tau_x6);
	ans(6,7) = 1/tau_x6;
	ans(6,8) = 0;
	ans(6,9) = 0;
	ans(6,10) = 0;
	ans(6,11) = 0;
	ans(6,12) = 0;
	ans(6,13) = 0;
	ans(6,14) = 0;
	ans(6,15) = 0;
	ans(7,0) = 0;
	ans(7,1) = 0;
	ans(7,2) = 0;
	ans(7,3) = 0;
	ans(7,4) = 0;
	ans(7,5) = 0;
	ans(7,6) = 1/tau_x6;
	ans(7,7) = -(1/tau_x6 + 1/tau_x7);
	ans(7,8) = 1/tau_x7;
	ans(7,9) = 0;
	ans(7,10) = 0;
	ans(7,11) = 0;
	ans(7,12) = 0;
	ans(7,13) = 0;
	ans(7,14) = 0;
	ans(7,15) = 0;
	ans(8,0) = 0;
	ans(8,1) = 0;
	ans(8,2) = 0;
	ans(8,3) = 0;
	ans(8,4) = 0;
	ans(8,5) = 0;
	ans(8,6) = 0;
	ans(8,7) = 1/tau_x7;
	ans(8,8) = -(1/tau_x7 + 1/tau_x8);
	ans(8,9) = 1/tau_x8;
	ans(8,10) = 0;
	ans(8,11) = 0;
	ans(8,12) = 0;
	ans(8,13) = 0;
	ans(8,14) = 0;
	ans(8,15) = 0;
	ans(9,0) = 0;
	ans(9,1) = 0;
	ans(9,2) = 0;
	ans(9,3) = 0;
	ans(9,4) = 0;
	ans(9,5) = 0;
	ans(9,6) = 0;
	ans(9,7) = 0;
	ans(9,8) = 1/tau_x8;
	ans(9,9) = -(1/tau_x8 + 1/tau_x9);
	ans(9,10) = 1/tau_x9;
	ans(9,11) = 0;
	ans(9,12) = 0;
	ans(9,13) = 0;
	ans(9,14) = 0;
	ans(9,15) = 0;
	ans(10,0) = 0;
	ans(10,1) = 0;
	ans(10,2) = 0;
	ans(10,3) = 0;
	ans(10,4) = 0;
	ans(10,5) = 0;
	ans(10,6) = 0;
	ans(10,7) = 0;
	ans(10,8) = 0;
	ans(10,9) = 1/tau_x9;
	ans(10,10) = -(1/tau_x10 + 1/tau_x9);
	ans(10,11) = 1/tau_x10;
	ans(10,12) = 0;
	ans(10,13) = 0;
	ans(10,14) = 0;
	ans(10,15) = 0;
	ans(11,0) = 0;
	ans(11,1) = 0;
	ans(11,2) = 0;
	ans(11,3) = 0;
	ans(11,4) = 0;
	ans(11,5) = 0;
	ans(11,6) = 0;
	ans(11,7) = 0;
	ans(11,8) = 0;
	ans(11,9) = 0;
	ans(11,10) = 1/tau_x10;
	ans(11,11) = -(1/tau_x10 + 1/tau_x11);
	ans(11,12) = 1/tau_x11;
	ans(11,13) = 0;
	ans(11,14) = 0;
	ans(11,15) = 0;
	ans(12,0) = 0;
	ans(12,1) = 0;
	ans(12,2) = 0;
	ans(12,3) = 0;
	ans(12,4) = 0;
	ans(12,5) = 0;
	ans(12,6) = 0;
	ans(12,7) = 0;
	ans(12,8) = 0;
	ans(12,9) = 0;
	ans(12,10) = 0;
	ans(12,11) = 1/tau_x11;
	ans(12,12) = -(1/tau_x11 + 1/tau_x12);
	ans(12,13) = 1/tau_x12;
	ans(12,14) = 0;
	ans(12,15) = 0;
	ans(13,0) = 0;
	ans(13,1) = 0;
	ans(13,2) = 0;
	ans(13,3) = 0;
	ans(13,4) = 0;
	ans(13,5) = 0;
	ans(13,6) = 0;
	ans(13,7) = 0;
	ans(13,8) = 0;
	ans(13,9) = 0;
	ans(13,10) = 0;
	ans(13,11) = 0;
	ans(13,12) = 1/tau_x12;
	ans(13,13) = -(1/tau_x12 + 1/tau_x13);
	ans(13,14) = 1/tau_x13;
	ans(13,15) = 0;
	ans(14,0) = 0;
	ans(14,1) = 0;
	ans(14,2) = 0;
	ans(14,3) = 0;
	ans(14,4) = 0;
	ans(14,5) = 0;
	ans(14,6) = 0;
	ans(14,7) = 0;
	ans(14,8) = 0;
	ans(14,9) = 0;
	ans(14,10) = 0;
	ans(14,11) = 0;
	ans(14,12) = 0;
	ans(14,13) = 1/tau_x13;
	ans(14,14) = -(1/tau_x13 + 1/tau_x14);
	ans(14,15) = 1/tau_x14;
	ans(15,0) = 0;
	ans(15,1) = 0;
	ans(15,2) = 0;
	ans(15,3) = 0;
	ans(15,4) = 0;
	ans(15,5) = 0;
	ans(15,6) = 0;
	ans(15,7) = 0;
	ans(15,8) = 0;
	ans(15,9) = 0;
	ans(15,10) = 0;
	ans(15,11) = 0;
	ans(15,12) = 0;
	ans(15,13) = 0;
	ans(15,14) = 1/tau_x14;
	ans(15,15) = -(1/tau_x14);
	return ans;
}

//////////// diffusion function ///////////
template<class Type>
matrix<Type> g__(Type sigma_x0, Type sigma_x1, Type sigma_x10, Type sigma_x11, Type sigma_x12, Type sigma_x13, Type sigma_x14, Type sigma_x15, Type sigma_x2, Type sigma_x3, Type sigma_x4, Type sigma_x5, Type sigma_x6, Type sigma_x7, Type sigma_x8, Type sigma_x9){
	matrix<Type> ans(16,1);
	ans(0,0) = sigma_x0;
	ans(1,0) = sigma_x1;
	ans(2,0) = sigma_x2;
	ans(3,0) = sigma_x3;
	ans(4,0) = sigma_x4;
	ans(5,0) = sigma_x5;
	ans(6,0) = sigma_x6;
	ans(7,0) = sigma_x7;
	ans(8,0) = sigma_x8;
	ans(9,0) = sigma_x9;
	ans(10,0) = sigma_x10;
	ans(11,0) = sigma_x11;
	ans(12,0) = sigma_x12;
	ans(13,0) = sigma_x13;
	ans(14,0) = sigma_x14;
	ans(15,0) = sigma_x15;
	return ans;
}

//////////// observation function ///////////
template<class Type>
vector<Type> h__(Type x0, Type x1, Type x10, Type x11, Type x12, Type x13, Type x14, Type x15, Type x2, Type x3, Type x4, Type x5, Type x6, Type x7, Type x8, Type x9){
	vector<Type> ans(16);
	ans(0) = x0;
	ans(1) = x1;
	ans(2) = x2;
	ans(3) = x3;
	ans(4) = x4;
	ans(5) = x5;
	ans(6) = x6;
	ans(7) = x7;
	ans(8) = x8;
	ans(9) = x9;
	ans(10) = x10;
	ans(11) = x11;
	ans(12) = x12;
	ans(13) = x13;
	ans(14) = x14;
	ans(15) = x15;
	return ans;
}

//////////// jacobian of observation function ///////////
template<class Type>
matrix<Type> dhdx__(Type state_void_filler){
	matrix<Type> ans(16,16);
	ans(0,0) = 1;
	ans(0,1) = 0;
	ans(0,2) = 0;
	ans(0,3) = 0;
	ans(0,4) = 0;
	ans(0,5) = 0;
	ans(0,6) = 0;
	ans(0,7) = 0;
	ans(0,8) = 0;
	ans(0,9) = 0;
	ans(0,10) = 0;
	ans(0,11) = 0;
	ans(0,12) = 0;
	ans(0,13) = 0;
	ans(0,14) = 0;
	ans(0,15) = 0;
	ans(1,0) = 0;
	ans(1,1) = 1;
	ans(1,2) = 0;
	ans(1,3) = 0;
	ans(1,4) = 0;
	ans(1,5) = 0;
	ans(1,6) = 0;
	ans(1,7) = 0;
	ans(1,8) = 0;
	ans(1,9) = 0;
	ans(1,10) = 0;
	ans(1,11) = 0;
	ans(1,12) = 0;
	ans(1,13) = 0;
	ans(1,14) = 0;
	ans(1,15) = 0;
	ans(2,0) = 0;
	ans(2,1) = 0;
	ans(2,2) = 1;
	ans(2,3) = 0;
	ans(2,4) = 0;
	ans(2,5) = 0;
	ans(2,6) = 0;
	ans(2,7) = 0;
	ans(2,8) = 0;
	ans(2,9) = 0;
	ans(2,10) = 0;
	ans(2,11) = 0;
	ans(2,12) = 0;
	ans(2,13) = 0;
	ans(2,14) = 0;
	ans(2,15) = 0;
	ans(3,0) = 0;
	ans(3,1) = 0;
	ans(3,2) = 0;
	ans(3,3) = 1;
	ans(3,4) = 0;
	ans(3,5) = 0;
	ans(3,6) = 0;
	ans(3,7) = 0;
	ans(3,8) = 0;
	ans(3,9) = 0;
	ans(3,10) = 0;
	ans(3,11) = 0;
	ans(3,12) = 0;
	ans(3,13) = 0;
	ans(3,14) = 0;
	ans(3,15) = 0;
	ans(4,0) = 0;
	ans(4,1) = 0;
	ans(4,2) = 0;
	ans(4,3) = 0;
	ans(4,4) = 1;
	ans(4,5) = 0;
	ans(4,6) = 0;
	ans(4,7) = 0;
	ans(4,8) = 0;
	ans(4,9) = 0;
	ans(4,10) = 0;
	ans(4,11) = 0;
	ans(4,12) = 0;
	ans(4,13) = 0;
	ans(4,14) = 0;
	ans(4,15) = 0;
	ans(5,0) = 0;
	ans(5,1) = 0;
	ans(5,2) = 0;
	ans(5,3) = 0;
	ans(5,4) = 0;
	ans(5,5) = 1;
	ans(5,6) = 0;
	ans(5,7) = 0;
	ans(5,8) = 0;
	ans(5,9) = 0;
	ans(5,10) = 0;
	ans(5,11) = 0;
	ans(5,12) = 0;
	ans(5,13) = 0;
	ans(5,14) = 0;
	ans(5,15) = 0;
	ans(6,0) = 0;
	ans(6,1) = 0;
	ans(6,2) = 0;
	ans(6,3) = 0;
	ans(6,4) = 0;
	ans(6,5) = 0;
	ans(6,6) = 1;
	ans(6,7) = 0;
	ans(6,8) = 0;
	ans(6,9) = 0;
	ans(6,10) = 0;
	ans(6,11) = 0;
	ans(6,12) = 0;
	ans(6,13) = 0;
	ans(6,14) = 0;
	ans(6,15) = 0;
	ans(7,0) = 0;
	ans(7,1) = 0;
	ans(7,2) = 0;
	ans(7,3) = 0;
	ans(7,4) = 0;
	ans(7,5) = 0;
	ans(7,6) = 0;
	ans(7,7) = 1;
	ans(7,8) = 0;
	ans(7,9) = 0;
	ans(7,10) = 0;
	ans(7,11) = 0;
	ans(7,12) = 0;
	ans(7,13) = 0;
	ans(7,14) = 0;
	ans(7,15) = 0;
	ans(8,0) = 0;
	ans(8,1) = 0;
	ans(8,2) = 0;
	ans(8,3) = 0;
	ans(8,4) = 0;
	ans(8,5) = 0;
	ans(8,6) = 0;
	ans(8,7) = 0;
	ans(8,8) = 1;
	ans(8,9) = 0;
	ans(8,10) = 0;
	ans(8,11) = 0;
	ans(8,12) = 0;
	ans(8,13) = 0;
	ans(8,14) = 0;
	ans(8,15) = 0;
	ans(9,0) = 0;
	ans(9,1) = 0;
	ans(9,2) = 0;
	ans(9,3) = 0;
	ans(9,4) = 0;
	ans(9,5) = 0;
	ans(9,6) = 0;
	ans(9,7) = 0;
	ans(9,8) = 0;
	ans(9,9) = 1;
	ans(9,10) = 0;
	ans(9,11) = 0;
	ans(9,12) = 0;
	ans(9,13) = 0;
	ans(9,14) = 0;
	ans(9,15) = 0;
	ans(10,0) = 0;
	ans(10,1) = 0;
	ans(10,2) = 0;
	ans(10,3) = 0;
	ans(10,4) = 0;
	ans(10,5) = 0;
	ans(10,6) = 0;
	ans(10,7) = 0;
	ans(10,8) = 0;
	ans(10,9) = 0;
	ans(10,10) = 1;
	ans(10,11) = 0;
	ans(10,12) = 0;
	ans(10,13) = 0;
	ans(10,14) = 0;
	ans(10,15) = 0;
	ans(11,0) = 0;
	ans(11,1) = 0;
	ans(11,2) = 0;
	ans(11,3) = 0;
	ans(11,4) = 0;
	ans(11,5) = 0;
	ans(11,6) = 0;
	ans(11,7) = 0;
	ans(11,8) = 0;
	ans(11,9) = 0;
	ans(11,10) = 0;
	ans(11,11) = 1;
	ans(11,12) = 0;
	ans(11,13) = 0;
	ans(11,14) = 0;
	ans(11,15) = 0;
	ans(12,0) = 0;
	ans(12,1) = 0;
	ans(12,2) = 0;
	ans(12,3) = 0;
	ans(12,4) = 0;
	ans(12,5) = 0;
	ans(12,6) = 0;
	ans(12,7) = 0;
	ans(12,8) = 0;
	ans(12,9) = 0;
	ans(12,10) = 0;
	ans(12,11) = 0;
	ans(12,12) = 1;
	ans(12,13) = 0;
	ans(12,14) = 0;
	ans(12,15) = 0;
	ans(13,0) = 0;
	ans(13,1) = 0;
	ans(13,2) = 0;
	ans(13,3) = 0;
	ans(13,4) = 0;
	ans(13,5) = 0;
	ans(13,6) = 0;
	ans(13,7) = 0;
	ans(13,8) = 0;
	ans(13,9) = 0;
	ans(13,10) = 0;
	ans(13,11) = 0;
	ans(13,12) = 0;
	ans(13,13) = 1;
	ans(13,14) = 0;
	ans(13,15) = 0;
	ans(14,0) = 0;
	ans(14,1) = 0;
	ans(14,2) = 0;
	ans(14,3) = 0;
	ans(14,4) = 0;
	ans(14,5) = 0;
	ans(14,6) = 0;
	ans(14,7) = 0;
	ans(14,8) = 0;
	ans(14,9) = 0;
	ans(14,10) = 0;
	ans(14,11) = 0;
	ans(14,12) = 0;
	ans(14,13) = 0;
	ans(14,14) = 1;
	ans(14,15) = 0;
	ans(15,0) = 0;
	ans(15,1) = 0;
	ans(15,2) = 0;
	ans(15,3) = 0;
	ans(15,4) = 0;
	ans(15,5) = 0;
	ans(15,6) = 0;
	ans(15,7) = 0;
	ans(15,8) = 0;
	ans(15,9) = 0;
	ans(15,10) = 0;
	ans(15,11) = 0;
	ans(15,12) = 0;
	ans(15,13) = 0;
	ans(15,14) = 0;
	ans(15,15) = 1;
	return ans;
}

//////////// UKF sigma points observation function ///////////
template<class Type>
struct ode_integration {
	vector<Type> X_next;
	matrix<Type> P_next;
	ode_integration(vector<Type> x0__, matrix<Type> p0__, Type t, Type dt__, int algo__, Type sigma_x0, Type sigma_x1, Type sigma_x10, Type sigma_x11, Type sigma_x12, Type sigma_x13, Type sigma_x14, Type sigma_x15, Type sigma_x2, Type sigma_x3, Type sigma_x4, Type sigma_x5, Type sigma_x6, Type sigma_x7, Type sigma_x8, Type sigma_x9, Type tau_x0, Type tau_x1, Type tau_x10, Type tau_x11, Type tau_x12, Type tau_x13, Type tau_x14, Type tau_x2, Type tau_x3, Type tau_x4, Type tau_x5, Type tau_x6, Type tau_x7, Type tau_x8, Type tau_x9){
		if(algo__==1){
			/*Forward Euler*/
			X_next = x0__ + f__(tau_x0__, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9)) * dt__;
			P_next = p0__ + (dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9)*p0__ + p0__*dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9).transpose() + g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9)*g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9).transpose()) * dt__;
		} else if (algo__==2){
			/*4th Order Runge-Kutta 4th*/
			vector<Type> X0__ = x0__;
			matrix<Type> P0__ = p0__;
			/**/
			vector<Type> k1__,k2__,k3__,k4__;
			matrix<Type> a1__,a2__,a3__,a4__;
			/*SOLVE ODE*/
			/*step 1*/
			k1__ = dt__ * f__(tau_x0__, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
			a1__ = dt__ * (dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9)*p0__ + p0__*dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9).transpose() + g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9)*g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9).transpose());
			/*step 2*/
			t = t + 0.5 * dt__;
			x0__ = X0__ + 0.5 * k1__;
			p0__ = P0__ + 0.5 * a1__;
      k2__ = dt__ * f__(tau_x0__, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
      a2__ = dt__ * (dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9)*p0__ + p0__*dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9).transpose() + g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9)*g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9).transpose());
			/*step 3*/
			x0__ = X0__ + 0.5 * k2__;
			p0__ = P0__ + 0.5 * a2__;
      k3__ = dt__ * f__(tau_x0__, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
      a3__ = dt__ * (dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9)*p0__ + p0__*dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9).transpose() + g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9)*g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9).transpose());
			/*step 4*/
			t = t + 0.5 * dt__;
			x0__ = X0__ + k3__;
			p0__ = P0__ + a3__;
      k4__ = dt__ * f__(tau_x0__, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
      a4__ = dt__ * (dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9)*p0__ + p0__*dfdx__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9).transpose() + g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9)*g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9).transpose());
			/*ODE UPDATE*/
			X_next = X0__ + (k1__ + 2.0*k2__ + 2.0*k3__ + k4__)/6.0;
			P_next = P0__ + (a1__ + 2.0*a2__ + 2.0*a3__ + a4__)/6.0;
		} else {
			/*nothing*/
		}
	}
};
template<class Type>
matrix<Type> construct_H__(matrix<Type> Xsp, Type void_filler){
  int nn = 33;
  int m = 16;
  matrix<Type> H(m,nn);
  vector<Type> x0__;
  for(int i=0;i<nn;i++){
  x0__ = Xsp.col(i);
  H.col(i) = h__(x0__(0),x0__(1),x0__(10),x0__(11),x0__(12),x0__(13),x0__(14),x0__(15),x0__(2),x0__(3),x0__(4),x0__(5),x0__(6),x0__(7),x0__(8),x0__(9));
  }
return H;
}

//////////// observation variance matrix function ///////////
template<class Type>
matrix<Type> obsvarFun__(Type sigma_X0, Type sigma_X1, Type sigma_X10, Type sigma_X11, Type sigma_X12, Type sigma_X13, Type sigma_X14, Type sigma_X15, Type sigma_X2, Type sigma_X3, Type sigma_X4, Type sigma_X5, Type sigma_X6, Type sigma_X7, Type sigma_X8, Type sigma_X9){
	matrix<Type> V(16,16);
	V(0,0) = sigma_X0;
	V(1,1) = sigma_X1;
	V(2,2) = sigma_X2;
	V(3,3) = sigma_X3;
	V(4,4) = sigma_X4;
	V(5,5) = sigma_X5;
	V(6,6) = sigma_X6;
	V(7,7) = sigma_X7;
	V(8,8) = sigma_X8;
	V(9,9) = sigma_X9;
	V(10,10) = sigma_X10;
	V(11,11) = sigma_X11;
	V(12,12) = sigma_X12;
	V(13,13) = sigma_X13;
	V(14,14) = sigma_X14;
	V(15,15) = sigma_X15;
	return V;
}

//////////// observation variance vector function ///////////
template<class Type>
vector<Type> obsvarFun_tmb__(Type sigma_X0, Type sigma_X1, Type sigma_X10, Type sigma_X11, Type sigma_X12, Type sigma_X13, Type sigma_X14, Type sigma_X15, Type sigma_X2, Type sigma_X3, Type sigma_X4, Type sigma_X5, Type sigma_X6, Type sigma_X7, Type sigma_X8, Type sigma_X9){
	vector<Type> ans(16);
	ans(0) = sigma_X0;
	ans(1) = sigma_X1;
	ans(2) = sigma_X2;
	ans(3) = sigma_X3;
	ans(4) = sigma_X4;
	ans(5) = sigma_X5;
	ans(6) = sigma_X6;
	ans(7) = sigma_X7;
	ans(8) = sigma_X8;
	ans(9) = sigma_X9;
	ans(10) = sigma_X10;
	ans(11) = sigma_X11;
	ans(12) = sigma_X12;
	ans(13) = sigma_X13;
	ans(14) = sigma_X14;
	ans(15) = sigma_X15;
	return ans;
}

//////////// UKF sigma points observation variance matrix function ///////////
template<class Type>
matrix<Type> obsvarFun_usingXsp__(matrix<Type> Xsp, Type sigma_X0, Type sigma_X1, Type sigma_X10, Type sigma_X11, Type sigma_X12, Type sigma_X13, Type sigma_X14, Type sigma_X15, Type sigma_X2, Type sigma_X3, Type sigma_X4, Type sigma_X5, Type sigma_X6, Type sigma_X7, Type sigma_X8, Type sigma_X9){
  matrix<Type> V;
  vector<Type> x0__ = Xsp.col(0);
  V = obsvarFun__(sigma_X0, sigma_X1, sigma_X10, sigma_X11, sigma_X12, sigma_X13, sigma_X14, sigma_X15, sigma_X2, sigma_X3, sigma_X4, sigma_X5, sigma_X6, sigma_X7, sigma_X8, sigma_X9);
  return V;
  }

//////////// objective function ///////////
template<class Type>
Type objective_function<Type>::operator() ()
{
	 DATA_INTEGER(estMethod__);
	 DATA_INTEGER(pred__);
	 DATA_INTEGER(algo__);
	 Type nll__ = 0;

//////////// EKF METHOD ///////////
	 if(estMethod__ == 1){

	 //////////// Estimation //////////////

	 //////////// Estimation //////////////
	 if(pred__ == 0){

//// observations ////
	 DATA_VECTOR(X0);
	 DATA_VECTOR(X1);
	 DATA_VECTOR(X2);
	 DATA_VECTOR(X3);
	 DATA_VECTOR(X4);
	 DATA_VECTOR(X5);
	 DATA_VECTOR(X6);
	 DATA_VECTOR(X7);
	 DATA_VECTOR(X8);
	 DATA_VECTOR(X9);
	 DATA_VECTOR(X10);
	 DATA_VECTOR(X11);
	 DATA_VECTOR(X12);
	 DATA_VECTOR(X13);
	 DATA_VECTOR(X14);
	 DATA_VECTOR(X15);

//// inputs ////
	 DATA_VECTOR(t);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// loss parameters ////
	 DATA_VECTOR(tukey_pars__);
	 DATA_INTEGER(which_loss__);
	 DATA_SCALAR(loss_c_value__);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(tau_x0);
	 PARAMETER(tau_x1);
	 PARAMETER(tau_x2);
	 PARAMETER(tau_x3);
	 PARAMETER(tau_x4);
	 PARAMETER(tau_x5);
	 PARAMETER(tau_x6);
	 PARAMETER(tau_x7);
	 PARAMETER(tau_x8);
	 PARAMETER(tau_x9);
	 PARAMETER(tau_x10);
	 PARAMETER(tau_x11);
	 PARAMETER(tau_x12);
	 PARAMETER(tau_x13);
	 PARAMETER(tau_x14);
	 PARAMETER(sigma_x0);
	 PARAMETER(sigma_x1);
	 PARAMETER(sigma_x2);
	 PARAMETER(sigma_x3);
	 PARAMETER(sigma_x4);
	 PARAMETER(sigma_x5);
	 PARAMETER(sigma_x6);
	 PARAMETER(sigma_x7);
	 PARAMETER(sigma_x8);
	 PARAMETER(sigma_x9);
	 PARAMETER(sigma_x10);
	 PARAMETER(sigma_x11);
	 PARAMETER(sigma_x12);
	 PARAMETER(sigma_x13);
	 PARAMETER(sigma_x14);
	 PARAMETER(sigma_x15);
	 PARAMETER(sigma_X0);
	 PARAMETER(sigma_X1);
	 PARAMETER(sigma_X2);
	 PARAMETER(sigma_X3);
	 PARAMETER(sigma_X4);
	 PARAMETER(sigma_X5);
	 PARAMETER(sigma_X6);
	 PARAMETER(sigma_X7);
	 PARAMETER(sigma_X8);
	 PARAMETER(sigma_X9);
	 PARAMETER(sigma_X10);
	 PARAMETER(sigma_X11);
	 PARAMETER(sigma_X12);
	 PARAMETER(sigma_X13);
	 PARAMETER(sigma_X14);
	 PARAMETER(sigma_X15);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);

//////////// storage variables ///////////
	 vector<vector<Type>> xPrior(t.size());
	 vector<matrix<Type>> pPrior(t.size());
	 vector<vector<Type>> xPost(t.size());
	 vector<matrix<Type>> pPost(t.size());
	 vector<vector<Type>> Innovation(t.size());
	 vector<matrix<Type>> InnovationCovariance(t.size());

//////////// set initial value ///////////
	 vector<Type> x0__ = X0__;
	 matrix<Type> p0__ = P0__;
	 xPrior(0) = X0__;
	 xPost(0) = X0__;
	 pPrior(0) = P0__;
	 pPost(0) = P0__;

	 //////////// initialize variables ///////////
	 int s__;
	 Type half_log2PI = Type(0.5)*log(2*M_PI);
	 vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;
	 matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;

	 //////////// identity matrix ///////////
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<t.size()-1 ; i++){

		 //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
		 for(int j=0 ; j<N__(i) ; j++){
			 ode_integration<Type> odelist = {x0__, p0__, t(i)+j*dt__(i), dt__(i), algo__, sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9, tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9};
			 x0__ = odelist.X_next;
			 p0__ = odelist.P_next;
		 }
		 xPrior(i+1) = x0__;
		 pPrior(i+1) = p0__;

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << X0(i+1), X1(i+1), X2(i+1), X3(i+1), X4(i+1), X5(i+1), X6(i+1), X7(i+1), X8(i+1), X9(i+1), X10(i+1), X11(i+1), X12(i+1), X13(i+1), X14(i+1), X15(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(sigma_X0 ,sigma_X1 ,sigma_X10 ,sigma_X11 ,sigma_X12 ,sigma_X13 ,sigma_X14 ,sigma_X15 ,sigma_X2 ,sigma_X3 ,sigma_X4 ,sigma_X5 ,sigma_X6 ,sigma_X7 ,sigma_X8 ,sigma_X9) * E__.transpose();
			 R__  = C__ * p0__ * C__.transpose() + V__;
			 Ri__ = R__.inverse();
			 K__ 	= p0__ * C__.transpose() * Ri__;
			 x0__ = x0__ + K__*e__;
			 p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();
			 nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction__((e__*(Ri__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + half_log2PI * asDouble(s__);
			 Innovation(i+1) = e__;
			 InnovationCovariance(i+1) = R__;
		 }
		 xPost(i+1) = x0__;
		 pPost(i+1) = p0__;
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(47);
		 vector<Type> map_pars__;
		 parvec__ << tau_x0, tau_x1, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, sigma_x0, sigma_x1, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_X0, sigma_X1, sigma_X2, sigma_X3, sigma_X4, sigma_X5, sigma_X6, sigma_X7, sigma_X8, sigma_X9, sigma_X10, sigma_X11, sigma_X12, sigma_X13, sigma_X14, sigma_X15;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }

	 //////////// Return/Report //////////////
	 REPORT(Innovation);
	 REPORT(InnovationCovariance);
	 REPORT(xPrior);
	 REPORT(xPost);
	 REPORT(pPrior);
	 REPORT(pPost);

	 //////////// Prediction //////////////

	 //////////// Prediction //////////////
	 } else if(pred__ == 1){

//// observations ////
	 DATA_VECTOR(X0);
	 DATA_VECTOR(X1);
	 DATA_VECTOR(X2);
	 DATA_VECTOR(X3);
	 DATA_VECTOR(X4);
	 DATA_VECTOR(X5);
	 DATA_VECTOR(X6);
	 DATA_VECTOR(X7);
	 DATA_VECTOR(X8);
	 DATA_VECTOR(X9);
	 DATA_VECTOR(X10);
	 DATA_VECTOR(X11);
	 DATA_VECTOR(X12);
	 DATA_VECTOR(X13);
	 DATA_VECTOR(X14);
	 DATA_VECTOR(X15);

//// inputs ////
	 DATA_VECTOR(t);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// parameters ////
	 PARAMETER(tau_x0);
	 PARAMETER(tau_x1);
	 PARAMETER(tau_x2);
	 PARAMETER(tau_x3);
	 PARAMETER(tau_x4);
	 PARAMETER(tau_x5);
	 PARAMETER(tau_x6);
	 PARAMETER(tau_x7);
	 PARAMETER(tau_x8);
	 PARAMETER(tau_x9);
	 PARAMETER(tau_x10);
	 PARAMETER(tau_x11);
	 PARAMETER(tau_x12);
	 PARAMETER(tau_x13);
	 PARAMETER(tau_x14);
	 PARAMETER(sigma_x0);
	 PARAMETER(sigma_x1);
	 PARAMETER(sigma_x2);
	 PARAMETER(sigma_x3);
	 PARAMETER(sigma_x4);
	 PARAMETER(sigma_x5);
	 PARAMETER(sigma_x6);
	 PARAMETER(sigma_x7);
	 PARAMETER(sigma_x8);
	 PARAMETER(sigma_x9);
	 PARAMETER(sigma_x10);
	 PARAMETER(sigma_x11);
	 PARAMETER(sigma_x12);
	 PARAMETER(sigma_x13);
	 PARAMETER(sigma_x14);
	 PARAMETER(sigma_x15);
	 PARAMETER(sigma_X0);
	 PARAMETER(sigma_X1);
	 PARAMETER(sigma_X2);
	 PARAMETER(sigma_X3);
	 PARAMETER(sigma_X4);
	 PARAMETER(sigma_X5);
	 PARAMETER(sigma_X6);
	 PARAMETER(sigma_X7);
	 PARAMETER(sigma_X8);
	 PARAMETER(sigma_X9);
	 PARAMETER(sigma_X10);
	 PARAMETER(sigma_X11);
	 PARAMETER(sigma_X12);
	 PARAMETER(sigma_X13);
	 PARAMETER(sigma_X14);
	 PARAMETER(sigma_X15);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);
	 DATA_INTEGER(last_pred_index);
	 DATA_INTEGER(k_step_ahead);
	 vector<matrix<Type>> xk__(last_pred_index);
	 vector<matrix<Type>> pk__(last_pred_index);
	 xk__.fill(matrix<Type>(k_step_ahead+1,n__));
	 pk__.fill(matrix<Type>(k_step_ahead+1,n__*n__));
	 matrix<Type> xk_temp__(k_step_ahead+1,n__);
	 array<Type> pk_temp__(n__,n__,k_step_ahead+1);

//////////// set initial value ///////////
	 vector<Type> x0__ = X0__;
	 matrix<Type> p0__ = P0__;

	 //////////// initialize variables ///////////
	 int s__;
	 vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;
	 matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;

	 //////////// identity matrix ///////////
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<last_pred_index ; i++){

		 xk_temp__.row(0) = x0__;
		 pk_temp__.col(0) = p0__;

	 //////////// K-STEP-AHEAD LOOP ///////////
	 for(int k=0 ; k < k_step_ahead ; k++){

		 //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
		 for(int j=0 ; j<N__(i+k) ; j++){
			 ode_integration<Type> odelist = {x0__, p0__, t(i+k)+j*dt__(i+k), dt__(i+k), algo__, sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9, tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9};
			 x0__ = odelist.X_next;
			 p0__ = odelist.P_next;
		 }

			 //////////// save k-step-ahead prediction ///////////
			 xk_temp__.row(k+1) = x0__;
			 pk_temp__.col(k+1) = p0__;
		 }

		 //////////// save all 0 to k step-ahead predictions ///////////
		 xk__(i) = xk_temp__;
		 for(int kk=0 ; kk < k_step_ahead+1 ; kk++){
			 pk__(i).row(kk) = vector<Type>(pk_temp__.col(kk).transpose());
		 }

		 //////////// rewrite x0 and p0 to one-step predictions ///////////
		 x0__ = xk_temp__.row(1);
		 p0__ = pk_temp__.col(1).matrix();

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << X0(i+1), X1(i+1), X2(i+1), X3(i+1), X4(i+1), X5(i+1), X6(i+1), X7(i+1), X8(i+1), X9(i+1), X10(i+1), X11(i+1), X12(i+1), X13(i+1), X14(i+1), X15(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0), x0__(1), x0__(10), x0__(11), x0__(12), x0__(13), x0__(14), x0__(15), x0__(2), x0__(3), x0__(4), x0__(5), x0__(6), x0__(7), x0__(8), x0__(9));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(sigma_X0 ,sigma_X1 ,sigma_X10 ,sigma_X11 ,sigma_X12 ,sigma_X13 ,sigma_X14 ,sigma_X15 ,sigma_X2 ,sigma_X3 ,sigma_X4 ,sigma_X5 ,sigma_X6 ,sigma_X7 ,sigma_X8 ,sigma_X9) * E__.transpose();
			 R__  = C__ * p0__ * C__.transpose() + V__;
			 Ri__ = R__.inverse();
			 K__ 	= p0__ * C__.transpose() * Ri__;
			 x0__ = x0__ + K__*e__;
			 p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();
		 }
	 }

	 //////////// Return/Report //////////////
	 REPORT(xk__);
	 REPORT(pk__);
	 }

//////////// TMB METHOD ///////////
} else if(estMethod__ == 3) {

//// observations ////
	 DATA_VECTOR(X0);
	 DATA_VECTOR(X1);
	 DATA_VECTOR(X2);
	 DATA_VECTOR(X3);
	 DATA_VECTOR(X4);
	 DATA_VECTOR(X5);
	 DATA_VECTOR(X6);
	 DATA_VECTOR(X7);
	 DATA_VECTOR(X8);
	 DATA_VECTOR(X9);
	 DATA_VECTOR(X10);
	 DATA_VECTOR(X11);
	 DATA_VECTOR(X12);
	 DATA_VECTOR(X13);
	 DATA_VECTOR(X14);
	 DATA_VECTOR(X15);

//// inputs ////
	 DATA_VECTOR(t);

//// state random effect vectors ////
	 PARAMETER_VECTOR(x0);
	 PARAMETER_VECTOR(x1);
	 PARAMETER_VECTOR(x2);
	 PARAMETER_VECTOR(x3);
	 PARAMETER_VECTOR(x4);
	 PARAMETER_VECTOR(x5);
	 PARAMETER_VECTOR(x6);
	 PARAMETER_VECTOR(x7);
	 PARAMETER_VECTOR(x8);
	 PARAMETER_VECTOR(x9);
	 PARAMETER_VECTOR(x10);
	 PARAMETER_VECTOR(x11);
	 PARAMETER_VECTOR(x12);
	 PARAMETER_VECTOR(x13);
	 PARAMETER_VECTOR(x14);
	 PARAMETER_VECTOR(x15);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);

//// time-step ////
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);
	 DATA_IVECTOR(Nc__);

//// iobs vectors ////
	 DATA_IVECTOR(iobs_X0);
	 DATA_IVECTOR(iobs_X1);
	 DATA_IVECTOR(iobs_X2);
	 DATA_IVECTOR(iobs_X3);
	 DATA_IVECTOR(iobs_X4);
	 DATA_IVECTOR(iobs_X5);
	 DATA_IVECTOR(iobs_X6);
	 DATA_IVECTOR(iobs_X7);
	 DATA_IVECTOR(iobs_X8);
	 DATA_IVECTOR(iobs_X9);
	 DATA_IVECTOR(iobs_X10);
	 DATA_IVECTOR(iobs_X11);
	 DATA_IVECTOR(iobs_X12);
	 DATA_IVECTOR(iobs_X13);
	 DATA_IVECTOR(iobs_X14);
	 DATA_IVECTOR(iobs_X15);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(tau_x0);
	 PARAMETER(tau_x1);
	 PARAMETER(tau_x2);
	 PARAMETER(tau_x3);
	 PARAMETER(tau_x4);
	 PARAMETER(tau_x5);
	 PARAMETER(tau_x6);
	 PARAMETER(tau_x7);
	 PARAMETER(tau_x8);
	 PARAMETER(tau_x9);
	 PARAMETER(tau_x10);
	 PARAMETER(tau_x11);
	 PARAMETER(tau_x12);
	 PARAMETER(tau_x13);
	 PARAMETER(tau_x14);
	 PARAMETER(sigma_x0);
	 PARAMETER(sigma_x1);
	 PARAMETER(sigma_x2);
	 PARAMETER(sigma_x3);
	 PARAMETER(sigma_x4);
	 PARAMETER(sigma_x5);
	 PARAMETER(sigma_x6);
	 PARAMETER(sigma_x7);
	 PARAMETER(sigma_x8);
	 PARAMETER(sigma_x9);
	 PARAMETER(sigma_x10);
	 PARAMETER(sigma_x11);
	 PARAMETER(sigma_x12);
	 PARAMETER(sigma_x13);
	 PARAMETER(sigma_x14);
	 PARAMETER(sigma_x15);
	 PARAMETER(sigma_X0);
	 PARAMETER(sigma_X1);
	 PARAMETER(sigma_X2);
	 PARAMETER(sigma_X3);
	 PARAMETER(sigma_X4);
	 PARAMETER(sigma_X5);
	 PARAMETER(sigma_X6);
	 PARAMETER(sigma_X7);
	 PARAMETER(sigma_X8);
	 PARAMETER(sigma_X9);
	 PARAMETER(sigma_X10);
	 PARAMETER(sigma_X11);
	 PARAMETER(sigma_X12);
	 PARAMETER(sigma_X13);
	 PARAMETER(sigma_X14);
	 PARAMETER(sigma_X15);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);
	 DATA_INTEGER(ng__);

	 //////////// initialize variables ///////////
	 vector<Type> F__(n__);
	 vector<Type> Xi__(n__);
	 vector<Type> Xip1__(n__);
	 vector<Type> Z__(n__);
	 matrix<Type> G__(n__,ng__);
	 matrix<Type> V__(n__,n__);
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();
	 I__ *= 1e-8;

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<t.size()-1 ; i++){
		 for(int j=0 ; j<N__(i) ; j++){
			 F__ = f__(tau_x0, tau_x1, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, x0(Nc__(i)+j), x1(Nc__(i)+j), x10(Nc__(i)+j), x11(Nc__(i)+j), x12(Nc__(i)+j), x13(Nc__(i)+j), x14(Nc__(i)+j), x15(Nc__(i)+j), x2(Nc__(i)+j), x3(Nc__(i)+j), x4(Nc__(i)+j), x5(Nc__(i)+j), x6(Nc__(i)+j), x7(Nc__(i)+j), x8(Nc__(i)+j), x9(Nc__(i)+j));
			 G__ = g__(sigma_x0, sigma_x1, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9);
			 Xi__ << x0(Nc__(i)+j), x1(Nc__(i)+j), x2(Nc__(i)+j), x3(Nc__(i)+j), x4(Nc__(i)+j), x5(Nc__(i)+j), x6(Nc__(i)+j), x7(Nc__(i)+j), x8(Nc__(i)+j), x9(Nc__(i)+j), x10(Nc__(i)+j), x11(Nc__(i)+j), x12(Nc__(i)+j), x13(Nc__(i)+j), x14(Nc__(i)+j), x15(Nc__(i)+j);
			 Xip1__ << x0(Nc__(i)+j+1), x1(Nc__(i)+j+1), x2(Nc__(i)+j+1), x3(Nc__(i)+j+1), x4(Nc__(i)+j+1), x5(Nc__(i)+j+1), x6(Nc__(i)+j+1), x7(Nc__(i)+j+1), x8(Nc__(i)+j+1), x9(Nc__(i)+j+1), x10(Nc__(i)+j+1), x11(Nc__(i)+j+1), x12(Nc__(i)+j+1), x13(Nc__(i)+j+1), x14(Nc__(i)+j+1), x15(Nc__(i)+j+1);
			 Z__ = Xip1__ - Xi__ - F__ * dt__(i);
			 V__ = (G__ * G__.transpose() + I__) * dt__(i);
			 nll__ += MVNORM(V__)(Z__);
		 }
	 }

		 //////////// DATA-UPDATE ///////////
	 matrix<Type> varDiag__(m__,t.size());
	 matrix<Type> varFun__(m__,t.size());
	 for(int i=0 ; i<t.size() ; i++){
		 varDiag__.col(i) = obsvarFun_tmb__(sigma_X0, sigma_X1, sigma_X10, sigma_X11, sigma_X12, sigma_X13, sigma_X14, sigma_X15, sigma_X2, sigma_X3, sigma_X4, sigma_X5, sigma_X6, sigma_X7, sigma_X8, sigma_X9);
		 varFun__.col(i) = h__(x0(Nc__(i)), x1(Nc__(i)), x10(Nc__(i)), x11(Nc__(i)), x12(Nc__(i)), x13(Nc__(i)), x14(Nc__(i)), x15(Nc__(i)), x2(Nc__(i)), x3(Nc__(i)), x4(Nc__(i)), x5(Nc__(i)), x6(Nc__(i)), x7(Nc__(i)), x8(Nc__(i)), x9(Nc__(i)));
	 }
	 for(int i=0 ; i<iobs_X0.size() ; i++){
		 int j = iobs_X0(i);
		 nll__ -= dnorm(X0(j),varFun__.col(j)(0),sqrt(varDiag__.col(j)(0)),true);
	 }
	 for(int i=0 ; i<iobs_X1.size() ; i++){
		 int j = iobs_X1(i);
		 nll__ -= dnorm(X1(j),varFun__.col(j)(1),sqrt(varDiag__.col(j)(1)),true);
	 }
	 for(int i=0 ; i<iobs_X2.size() ; i++){
		 int j = iobs_X2(i);
		 nll__ -= dnorm(X2(j),varFun__.col(j)(2),sqrt(varDiag__.col(j)(2)),true);
	 }
	 for(int i=0 ; i<iobs_X3.size() ; i++){
		 int j = iobs_X3(i);
		 nll__ -= dnorm(X3(j),varFun__.col(j)(3),sqrt(varDiag__.col(j)(3)),true);
	 }
	 for(int i=0 ; i<iobs_X4.size() ; i++){
		 int j = iobs_X4(i);
		 nll__ -= dnorm(X4(j),varFun__.col(j)(4),sqrt(varDiag__.col(j)(4)),true);
	 }
	 for(int i=0 ; i<iobs_X5.size() ; i++){
		 int j = iobs_X5(i);
		 nll__ -= dnorm(X5(j),varFun__.col(j)(5),sqrt(varDiag__.col(j)(5)),true);
	 }
	 for(int i=0 ; i<iobs_X6.size() ; i++){
		 int j = iobs_X6(i);
		 nll__ -= dnorm(X6(j),varFun__.col(j)(6),sqrt(varDiag__.col(j)(6)),true);
	 }
	 for(int i=0 ; i<iobs_X7.size() ; i++){
		 int j = iobs_X7(i);
		 nll__ -= dnorm(X7(j),varFun__.col(j)(7),sqrt(varDiag__.col(j)(7)),true);
	 }
	 for(int i=0 ; i<iobs_X8.size() ; i++){
		 int j = iobs_X8(i);
		 nll__ -= dnorm(X8(j),varFun__.col(j)(8),sqrt(varDiag__.col(j)(8)),true);
	 }
	 for(int i=0 ; i<iobs_X9.size() ; i++){
		 int j = iobs_X9(i);
		 nll__ -= dnorm(X9(j),varFun__.col(j)(9),sqrt(varDiag__.col(j)(9)),true);
	 }
	 for(int i=0 ; i<iobs_X10.size() ; i++){
		 int j = iobs_X10(i);
		 nll__ -= dnorm(X10(j),varFun__.col(j)(10),sqrt(varDiag__.col(j)(10)),true);
	 }
	 for(int i=0 ; i<iobs_X11.size() ; i++){
		 int j = iobs_X11(i);
		 nll__ -= dnorm(X11(j),varFun__.col(j)(11),sqrt(varDiag__.col(j)(11)),true);
	 }
	 for(int i=0 ; i<iobs_X12.size() ; i++){
		 int j = iobs_X12(i);
		 nll__ -= dnorm(X12(j),varFun__.col(j)(12),sqrt(varDiag__.col(j)(12)),true);
	 }
	 for(int i=0 ; i<iobs_X13.size() ; i++){
		 int j = iobs_X13(i);
		 nll__ -= dnorm(X13(j),varFun__.col(j)(13),sqrt(varDiag__.col(j)(13)),true);
	 }
	 for(int i=0 ; i<iobs_X14.size() ; i++){
		 int j = iobs_X14(i);
		 nll__ -= dnorm(X14(j),varFun__.col(j)(14),sqrt(varDiag__.col(j)(14)),true);
	 }
	 for(int i=0 ; i<iobs_X15.size() ; i++){
		 int j = iobs_X15(i);
		 nll__ -= dnorm(X15(j),varFun__.col(j)(15),sqrt(varDiag__.col(j)(15)),true);
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(47);
		 vector<Type> map_pars__;
		 parvec__ << tau_x0, tau_x1, tau_x2, tau_x3, tau_x4, tau_x5, tau_x6, tau_x7, tau_x8, tau_x9, tau_x10, tau_x11, tau_x12, tau_x13, tau_x14, sigma_x0, sigma_x1, sigma_x2, sigma_x3, sigma_x4, sigma_x5, sigma_x6, sigma_x7, sigma_x8, sigma_x9, sigma_x10, sigma_x11, sigma_x12, sigma_x13, sigma_x14, sigma_x15, sigma_X0, sigma_X1, sigma_X2, sigma_X3, sigma_X4, sigma_X5, sigma_X6, sigma_X7, sigma_X8, sigma_X9, sigma_X10, sigma_X11, sigma_X12, sigma_X13, sigma_X14, sigma_X15;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }
}
return nll__;
}
