
#include "model_inference_header.h"



using namespace std;

int main (){

int n_iter = 2000000; //number of iterations for MCMC


para_key para_true;

para_aux para_other;

epi_struct epi_final;

// grid_struct grid_data;

grid_lines_struct grid_lines;

vector<int> index;

IO_simpara (para_true, para_other);  //Importing aux/key parameters 

vector < vector<double> > coordinate(para_other.n,vector<double>(2));

vector < vector<double> > pop_grid (para_other.n_row_grid,vector<double>(para_other.n_col_grid));

IO_simdata(para_true, para_other, coordinate, epi_final,index, pop_grid, grid_lines); //Importing simulated/real data 


/*----------------------------*/

ifstream myfile_in;
ofstream myfile_out; 

vector<int> xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR, xi_null_onset, xi_onset, xi_onset_2; // indices sets indicating the individuals stay in S OR have gone through the other classes (E OR I OR R), and individuals hve gone through E but not I (EnI) and I but not R (InR)


xi_U.reserve(para_other.n); //dynamically updated if necessary
xi_E.reserve(para_other.n);
xi_E_minus.reserve(para_other.n);
xi_I.reserve(para_other.n);
xi_R.reserve(para_other.n);
xi_EnI.reserve(para_other.n);
xi_InR.reserve(para_other.n);
xi_null_onset.reserve(para_other.n); // set of cases without trustworthy onset times (lab tested postive)
xi_onset.reserve(para_other.n); // set of cases with trustworthy onset times (lab tested postive)
xi_onset_2.reserve(para_other.n); // set of cases with trustworthy onset times (lab tested postive); but need to imputed within 1 day range



for (int i=0; i<=(para_other.n-1);i++){
if (epi_final.t_e.at(i)==para_other.unassigned_time) xi_U.push_back(i);
if (epi_final.t_e.at(i)!=para_other.unassigned_time) xi_E.push_back(i);
if (epi_final.t_i.at(i)!=para_other.unassigned_time) xi_I.push_back(i);
if (epi_final.t_r.at(i)!=para_other.unassigned_time & epi_final.t_i.at(i)!=para_other.unassigned_time) xi_R.push_back(i);
if (epi_final.onset_type.at(i)==0) xi_null_onset.push_back(i);
if (epi_final.onset_type.at(i)==1) xi_onset.push_back(i);
if (epi_final.onset_type.at(i)==2) xi_onset_2.push_back(i);

}

int num_null_onset = xi_null_onset.size();
int num_onset = xi_onset.size();
int num_onset_2 = xi_onset_2.size();


xi_E_minus = xi_E;
for (int i=0; i<= (int)(index.size()-1);i++){
xi_E_minus.erase(find(xi_E_minus.begin(),xi_E_minus.end(),index.at(i)));
} // E set excluding index

xi_EnI = xi_E;
for (int i=0; i<= (int)(xi_I.size()-1);i++){
xi_EnI.erase(find(xi_EnI.begin(),xi_EnI.end(),xi_I.at(i)));
} // E set excluding I

xi_InR = xi_I;
for (int i=0; i<= (int)(xi_R.size()-1);i++){
xi_InR.erase(find(xi_InR.begin(),xi_InR.end(),xi_R.at(i)));
} // I set excluding R

myfile_out.open((string(path4)+string("xi_null_onset.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_null_onset.empty()!=1){
for (int i=0; i<=((int)xi_null_onset.size()-1);i++){
myfile_out << xi_null_onset.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_null_onset.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_onset_2.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_onset_2.empty()!=1){
for (int i=0; i<=((int)xi_onset_2.size()-1);i++){
myfile_out << xi_onset_2.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_onset_2.size();
myfile_out.close();


myfile_out.open((string(path4)+string("xi_onset.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_onset.empty()!=1){
for (int i=0; i<=((int)xi_onset.size()-1);i++){
myfile_out << xi_onset.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_onset.size();
myfile_out.close();


/*----------------------------*/

lh_SQUARE lh_square; // the object contains the information of the likelihood contribution of each individual, and it is dynamic and changed during MCMC sampling

lh_square.f_U.assign(para_other.n,1.0);
lh_square.q_T.assign(para_other.n,0.0);
lh_square.kt_sum_U.assign(para_other.n,0.0);
lh_square.f_E.assign(para_other.n,1.0);
lh_square.g_E.assign(para_other.n,1.0);
lh_square.h_E.assign(para_other.n,1.0);
lh_square.k_sum_E.assign(para_other.n,0.0);
lh_square.q_E.assign(para_other.n,0.0);
lh_square.kt_sum_E.assign(para_other.n,0.0);
lh_square.f_I.assign(para_other.n,1.0);
lh_square.f_R.assign(para_other.n,1.0);
lh_square.f_EnI.assign(para_other.n,1.0);
lh_square.f_InR.assign(para_other.n,1.0);

lh_square.f_Grid.assign(para_other.n,1.0);
lh_square.f_Arc.assign(para_other.n,1.0);

lh_square.f_Surv=1.0;

lh_square.f_RR.assign(para_other.n,1.0);
lh_square.p_RorRR=1.0;



/*----------------------------*/


// vector < vector<double> > kernel_mat(para_other.n, vector<double>(para_other.n)), delta_mat(para_other.n, vector<double>(para_other.n));; // a dynamic matrix contain the "kernel distance" between each individual (note:named distance_mat in simulation code)

vector <double> norm_const(para_other.n);


//--
// vector < vector<double> > distance_mat(para_other.n, vector<double>(para_other.n));

// for (int i=0;i<=(para_other.n-1);i++) {
//  for (int j=0;j<=(para_other.n-1);j++) {
//  if (i==j) distance_mat[i][j]=0.0;
//  if (i<j) distance_mat[i][j] = func_distance (coordinate[i][0],coordinate[i][1],coordinate[j][0],coordinate[j][1]);
//  if (i>j) distance_mat[i][j]=distance_mat[j][i];
//  }
// }

//---
// myfile_out.open((string(path4)+string("kernel_matrix_before.txt")).c_str(),ios::app);
// for (int i=0;i<=(para_other.n-1);i++) {
// for (int j=0;j<=(para_other.n-1);j++) {
// if (j<(para_other.n-1)) myfile_out << kernel_mat[i][j] << ",";
// if (j==(para_other.n-1)) myfile_out << kernel_mat[i][j] << " " << endl;
// }
// }
// myfile_out.close();

vector<double> t_e = epi_final.t_e, t_i=epi_final.t_i, t_r=epi_final.t_r; 

vector<int> infected_source = epi_final.infected_source;

// vector<int> gp_stb= epi_final.gp_stb;  
vector<double> stb = epi_final.stb;

myfile_out.open((string(path4)+string("stb_true.txt")).c_str(),ios::out);
for (int i=0; i<=((int)stb.size()-1);i++){
myfile_out  << stb.at(i)  << endl;
}
myfile_out.close();

//para_key para_current = para_true; 






// FUNC func;

// func.set_para(para_true, para_other, coordinate, xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR,t_e, t_i, t_r, stb, index);
// func.initialize_kernel_mat(kernel_mat, norm_const); // compute the kernel matrix
// func.initialize_delta_mat(delta_mat); // compute the kernel matrix

// func.initialize_lh_square(lh_square, kernel_mat, delta_mat, norm_const); //initialize lh_square

// double log_lh_true = log_lh_func (lh_square, para_other.n); // the log-likelihood value when using true values of the parameters


// myfile_out.open((string(path4)+string("kernel_matrix_after.txt")).c_str(),ios::out);
// for (int i=0;i<=(para_other.n-1);i++) {
// for (int j=0;j<=(para_other.n-1);j++) {
// if (j<(para_other.n-1)) myfile_out << kernel_mat[i][j] << ",";
// if (j==(para_other.n-1)) myfile_out << kernel_mat[i][j] << " " << endl;
// }
// }
// myfile_out.close();
/*
myfile_out.open((string(path4)+string("lh_square_q_E_Minus_epi_final_q.txt")).c_str(),ios::out);
//for (int i=0;i<=(para_other.n-1);i++) {
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.q_E.at(xi_E.at(i)) - epi_final.q.at(xi_E.at(i)) << "," <<endl;
}// it should give 0 for all entries, the small deviation seems to do with rounding error in simulation
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_f_U.txt")).c_str(),ios::out);

for (int i=0;i<=(para_other.n-1);i++) {
myfile_out << lh_square.f_U.at(i) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_f_E.txt")).c_str(),ios::out);

for (int i=0;i<=(para_other.n-1);i++) {
myfile_out << lh_square.f_E.at(i) <<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_g_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.g_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_q_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.q_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_h_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.h_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_k_sum_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.k_sum_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_kt_sum_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.kt_sum_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();*/

// myfile_out.open((string(path4)+string("log_value_true.txt")).c_str(),ios::out);
// myfile_out << log_lh_true <<endl;
// myfile_out.close();

/*----------------------------*/

//lh_plot(para_true,para_other, coordinate, xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR, t_e, t_i, t_r, index); // this generates likelihood values at different values of parameters


/*-----------------------------------------------------Start of MCMC sampling------------------------------------------*/

para_key para_current = para_true; 
vector<int> xi_I_current = xi_I;
vector<int> xi_U_current = xi_U;
vector<int> xi_E_current = xi_E;
vector<int> xi_E_minus_current = xi_E_minus;
vector<int> xi_R_current= xi_R;
vector<int> xi_EnI_current = xi_EnI;
vector <int> xi_InR_current =xi_InR;
vector<double> t_e_current =t_e;
vector<double>t_i_current = t_i;
vector<double>t_r_current = t_r;
vector<int>index_current = index;

vector<double> t_obs = t_i; // note that only the entries for valid onset would be used in t_i_update_with_obs

vector<double> stb_current = stb;
// vector<int> gp_stb_current = gp_stb;

//----------------------------------------//
// stb_current.assign(para_other.n,1.0); // when ignoring the stb
//------------------------------------//

vector<int> infected_source_current = infected_source;

vector<double> infected_distance_current(para_other.n);
infected_distance_current.assign(para_other.n,-99);

for (int j=0; j<=((int)xi_E_minus_current.size()-1);j++){

	int source = infected_source_current.at(xi_E_minus_current.at(j));

	if(source!=9999){
		infected_distance_current.at(xi_E_minus_current.at(j)) = func_distance_for_main (coordinate[source][0],coordinate[source][1],coordinate[xi_E_minus_current.at(j)][0],coordinate[xi_E_minus_current.at(j)][1]);
	}
}

//--


lh_SQUARE lh_square_current ;

// vector < vector<double> > kernel_mat_current(para_other.n, vector<double>(para_other.n)), delta_mat_current(para_other.n, vector<double>(para_other.n)); // a dynamic matrix contain the "kernel distance" 
vector < vector<double> > kernel_mat_current,  delta_mat_current; // will not be used

vector <double> norm_const_current(para_other.n);


// para_current.stb_1 = 1.0; //initialization of parameter to be estimated
// para_current.stb_2 = 1.0; //initialization of parameter to be estimated

para_current.alpha = 0.01; //initialization of parameter to be estimated
// para_current.beta = 0.2 ; //initialization of parameter to be estimated
para_current.mu_lat = 7; //mean_IH
para_current.var_lat = 1; //initialization of parameter to be estimated
para_current.c = 10; //initialization of parameter to be estimated
para_current.d = 1; //initialization of parameter to be estimated
para_current.k_1 = 0.1; //initialization of parameter to be estimated
para_current.k_2 = 2; //initialization of parameter to be estimated

para_current.mean_IR = 10 ; // mean IR

para_current.beta_age.resize(5); //initialization of betas for different age gp
para_current.beta_age.at(0) = 0.03; // 0.01
para_current.beta_age.at(1) = 0.01; // 0.01

// para_current.beta_age.at(2) = 0.01;
// para_current.beta_age.at(3) = 0.02;
// para_current.beta_age.at(4) = 0.5;


para_current.gamma_age.resize(5); //initialization of gammas for different age gp
para_current.gamma_age.at(0) = 1;//
para_current.gamma_age.at(1) = 2; //

// para_current.gamma_age.resize(5); //initialization of gammas for different age gp
// para_current.gamma_age.at(0) = 1;//set as 1
// para_current.gamma_age.at(1) = 1.5;
// para_current.gamma_age.at(2) = 2;
// para_current.gamma_age.at(3) = 1;
// para_current.gamma_age.at(4) = 3;

para_current.omega=0.043;

para_current.eta=1;

// https://www.thelancet.com/action/showFullTableHTML?isHtml=true&tableId=tbl3&pii=S1473-3099%2820%2930243-7

para_current.cRR = 14.0; //initialization of parameter to be estimated

para_current.p_hosp_age.resize(5); // (known) initialization of hosp rates for different age gp
para_current.p_hosp_age.at(0) = 0.06 ; // 0-60
para_current.p_hosp_age.at(1) = 0.17; // >60

////---- intialization of xi_E and xi_EnI and xi_U ----////

xi_E_current = xi_I_current; // individuals gone through I (assumed known here) would be initialized as infected
xi_EnI_current.clear();
xi_U_current.clear();
for (int i=0; i<= (int)(para_other.n-1);i++){
if(find(xi_E_current.begin(),xi_E_current.end(),i)==xi_E_current.end()){
xi_U_current.push_back(i);
t_e_current.at(i)=para_other.unassigned_time;
}
}


myfile_out.open((string(path4)+string("xi_E_initial.txt")).c_str(),ios::out);
myfile_out << "k"  << endl;
if (xi_E_current.empty()!=1){
for (int i=0; i<=((int)xi_E_current.size()-1);i++){
myfile_out << xi_E_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_E_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_I_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_I_current.empty()!=1){
for (int i=0; i<=((int)xi_I_current.size()-1);i++){
myfile_out << xi_I_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_I_current.size();
myfile_out.close();



myfile_out.open((string(path4)+string("xi_EnI_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_EnI_current.empty()!=1){
for (int i=0; i<=((int)xi_EnI_current.size()-1);i++){
myfile_out << xi_EnI_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_EnI_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_U_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_U_current.empty()!=1){
for (int i=0; i<=((int)xi_U_current.size()-1);i++){
myfile_out << xi_U_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_U_current.size();
myfile_out.close();


//////------------------------------------------------------------------------//////

lh_square_current.f_U.assign(para_other.n,1.0);
lh_square_current.q_T.assign(para_other.n,0.0);
lh_square_current.kt_sum_U.assign(para_other.n,0.0);
lh_square_current.f_E.assign(para_other.n,1.0);
lh_square_current.g_E.assign(para_other.n,1.0);
lh_square_current.h_E.assign(para_other.n,1.0);
lh_square_current.k_sum_E.assign(para_other.n,0.0);
lh_square_current.q_E.assign(para_other.n,0.0);
lh_square_current.kt_sum_E.assign(para_other.n,0.0);
lh_square_current.f_I.assign(para_other.n,1.0);
lh_square_current.f_R.assign(para_other.n,1.0);
lh_square_current.f_EnI.assign(para_other.n,1.0);
lh_square_current.f_InR.assign(para_other.n,1.0);

lh_square_current.f_Grid.assign(para_other.n,1.0);
lh_square_current.f_Arc.assign(para_other.n,1.0);

lh_square_current.f_Surv=1.0;

lh_square_current.f_RR.assign(para_other.n,1.0);
lh_square_current.p_RorRR =1.0;

FUNC func_mcmc;

func_mcmc.set_para(para_current, para_other, coordinate, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_R_current, xi_EnI_current, xi_InR_current, t_e_current, t_i_current, t_r_current,epi_final.pop, stb_current, index_current, epi_final.age_gp,infected_source_current , pop_grid, grid_lines, epi_final.spm_date,epi_final.onset_type,epi_final.recover_type);




func_mcmc.initialize_lh_square(lh_square_current, kernel_mat_current, delta_mat_current,norm_const_current); //initialize lh_square

myfile_out.open((string(path4)+string("initial_norm_const.txt")).c_str(),ios::app);
for (int i=0;i<=((int)para_other.n-1);i++){
myfile_out<<norm_const_current.at(i)<<endl;
}
myfile_out.close();

double log_lh_current = log_lh_func (lh_square_current, para_other.n); // initialization of log-likelihood value

myfile_out.open((string(path4)+string("initial_lh.txt")).c_str(),ios::app);
myfile_out<<log_lh_current<<endl; //NOTE: this must be defined, otherwise has to re-initialize some components
myfile_out.close();

//---
myfile_out.open((string(path4)+string("initial_f_I.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_I_current.size()-1);i++){
myfile_out<<lh_square_current.f_I.at(xi_I_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_EnI.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_EnI_current.size()-1);i++){
myfile_out<<lh_square_current.f_EnI.at(xi_EnI_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_kt_sum_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.kt_sum_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_g_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.g_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();myfile_out.open((string(path4)+string("initial_h_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.h_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_q_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.q_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.f_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_U.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_U_current.size()-1);i++){
myfile_out<<lh_square_current.f_U.at(xi_U_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_InR.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_InR_current.size()-1);i++){
myfile_out<<lh_square_current.f_InR.at(xi_InR_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_R.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_R_current.size()-1);i++){
myfile_out<<lh_square_current.f_R.at(xi_R_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_RR.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_R_current.size()-1);i++){
myfile_out<<lh_square_current.f_RR.at(xi_R_current.at(i))<<endl;
}
myfile_out.close();


myfile_out.open((string(path4)+string("initial_f_Grid.txt")).c_str(),ios::app);
for (int i=0;i<=(para_other.n-1);i++){
myfile_out<<lh_square_current.f_Grid.at(i)<<endl;
}
myfile_out.close();


myfile_out.open((string(path4)+string("initial_f_Arc.txt")).c_str(),ios::app);
for (int i=0;i<=(para_other.n-1);i++){
myfile_out<<lh_square_current.f_Arc.at(i)<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_Surv.txt")).c_str(),ios::app);
myfile_out<<lh_square_current.f_Surv<<endl;
myfile_out.close();

myfile_out.open((string(path4)+string("initial_p_RorRR.txt")).c_str(),ios::app);
myfile_out<<lh_square_current.p_RorRR<<endl;
myfile_out.close();

//----
// myfile_out.open((string(path4)+string("t_r.txt")).c_str(),ios::app);
// for (int i=0;i<=((int)xi_R_current.size()-1);i++){
// myfile_out<<t_r_current.at(xi_R_current.at(i))<<endl;
// }
// myfile_out.close();

myfile_out.open((string(path4)+string("xi_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_current.size()-1);i++){
myfile_out<<xi_E_current.at(i)<<endl;
}
myfile_out.close();


/*--------------------*/
ifstream myfile1_in, myfile2_in, myfile3_in,myfile4_in, myfile5_in,myfile6_in,myfile7_in,myfile8_in;
ofstream myfile1_out, myfile2_out, myfile3_out,myfile4_out ,myfile5_out,myfile6_out,myfile7_out,myfile8_out; 



mcmc_UPDATE mcmc_update;
mcmc_update.set_para (para_other, coordinate, epi_final.age_gp, epi_final.pop, pop_grid, grid_lines, epi_final.spm_date,epi_final.onset_type,epi_final.recover_type);



myfile1_out.open((string(path4)+string("parameters_current.txt")).c_str(),ios::app);


myfile1_out << "alpha" << "," << "beta_1" << ","<< "beta_2" << "," << "beta_3" << ","<< "beta_4" << ","<< "mu_lat" << "," << "var_lat" << "," << "c" << "," << "d" << "," << "k_1" << "," << "k_2" << "," << "omega" << ","<<  "gamma_1" << ","<< "gamma_2" << "," << "gamma_3" << ","<< "gamma_4" <<endl;


myfile2_out.open((string(path4)+string("log_lh_check.txt")).c_str(),ios::app);
myfile3_out.open((string(path4)+string("log_plh_check.txt")).c_str(),ios::app);

myfile4_out.open((string(path4)+string("index_current.txt")).c_str(),ios::app);

myfile_out.open((string(path4)+string("t_i_mcmc.txt")).c_str(),ios::app);

myfile8_out.open((string(path4)+string("t_r_mcmc.txt")).c_str(),ios::app);


myfile5_out.open((string(path4)+string("t_e_mcmc.txt")).c_str(),ios::app);

myfile6_out.open((string(path4)+string("infected_source_current.txt")).c_str(),ios::app);

myfile7_out.open((string(path4)+string("infected_distance_current.txt")).c_str(),ios::app);

int freq_t_e_up_1 = 1 ; // frequency to translate an infection time
int freq_t_e_up_2 = 1 ; // frequency to update a source


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,0); // set a seed


for (int i=0;i<=(n_iter-1);i++){

mcmc_update.alpha_omega_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);


// mcmc_update.alpha_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);

// mcmc_update.omega_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);



mcmc_update.beta_1_2_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);

// mcmc_update.beta_1_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);
// mcmc_update.beta_2_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,i);




mcmc_update.mu_lat_update(lh_square_current, log_lh_current,xi_U_current, xi_I_current, xi_EnI_current, t_i_current, t_e_current,t_r_current,  index_current, para_current, infected_source_current,i);


mcmc_update.c_update(lh_square_current,log_lh_current, xi_U_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, index_current, para_current, infected_source_current,i);

mcmc_update.k_1_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,i);



for (int j=0; j<=(freq_t_e_up_1 -1);j++){
// for (int j=0; j<=(xi_E_minus_current.size() -1);j++){

	int subject_proposed_1 = xi_E_minus_current.at(gsl_rng_uniform_int (r_c, xi_E_minus_current.size()));
	int subject_proposed_2 = xi_E_minus_current.at(gsl_rng_uniform_int (r_c, xi_E_minus_current.size()));
	int subject_proposed_3 = xi_E_minus_current.at(gsl_rng_uniform_int (r_c, xi_E_minus_current.size()));

	mcmc_update.t_i_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_EnI_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, t_e_current, index_current, para_current, stb_current,  norm_const_current,infected_source_current, subject_proposed_1);
	// mcmc_update.t_i_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_EnI_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, t_e_current, index_current, para_current, stb_current,  norm_const_current,infected_source_current, xi_E_minus_current.at(j));

	mcmc_update.t_i_xx_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_EnI_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, t_e_current, index_current, para_current, stb_current,  norm_const_current,infected_source_current, subject_proposed_2);


	mcmc_update.t_r_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_EnI_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, t_e_current, index_current, para_current, stb_current,  norm_const_current,infected_source_current, subject_proposed_3);


}


// for (int j=0; j<=(freq_t_e_up_2 -1);j++){
for (int j=0; j<=(xi_E_minus_current.size() -1);j++){

	int subject_proposed = xi_E_minus_current.at(j);
	// int subject_proposed = xi_E_minus_current.at(gsl_rng_uniform_int (r_c, xi_E_minus_current.size()));

	mcmc_update.source_t_e_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,infected_distance_current, subject_proposed);

}



myfile1_out << para_current.alpha << "," << para_current.beta_age.at(0) << ","<< para_current.beta_age.at(1) << ","<< para_current.beta_age.at(2) << ","<< para_current.beta_age.at(3) <<  "," << para_current.mu_lat << "," << para_current.var_lat << "," << para_current.c << "," << para_current.d << "," << para_current.k_1 << "," <<  para_current.k_2 << "," <<  para_current.omega << "," << para_current.gamma_age.at(0) << ","<< para_current.gamma_age.at(1) << ","<< para_current.gamma_age.at(2) << ","<< para_current.gamma_age.at(3) << endl;


for (int i_index=0; i_index<=((int)index_current.size()-1);i_index++){
myfile4_out << index_current.at(i_index)<<endl;
}



//---------------//

div_t div_iter;
div_iter = div (i,500);

if ((i>= 500)  & (div_iter.rem==0) & (xi_E_minus_current.empty()==0)){



// output t_i//

	for (int ii=0; ii<= (int) (t_i_current.size()-1);ii++){
		if (ii<(t_i_current.size()-1)) myfile_out  << t_i_current.at(ii) << ",";
		if (ii==(t_i_current.size()-1)) myfile_out  << t_i_current.at(ii) << endl;

	}



	for (int ii=0; ii<= (int) (t_e_current.size()-1);ii++){
		if (ii<(t_e_current.size()-1)) myfile5_out  << t_e_current.at(ii) << ",";
		if (ii==(t_e_current.size()-1)) myfile5_out  << t_e_current.at(ii) << endl;

	}


	for (int ii=0; ii<= (int) (t_r_current.size()-1);ii++){
		if (ii<(t_r_current.size()-1)) myfile8_out  << t_r_current.at(ii) << ",";
		if (ii==(t_r_current.size()-1)) myfile8_out  << t_r_current.at(ii) << endl;

	}


}

//-------------------------------//

div_t div_iter_lh;
div_iter_lh = div (i,100);

if (div_iter_lh.rem==0){

lh_SQUARE lh_square_check;

lh_square_check.f_U.assign(para_other.n,1.0);
lh_square_check.q_T.assign(para_other.n,0.0);
lh_square_check.kt_sum_U.assign(para_other.n,0.0);
lh_square_check.f_E.assign(para_other.n,1.0);
lh_square_check.g_E.assign(para_other.n,1.0);
lh_square_check.h_E.assign(para_other.n,1.0);
lh_square_check.k_sum_E.assign(para_other.n,0.0);
lh_square_check.q_E.assign(para_other.n,0.0);
lh_square_check.kt_sum_E.assign(para_other.n,0.0);
lh_square_check.f_I.assign(para_other.n,1.0);
lh_square_check.f_R.assign(para_other.n,1.0);
lh_square_check.f_EnI.assign(para_other.n,1.0);
lh_square_check.f_InR.assign(para_other.n,1.0);

lh_square_check.f_RR.assign(para_other.n,1.0);
lh_square_check.p_RorRR = 1.0;


lh_square_check.f_Grid.assign(para_other.n,1.0);
lh_square_check.f_Arc.assign(para_other.n,1.0);

lh_square_check.f_Surv = 1.0;

vector < vector<double> > kernel_check = kernel_mat_current;
vector < vector<double> > delta_mat_check = delta_mat_current;
vector<double> norm_const_check =  norm_const_current;

vector<int> infected_source_check =  infected_source_current;

FUNC func_check;

func_check.set_para(para_current, para_other, coordinate, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_R_current, xi_EnI_current, xi_InR_current, t_e_current, t_i_current, t_r_current, epi_final.pop, stb_current, index_current, epi_final.age_gp, infected_source_check, pop_grid, grid_lines, epi_final.spm_date,epi_final.onset_type,epi_final.recover_type);
// func_check.initialize_kernel_mat(kernel_check, norm_const_check); // initialize the kernel matrix
// func_check.initialize_delta_mat(delta_mat_check); // initialize the kernel matrix
func_check.initialize_lh_square(lh_square_check, kernel_check, delta_mat_check, norm_const_check); //initialize lh_square


double log_lh_check = log_lh_func (lh_square_check, para_other.n); // initialization of log-likelihood value

myfile2_out <<  log_lh_current << "," <<   log_lh_check << endl;



double log_plh_current=0.0;
double log_plh_check=0.0;
// log_plh_current = log_plh_current + log(lh_square_current.f_U.at(ik)*lh_square_current.f_I.at(ik)*lh_square_current.f_R.at(ik)*lh_square_current.f_InR.at(ik));
// log_plh_check = log_plh_check + log(lh_square_check.f_U.at(ik)*lh_square_check.f_I.at(ik)*lh_square_check.f_R.at(ik)*lh_square_check.f_InR.at(ik));
log_plh_current = log_plh_current + log(lh_square_current.f_Surv);
log_plh_check = log_plh_check + log(lh_square_check.f_Surv);

myfile3_out << log_plh_current << "," <<   log_plh_check << endl;



for (int js=0;js<=(para_other.n-1);js++){
int rem = (js+1)%para_other.n;
if ((rem!=0) | (js==0)) myfile6_out<< infected_source_current.at(js) << ",";
if ((rem==0) & (js!=0)) myfile6_out <<  infected_source_current.at(js) << " " << endl;
}

for (int js=0;js<=(para_other.n-1);js++){
int rem = (js+1)%para_other.n;
if ((rem!=0) | (js==0)) myfile7_out<< infected_distance_current.at(js) << ",";
if ((rem==0) & (js!=0)) myfile7_out <<  infected_distance_current.at(js) << " " << endl;
}





}


//------------------

} // end of MCMC loop


gsl_rng_free(r_c);


myfile1_out.close();
myfile2_out.close();
myfile3_out.close();
myfile4_out.close();
myfile_out.close();
myfile5_out.close();
myfile6_out.close();
myfile7_out.close();
myfile8_out.close();


myfile1_out.open((string(path4)+string("xi_E_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U"  <<","<< "f_E" << endl;
if (xi_E_current.empty()!=1){
for (int i=0; i<=((int)xi_E_current.size()-1);i++){
myfile1_out << xi_E_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_E_current.at(i)) <<","<< lh_square_current.f_U.at(xi_E_current.at(i)) <<","<< lh_square_current.f_E.at(xi_E_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_E_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_E_minus_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_E_minus_current.empty()!=1){
for (int i=0; i<=((int)xi_E_minus_current.size()-1);i++){
myfile1_out << xi_E_minus_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_E_minus_current.at(i)) <<","<< lh_square_current.f_U.at(xi_E_minus_current.at(i)) <<","<< lh_square_current.f_E.at(xi_E_minus_current.at(i))<< endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_E_minus_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_EnI_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_EnI_current.empty()!=1){
for (int i=0; i<=((int)xi_EnI_current.size()-1);i++){
myfile1_out << xi_EnI_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_EnI_current.at(i)) <<","<< lh_square_current.f_U.at(xi_EnI_current.at(i)) <<","<< lh_square_current.f_E.at(xi_EnI_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_EnI_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_U_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_U_current.empty()!=1){
for (int i=0; i<=((int)xi_U_current.size()-1);i++){
myfile1_out << xi_U_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_U_current.at(i)) <<","<< lh_square_current.f_U.at(xi_U_current.at(i)) <<","<< lh_square_current.f_E.at(xi_U_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_U_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_I_final.txt")).c_str(),ios::out);
myfile1_out << "k" << endl;
if (xi_I_current.empty()!=1){
for (int i=0; i<=((int)xi_I_current.size()-1);i++){
myfile1_out << xi_I_current.at(i) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_I_current.size();
myfile1_out.close();



myfile1_out.open((string(path4)+string("xi_InR_final.txt")).c_str(),ios::out);
myfile1_out << "k" << endl;
if (xi_InR_current.empty()!=1){
for (int i=0; i<=((int)xi_InR_current.size()-1);i++){
myfile1_out << xi_InR_current.at(i) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_InR_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("index_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "t_e"<< "," <<  "min_t_e"<<endl;
if (index_current.empty()!=1){
for (int i=0; i<=((int)index_current.size()-1);i++){
myfile1_out << index_current.at(i) << "," << t_e_current.at(index_current.at(0)) << ","<< *min_element(t_e_current.begin(),t_e_current.end())<<endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << index_current.size();
myfile1_out.close();

//------------------

return(0);
}


