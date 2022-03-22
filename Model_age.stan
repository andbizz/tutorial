functions{
  real[] SIR_age_equations(real t,      
             real[] y,     
             real[] theta, 
             real[] x_r, 
             int[] x_i)  
{ 
real S1;
real I1;
real R1;
real S2;
real I2;
real R2;

int N1;
int N2;
 
real C11;
real C21;
real C12;
real C22;
 
real dS1_dt;
real dI1_dt;
real dR1_dt;
real dS2_dt;
real dI2_dt;
real dR2_dt;

real beta;
real gamma1;
real gamma2;

real lambda1;
real lambda2;

S1 = y[1];
I1 = y[2];
R1 = y[3];
S2 = y[4];
I2 = y[5];
R2 = y[6];

N1 = x_i[1];
N2 = x_i[2];

C11 = x_r[1];
C12 = x_r[2];
C21 = x_r[3];
C22 = x_r[4];

beta = theta[1];
gamma1 = theta[2];
gamma2 = theta[2];  
  
lambda1 = beta * (C11 * I1 / N1 + C12 * I2 / N2);
lambda2 = beta * (C21 * I1 / N1 + C22 * I2 / N2);    
    
  dS1_dt = - lambda1 * S1;
  dI1_dt =   lambda1 * S1 - gamma1 * I1;
  dR1_dt =                          gamma1 * I1;
  dS2_dt = - lambda2 * S2;
  dI2_dt =   lambda2 * S2 - gamma2 * I2;
  dR2_dt =                          gamma2 * I2;  
  
return{dS1_dt, dI1_dt, dR1_dt, dS2_dt, dI2_dt, dR2_dt};
}
}            
             
data {
       int  time_simulation; 
       int  number_data;
       real initial_time;    
       real time_sequence[time_simulation];   
       int  N1;
       int  N2;
       real contact_matrix[2,2];
       
       real rel_tol;
       real abs_tol;
       real max_num_steps;
       
       int Y1_data[number_data];
       int Y2_data[number_data];
       
}             
             
transformed data {
 
    real C11 = contact_matrix[1,1];
    real C12 = contact_matrix[1,2];    
    real C21 = contact_matrix[2,1];    
    real C22 = contact_matrix[2,2];
    
    real x_r[4] = {C11, C12, C21, C22} ;
    int  x_i[2] = {N1, N2};
    
}         

parameters {

real<lower=0, upper=1> beta;
real<lower=0, upper=1> gamma1;
real<lower=0, upper=1> gamma2;
  
real<lower=1, upper=10> I10;
real<lower=1, upper=10> I20;
  
real<lower=0> d_inv;
  }
             
transformed parameters{

  real SIR_age_ode[time_simulation, 6]; 
  real Initial_Conditions[6];
  real theta[3];
  
  real  I1_incidence[time_simulation];
  real  I2_incidence[time_simulation];
  
  real  d;
  
  Initial_Conditions = {N1-I10, I10, 0,N2-I20, I20, 0 };
    
  theta[1] = beta;
  theta[2] = gamma1;
  theta[3] = gamma2; 
    
  SIR_age_ode = integrate_ode_rk45(SIR_age_equations, Initial_Conditions, initial_time, time_sequence, theta, x_r, x_i, rel_tol , abs_tol ,max_num_steps );
 
 for (i in 1 : time_simulation){
    I1_incidence[i]   = beta * (C11 * SIR_age_ode[i,2]/N1 + C12 * SIR_age_ode[i,5]/N2) * SIR_age_ode[i,1];
    I2_incidence[i]   = beta * (C21 * SIR_age_ode[i,2]/N1 + C22 * SIR_age_ode[i,5]/N2) * SIR_age_ode[i,4];
                               }
    d = 1/d_inv;
}             
             
model {
   beta ~ normal(0.3, 0.1);
   gamma1 ~ normal(0.1,0.5);
   gamma2 ~ normal(0.1,0.5);
   
   I10 ~ normal(2, 10);  
   I20 ~ normal(2, 10);
   d_inv ~ exponential(5);
   
     for (i in 1 : number_data){
     target += neg_binomial_2_lpmf(Y1_data[i]   | I1_incidence[i] + 1e-4, d);
     target += neg_binomial_2_lpmf(Y2_data[i]   | I2_incidence[i] + 1e-4, d);
     }
   
}

generated quantities {
real prediction1_SIR[time_simulation];
real prediction2_SIR[time_simulation];

for (i in 1 : time_simulation){
prediction1_SIR[i]  =  neg_binomial_2_rng(I1_incidence[i]+1e-5,d);
prediction2_SIR[i]  =  neg_binomial_2_rng(I2_incidence[i]+1e-5,d);
}
}

