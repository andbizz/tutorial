functions{
  real[] SIR_equations(real t,      
             real[] y,     
             real[] theta, 
             real[] x_r, 
             int[] x_i)  
{ 
real S;
real I;
real R;

int N;           
 
real dS_dt;
real dI_dt;
real dR_dt;

real beta;
real gamma;

int t_seed;

S = y[1];
I = y[2];
R = y[3];

N = x_i[1];

beta = theta[1];
gamma = theta[2];
  
  
  dS_dt = - beta * I/N * S             ;
  dI_dt =   beta * I/N * S - gamma * I ;
  dR_dt =                    gamma * I ;
  
return{dS_dt,   dI_dt,   dR_dt};
}
}            
             
data {
       int  time_simulation; 
       int  number_data;
       real initial_time;    
       real time_sequence[time_simulation];   
       int  N;
       
       real rel_tol;
       real abs_tol;
       real max_num_steps;
       
       int Y_data[number_data];  
}             
             
transformed data {

    real x_r[0] ;
    int x_i[1] = {N};

}         

parameters {

real<lower=0, upper=1> beta;
real<lower=0, upper=1> gamma;
  
real<lower=1, upper=10> I0;
  
real<lower=0> d_inv;
  }
             
transformed parameters{

  real SIR_ode[time_simulation, 3]; 
  real Initial_Conditions[3];
  real theta[2];
  
  real  I_incidence[time_simulation];
  
  real  d;
  
  Initial_Conditions = {N-I0, I0, 0};
    
  theta[1] = beta;
  theta[2] = gamma;
    
  SIR_ode = integrate_ode_rk45(SIR_equations, Initial_Conditions, initial_time, time_sequence, theta, x_r, x_i, rel_tol , abs_tol ,max_num_steps );
 
 for (i in 1 : time_simulation){
      
    I_incidence[i]   = beta * SIR_ode[i,2]/N * SIR_ode[i,1];
                               }
    d = 1/d_inv;
}             
             
model {
   beta ~ normal(0.3, 0.1);
   gamma ~ normal(0.15,0.1);
   
   I0 ~ normal(2, 10);  
   d_inv ~ exponential(5);
   
     for (i in 1 : number_data){
     target += neg_binomial_2_lpmf(Y_data[i]   | I_incidence[i] + 1e-5, d);}
   
}

generated quantities {
real prediction_SIR[time_simulation];
real R0;

for (i in 1 : time_simulation){
  
prediction_SIR[i]  =  neg_binomial_2_rng(I_incidence[i]+1e-5,d);}
R0 = beta/gamma;
}

