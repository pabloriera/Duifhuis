void rk4(double X[], int nX, double t, double dt,
  void (*derivsRK)(int nX, double X[], double deriv[], double t ) ) {
// Runge-Kutta integrator (4th order)
// Inputs
//   X          Current value of dependent variable
//   nX         Number of elements in dependent variable X
//   t          Independent variable (usually time)
//   dt        Step size (usually time step)
//   derivsRK   Right hand side of the ODE; derivsRK is the
//              name of the function which returns dX/dt
//              Calling format derivsRK(X,t,param,dXdt).
//   param      EXtra parameters passed to derivsRK
// Output
//   X          New value of X after a step of size dt

  double *F1, *F2, *F3, *F4, *Xtemp;
  Xtemp = new double [nX];
  F1 = new double [nX];  F2 = new double [nX];
  F3 = new double [nX];  F4 = new double [nX];


  //* Evaluate F1 = f(X,t).
  (*derivsRK)( nX, X, F1, t);

  //* Evaluate F2 = f( X+dt*F1/2, t+dt/2 ).
  double half_dt = 0.5*dt;
  double t_half = t + half_dt;
  int i;
  for( i=0; i<nX; i++ )
    Xtemp[i] = X[i] + half_dt*F1[i];
  (*derivsRK)( nX, Xtemp, F2, t_half);

  //* Evaluate F3 = f( X+dt*F2/2, t+dt/2 ).
  for( i=0; i<nX; i++ )
    Xtemp[i] = X[i] + half_dt*F2[i];
  (*derivsRK)(nX,  Xtemp, F3, t_half );

  //* Evaluate F4 = f( X+dt*F3, t+dt ).
  double t_full = t + dt;
  for( i=0; i<nX; i++ )
    Xtemp[i] = X[i] + dt*F3[i];
  (*derivsRK)(nX, Xtemp, F4, t_full );

  //* Return X(t+dt) computed from fourth-order R-K.
  for( i=0; i<nX; i++ )
    X[i] += dt/6.*(F1[i] + F4[i] + 2.*(F2[i]+F3[i]));

  delete [] Xtemp;
  delete [] F1;
  delete [] F2;
  delete [] F3;
  delete [] F4;

}
