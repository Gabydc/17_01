function[px]=deflatevect(z,ei,a,x)
px=z'*x;
px=ei*px;
px=z*px;
px=x-a*px;

