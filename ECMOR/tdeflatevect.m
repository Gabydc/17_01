function[px]=tdeflatevect(z,ei,a,x)
px=a'*x;
px=z'*px;
px=ei*px;
px=z*px;
px=x-px;

