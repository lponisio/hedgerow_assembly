function gi=distLin(kp)

kp_min=min(kp);
kp_max=max(kp);
gi_min=1;% 0.1
gi_max=0.3;% 1
m=(gi_max-gi_min)/(kp_max-kp_min);

gi=m*kp - m*kp_min + gi_min;

end