function ci=interv_conf(datos,dimension)
% Calcula intervalo de confianza con distribucion T?
% dimension=1, calcula el ci de los vectores columna.
[m n]=size(datos);
tmp=isnan(datos);
nan_cols=sum(tmp,2);
nan_rows=sum(tmp);
if dimension==1
    ci=1.96*nanstd(datos,0,1)./sqrt(m-nan_rows);
else
    ci=1.96*nanstd(datos,0,2)./sqrt(n-nan_cols);
end