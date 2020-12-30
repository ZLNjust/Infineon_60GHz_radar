function per_DBF=DBF(A_Calibration,P_Calibration,Wr,per_data_matrix_antenna)
    
    Phase1 = exp(1i*2*pi*d/lamda*[0:Ant_Num-1]'*cos(pi-theta0/180*pi));
    W1= P_Calibration.*Wr.*Phase1.*A_Calibration;
    Rx_DBF_h_3 = W1' * per_data_matrix_antenna([2,3,1],:);
    DBF_h_data_matrix_3=reshape(Rx_DBF_h_3',NL,[])';
    per_DBF = DBF_h_data_matrix_3;
    %per_DBF_3=per_DBF(~any(isnan(per_DBF),2),:);
    
en