%close all
if numberhw > 0
figure(1)
[x1,k1] = plot_hw_kali(hw1(:,2),hw1(:,1),1);
tempkorr1 = (215-mean(hw1(:,3)))/ (215 - temperature);
k1_corr = k1.*tempkorr1;
if numberhw > 1
figure(2)
[x2,k2] = plot_hw_kali(hw2(:,2),hw2(:,1),2);
tempkorr2 = (215-mean(hw2(:,3)))/ (215 - temperature);
k2_corr = k2.*tempkorr2;
if numberhw > 2
figure(3)
[x3,k3] = plot_hw_kali(hw3(:,2),hw3(:,1),3);
tempkorr3 = (215-mean(hw3(:,3)))/ (215 - temperature);
k3_corr = k3.*tempkorr3;
if numberhw > 3 
figure(4)
[x4,k4] = plot_hw_kali(hw4(:,2),hw4(:,1),4);
tempkorr4 = (215-mean(hw4(:,3)))/ (215 - temperature);
k4_corr = k4.*tempkorr4;
if numberhw > 4
    error('too many hot wires. Choose between 0 and 4')
end
end
end
end
end

