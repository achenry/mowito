function [x,k] = plot_hw_kali(HW,u,which_HW)
k=polyfit(HW,u,4);
x=linspace(0,max(HW),100);
fit=k(1)*x.^4 + k(2)*x.^3 + k(3)*x.^2 + k(4)*x.^1 + k(5)*x.^0;

% figure(which_HW)
plot(x,fit,'b-');hold on;grid on
plot(HW,u,'ro');
title(['HW # ',num2str(which_HW),])
end

