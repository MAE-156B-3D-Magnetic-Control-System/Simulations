%% BB Magnet Simulation

syms Br
Br = 1;                                %% residual flux
a = [1 2];                                 %% cube side length
syms z   %% z is the distance from the magnet

fun = @(z) B(a, Br, z);             %% b field outside magnet

%% To get rid of the warnings:

figure
hold on
xlabel('z/a (distance above the magnet)') %% normalized distance to magnet
ylabel('B', 'rotation', 0, 'VerticalAlignment', 'bottom')

fplot(fun, [0 3]);

xticks([0 1/2 1 3/2 2 3])
xticklabels({'0', '1/2', '1', '3/2', '2', '3'})

yticks([Br/2 Br])
yticklabels({'Br'})


%% Functions


function B = B(a,Br, z)

    t1 = atan((a).^2./(2*z*sqrt((2*z).^2+2*(a).^2)));      
    t2 = atan((a).^2./(2*(a+z).*sqrt(4*(a+z).^2+2.*(a).^2)));   
    
    B = Br/pi*(t1-t2);                                 %% B flux

end
