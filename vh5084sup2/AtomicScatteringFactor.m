% Created by Celine Lichtensteiger
% Calculates the different atomic scattering factors used in this program,
% based on the International tables for crystallography Volume C
% http://www.isis.rl.ac.uk/ISISPublic/reference/Xray_scatfac.htm, 
% http://wwwisis2.isis.rl.ac.uk/reference/Xray_scatfac.htm
%*********************************
function[f]=AtomicScatteringFactor(entry,QQ)

switch entry
    case 'Al'
        a1=6.42020; b1=3.03870; a2=1.90020; b2=0.742600; a3=1.59360; b3=31.5472; a4=1.96460; b4=85.0886; c=1.11510;
    case 'Ba'
        a1=20.3361; b1=3.216; a2=19.297; b2=0.2756; a3=10.888; b3=20.2073; a4=2.6959; b4=167.202; c=2.7731;
    case 'Bi'
        a1=33.3689; b1=0.70400; a2=12.9510; b2=2.92380; a3=16.5877; b3=8.79370; a4=6.46920; b4=48.0093; c=13.5782;
    case 'Ca'
        a1=8.62660; b1=10.4421; a2=7.38730; b2=0.659900; a3=1.58990; b3=85.7484; a4=1.02110; b4=178.437; c=1.37510;
    case 'Cu'
        a1=13.3380; b1=3.58280; a2=7.16760; b2=0.247000; a3=5.61580; b3=11.3966; a4=1.67350; b4=64.8126; c=1.19100;
    case 'Dy'
        a1=26.5070; b1=2.18020; a2=17.6383; b2=0.202172; a3=14.5596; b3=12.1899; a4=2.96577; b4=111.874; c=4.29728;
     case 'Fe'
        a1=11.7695; b1=4.76110; a2=7.35730; b2=0.307200; a3=3.52220; b3=15.3535; a4=2.30450; b4=76.8805; c=1.03690;
    case 'Ga'
        a1=15.2354; b1=3.06690; a2=6.70060; b2=0.241200; a3=4.35910; b3=10.7805; a4=2.96230; b4=61.4135; c=1.71890;
    case 'Gd'
        a1=25.0709; b1=2.25341; a2=19.0798; b2=0.181951; a3=13.8518; b3=12.9331; a4=3.54545; b4=101.398; c=2.41960;
    case 'K'
       a1=8.21860; b1=12.7949; a2=7.4398; b2=0.7748; a3=1.0519; b3=213.187; a4=0.8659; b4=41.6841; c=1.4228;
    case 'La'
        a1=20.5780; b1=2.94817; a2=19.5990; b2=0.244475; a3=11.3727; b3=18.7726; a4=3.28719; b4=133.124; c=2.14678;
    case 'Mn'
        a1=11.2819; b1=5.34090; a2=7.35730; b2=0.343200; a3=3.01930; b3=17.8674; a4=2.24410; b4=83.7543; c=1.08960;
    case 'Nd'
        a1=22.6845; b1=2.66248; a2=19.6847; b2=0.210628; a3=12.7740; b3=15.8850; a4=2.85137; b4=137.903; c=1.98486;
    case 'Ni'
        a1=12.8376; b1=3.87850; a2=7.29200; b2=0.256500; a3=4.44380; b3=12.1763; a4=2.38000; b4=66.3421; c=1.03410;
    case 'O'
        a1=3.04850; b1=13.2771; a2=2.28680; b2=5.70110; a3=1.54630; b3=0.323900; a4=0.867000; b4=32.9089; c=0.250800;     
    case 'Pb'
        a1=31.0617; b1=0.69020; a2=13.0637; b2=2.35760; a3=18.4420; b3=8.61800; a4=5.96960; b4=47.2579; c=13.4118;
    case 'Pb2+'
        a1=21.7886; b1=1.33660 ; a2=19.5682; b2=0.488383; a3=19.1406; b3=6.77270; a4=7.01107; b4=23.8132; c=12.4734; 
    case 'Pr'
        a1=22.0440; b1= 2.77393; a2=19.6697; b2=0.222087; a3=12.3856; b3=16.7669; a4=2.82428; b4=143.644; c=2.05830;
    case 'Ru'
        a1=19.2674; b1=0.80852; a2=12.9182; b2=8.43467; a3=4.86337; a3=4.86337; b3=24.7997; a4=1.56756; b4=94.2928; c=5.37874;
    case 'Sc'
        a1=9.18900; b1=9.02130; a2=7.36790; b2=0.572900; a3=1.64090; b3=136.108; a4=1.46800; b4=51.3531; c=1.33290; 
    case 'Si'
        a1=6.29150; b1=2.43860; a2=3.03530; b2=32.3337; a3=1.98910; b3=0.678500; a4=1.54100; b4=81.6937; c=1.14070;  
    case 'Sm'
        a1=24.0042; b1=2.47274; a2=19.4258; b2=0.196451; a3=13.4396; b3=14.3996; a4=2.89604; b4=128.007; c=2.20963; 
    case 'Sr'
        a1=17.5663; b1=1.55640; a2=9.81840; b2=14.0988; a3=5.42200; b3=0.166400; a4=2.66940; b4=132.376; c=2.50640;
    case 'Ta'
        a1=29.2024; b1=1.77333; a2=15.2293; b2=9.37046; a3=14.5135; b3=0.295977; a4=4.76492; b4=63.3644; c=9.24354;
    case 'Tb'
        a1=25.8976; b1=2.24256; a2=18.2185; b2=0.196143; a3=14.3167; b3=12.6648; a4=2.95354; b4=115.362; c=3.58024;
    case 'Ti'
        a1=9.75950; b1=7.85080; a2=7.35580; b2=0.5000000; a3=1.69910; b3=35.6338; a4=1.90210; b4=116.105; c=1.28070;   
    case 'Ti4+'
        a1=19.5114; b1=0.17884; a2=8.23473; b2=6.67018; a3=2.01341; b3=-0.29263; a4=1.52080; b4=12.9464; c=-13.280;   
    case 'V'
        a1=10.2971; b1=6.86570; a2=7.35110; b2=0.438500; a3=2.07030; b3=26.8938; a4=2.05710; b4=102.478; c=1.21990;
    case 'Y'
        a1=17.7760; b1=1.40290; a2=10.2946; b2=12.8006; a3=5.72629; b3=0.125599; a4=3.26588; b4=104.354; c=1.91213; 
    case 'Zr'
        a1=17.8765; b1=1.27618; a2=10.9480; b2=11.9160; a3=5.41732; b3=0.117622; a4=3.65721; b4=87.6627; c=2.06929; 
end;

f=a1.*exp(-b1.*QQ./16./pi^2)+a2.*exp(-b2.*QQ./16./pi^2)+a3.*exp(-b3.*QQ./16./pi^2)+a4.*exp(-b4.*QQ./16./pi^2)+c;exp(-b3.*QQ./16./pi^2)+a4.*exp(-b4.*QQ./16./pi^2)+c;