% Rhabdom params
rhab_z1 = [71164;10142;41887];
rhab_z2 = [70360;21978;35614];

rhab_z = norm(rhab_z1-rhab_z2)

rhab_r1 = [70392;10277;41834];
rhab_r2 = [72065;10804;42260];

rhab_r = norm(rhab_r1-rhab_r2)/2

con_r11 = [69613;7103;43529];
con_r12 = [73822;7671;44099];

con_rmax = norm(con_r11-con_r12)/2

con_r21 = [71277;9499;42079];
con_r22 = [72176;10047;42336];
con_rmin = rhab_r*0.8

con_z = norm(con_r21-con_r11)

con_rhab_sep = con_rmin/3;

con_r31 = [71924;6680;42344];
con_r32 = [72896;7164;42852];

con_rin = norm(con_r31-con_r32)/2;

rhab_pts = [0,0;rhab_r/10,0; rhab_r,rhab_z; 0,rhab_z];
con_z1 = rhab_z + con_rhab_sep;
con_z2 = con_z1 + 4*con_z/5;
con_z3 = con_z1 + con_z;
con_rtmp = 0.75*con_rmax + 0.25*con_rmin;
con_pts = [0,con_z1;con_rmin,con_z1;con_rmax,con_z2;
      con_rtmp,con_z3;
      con_rin, con_z2;
      0,con_z2];

con_lens_sep = con_rmin/6;

lens_z1 = con_z2 + con_lens_sep;
lens_z3 = con_z2 + con_lens_sep + 2/5*con_z;

slope = (con_z3 - con_z2)/(con_rtmp - con_rmax);

lens_z2 = con_z3 - slope/sqrt(slope^2+1)*con_lens_sep;
lens_r2 = con_rtmp - 1.0/sqrt(slope^2+1)*con_lens_sep;
lens_r1 = con_rin -1.0/sqrt(slope^2+1)*con_lens_sep;
lens_pts = [0, lens_z1; lens_r1, lens_z1; lens_r2,lens_z2; lens_r1, lens_z3; 0, lens_z3];

figure(1)
clf()
plot(rhab_pts(:,1),rhab_pts(:,2),'k.-','MarkerSize',20), hold on;
plot(con_pts(:,1),con_pts(:,2),'r.-','MarkerSize',20), hold on;
plot(lens_pts(:,1),lens_pts(:,2),'b.-','MarkerSize',20)
axis equal

x = [con_rtmp;con_rmax;lens_r2];
y = [con_z3;con_z2;lens_z2];
writematrix(rhab_pts,'rhabdom_axissym_pts.dat','Delimiter',' ');
writematrix(con_pts,'cone_axissym_pts.dat','Delimiter',' ');
writematrix(lens_pts,'lens_axissym_pts.dat','Delimiter',' ');



ys = 1000:1000:17000;
xs = zeros(size(ys));
y1s = zeros(size(ys));
y2s = zeros(size(ys));

x1s = zeros(size(ys));
x2s = zeros(size(ys));


y1s(1:13) = rhab_pts(2,2);
y2s(1:13) = rhab_pts(3,2);
x1s(1:13) = rhab_pts(2,1);
x2s(1:13) = rhab_pts(3,1);



y1s(14:16) = con_pts(2,2);
y2s(14:16) = con_pts(3,2);
x1s(14:16) = con_pts(2,1);
x2s(14:16) = con_pts(3,1);


y1s(17) = lens_pts(2,2);
y2s(17) = lens_pts(3,2);
x1s(17) = lens_pts(2,1);
x2s(17) = lens_pts(3,1);


xs = (x2s.*(ys-y1s)-x1s.*(ys-y2s))./(y2s-y1s);
plot(xs,ys,'m.','MarkerSize',20)
