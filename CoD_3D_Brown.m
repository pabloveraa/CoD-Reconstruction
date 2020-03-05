function CoD_3D_Brown()
%
% Simulation for computing the reconstruction error that results from
% an inaccurate estimation of the camera lens center of distortion using
% the Brown model with 2 radial distortion coefficients
%
create_figure();  %figure for the interface
fun_defaults();  %write defaults in text boxes

%========================================================
function create_figure()
%create the main figure for the interface
main_fig = figure('name','fig_main','units','normalized','position', ...
    [0.05 0.25 0.8 0.6],'menubar','none','numbertitle','off');
set_frame_world(main_fig);
set_frame_intrinsics(main_fig);
set_frame_show_plots(main_fig);
set_frame_view1(main_fig);
set_frame_view2(main_fig);
set_frame_error(main_fig);

%=================================================
function set_frame_world(main_fig)
global r
p = uipanel(main_fig,'title','3D world points','units','normalized', ...
    'position',[0.03 0.51 0.44 0.47],'fontsize',12);
x = [0.05 0.18 0.41 0.64 0.77];
y = [0.07 0.21 0.37 0.54 0.70 0.84];
w = 0.09;
h = 0.12;
fsz = 11;  %font size

%coordinate limits
uicontrol(p,'style','text','units','normalized','position',[x(1) y(5) w h], ...
    'string','X:','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(1) y(4) w h], ...
    'string','Y:','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(1) y(3) w h], ...
    'string','Z:','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(2) y(6) 2*w h], ...
    'string','min','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(3) y(6) 2*w h], ...
    'string','max','fontsize',fsz);
r.Xmin = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(5) 2*w h],'fontsize',fsz);
r.Xmax = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(5) 2*w h],'fontsize',fsz);
r.Ymin = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(4) 2*w h],'fontsize',fsz);
r.Ymax = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(4) 2*w h],'fontsize',fsz);
r.Zmin = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(3) 2*w h],'fontsize',fsz);
r.Zmax = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(3) 2*w h],'fontsize',fsz);

%angles
uicontrol(p,'style','text','units','normalized','position',[x(4) y(6) 3*w h], ...
    'string','angles (degrees)','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(4) y(5) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\alpha_x$:','interpreter','latex','fontsize',13);
annotation(p,'textbox','units','normalized','position',[x(4) y(4) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\alpha_y$:','interpreter','latex','fontsize',13);
annotation(p,'textbox','units','normalized','position',[x(4) y(3) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\alpha_z$:','interpreter','latex','fontsize',13);
r.alpha_x = uicontrol(p,'style','edit','units','normalized','position',[x(5) y(5) 2*w h],'fontsize',fsz);
r.alpha_y = uicontrol(p,'style','edit','units','normalized','position',[x(5) y(4) 2*w h],'fontsize',fsz);
r.alpha_z = uicontrol(p,'style','edit','units','normalized','position',[x(5) y(3) 2*w h],'fontsize',fsz);

%metric units
uicontrol(p,'style','text','units','normalized','position',[x(1) y(2) w h], ...
    'string','units:','fontsize',fsz);
bg = uibuttongroup(p,'visible','on','units','normalized','position',[x(2) y(2) 2.4*w 1.2*h]);
r.units_m = uicontrol(bg,'style','radiobutton','units','normalized', 'position',[0.05 0.05 0.4 0.9], ...
    'string','m','fontsize',fsz);
r.units_cm = uicontrol(bg,'style','radiobutton','units','normalized', 'position',[0.5 0.05 0.45 0.9], ...
    'string','cm','fontsize',fsz);

%number of points and random seed generator
uicontrol(p,'style','text','units','normalized','position',[x(4) y(2) 3*w h], ...
    'string','number of points:','fontsize',fsz);
r.num_points = uicontrol(p,'style','edit','units','normalized','position',[x(4) y(1) 2*w h],'fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(1) y(1) 3*w h], ...
    'string','random seed:','fontsize',fsz);
r.rand_seed = uicontrol(p,'style','popupmenu','units','normalized','position',[x(2)+2*w y(1) w h], ...
    'string','1|2|3|4|5|6|7|8|9|10','fontsize',fsz);

%=================================================
function set_frame_intrinsics(main_fig)
global r
p = uipanel(main_fig,'title','intrinsic parameters','units','normalized', ...
    'position',[0.50 0.51 0.33 0.47],'fontsize',12);
x = [0.06 0.24 0.53 0.71];
y = [0.09 0.25 0.45 0.61 0.80];
w = 0.12;
h = 0.13;
fsz = 11;  %font size

%focal length
uicontrol(p,'style','text','units','normalized','position',[x(2) y(5) 4*w h], ...
    'string','focal length (f) (pixels):','fontsize',fsz);
r.focal_length = uicontrol(p,'style','edit','units','normalized','position',[x(4) y(5) 2*w h],'fontsize',fsz);

%center of distortion
uicontrol(p,'style','text','units','normalized','position',[x(1) y(4) 6*w h], ...
    'string','center of distortion (CoD) (pixels):','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(3) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$x_0$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(3) y(3) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$y_0$:','interpreter','latex','fontsize',14);
r.x_0 = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(3) 2*w h],'fontsize',fsz);
r.y_0 = uicontrol(p,'style','edit','units','normalized','position',[x(4) y(3) 2*w h],'fontsize',fsz);

%image resolution
uicontrol(p,'style','text','units','normalized','position',[x(1) y(2) 4*w h], ...
    'string','image resolution (pixels):','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(1) y(1) 1.5*w h], ...
    'string','width:','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(3) y(1) 1.5*w h], ...
    'string','height:','fontsize',fsz);
r.cols = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(1) 2*w h],'fontsize',fsz);
r.rows = uicontrol(p,'style','edit','units','normalized','position',[x(4) y(1) 2*w h],'fontsize',fsz);

%=======================================================
function set_frame_show_plots(main_fig)
p = uipanel(main_fig,'title','show plots','units','normalized', ...
    'position',[0.86 0.51 0.11 0.47],'fontsize',12);
x = 0.05;
y = [0.04 0.23 0.42 0.62  0.81];
w = 0.9;
h = 0.15;
fsz = 11;  %font size
uicontrol(p,'style','pushbutton','units','normalized','position',[x y(5) w h], ...
    'string','3D plot','fontsize',fsz,'callback',@fun_3d_plot);
uicontrol(p,'style','pushbutton','units','normalized','position',[x y(4) w h], ...
    'string','2D plot','fontsize',fsz,'callback',@fun_2d_plot);
uicontrol(p,'style','pushbutton','units','normalized','position',[x y(3) w h], ...
    'string','error plot','fontsize',fsz,'callback',@fun_error_plot);
uicontrol(p,'style','pushbutton','units','normalized','position',[x y(2) w h], ...
    'string','defaults','fontsize',fsz,'callback',@fun_defaults);
uicontrol(p,'style','pushbutton','units','normalized','position',[x y(1) w h], ...
    'string','exit','fontsize',fsz,'callback',{@fun_exit,main_fig});

%=======================================================
function set_frame_view1(main_fig)
global r
p = uipanel(main_fig,'title','view 1 (reference frame)','units','normalized', ...
    'position',[0.03 0.02 0.33 0.47],'fontsize',12);
x = [0.06 0.38 0.69];
y = [0.08 0.22 0.37 0.54 0.68 0.83];
w = 0.25;
h = 0.12;
fsz = 11;  %font size

%rotation angles
uicontrol(p,'style','text','units','normalized','position',[x(1) y(6) 3*w h], ...
    'string','rotation angles (degrees)','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_x$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(2) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_y$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(3) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_z$:','interpreter','latex','fontsize',14);
r.theta_x1 = uicontrol(p,'style','edit','units','normalized','position',[x(1) y(4) w h], ...
    'fontsize',fsz,'enable','off');
r.theta_y1 = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(4) w h], ...
    'fontsize',fsz,'enable','off');
r.theta_z1 = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(4) w h], ...
    'fontsize',fsz,'enable','off');

%translation vector
uicontrol(p,'style','text','units','normalized','position',[x(1) y(3) 2*w h], ...
    'string','traslation vector','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_x$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(2) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_y$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(3) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_z$:','interpreter','latex','fontsize',14);
r.tx1 = uicontrol(p,'style','edit','units','normalized','position',[x(1) y(1) w h], ...
    'fontsize',fsz,'enable','off');
r.ty1 =uicontrol(p,'style','edit','units','normalized','position',[x(2) y(1) w h], ...
    'fontsize',fsz,'enable','off');
r.tz1 = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(1) w h], ...
    'fontsize',fsz,'enable','off');

%=======================================================
function set_frame_view2(main_fig)
global r
p = uipanel(main_fig,'title','view 2','units','normalized', ...
    'position',[0.38 0.02 0.33 0.47],'fontsize',12);
x = [0.06 0.38 0.69];
y = [0.08 0.22 0.37 0.54 0.68 0.83];
w = 0.25;
h = 0.12;
fsz = 11;  %font size

%rotation angles
uicontrol(p,'style','text','units','normalized','position',[x(1) y(6) 3*w h], ...
    'string','rotation angles (degrees)','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_x$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(2) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_y$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(3) y(5) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\theta_z$:','interpreter','latex','fontsize',14);
r.theta_x2 = uicontrol(p,'style','edit','units','normalized','position',[x(1) y(4) w h],'fontsize',fsz);
r.theta_y2 = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(4) w h],'fontsize',fsz);
r.theta_z2 = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(4) w h],'fontsize',fsz);

%translation vector
uicontrol(p,'style','text','units','normalized','position',[x(1) y(3) 2*w h], ...
    'string','traslation vector','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_x$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(2) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_y$:','interpreter','latex','fontsize',14);
annotation(p,'textbox','units','normalized','position',[x(3) y(2) 0.5*w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$t_z$:','interpreter','latex','fontsize',14);
r.tx2 = uicontrol(p,'style','edit','units','normalized','position',[x(1) y(1) w h],'fontsize',fsz);
r.ty2 = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(1) w h],'fontsize',fsz);
r.tz2 = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(1) w h],'fontsize',fsz);

%=================================================
function set_frame_error(main_fig)
global r
p = uipanel(main_fig,'title','reconstruction error','units','normalized', ...
    'position',[0.73 0.02 0.24 0.47],'fontsize',12);
x = [0.025 0.175 0.45 0.725];
y = [0.058 0.207 0.356 0.506 0.678 0.828];
w = 0.125;
h = 0.115;
fsz = 11;  %font size

%CoD displacement
uicontrol(p,'style','text','units','normalized','position',[x(1) y(6) 4*w h], ...
    'string','CoD displacement','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(5) w h],'edgecolor',[0.8 0.8 0.8], ... 
    'string','$\delta$:','interpreter','latex','fontsize',14);
r.delta_cod = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(5) 2*w h],'fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(3) y(5) w h], ...
    'string','pixels','fontsize',fsz);

%distortion coefficients
uicontrol(p,'style','text','units','normalized','position',[x(1) y(4) 6*w h], ...
    'string','distortion coefficients','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(2) y(3) 2*w h], ...
    'string','min','fontsize',fsz);
uicontrol(p,'style','text','units','normalized','position',[x(3) y(3) 2*w h], ...
    'string','max','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(2) w h],'edgecolor',[0.8 0.8 0.8], ...
    'string','$\kappa_1$:','interpreter','latex','fontsize',14);
r.kappa1_min = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(2) 2*w h],'fontsize',fsz);
r.kappa1_max = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(2) 2*w h],'fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(4) y(2) 2*w h],'edgecolor',[0.8 0.8 0.8], ...
    'string','pixels$^{-2}$','interpreter','latex','fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(1) y(1) w h],'edgecolor',[0.8 0.8 0.8], ...
    'string','$\kappa_2$:','interpreter','latex','fontsize',14);
r.kappa2_min = uicontrol(p,'style','edit','units','normalized','position',[x(2) y(1) 2*w h],'fontsize',fsz);
r.kappa2_max = uicontrol(p,'style','edit','units','normalized','position',[x(3) y(1) 2*w h],'fontsize',fsz);
annotation(p,'textbox','units','normalized','position',[x(4) y(1) 2*w h],'edgecolor',[0.8 0.8 0.8], ...
    'string','pixels$^{-4}$','interpreter','latex','fontsize',fsz);

%===============================================
function fun_3d_plot(~,~)
global r
%compute the 3D world points
[X,Y,Z] = compute_world_points();
%camera view angles and translation
[rg,t] = camera_poses();
%plot the points and the camera poses
figure('name','world_points','numbertitle','off');
plot3(X,Y,Z,'b.'); hold on; grid on;
set(gca,'fontsize',12);
if( r.units_m.Value==1 )
    units = '(m)';
else
    units = '(cm)';
end
xlabel(strcat('$X ',units,'$'),'interpreter','latex','fontsize',16);
ylabel(strcat('$Y ',units,'$'),'interpreter','latex','fontsize',16);
zlabel(strcat('$Z ',units,'$'),'interpreter','latex','fontsize',16);

R1 = rotation_matrix(rg{1}(1),rg{1}(2),rg{1}(3));
R2 = rotation_matrix(rg{2}(1),rg{2}(2),rg{2}(3));
plotCamera('Location',t{1},'Orientation',R1,'Size',6);
plotCamera('Location',t{2},'Orientation',R2,'Size',6);

%===============================================
function fun_2d_plot(~,~)
global r
%compute the 3D world points
[X,Y,Z] = compute_world_points();
%camera view angles and translation
[rg,t] = camera_poses();
[x1,y1,indz1] = ideal_point_projection(X,Y,Z,rg{1},t{1});
[x2,y2,indz2] = ideal_point_projection(X,Y,Z,rg{2},t{2});
ind = (indz1==1) & (indz2==1);

%image resolution
cols = str2double(r.cols.String);
rows = str2double(r.rows.String);

hf = figure('name','image_points','units','normalized','position',[0.12 0.08 0.34 0.78],'numbertitle','off');
p = uipanel(hf,'title','camera view 1 - ideal points','units','normalized','position',[0.05 0.52 0.9 0.46]);
ax = axes(p);
plot(ax,x1(ind),y1(ind),'b.'); grid on;
axis([1 cols 1 rows]);
xlabel('x (pixels)');
ylabel('y (pixels)');

p = uipanel(hf,'title','camera view 2 - ideal points','units','normalized','position',[0.05 0.03 0.9 0.46]);
ax = axes(p);
plot(ax,x2(ind),y2(ind),'b.'); grid on;
axis([1 cols 1 rows]);
xlabel('x (pixels)');
ylabel('y (pixels)');

%==================================================
function fun_error_plot(~,~)
global r
%compute the 3D world points
[X,Y,Z] = compute_world_points();

%camera view angles and translation
[rg,t] = camera_poses();
%(x1,y1) and (x2,y2): points projected to the images 1 and 2 respectively
%points infront of the camera: indz1=1 and indz2=1
%points behind the camera: indz1=0 and indz2=0
[x1,y1,indz1] = ideal_point_projection(X,Y,Z,rg{1},t{1});
[x2,y2,indz2] = ideal_point_projection(X,Y,Z,rg{2},t{2});
indz = (indz1==1) & (indz2==1);
ptu1 = [x1(indz), y1(indz)];
ptu2 = [x2(indz), y2(indz)];
P = [X(indz), Y(indz), Z(indz)];

%find best control point
x0 = str2double(r.x_0.String);
y0 = str2double(r.y_0.String);
r1 = sqrt((x1(indz)-x0).^2 + (y1(indz)-y0).^2);
r2 = sqrt((x2(indz)-x0).^2 + (y2(indz)-y0).^2);
rs = r1 + r2;
cp = find(rs==min(rs));
cp = cp(1);

%image resolution
cols = str2double(r.cols.String);
rows = str2double(r.rows.String);
ft = str2double(r.focal_length.String);  %focal length

nd = 5;
kappa1_min = str2double(r.kappa1_min.String);
kappa1_max = str2double(r.kappa1_max.String);
kappa2_min = str2double(r.kappa2_min.String);
kappa2_max = str2double(r.kappa2_max.String);
d_kappa1 = (kappa1_max - kappa1_min)/(nd-1);
d_kappa2 = (kappa2_max - kappa2_min)/(nd-1);
[kappa1,kappa2] = meshgrid(kappa1_min:d_kappa1:kappa1_max,kappa2_min:d_kappa2:kappa2_max);

delta = str2double(r.delta_cod.String);
theta = (0:45:315)';
e_mean = zeros(nd,nd,2);
e_std = zeros(nd,nd,2);
nc_min = inf;  %minimum number of correspondences
wb = waitbar(0,'computing reconstruction error');
for i=1:nd
    for j=1:nd
        ptd1 = apply_distortion(ptu1,kappa1(i,j),kappa2(i,j),x0,y0);
        ptd2 = apply_distortion(ptu2,kappa1(i,j),kappa2(i,j),x0,y0);
        %points inside the images
        ind = (ptd1(:,1)>=1 & ptd1(:,1)<=cols) & (ptd1(:,2)>=1 & ptd1(:,2)<=rows) ...
            & (ptd2(:,1)>=1 & ptd2(:,1)<=cols) & (ptd2(:,2)>=1 & ptd2(:,2)<=rows);
        nc_min = min(nc_min, sum(ind));
        if( sum(ind)>=8 )
            epsilon = reconstruction_error(ptd1,ptd2,x0,y0,delta,theta,kappa1(i,j),kappa2(i,j),P,cp,ft,ind);
            e_mean(i,j,:) = mean(epsilon);
            e_std(i,j,:) = std(epsilon);
        end
        waitbar((i-1)/nd + j/nd^2,wb,'computing reconstruction error');
    end
end
close(wb);
plot_error(e_mean,e_std,kappa1,kappa2);
if( nc_min<8 )
    msgbox('Image correspondences are insufficient. Results might be inaccurate.','error','error');
end

%====================================================
function fun_defaults(~,~)
%write defaults in text boxes
% pr: intrinsic parameters, rg: rotation angles, t: translation vectors
% wp: parameters to generate 3D world points
global r
pr = load_intrinsic_parameters();
[rg,t] = load_extrinsic_parameters();
wp = load_world_parameters();

%3D world points frame
r.Xmin.String = sprintf('%0.1f',wp.axis(1));
r.Xmax.String = sprintf('%0.1f',wp.axis(2));
r.Ymin.String = sprintf('%0.1f',wp.axis(3));
r.Ymax.String = sprintf('%0.1f',wp.axis(4));
r.Zmin.String = sprintf('%0.1f',wp.axis(5));
r.Zmax.String = sprintf('%0.1f',wp.axis(6));
r.alpha_x.String = sprintf('%d',wp.angles(1));
r.alpha_y.String = sprintf('%d',wp.angles(2));
r.alpha_z.String = sprintf('%d',wp.angles(3));
r.num_points.String = sprintf('%d',wp.num_points);
r.rand_seed.Value = wp.rseed;
if( wp.units==1 )
    r.units_m.Value = 1;
else
    r.units_cm.Value = 1;
end

%instrinsic parameters frame
r.focal_length.String = sprintf('%d',pr.ft);
r.x_0.String = sprintf('%0.1f',pr.x0);
r.y_0.String = sprintf('%0.1f',pr.y0);
r.cols.String = sprintf('%d',pr.cols);
r.rows.String = sprintf('%d',pr.rows);

%view 1 frame
r.theta_x1.String = sprintf('%0.1f',rg{1}(1));
r.theta_y1.String = sprintf('%0.1f',rg{1}(2));
r.theta_z1.String = sprintf('%0.1f',rg{1}(3));
r.tx1.String = sprintf('%0.1f',t{1}(1));
r.ty1.String = sprintf('%0.1f',t{1}(2));
r.tz1.String = sprintf('%0.1f',t{1}(3));

%view 2 frame
r.theta_x2.String = sprintf('%0.1f',rg{2}(1));
r.theta_y2.String = sprintf('%0.1f',rg{2}(2));
r.theta_z2.String = sprintf('%0.1f',rg{2}(3));
r.tx2.String = sprintf('%0.1f',t{2}(1));
r.ty2.String = sprintf('%0.1f',t{2}(2));
r.tz2.String = sprintf('%0.1f',t{2}(3));

%error frame
r.delta_cod.String = sprintf('%d',pr.delta);
r.kappa1_min.String = sprintf('%0.1e',pr.kappa1(1));
r.kappa1_max.String = sprintf('%0.1e',pr.kappa1(2));
r.kappa2_min.String = sprintf('%0.1e',pr.kappa2(1));
r.kappa2_max.String = sprintf('%0.1e',pr.kappa2(2));

%==================================================
function fun_exit(~,~,main_fig)
%close the interface main figure
close(main_fig);

%==================================================
function plot_error(e_mean,e_std,kappa1,kappa2)
global r
figure('name','reconstruction_error','numbertitle','off');
grid on;
hold on;
surface(kappa1,kappa2,e_mean(:,:,1),'EdgeColor','b','FaceColor','none');
surface(kappa1,kappa2,e_mean(:,:,2),'EdgeColor','r','FaceColor','none');
surface(kappa1,kappa2,e_mean(:,:,1)-e_std(:,:,1),'EdgeColor',[0.5 0.5 1], ...
    'FaceColor','none','LineStyle','--','HandleVisibility','off');
surface(kappa1,kappa2,e_mean(:,:,1)+e_std(:,:,1),'EdgeColor',[0.5 0.5 1], ...
    'FaceColor','none','LineStyle','--','HandleVisibility','off');
surface(kappa1,kappa2,e_mean(:,:,2)-e_std(:,:,2),'EdgeColor',[1 0.5 0.5], ...
    'FaceColor','none','LineStyle','--','HandleVisibility','off');
surface(kappa1,kappa2,e_mean(:,:,2)+e_std(:,:,2),'EdgeColor',[1 0.5 0.5], ...
    'FaceColor','none','LineStyle','--','HandleVisibility','off');
set(gca,'fontsize',12);
if( r.units_m.Value==1 )
    units = '(m)';
else
    units = '(cm)';
end
xlabel('$\kappa_1$ (pixels$^{-2}$)','interpreter','latex','fontsize',16);
ylabel('$\kappa_2$ (pixels$^{-4}$)','interpreter','latex','fontsize',16);
zlabel(strcat('$\epsilon ',units,'$'),'interpreter','latex','fontsize',16);
legend('case A: principal point = CoD','case B: principal point fixed');
view(-30,30);

%========================================================================
function [X,Y,Z] = compute_world_points()
global r
num_points = str2double(r.num_points.String);
w_axis = [str2double(r.Xmin.String), str2double(r.Xmax.String), str2double(r.Ymin.String), ...
    str2double(r.Ymax.String), str2double(r.Zmin.String), str2double(r.Zmax.String)];
angles = [str2double(r.alpha_x.String), str2double(r.alpha_y.String), str2double(r.alpha_z.String)];
rseed = r.rand_seed.Value;
rng(rseed);
Xw = w_axis(1) + (w_axis(2) - w_axis(1))*rand(num_points,1);
Yw = w_axis(3) + (w_axis(4) - w_axis(3))*rand(num_points,1);
Zw = w_axis(5) + (w_axis(6) - w_axis(5))*rand(num_points,1);
R_mat = rotation_matrix(angles(1),angles(2),angles(3));
X = [Xw Yw Zw]*R_mat(1,:)';
Y = [Xw Yw Zw]*R_mat(2,:)';
Z = [Xw Yw Zw]*R_mat(3,:)';

%=================================================
function [rg,t] = camera_poses()
%camera view angles and translation
global r
rg{1} = [str2double(r.theta_x1.String); str2double(r.theta_y1.String); str2double(r.theta_z1.String)];
t{1} = [str2double(r.tx1.String); str2double(r.ty1.String); str2double(r.tz1.String)];
rg{2} = [str2double(r.theta_x2.String); str2double(r.theta_y2.String); str2double(r.theta_z2.String)];
t{2} = [str2double(r.tx2.String); str2double(r.ty2.String); str2double(r.tz2.String)];

%=====================================================================================
function epsilon = reconstruction_error(ptd1,ptd2,x0,y0,delta,theta,kappa1,kappa2,P,cp,ft,ind)
%compute the reconstruction error for the cases A and B
%case A: the principal point is the CoD
%case B: the principal point is fixed at (x0,y0)
n = length(theta);
epsilon = zeros(n,2);
for i=1:n
    %CoD displacement
    xd = x0 + delta*cosd(theta(i));
    yd = y0 + delta*sind(theta(i));
    ptu1 = remove_distortion(ptd1,xd,yd,kappa1,kappa2);
    ptu2 = remove_distortion(ptd2,xd,yd,kappa1,kappa2);
    F = estimateFundamentalMatrix(ptu1(ind,:),ptu2(ind,:),'Method','Norm8Point');
    for j=1:2
        switch j  %[xp, yp] principal point for cases A (1) and B (2)
            case 1;  xp = xd;  yp = yd;
            case 2;  xp = x0;  yp = y0;
        end
        %recover camera poses
        R{1} = [1 0 0; 0 1 0; 0 0 1];
        t{1} = [0; 0; 0];
        K = [ft, 0, xp; 0, ft, yp; 0, 0, 1];
        camPars = cameraParameters('IntrinsicMatrix',K');
        [Rrel,trel] = relativeCameraPose(F,camPars,ptu1(ind,:),ptu2(ind,:));
        R{2} = Rrel;
        t{2} = -Rrel*trel';
        Q = point_reconstruction(ptu1,ptu2,K,R,t);
        scale = (P(cp,1)*Q(cp,1) + P(cp,2)*Q(cp,2) + P(cp,3)*Q(cp,3)) / ...
            (Q(cp,1)^2 + Q(cp,2)^2 + Q(cp,3)^2);
        Q = scale * Q;
        epsilon(i,j) = sqrt(mean((Q(ind,1) - P(ind,1)).^2 + ...
            (Q(ind,2) - P(ind,2)).^2 + (Q(ind,3) - P(ind,3)).^2));
    end
end

%=============================================================
function Q = point_reconstruction(ptu1,ptu2,K,R,t)
xp = K(1,3);
yp = K(2,3);
ft = K(1,1);
x1 = (ptu1(:,1) - xp)/ft;
y1 = (ptu1(:,2) - yp)/ft;
x2 = (ptu2(:,1) - xp)/ft;
y2 = (ptu2(:,2) - yp)/ft;
n = length(x1);
Q = zeros(n,3);
for i=1:n
    if( not(isnan(x1(i))) && not(isnan(x2(i))) )
        M = [R{1}(1,1)-x1(i)*R{1}(3,1), R{1}(1,2)-x1(i)*R{1}(3,2), R{1}(1,3)-x1(i)*R{1}(3,3); ...
            R{1}(2,1)-y1(i)*R{1}(3,1), R{1}(2,2)-y1(i)*R{1}(3,2), R{1}(2,3)-y1(i)*R{1}(3,3); ...
            R{2}(1,1)-x2(i)*R{2}(3,1), R{2}(1,2)-x2(i)*R{2}(3,2), R{2}(1,3)-x2(i)*R{2}(3,3); ...
            R{2}(2,1)-y2(i)*R{2}(3,1), R{2}(2,2)-y2(i)*R{2}(3,2), R{2}(2,3)-y2(i)*R{2}(3,3)];
        N = [x1(i)*t{1}(3) - t{1}(1); y1(i)*t{1}(3) - t{1}(2); ...
            x2(i)*t{2}(3) - t{2}(1); y2(i)*t{2}(3) - t{2}(2)];
        s = (M'*M)\(M'*N);
        Q(i,:) = s';
    end
end

%============================================================
function R_mat = rotation_matrix(thx,thy,thz)
Rx = [1 0 0; 0 cosd(thx) -sind(thx); 0 sind(thx) cosd(thx)];
Ry = [cosd(thy) 0 sind(thy); 0 1 0; -sind(thy) 0 cosd(thy)];
Rz = [cosd(thz) -sind(thz) 0; sind(thz) cosd(thz) 0; 0 0 1];
R_mat = Rx * Ry * Rz;

%========================================================
function [x,y,indz] = ideal_point_projection(X,Y,Z,rg,t)
global r
R_mat = rotation_matrix(rg(1),rg(2),rg(3));
n = length(X);
P = R_mat*[X Y Z]' + repmat(t,1,n);
f = str2double(r.focal_length.String);
x0 = str2double(r.x_0.String);
y0 = str2double(r.y_0.String);
x = f * P(1,:)'./P(3,:)' + x0;
y = f * P(2,:)'./P(3,:)' + y0;
indz = (P(3,:)'>0);

%=============================================================
function ptd = apply_distortion(ptu,kappa1,kappa2,x0,y0)
r2 = (ptu(:,1)-x0).^2 + (ptu(:,2)-y0).^2;
f_dist = 1 + kappa1.*r2 + kappa2.*r2.^2;
xd = (ptu(:,1)-x0).*f_dist + x0;
yd = (ptu(:,2)-y0).*f_dist + y0;
ptd = [xd, yd];

%============================================================
function ptu = remove_distortion(ptd,xd,yd,kappa1,kappa2)
%remove distortion using Newton-Raphson method
%f(r) = r + kappa1*r^3 + kappa2*r^5 - rd = 0
%f'(r) = 1 + 3*kappa1*r^2 + 5*kappa2*r^4
rd = sqrt((ptd(:,1)-xd).^2 + (ptd(:,2)-yd).^2);
r = rd;
ftol = 1e-3;
while(1)
    f_r = r + kappa1*r.^3 + kappa2*r.^5 - rd;
    fd_r = 1 + 3*kappa1*r.^2 + 5*kappa2*r.^4;
    r = r - f_r./fd_r;
    if( max(abs(f_r))<=ftol )
        break;
    end
end
xu = (ptd(:,1)-xd).*r./rd + xd;
yu = (ptd(:,2)-yd).*r./rd + yd;
ptu = [xu, yu];
ind = (abs(rd)<1e-6);
ptu(ind,:) = ptd(ind,:);

%=========================================================
function pr = load_intrinsic_parameters()
pr.ft = 1750;
pr.x0 = 2012.6;
pr.y0 = 1544.6;
pr.cols = 4000;
pr.rows = 3000;
pr.delta = 10;
pr.kappa1 = [-6e-8, 0];
pr.kappa2 = [0, 1e-14];

%=====================================================
function [rg,t] = load_extrinsic_parameters()
rg{1} = [0; 0; 0];
t{1} = [0; 0; 0];
rg{2} = [-1; 2; -5];
t{2} = [-8; -4; -5];

%==================================================
function wp = load_world_parameters()
wp.num_points = 1000;
wp.axis = [-90 90 -70 70 95 105];
wp.angles = [0, 0, 0];
wp.units = 1;
wp.rseed = 1;