%============ check station and monitor settings==============
station_file='station_japan.txt';
station_file_out='tmp_orgsetting.txt';
expected_moni_center=[135.5 35.1];
delete_stn_id=[];
%{
station_file: the station file you plan to use, but you need to select proper
              stations to adapt the monitoring target area.
station_file_out: the selected stations and monitoring center, latter, you
                  may use this file to construct data_info.json file.
expected_moni_center: set the longitude and latitude of the monitoring
                      center you expected, and you should adjust it again
                      until the mean xy coorinates of the selected stations 
                     ("mean center") are close to the "expected center".
delete_stn_id: you can delete some stations to adapt the monitoring setting
               according to your plots.
explain plots: the black and red rectangles are "Expected range" and "Mean range"
               for monitoring areas. "Expected range" is centered at the
               "Expected center", while "Mean range" is centered at the
               "Mean center". Blue triangles are the stations not used, and
               Red triangles are the stations selected for monitoring. The
               stations you delected in "delete_stn_id" and out of the
               "Expected range" are displayed in blue color. The red
               triangles are you final monitoring stations.
               The inside rectangles are the monitor ranges of earthquake
               location, and the outside rectangles are ranges of station
               locations.

%}
%==============================================================
[stns, stnnams] = readStationFile(station_file);
center=expected_moni_center;

%tmp=lonlat2km(center(1),center(2),center(1)+1,center(2));
[a,b,c]=distaz(center(1),center(2),center(1)+1,center(2));tmp=a*111.19;
scale=[tmp,111.19];scale1=scale;
stn_rg=[(center(1)*scale(1)-41)/scale(1),(center(1)*scale(1)+41)/scale(1), ...
        (center(2)*scale(2)-50)/scale(2),(center(2)*scale(2)+50)/scale(2)];
src_rg=[(center(1)*scale(1)-25)/scale(1),(center(1)*scale(1)+25)/scale(1), ...
        (center(2)*scale(2)-50)/scale(2),(center(2)*scale(2)+50)/scale(2)];
stn_rg1=stn_rg;src_rg1=src_rg;
figure()
center_plot=plot(center(1),center(2),'kp','MarkerFaceColor', 'black');hold on

stns_plot=plot(stns(:,1),stns(:,2),'b^','MarkerFaceColor', 'blue');
for i=1:length(stns(:,1))
    text(stns(i,1),stns(i,2), num2str(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
stns(delete_stn_id,:)=-999;
for i=1:length(stns(:,1))
    [a,b,c]=distaz(center(1),center(2),center1(1),center1(2));dist(i)=a*111.19;
end
tmp=sort(dist);tmp1=tmp(12);
for i=1:length(stns(:,1))
    if dist(i)>tmp1+0.1E-23
        stns(i,:)=-999;
    end
end

[stn1,stn1_id]=inbox(stns,stn_rg);
stn1_plot=plot(stn1(:,1),stn1(:,2),'r^','MarkerFaceColor', 'red');
stn_rg_plot0=patch([stn_rg(1),stn_rg(1),stn_rg(2),stn_rg(2)],[stn_rg(3),stn_rg(4),stn_rg(4),stn_rg(3)], 'k', 'LineWidth', 2, 'FaceColor', 'none');
patch([src_rg(1),src_rg(1),src_rg(2),src_rg(2)],[src_rg(3),src_rg(4),src_rg(4),src_rg(3)], 'k', 'LineWidth', 2, 'FaceColor', 'none')
center1=[mean(stn1(:,1)),mean(stn1(:,2))];
%tmp=lonlat2km(center1(1),center1(2),center1(1)+1,center1(2));
[a,b,c]=distaz(center1(1),center1(2),center1(1)+1,center1(2));tmp=a*111.19;
scale=[tmp,111.19];
stn_rg=[(center1(1)*scale(1)-41)/scale(1),(center1(1)*scale(1)+41)/scale(1), ...
        (center1(2)*scale(2)-50)/scale(2),(center1(2)*scale(2)+50)/scale(2)];
src_rg=[(center1(1)*scale(1)-25)/scale(1),(center1(1)*scale(1)+25)/scale(1), ...
        (center1(2)*scale(2)-50)/scale(2),(center1(2)*scale(2)+50)/scale(2)];
stn_rg_plot=patch([stn_rg(1),stn_rg(1),stn_rg(2),stn_rg(2)],[stn_rg(3),stn_rg(4),stn_rg(4),stn_rg(3)], 'r', 'LineWidth', 1, 'FaceColor', 'none','EdgeColor','r');
patch([src_rg(1),src_rg(1),src_rg(2),src_rg(2)],[src_rg(3),src_rg(4),src_rg(4),src_rg(3)], 'r', 'LineWidth', 1, 'FaceColor', 'none','EdgeColor','r')
%plot(center(1),center(2),'k.')
center_plot1=plot(center1(1),center1(2),'rp','MarkerFaceColor', 'red');
legend([center_plot,center_plot1,stns_plot,stn1_plot,stn_rg_plot0,stn_rg_plot],'Expected center','Mean center', ...
'Not used sta.','Selected sta.','Expected range','Mean range','Location','EastOutside')
set(gcf,'position',[100 100 800 600]);
pos=axis;
%center_err=lonlat2km(center(1),center(2),center1(1),center1(2));
[center_err,a,b]=distaz(center(1),center(2),center1(1),center1(2));center_err=center_err*111.19;
text(pos(2)+0.05,pos(4)-0.1,{'The diffrence between the','"expected center" ','and "mean center" is ', [num2str(center_err) ' km']}, ...
      'fontsize',12)
fid = fopen(station_file_out, 'w');
%fprintf(fid, '%s %f %f %s %f %f %f %f\n', 'stn_mean_moni_center:',center1(1),center1(2),'expected_moni_center:',center(1),center(2),scale1(1),scale1(2));
fprintf(fid, '%s %f %f %f %f\n', 'expected_moni_center:',center(1),center(2),scale1(1),scale1(2));
fprintf(fid, '%f %f %f %f\n', stn_rg1(1),stn_rg1(2),stn_rg1(3),stn_rg1(4));
fprintf(fid, '%f %f %f %f\n', src_rg1(1),src_rg1(2),src_rg1(3),src_rg1(4));
fprintf(fid, '%s %f %f %f %f\n', 'stn_mean_moni_center:',center1(1),center1(2),scale(1),scale(2));
fprintf(fid, '%f %f %f %f\n', stn_rg(1),stn_rg(2),stn_rg(3),stn_rg(4));
fprintf(fid, '%f %f %f %f\n', src_rg(1),src_rg(2),src_rg(3),src_rg(4));
for i = 1:length(stn1_id)
    fprintf(fid, '%d %f %f %s\n', stn1_id(i), stn1(i,1),stn1(i,2),stnnams{stn1_id(i)});
end
fclose(fid);
