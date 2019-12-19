clc; clear all; close all;

match_threshold = 3; % 2 pixels

ADTdir=['../../allsat_sla/'];

jetfiles = dir('../jet_long_revise*.mat');
lat_jet = ncread(['../jet_detection_Southern_Ocean_1993','.nc'],'lat');
lon_jet = ncread(['../jet_detection_Southern_Ocean_1993','.nc'],'lon');

Aeddyfiles = dir('../../allsat_sla/eddy_find/all2/anticyclonic/*.mat');
Ceddyfiles = dir('../../allsat_sla/eddy_find/all2/cyclonic/*.mat');

load('G:\LIURAN_SCIENCE\allsat_sla\eddy_find\ACC_SLA_eddyTracks.mat');

day_index = [];h = 1;
day_index(1,1) = str2num(Aeddyfiles(1).name(end-11:end-4));
day_index(1,2) = str2num(Aeddyfiles(1).name(end-11:end-8));
day_index(1,3) = 1;
for d= 2:length(Aeddyfiles)
    day_index(d,1) = str2num(Aeddyfiles(d).name(end-11:end-4));
    day_index(d,2) = str2num(Aeddyfiles(d).name(end-11:end-8));
    
    if day_index(d,2) == day_index(d-1,2)
        h = h+1;
    else
        h = 1;
    end
    day_index(d,3) = h;
end

% % select eastward eddy track
% tic
% 
% for i = 1:length(atracks)
%     track_cache = atracks{i};
%     lon_cache = track_cache(:,2);
%     lon_cache(lon_cache<0) = lon_cache(lon_cache<0) +360;
%     diff_cache = diff(lon_cache);
%     it_cross360e = find(diff_cache < -350);
%     it_cross360w = find(diff_cache > 350);  
%     if length(it_cross360e>0)
%         diff_cache(it_cross360e) = diff_cache(it_cross360e) + 360;
% %         lon_cache(it_cross360+1:end) = lon_cache(it_cross360+1:end) + 360;
%     end
%     if length(it_cross360w>0)
%         diff_cache(it_cross360w) = diff_cache(it_cross360w) - 360;
%     end
%     alife(i) = length(lon_cache);
%     adlon(i) = sum(diff_cache);
% end
% 
% it_eastward = find(adlon >= 0 & alife>=28);
% atrack_eastward = atracks(it_eastward);
% 
% toc
% %disp(['运行时间: ',num2str(toc)]);
% tic
% 
% for i = 1:length(ctracks)
%     track_cache = ctracks{i};
%     lon_cache = track_cache(:,2);
%     lon_cache(lon_cache<0) = lon_cache(lon_cache<0) +360;
%     diff_cache = diff(lon_cache);
%     it_cross360e = find(diff_cache < -350);
%     it_cross360w = find(diff_cache > 350);
%     if length(it_cross360e>0)
%         diff_cache(it_cross360e) = diff_cache(it_cross360e) + 360;
%         %         lon_cache(it_cross360+1:end) = lon_cache(it_cross360+1:end) + 360;
%     end
%     if length(it_cross360w>0)
%         diff_cache(it_cross360w) = diff_cache(it_cross360w) - 360;
%     end
%     clife(i) = length(lon_cache);
%     cdlon(i) = sum(diff_cache);
% end
% 
% it_eastward = find(cdlon >= 0 & clife>=28);
% ctrack_eastward = ctracks(it_eastward);
% toc

% dd = [30:10:300];
% dd_center = dd(1:end-1) + (dd(2:end)-dd(1:end-1))/2;
% ll = [-10:0.5:10];
% ll_center = ll(1:end-1) + (ll(2:end)-ll(1:end-1))/2;
% dlon_life_box = zeros(length(dd_center),length(ll_center));
% for i = 1:length(dd_center)
%     for j = 1:length(ll_center)
%         it_abox = find(alife>dd_center(i)-5 & alife<=dd_center(i)+5 &...
%             adlon > ll(j)-0.25 & adlon<=ll(j)+0.25 );
%         it_cbox = find(clife>dd_center(i)-5 & clife<=dd_center(i)+5 &...
%             cdlon > ll(j)-0.25 & cdlon<=ll(j)+0.25 );
%         dlon_life_box(i,j,1) = length(it_abox);;
%         dlon_life_box(i,j,2) = length(it_cbox);;
%         dlon_life_box(i,j,3) = length(it_abox) + length(it_cbox);;
%     end
% end
% % dlon_life_box(1,1) = length( find(alife>0 & alife<=dd(1) &...
% %                                 adlon <= ll(1)) );
% % dlon_life_box(1,1) = length( find(alife>0 & alife<=dd(1) &...
% %                                 adlon <= ll(1)) );
% dlon_life_box(dlon_life_box==0) = nan;

% figure
% subplot(2,2,1)
% contourf(dd_center,ll_center,log10(squeeze(dlon_life_box(:,:,1))'),[0:0.2:4])
% caxis([0, 4])
% colorbar
% grid on;
% hold on;
% plot([min(dd_center),max(dd_center)],[0,0],'w-','linewidth',1.5)
% subplot(2,2,2)
% contourf(dd_center,ll_center,log10(squeeze(dlon_life_box(:,:,2))'),[0:0.2:4])
% caxis([0, 4])
% colorbar
% grid on;
% hold on;
% plot([min(dd_center),max(dd_center)],[0,0],'w-','linewidth',1.5)
% subplot(2,2,3)
% contourf(dd_center,ll_center,log10(squeeze(dlon_life_box(:,:,3))'),[0:0.2:4])
% caxis([0, 4])
% colorbar
% grid on;
% hold on;
% plot([min(dd_center),max(dd_center)],[0,0],'w-','linewidth',1.5)
% print('-dpng','-r600','dlon_life_contourf')

%% eddy track to 2dim

% atracks_assign2dim = [];
% for t = 1:length(atrack_eastward)
%     track_cache = atrack_eastward{t};
%     track_cache(:,6) = 1;   %************ cyc information ************
%     track_cache(:,7) = t;   %************ track index information ************
%     atracks_assign2dim = [atracks_assign2dim; track_cache];
% end
% ctracks_assign2dim = [];
% for t = 1:length(ctrack_eastward)
%     track_cache = ctrack_eastward{t};
%     track_cache(:,6) = 1;   % cyc information
%     track_cache(:,7) = t;   % track index information
%     ctracks_assign2dim = [ctracks_assign2dim; track_cache];
% end
% save('tracks_28day_long2dim.mat','atracks_assign2dim','ctracks_assign2dim');

%% jet to 2dim
% for y = 1:24
%     load([jetfiles(y).folder,'/',jetfiles(y).name]);
% 
%     for d = 1:length(jet_long_year)
%         jet_day = [];
%         jet_cache = jet_long_year{d};
%         start = 1;
%         for i = 1:length(jet_cache)
%             cache = jet_cache{i};
%             len_cache = size(cache,1);
%             jet_day(start:start+len_cache-1,2:3) = cache;
%             jet_day(start:start+len_cache-1,1) = i;
%             start = start+len_cache;
%         end
%         it = find(day_index(:,2)==1992+y & day_index(:,3)==d);
%         if length(it)==0
%             break;
%         end
%         save(['./jet_strings_2dim/',num2str(day_index(it,1)),'.mat'],'jet_day');
%     end
% end

%% match
load('tracks_eastward28day_long2dim.mat')
match_day_aeddymask{size(atracks_assign2dim,1),9}=[]; match_day_ceddymask{size(ctracks_assign2dim,1),9}=[];%*****1,eddy_track_index,2,eddy_lat,3,eddy_lon,4,amp,5,cyc,6,day,7,eddy index in one day
%***********************8,match jet_index cell,9,match jet_points cell,10, match jet_points vel cell
match_day_aeddymask_jet{size(atracks_assign2dim,1),12}=[]; match_day_ceddymask_jet{size(ctracks_assign2dim,1),12}=[];%*****1,eddy_track_index,2,eddy_lat,3,eddy_lon,4,amp,5,cyc,6,day,7,eddy index in one day

match_day_jetmask = []; %*****
%*****************************1,Aeddyindex in tracks2dim(none=0),2,-1*dist_jet2center(在涡旋内部)/0（涡旋边界上）/value(涡旋外多少km),
%*****************************3,Aeddy amp, 4, Aarea,5 lat,6,lon,7,cyc

for d = 1:length(Aeddyfiles)
    d
    load(['./jet_strings_2dim/',num2str(day_index(d,1)),'.mat']);
    adtfilesall = dir([ADTdir,num2str(day_index(d,2)),'/*.nc']);
    u = ncread([adtfilesall(day_index(d,3)).folder,'/',adtfilesall(day_index(d,3)).name],'ugos',[1,1,1],[1440,240,1]);
    v = ncread([adtfilesall(day_index(d,3)).folder,'/',adtfilesall(day_index(d,3)).name],'vgos',[1,1,1],[1440,240,1]);
    
    clear match_day_jetmask;
    match_day_jetmask{size(jet_day,1),14}=[];
    jet_oneday = jet_day;
    jet_oneday(:,2) = lon_jet(jet_day(:,2));
    jet_oneday(:,3) = lat_jet(jet_day(:,3));
    
    % anticyclonic
    it_Aeddy_readymatch = find(atracks_assign2dim(:,3)==d);
    Aeddy_readymatch = atracks_assign2dim(it_Aeddy_readymatch,:);
    
    Aeddys = load([Aeddyfiles(d).folder,'/',Aeddyfiles(d).name]);
    Aeddys = Aeddys.eddies;
    
    for i = 1:size(Aeddy_readymatch,1)
        if Aeddy_readymatch(i,5) == 0  % real eddy
            eddy_cache = Aeddys(Aeddy_readymatch(i,4));
            eddy_center_lon = eddy_cache.Lon;
            if eddy_center_lon<0
                eddy_center_lon = eddy_center_lon +360;
            end
            eddy_center_lat = eddy_cache.Lat;
            
            eddypoint = eddy_cache.Stats.PixelIdxList;
            eddypoint_lonIndex=zeros(length(eddypoint),1);
            eddypoint_latIndex=zeros(length(eddypoint),1);
            for j=1:length(eddypoint)
                eddypoint_lonIndex(j,1)=floor(eddypoint(j)/240)+1; %lon
                eddypoint_latIndex(j,1)=mod(eddypoint(j),240);   %lat
            end
            eddypoint_lon = lon_jet(eddypoint_lonIndex);  %***** 0-360 ******
            eddypoint_lat = lat_jet(eddypoint_latIndex);
            
            if abs(max(eddypoint_lon)-min(eddypoint_lon))>300
                eddypoint_lon_360 = eddypoint_lon;
                eddypoint_lon_360(eddypoint_lon_360>300) = eddypoint_lon_360(eddypoint_lon_360>300)-360;
                shp = alphaShape(eddypoint_lon_360,eddypoint_lat);
                bf = boundaryFacets(shp);    
                plot(shp)
            else
                shp = alphaShape(eddypoint_lon,eddypoint_lat);
                bf = boundaryFacets(shp);
            end
            if length(bf)>0
                edge_lon = eddypoint_lon(bf(:,1));
                edge_lat = eddypoint_lat(bf(:,1));
                edge_x = eddypoint_lonIndex(bf(:,1));
                edge_y = eddypoint_latIndex(bf(:,1));
                
                left_boundry = min(edge_lon)-1;
                right_boundry = max(edge_lon)+1;
                if abs(left_boundry-right_boundry)>300
                    left_boundry = min(edge_lon(find(edge_lon>300)))-1;
                    right_boundry = max(edge_lon(find(edge_lon<50)))+1;
                    it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                        jet_oneday(:,3) <= max(edge_lat)+1 & ...
                        (jet_oneday(:,2) <= right_boundry | ...
                        jet_oneday(:,2) >= left_boundry));                        
                else
                    if left_boundry<0
                        left_boundry = min(edge_lon)-1+360;
                        right_boundry = max(edge_lon)+1;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            (jet_oneday(:,2) <= right_boundry | ...
                            jet_oneday(:,2) >= left_boundry));   
                    elseif right_boundry>360
                        left_boundry = min(edge_lon)-1;
                        right_boundry = max(edge_lon)+1-360;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            (jet_oneday(:,2) <= right_boundry | ...
                            jet_oneday(:,2) >= left_boundry));                          
                    else
                        left_boundry = min(edge_lon)-1;
                        right_boundry = max(edge_lon)+1;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            jet_oneday(:,2) <= right_boundry & ...
                            jet_oneday(:,2) >= left_boundry);                          
                    end
                end
                
                in = [];on=[];
                if length(it_jet_readymatch)>0
                    for j = 1:length(it_jet_readymatch)
                        lon_jetpoint_cache = jet_oneday(it_jet_readymatch(j),2);
                        lat_jetpoint_cache = jet_oneday(it_jet_readymatch(j),3);
                        if abs(left_boundry-right_boundry)>300
                            lon_jetpoint_cache(lon_jetpoint_cache>300)=lon_jetpoint_cache(lon_jetpoint_cache>300)-360;
                            edge_lon360 = edge_lon;
                            edge_lon360(edge_lon360>300) = edge_lon360(edge_lon360>300)-360;
                            [in(j), on(j) ]= inpolygon(lon_jetpoint_cache,lat_jetpoint_cache,...
                                                       edge_lon360,edge_lat);
                        else
                            [in(j), on(j) ]= inpolygon(lon_jetpoint_cache,lat_jetpoint_cache,...
                                                       edge_lon,edge_lat);
                        end
                    end
                    
                    it_on = find(on == 1);
                    if length(it_on)>0
%                         for n_it_on = 1:length(it_on)
%                             if length(match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),1}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),1} = Aeddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),2} = 0;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),3} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),4} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),5} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),6} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),7} = eddy_cache.Cyc;
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),1} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),1},Aeddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),2} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),2},0];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),3} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),3},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),4} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),4},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),5} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),5},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),6} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),6},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),7} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),7},eddy_cache.Cyc];
%                             end
%                             
%                         end
                        %match_day_aeddymask{it_Aeddy_readymatch(i),1} = it_jet_readymatch(it_on);                % jet numberth in day
                        %match_day_aeddymask{it_Aeddy_readymatch(i),2} = jet_oneday(it_jet_readymatch(it_on),1); % jet index
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),1} = single(jet_oneday(it_jet_readymatch(it_on),2)); % jet lon
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),2} = single(jet_oneday(it_jet_readymatch(it_on),3)); % jet lat
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),3} = single(diag(u(jet_day(it_jet_readymatch(it_on),2),jet_day(it_jet_readymatch(it_on),3)))); % jet u
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),4} = single(diag(v(jet_day(it_jet_readymatch(it_on),2),jet_day(it_jet_readymatch(it_on),3)))); % jet v                      
                        %match_day_aeddymask{it_Aeddy_readymatch(i),3} = 0;
                        
                    end
                    
                    it_in = find(in == 1 & on == 0);
                    if length(it_in)>0
                        dist_jet2center = [];
                        for n_in = 1:length(it_in)
                            dist_jet2center(n_in) = -1 * sw_dist([eddy_center_lat,jet_oneday(it_jet_readymatch(it_in(n_in)),3)],...
                                [eddy_center_lon,jet_oneday(it_jet_readymatch(it_in(n_in)),2)],'km');
                        end
%                         for n_it_in = 1:length(it_in)
%                             if length(match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),1}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),1} = Aeddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),2} = dist_jet2center(n_it_in);
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),3} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),4} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),5} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),6} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),7} = eddy_cache.Cyc;
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),1} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),1},Aeddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),2} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),2},dist_jet2center(n_it_in);];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),3} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),3},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),4} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),4},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),5} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),5},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),6} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),6},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),7} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),7},eddy_cache.Cyc];
%                             end
%                         end
%                         match_day_aeddymask{it_Aeddy_readymatch(i),3} = it_jet_readymatch(it_in);
%                         match_day_aeddymask{it_Aeddy_readymatch(i),4} = -1*dist_jet2center;
                        %match_day_aeddymask{it_Aeddy_readymatch(i),4} = it_jet_readymatch(it_in);                % jet numberth in day
                        %match_day_aeddymask{it_Aeddy_readymatch(i),5} = jet_oneday(it_jet_readymatch(it_in),1); % jet index
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),5} = single(jet_oneday(it_jet_readymatch(it_in),2)); % jet lon
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),6} = single(jet_oneday(it_jet_readymatch(it_in),3)); % jet lat
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),7} = single(diag(u(jet_day(it_jet_readymatch(it_in),2),jet_day(it_jet_readymatch(it_in),3)))); % jet u
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),8} = single(diag(v(jet_day(it_jet_readymatch(it_in),2),jet_day(it_jet_readymatch(it_in),3)))); % jet v                      
                        %match_day_aeddymask{it_Aeddy_readymatch(i),6} = -1*dist_jet2center;
                    end
                    
                    it_out = find(in == 0);
                    if length(it_out)>0
                        dist_jet2edge = [];dgrid_jet2edge = [];mindegree_index=[];
                        for n_out = 1:length(it_out)
                            dist_cacahe = []; dgrid_cache = [];
                            for n_edge = 1:length(edge_lat)
                                dist_cacahe(n_edge) = sw_dist([edge_lat(n_edge),jet_oneday(it_jet_readymatch(it_out(n_out)),3)],...
                                    [edge_lon(n_edge),jet_oneday(it_jet_readymatch(it_out(n_out)),2)],'km');
                                
                                dgrid_cache(n_edge) = ((edge_y(n_edge)-jet_day(it_jet_readymatch(it_out(n_out)),3))^2 +...
                                    (edge_x(n_edge)-jet_day(it_jet_readymatch(it_out(n_out)),2))^2)^0.5;
                                
                            end
                            dist_jet2edge(n_out) = min(dist_cacahe);
                            [dgrid_jet2edge(n_out),mindegree_index(n_out)]= min(dgrid_cache);
                            
                        end
                        it_outmatch = find(dgrid_jet2edge <= match_threshold); % 2 pixels for eddy-jet match
%                         for n_it_out = 1:length(it_outmatch)
%                             if length(match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),1}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),1} = Aeddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),2} = dist_jet2edge(it_outmatch(n_it_out));
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),3} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),4} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),5} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),6} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),7} = eddy_cache.Cyc;
%                                 
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),1} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),1},Aeddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),2} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),2},dist_jet2edge(it_outmatch(n_it_out))];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),3} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),3},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),4} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),4},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),5} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),5},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),6} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),6},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),7} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),7},eddy_cache.Cyc];
%                             end
%                         end
%                         match_day_aeddymask{it_Aeddy_readymatch(i),5} = it_jet_readymatch(it_out(it_outmatch));
%                         match_day_aeddymask{it_Aeddy_readymatch(i),6} = dist_jet2edge(it_outmatch);
                        %match_day_aeddymask{it_Aeddy_readymatch(i),7} = it_jet_readymatch(it_out(it_outmatch));                % jet numberth in day
                        %match_day_aeddymask{it_Aeddy_readymatch(i),8} = jet_oneday(it_jet_readymatch(it_out(it_outmatch)),1); % jet index
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),9} = single(jet_oneday(it_jet_readymatch(it_out(it_outmatch)),2)); % jet lon
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),10} = single(jet_oneday(it_jet_readymatch(it_out(it_outmatch)),3)); % jet lat
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),11} = single(diag(u(jet_day(it_jet_readymatch(it_out(it_outmatch)),2),jet_day(it_jet_readymatch(it_out(it_outmatch)),3)))); % jet u
                        match_day_aeddymask_jet{it_Aeddy_readymatch(i),12} = single(diag(v(jet_day(it_jet_readymatch(it_out(it_outmatch)),2),jet_day(it_jet_readymatch(it_out(it_outmatch)),3)))); % jet v                      
                        %match_day_aeddymask{it_Aeddy_readymatch(i),9} = dist_jet2edge(it_outmatch);
                    end
                    
                end
            end
            
        end % real eddy
        
    end
 
    
    % cyclonic
    it_Ceddy_readymatch = find(ctracks_assign2dim(:,3)==d);
    Ceddy_readymatch = ctracks_assign2dim(it_Ceddy_readymatch,:);
    
    Ceddys = load([Ceddyfiles(d).folder,'/',Ceddyfiles(d).name]);
    Ceddys = Ceddys.eddies;
    
    for i = 1:size(Ceddy_readymatch,1)
        if Ceddy_readymatch(i,5) == 0  % real eddy
            eddy_cache = Ceddys(Ceddy_readymatch(i,4));
            eddy_center_lon = eddy_cache.Lon;
            if eddy_center_lon<0
                eddy_center_lon = eddy_center_lon +360;
            end
            eddy_center_lat = eddy_cache.Lat;
            
            eddypoint = eddy_cache.Stats.PixelIdxList;
            eddypoint_lonIndex=zeros(length(eddypoint),1);
            eddypoint_latIndex=zeros(length(eddypoint),1);
            for j=1:length(eddypoint)
                eddypoint_lonIndex(j,1)=floor(eddypoint(j)/240)+1; %lon
                eddypoint_latIndex(j,1)=mod(eddypoint(j),240);   %lat
            end
            eddypoint_lon = lon_jet(eddypoint_lonIndex);  %***** 0-360 ******
            eddypoint_lat = lat_jet(eddypoint_latIndex);
            
            if abs(max(eddypoint_lon)-min(eddypoint_lon))>300
                eddypoint_lon_360 = eddypoint_lon;
                eddypoint_lon_360(eddypoint_lon_360>300) = eddypoint_lon_360(eddypoint_lon_360>300)-360;
                shp = alphaShape(eddypoint_lon_360,eddypoint_lat);
                bf = boundaryFacets(shp);    
                plot(shp)
            else
                shp = alphaShape(eddypoint_lon,eddypoint_lat);
                bf = boundaryFacets(shp);
            end
            if length(bf)>0                              %**** if not line shape eddy
                edge_lon = eddypoint_lon(bf(:,1));
                edge_lat = eddypoint_lat(bf(:,1));
                edge_x = eddypoint_lonIndex(bf(:,1));
                edge_y = eddypoint_latIndex(bf(:,1));
                

                left_boundry = min(edge_lon)-1;
                right_boundry = max(edge_lon)+1;
                if abs(left_boundry-right_boundry)>300
                    left_boundry = min(edge_lon(find(edge_lon>300)))-1;
                    right_boundry = max(edge_lon(find(edge_lon<50)))+1;
                    it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                        jet_oneday(:,3) <= max(edge_lat)+1 & ...
                        (jet_oneday(:,2) <= right_boundry | ...
                        jet_oneday(:,2) >= left_boundry));                        
                else
                    if left_boundry<0
                        left_boundry = min(edge_lon)-1+360;
                        right_boundry = max(edge_lon)+1;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            (jet_oneday(:,2) <= right_boundry | ...
                            jet_oneday(:,2) >= left_boundry));   
                    elseif right_boundry>360
                        left_boundry = min(edge_lon)-1;
                        right_boundry = max(edge_lon)+1-360;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            (jet_oneday(:,2) <= right_boundry | ...
                            jet_oneday(:,2) >= left_boundry));                          
                    else
                        left_boundry = min(edge_lon)-1;
                        right_boundry = max(edge_lon)+1;
                        it_jet_readymatch = find(jet_oneday(:,3) >= min(edge_lat)-1 &...
                            jet_oneday(:,3) <= max(edge_lat)+1 & ...
                            jet_oneday(:,2) <= right_boundry & ...
                            jet_oneday(:,2) >= left_boundry);                          
                    end
                end
                in = [];on=[];
                if length(it_jet_readymatch)>0
                    for j = 1:length(it_jet_readymatch)
                        lon_jetpoint_cache = jet_oneday(it_jet_readymatch(j),2);
                        lat_jetpoint_cache = jet_oneday(it_jet_readymatch(j),3);
                        if abs(left_boundry-right_boundry)>300
                            lon_jetpoint_cache(lon_jetpoint_cache>300)=lon_jetpoint_cache(lon_jetpoint_cache>300)-360;
                            edge_lon360 = edge_lon;
                            edge_lon360(edge_lon360>300) = edge_lon360(edge_lon360>300)-360;
                            [in(j), on(j) ]= inpolygon(lon_jetpoint_cache,lat_jetpoint_cache,...
                                                       edge_lon360,edge_lat);
                        else
                            [in(j), on(j) ]= inpolygon(lon_jetpoint_cache,lat_jetpoint_cache,...
                                                       edge_lon,edge_lat);
                        end
                    end
                    
                    it_on = find(on == 1);
                    if length(it_on)>0
%                         for n_it_on = 1:length(it_on)
%                             if length(match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),8}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),8} = Ceddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),9} = 0;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),10} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),11} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),12} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),13} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),14} = eddy_cache.Cyc;
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),8} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),8},Ceddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),9} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),9},0];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),10} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),10},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),11} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),11},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),12} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),12},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),13} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),13},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),14} = [match_day_jetmask{it_jet_readymatch(it_on(n_it_on)),14},eddy_cache.Cyc];
%                             end
%                             
%                         end
%                         match_day_ceddymask{it_Ceddy_readymatch(i),1} = it_jet_readymatch(it_on);
%                         match_day_ceddymask{it_Ceddy_readymatch(i),2} = 0;
                        %match_day_ceddymask{it_Ceddy_readymatch(i),1} = it_jet_readymatch(it_on);                % jet numberth in day
                        %match_day_ceddymask{it_Ceddy_readymatch(i),2} = jet_oneday(it_jet_readymatch(it_on),1); % jet index
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),1} = single(jet_oneday(it_jet_readymatch(it_on),2)); % jet lon
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),2} = single(jet_oneday(it_jet_readymatch(it_on),3)); % jet lat
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),3} = single(diag(u(jet_day(it_jet_readymatch(it_on),2),jet_day(it_jet_readymatch(it_on),3)))); % jet u
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),4} = single(diag(v(jet_day(it_jet_readymatch(it_on),2),jet_day(it_jet_readymatch(it_on),3)))); % jet v                      
                        %match_day_ceddymask{it_Ceddy_readymatch(i),3} = 0;
                    end
                    
                    it_in = find(in == 1 & on == 0);
                    if length(it_in)>0
                        dist_jet2center = [];
                        for n_in = 1:length(it_in)
                            dist_jet2center(n_in) = -1 * sw_dist([eddy_center_lat,jet_oneday(it_jet_readymatch(it_in(n_in)),3)],...
                                [eddy_center_lon,jet_oneday(it_jet_readymatch(it_in(n_in)),2)],'km');
                        end
%                         for n_it_in = 1:length(it_in)
%                             if length(match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),8}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),8} = Ceddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),9} = dist_jet2center(n_it_in);
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),10} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),11} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),12} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),13} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),14} = eddy_cache.Cyc;
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),8} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),8},Ceddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),9} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),9},dist_jet2center(n_it_in);];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),10} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),10},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),11} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),11},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),12} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),12},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),13} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),13},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),14} = [match_day_jetmask{it_jet_readymatch(it_in(n_it_in)),14},eddy_cache.Cyc];
%                             end
%                         end
%                         match_day_ceddymask{it_Ceddy_readymatch(i),3} = it_jet_readymatch(it_in);
%                         match_day_ceddymask{it_Ceddy_readymatch(i),4} = -1*dist_jet2center;
                        %match_day_ceddymask{it_Ceddy_readymatch(i),4} = it_jet_readymatch(it_in);                % jet numberth in day
                        %match_day_ceddymask{it_Ceddy_readymatch(i),5} = jet_oneday(it_jet_readymatch(it_in),1); % jet index
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),5} = single(jet_oneday(it_jet_readymatch(it_in),2)); % jet lon
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),6} = single(jet_oneday(it_jet_readymatch(it_in),3)); % jet lat
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),7} = single(diag(u(jet_day(it_jet_readymatch(it_in),2),jet_day(it_jet_readymatch(it_in),3)))); % jet u
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),8} = single(diag(v(jet_day(it_jet_readymatch(it_in),2),jet_day(it_jet_readymatch(it_in),3)))); % jet v                      
                        %match_day_ceddymask{it_Ceddy_readymatch(i),6} = -1*dist_jet2center;
                    end
                    
                    it_out = find(in == 0);
                    if length(it_out)>0
                        dist_jet2edge = [];dgrid_jet2edge = [];mindegree_index=[];
                        for n_out = 1:length(it_out)
                            dist_cacahe = []; dgrid_cache = [];
                            for n_edge = 1:length(edge_lat)
                                dist_cacahe(n_edge) = sw_dist([edge_lat(n_edge),jet_oneday(it_jet_readymatch(it_out(n_out)),3)],...
                                    [edge_lon(n_edge),jet_oneday(it_jet_readymatch(it_out(n_out)),2)],'km');
                                
                                dgrid_cache(n_edge) = ((edge_y(n_edge)-jet_day(it_jet_readymatch(it_out(n_out)),3))^2 +...
                                    (edge_x(n_edge)-jet_day(it_jet_readymatch(it_out(n_out)),2))^2)^0.5;
                                
                            end
                            dist_jet2edge(n_out) = min(dist_cacahe);
                            [dgrid_jet2edge(n_out),mindegree_index(n_out)]= min(dgrid_cache);
                            
                        end
                        it_outmatch = find(dgrid_jet2edge <= match_threshold);
%                         for n_it_out = 1:length(it_outmatch)
%                             if length(match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),8}) ==0
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),8} = Ceddy_readymatch(i,4);
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),9} = dist_jet2edge(it_outmatch(n_it_out));
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),10} = eddy_cache.Amplitude; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),11} = eddy_cache.SurfaceArea;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),12} = eddy_cache.Lat;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),13} = eddy_cache.Lon;
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),14} = eddy_cache.Cyc;
%                                 
%                             else
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),8} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),8},Ceddy_readymatch(i,4)];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),9} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),9},dist_jet2edge(it_outmatch(n_it_out))];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),10} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),10},eddy_cache.Amplitude]; %cm
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),11} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),11},eddy_cache.SurfaceArea];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),12} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),12},eddy_cache.Lat];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),13} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),13},eddy_cache.Lon];
%                                 match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),14} = [match_day_jetmask{it_jet_readymatch(it_outmatch(n_it_out)),14},eddy_cache.Cyc];
%                             end
%                         end
%                         match_day_ceddymask{it_Ceddy_readymatch(i),5} = it_jet_readymatch(it_out(it_outmatch));
%                         match_day_ceddymask{it_Ceddy_readymatch(i),6} = dist_jet2edge(it_outmatch);
                        %match_day_ceddymask{it_Ceddy_readymatch(i),7} = it_jet_readymatch(it_out(it_outmatch));                % jet numberth in day
                        %match_day_ceddymask{it_Ceddy_readymatch(i),8} = jet_oneday(it_jet_readymatch(it_out(it_outmatch)),1); % jet index
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),9} = single(jet_oneday(it_jet_readymatch(it_out(it_outmatch)),2)); % jet lon
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),10} = single(jet_oneday(it_jet_readymatch(it_out(it_outmatch)),3)); % jet lat
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),11} = single(diag(u(jet_day(it_jet_readymatch(it_out(it_outmatch)),2),jet_day(it_jet_readymatch(it_out(it_outmatch)),3)))); % jet u
                        match_day_ceddymask_jet{it_Ceddy_readymatch(i),12} = single(diag(v(jet_day(it_jet_readymatch(it_out(it_outmatch)),2),jet_day(it_jet_readymatch(it_out(it_outmatch)),3)))); % jet v                      
                        %match_day_ceddymask{it_Ceddy_readymatch(i),9} = dist_jet2edge(it_outmatch);                        
                    end
                end
                
            end
            
        end % real eddy
        
    end

%     save(['./jets_match_eddy/jets_tracks_withEddy_',num2str(day_index(d,1)),'.mat'],'match_day_jetmask','day_index');
    
end

% match_day_ceddymask_jetsinfo = match_day_ceddymask(:,[1,3,4,6,7,9]);
% match_day_aeddymask_jetsinfo = match_day_aeddymask(:,[1,3,4,6,7,9]);
% save(['./eddytracks_withJets.mat'],'match_day_aeddymask_jetsinfo','match_day_ceddymask_jetsinfo','day_index');
% 
% match_day_ceddymask_jetssery = match_day_ceddymask(:,[2,5,8]);
% match_day_aeddymask_jetssery = match_day_aeddymask(:,[2,5,8]);
% save(['./eddytracks_withJets_sery.mat'],'match_day_ceddymask_jetssery','match_day_aeddymask_jetssery','day_index');


match_day_ceddymask_jetslatlon = match_day_ceddymask_jet(:,[1,2,5,6,9,10]);
match_day_aeddymask_jetslatlon = match_day_aeddymask_jet(:,[1,2,5,6,9,10]);
save(['./eddytracks_withJets_latlon_3pixels.mat'],'match_day_ceddymask_jetslatlon','match_day_aeddymask_jetslatlon','day_index');

match_day_ceddymask_jetsVel = match_day_ceddymask_jet(:,[3,4,7,8,11,12]);
match_day_aeddymask_jetsVel = match_day_aeddymask_jet(:,[3,4,7,8,11,12]);
save(['./eddytracks_withJets_vel_3pixels.mat'],'match_day_ceddymask_jetsVel','match_day_aeddymask_jetsVel','day_index');





% match_day_ceddymask_jetssery = match_day_ceddymask(:,[2,5,8]);
% match_day_aeddymask_jetssery = match_day_aeddymask(:,[2,5,8]);
% save(['./eddytracks_withJetsvel.mat'],'match_day_ceddymask_jetssery','match_day_aeddymask_jetssery','day_index');

% match_day_ceddymask_jetsvel = match_day_ceddymask(:,[2,9,16]);
% match_day_aeddymask_jetsvel = match_day_aeddymask(:,[2,9,16]);
% save(['./eddytracks_withJets_sery.mat'],'match_day_ceddymask_jetssery','match_day_aeddymask_jetssery','day_index');


% match_day_ceddymask_on = match_day_ceddymask(:,1:7);
% match_day_ceddymask_in = match_day_ceddymask(:,8:14);
% match_day_ceddymask_out = match_day_ceddymask(:,15:21);
% clear match_day_ceddymask
% match_day_aeddymask_on = match_day_aeddymask(:,1:7);
% match_day_aeddymask_in = match_day_aeddymask(:,8:14);
% match_day_aeddymask_out = match_day_aeddymask(:,15:21);
% clear match_day_aeddymask
% 
% save(['./eddytracks_withJets_on.mat'],'match_day_ceddymask_on','match_day_aeddymask_on','day_index','-v7.3');
% save(['./eddytracks_withJets_in.mat'],'match_day_ceddymask_in','match_day_aeddymask_in','day_index','-v7.3');
% save(['./eddytracks_withJets_out.mat'],'match_day_ceddymask_out','match_day_aeddymask_out','day_index','-v7.3');

% match_day_ceddymask = match_day_ceddymask_on(:,[1,7,8,14,15,21]);
% match_day_aeddymask = match_day_aeddymask_on(:,[1,7,8,14,15,21]);
% save(['./eddytracks_withJets.mat'],'match_day_aeddymask','match_day_ceddymask','day_index');
% 
% 
% 
% match_day_ceddymask_on_index = match_day_ceddymask_on(:,[1:2,7]);
% match_day_aeddymask_on_index = match_day_aeddymask_on(:,[1:2,7]);
% save(['./eddytracks_withJets_on_index.mat'],'match_day_ceddymask_on_index','match_day_aeddymask_on_index','day_index','-v7.3');
% 

% 
% %% test
% adtfiles = dir([ADTdir,num2str(day_index(d,2)),'/*.nc']);
% sla = ncread([adtfiles(d).folder,'/',adtfiles(d).name],'sla',[1,1,1],[1440,240,1]);
% lon1480(1:40) = lon_jet(1401:1440)-360;
% lon1480(41:1480) = lon_jet;
% 
% sla1480(1:40,:) = sla(1401:1440,:);
% sla1480(41:1480,:) = sla;
% eddypoint_lon_360 = eddypoint_lon;
% eddypoint_lon_360(eddypoint_lon_360>300) = eddypoint_lon_360(eddypoint_lon_360>300) -360;
% edge_lon_360 = edge_lon;
% edge_lon_360(edge_lon_360>300) = edge_lon_360(edge_lon_360>300) - 360;
% jet_oneday_360 = jet_oneday;
% jet_oneday_360(jet_oneday_360(:,2)>300,2) = jet_oneday_360(jet_oneday_360(:,2)>300,2)-360 ;
% 
% 
% figure
% % m_proj('stereographic','lat',-90,'long',-70,'radius',60);
% m_proj('miller','lat',[min(eddypoint_lat)-2 max(eddypoint_lat)+2],...
%                 'long',[min(eddypoint_lon_360)-2 max(eddypoint_lon_360)+2]);
% m_grid('xtick',7,'tickdir','out','linest',':');
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% hold on
% m_contour(lon1480,lat_jet,sla1480',25)
% m_plot(eddypoint_lon_360,eddypoint_lat,'ro','markersize',5,'markerfacecolor','r')
% m_plot(edge_lon_360,edge_lat,'o-','markersize',4,'markerfacecolor','w')
% % m_plot(jet_oneday_360(it_jet_readymatch,2),jet_oneday_360(it_jet_readymatch,3),'o-','markersize',2,'markerfacecolor','w')
% m_plot(jet_oneday_360(it_jet_readymatch(it_out(it_outmatch)),2),jet_oneday_360(it_jet_readymatch(it_out(it_outmatch)),3),'>','markerfacecolor','k')
% m_plot(jet_oneday_360(it_jet_readymatch(it_on),2),jet_oneday_360(it_jet_readymatch(it_on),3),'s-','markersize',8,'markerfacecolor','g')
% m_plot(jet_oneday_360(it_jet_readymatch(it_in),2),jet_oneday_360(it_jet_readymatch(it_in),3),'p-','markersize',5,'markeredgecolor','g')
% 
% figure
% % m_proj('stereographic','lat',-90,'long',-70,'radius',60);
% m_proj('miller','lat',[min(eddypoint_lat)-2 max(eddypoint_lat)+2],...
%                 'long',[min(eddypoint_lon_360)-2 max(eddypoint_lon_360)+2]);
% m_grid('xtick',7,'tickdir','out','linest',':');
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% hold on
% m_contour(lat_jet,lat_jet,sla',25)
% m_plot(eddypoint_lon,eddypoint_lat,'ro','markersize',5,'markerfacecolor','r')
% m_plot(edge_lon,edge_lat,'o-','markersize',4,'markerfacecolor','w')
% % m_plot(jet_oneday_360(it_jet_readymatch,2),jet_oneday_360(it_jet_readymatch,3),'o-','markersize',2,'markerfacecolor','w')
% m_plot(jet_oneday(it_jet_readymatch(it_out(it_outmatch)),2),jet_oneday(it_jet_readymatch(it_out(it_outmatch)),3),'>','markerfacecolor','k')
% m_plot(jet_oneday(it_jet_readymatch(it_on),2),jet_oneday(it_jet_readymatch(it_on),3),'s-','markersize',8,'markerfacecolor','g')
% m_plot(jet_oneday(it_jet_readymatch(it_in),2),jet_oneday(it_jet_readymatch(it_in),3),'p-','markersize',5,'markeredgecolor','g')

% figure(2)
% shp = alphaShape(eddypoint_lon,eddypoint_lat);
% plot(shp)
