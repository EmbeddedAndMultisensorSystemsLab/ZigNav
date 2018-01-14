%   This is a DEMO of ZigNav, a wireless positioning system that uses
%   dynamically generated radiomaps to localize an object using received
%   signal strength measurements.

%   Some definitions
%   Anchors     :   wireless nodes with known location
%   mobile node :   wireless nodes to be localized

%   ZigNav uses cooperative measurements exchanged between wireless nodes
%   to construct a dense radio map and to use the generated maps to
%   localize mobile objects.

%   ZigNav uses a hybrid approach that combines gaussian process
%   regression models and pathloss models to generate radiomaps from the
%   cooperative measurements exchanged between wireless nodes.

%   Radiomap generation code is not given in this package. Only the
%   dynamically generated maps data are given. The DEMO show the 
%   localization process using k-nearest-neighbourhood (knn) algorithm.

%   This is a 20-minutes stationary testing scenario. Dynamic scenarios and 
%   integration with other sensors such as inertial measurement units 
%   and vision sensors are currently under development.

%   Copyright (C) 2018, Mohamed Atia, all rights reserved.
%   The software is given under GNU Lesser General Public License
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.

%   For the hybrid radiomap generation method, check the following paper
%   M.M. Atia, M. Korenberg, A. Noureldin, “Dynamic Online-Calibrated 
%   Radio Maps for Indoor Positioning In Wireless Local Area Networks”, 
%   IEEE Trans. on Mobile Computing, Vol.12, no.9, pp. 1774 - 1787, 2013

%   For commercial use, please contact mohamed.atia@carleton.ca

%--> settings
set(0,'DefaultFigureWindowStyle','docked');
clear;
close all;
clc;

%--> Experimental Setup
number_of_wireless_nodes = 5;
mobile_node_index = 1;

%--> Radiomap generation method
%0: Hybrid: Gaussian process regression and pathloss models
%1: Gaussian process regression only
%2: Pathloss models only
ModelingTechnique = 0;

%--> Settings for positioning algorithm
rss_var_threshold = 8;
knn_num_of_records = 7;

%--> Area dimension
fid = fopen('..\experiment_data\Area_dimension.txt');
fgets(fid);line = fgets(fid);
line_data = textscan(line,'%d %d %f %f','delimiter', ',');
SCALE = line_data{3}/line_data{4}; %pixels per meters
max_x = line_data{1}/SCALE;% convert max x pixel dimension to meters
max_y = line_data{2}/SCALE;% convert max y pixel dimension to meters
fclose(fid);

%--> Area width and height
W = max_x;
L = max_y;

%--> Radiomaps data structure suitable for knn search
EstimatedGlobalRadioMap_For_Search = zeros(W*L,2+number_of_wireless_nodes);
EstimatedGlobalRadioMapVariances_For_Search = zeros(W+L,2+number_of_wireless_nodes);
EstimatedGlobalRadioMapVariancesInMeters_For_Search = zeros(W+L,2+number_of_wireless_nodes);

%--> error histograms initialization
error_bins = -5:0.1:5;
x_error_vector(1) = 0.0;
y_error_vector(1) = 0.0;
x_error_histogram_vector = histcounts(x_error_vector,error_bins);
y_error_histogram_vector = histcounts(y_error_vector,error_bins);

%--> Create to display display wireless nodes location information
figure;

%error histograms subplots
subplot(1,2,1);
h_x_histogram = bar(error_bins(1:end-1),x_error_histogram_vector,'b');grid on;
xlabel('x-error histogram(m)');ylabel('%');
subplot(1,2,2);
h_y_histogram = bar(error_bins(1:end-1),y_error_histogram_vector,'b');grid on;
xlabel('y-error histogram(m)');ylabel('%');

%save plots handles to refresh the plots with new values
set(h_x_histogram,'XDataSource','error_bins(1:end-1)');
set(h_x_histogram,'YDataSource','x_error_histogram_vector');
set(h_y_histogram,'XDataSource','error_bins(1:end-1)');
set(h_y_histogram,'YDataSource','y_error_histogram_vector');

%Subplot to display the area with wireless nodes location
figure;
%subplot(2,2,[3 4]);
axis([-1 (1.2)*max_x -1 (1.2)*max_y]);grid on;
title(sprintf('Experiment Area:%dm x %dm',max_x,max_y),'Fontsize',8);
xlabel('x(m)');ylabel('y(m)');

%--> extract anchor locations
fid = fopen('..\experiment_data\AP_Pixel_Based_Locations.txt');

ap_anchor_locations{1}.number_of_ap_anchors_pairs = 0;
while ~feof(fid)
    line = fgets(fid);
    line_data = textscan(line,'%d %d %s %s','delimiter', ';');
    anchor = char(line_data{4});
    ap = char(line_data{3});
    x = line_data{1};
    y = line_data{2};
    anchor_idx = 0;
    for i = 1:ap_anchor_locations{1}.number_of_ap_anchors_pairs
        if strcmpi(anchor , ap_anchor_locations{i+1}.anchor_id) == 1
            anchor_idx = i+1;
            break;
        end
    end
    if anchor_idx == 0 %new anchor?
        ap_anchor_locations{1}.number_of_ap_anchors_pairs = ap_anchor_locations{1}.number_of_ap_anchors_pairs +1;
        anchor_idx = ap_anchor_locations{1}.number_of_ap_anchors_pairs+1;
        ap_anchor_locations{anchor_idx}.anchor_id = anchor;
        ap_anchor_locations{anchor_idx}.ap_mac = ap;
        ap_anchor_locations{anchor_idx}.x = double(x)/SCALE;
        ap_anchor_locations{anchor_idx}.y = double(max_y)-double(y)/SCALE;%to convert from image pixel axis system to normal Cartesian axis system
    end
    
    %--> set anchor location info
    info = sprintf('%s\r\n[%.0f,%.0f]',ap_anchor_locations{anchor_idx}.ap_mac,...
        ap_anchor_locations{anchor_idx}.x,ap_anchor_locations{anchor_idx}.y);
    
    %--> if this is the wireless node to be localized, record and display
    %some information such as ground-truth location.
    if anchor_idx == mobile_node_index + 1
        ap_anchor_locations{anchor_idx}.est_x = ap_anchor_locations{anchor_idx}.x;
        ap_anchor_locations{anchor_idx}.est_y = ap_anchor_locations{anchor_idx}.y;
        text(ap_anchor_locations{anchor_idx}.x , ap_anchor_locations{anchor_idx}.y-1 , info,'Color','red','FontSize',7,'color','r');
        hold on;
        plot(ap_anchor_locations{anchor_idx}.x,ap_anchor_locations{anchor_idx}.y,'*','MarkerSize',5);
        mobile_node_est_x = ap_anchor_locations{anchor_idx}.est_x;
        mobile_node_est_y = ap_anchor_locations{anchor_idx}.est_y;
        mobile_node_ref_x = mobile_node_est_x;
        mobile_node_ref_y = mobile_node_est_y;
        h_mobile_node_plot = plot(mobile_node_est_x,mobile_node_est_y,'O','MarkerSize',12,'color','r');
        set(h_mobile_node_plot,'XDataSource','mobile_node_est_x');
        set(h_mobile_node_plot,'YDataSource','mobile_node_est_y');
    else
        text(double(ap_anchor_locations{anchor_idx}.x) , double(ap_anchor_locations{anchor_idx}.y) , info,'Color','red','FontSize',7,'color','k');
    end
    
end
fclose(fid);

%--> fill online anchors data
online_measurements_table{1,1}.number_of_anchors = number_of_wireless_nodes;
online_measurements_table{1,1}.number_of_aps = number_of_wireless_nodes;

for i = 1:number_of_wireless_nodes
    online_measurements_table{1+i,1}.ap_mac = ap_anchor_locations{i+1}.ap_mac;
    online_measurements_table{1+i,1}.x = ap_anchor_locations{i+1}.x;
    online_measurements_table{1+i,1}.y = ap_anchor_locations{i+1}.y;
    online_measurements_table{1+i,1}.info_updated = 0;
    online_measurements_table{1,1+i}.anchor_id = ap_anchor_locations{i+1}.anchor_id;
    online_measurements_table{1,1+i}.x = ap_anchor_locations{i+1}.x;
    online_measurements_table{1,1+i}.y = ap_anchor_locations{i+1}.y;
end

%--> Initialize plots handles
for ap_idx = 2:number_of_wireless_nodes+1
    h_radiommap{ap_idx} = [];
end

% --> read cooperative rss measurements
fid_rssi_data = fopen('..\experiment_data\ZigNav_logs_Dec.6.2017.txt');

processing_cycle_iteration = 1;

while (~feof(fid_rssi_data))
    %--> read one line record of cooperative rss measurements
    data_line = fgets(fid_rssi_data);
    data_items = strsplit(data_line);
    internal_system_time = data_items{2};
    line_data{1} = data_items{2};
    
    %--> convert rss values into double
    for i = 5:68
        str_value = data_items{i};
        str_value = str_value(1:end-1);
        line_data{i-3} = str2double(str_value);
    end
    
    %--> pre_process the received cooperative rss values
    %--> for each wireless node, fill the rss measurement vector
    %--> the rss value 99 is reserved to indicates nonapplicable rss
    %--> e.g. if a node cannot see a certain node, its rss is set to '99'
    %--> 99 is then replaced with a very weak rss value (-110dBm)
    for i = 1:number_of_wireless_nodes+1
        for j = 1:number_of_wireless_nodes+1
            online_measurements_table{i+1,1}.info_updated = 0;
            RSS_value = line_data((i-1)*8+j+1);
            RSS_value = RSS_value{1};
            online_measurements_table{1+i,1+j}.RSS = RSS_value;
            if (online_measurements_table{1+i,1+j}.RSS == 99)
                if i == j
                    online_measurements_table{1+i,1+j}.RSS = -5;
                else
                    online_measurements_table{1+i,1+j}.RSS = -110;
                end
            end
        end
    end
    
    %--> load the dynamic radiomap file generated at this iteration
    if ModelingTechnique == 0
        radio_map_file = sprintf('..\\dynamic_radiomaps\\hybrid_radiomaps\\hybrid_radio_map_in_cycle_%d.mat',processing_cycle_iteration);
    elseif ModelingTechnique == 1
        radio_map_file = sprintf('..\\dynamic_radiomaps\\zero_mean_gpr__radiomaps\\gpr_only_radio_map_in_cycle_%d.mat',processing_cycle_iteration);
    elseif ModelingTechnique == 2
        radio_map_file = sprintf('..\\dynamic_radiomaps\\pathloss_model_radiomaps\\pm_only_radio_map_in_cycle_%d.mat',processing_cycle_iteration);
    else
        radio_map_file = sprintf('..\dynamic_radiomaps\hybrid_radiomaps\radio_map_in_cycle_%d.mat',processing_cycle_iteration);
    end
    load(radio_map_file);
    
    %--> Draw the radiomap of each anchor
    for ap_idx = 2:1+ online_measurements_table{1,1}.number_of_aps
        if ap_idx == mobile_node_index+1
            continue;
        end
        rss_map_gpr = EstimatedGlobalRadioMap(:,:,ap_idx-1);
        
        if processing_cycle_iteration == 1
            figure(ap_idx);hold on;
        end
        
        if  isempty(h_radiommap{ap_idx})
            h_radiommap{ap_idx} = surf(rss_map_gpr');colorbar('vert');xlabel('x');ylabel('y');view(-47,26);grid on;
            rss_map_gpr_t = rss_map_gpr';
            set(h_radiommap{ap_idx},'ZDataSource','rss_map_gpr_t');
            title(sprintf('Radiomap of anchor ID:%s, anchor location(m) : {%.0f,%.0f}',...
                online_measurements_table{ap_idx,1}.ap_mac,...
                online_measurements_table{ap_idx,1}.x, online_measurements_table{ap_idx,1}.y),'Fontsize',8);
            zlabel('rss(dBm)');
        else
            rss_map_gpr_t = rss_map_gpr';
            refreshdata(h_radiommap(ap_idx));
            Delay(.01);
        end
        
        %--> Prepare radiomap structures for suitable knn search
        Cell_Idx = 1;
        for x=1:max_x*1.5
            for y=1:max_y*1.5
                
                EstimatedGlobalRadioMap_For_Search(Cell_Idx,1) = x;
                EstimatedGlobalRadioMap_For_Search(Cell_Idx,2) = y;
                EstimatedGlobalRadioMap_For_Search(Cell_Idx,2+ap_idx-1) = EstimatedGlobalRadioMap(x,y,ap_idx-1);
                
                EstimatedGlobalRadioMapVariances_For_Search(Cell_Idx,1) = x;
                EstimatedGlobalRadioMapVariances_For_Search(Cell_Idx,2) = y;
                EstimatedGlobalRadioMapVariances_For_Search(Cell_Idx,2+ap_idx-1) = EstimatedGlobalRadioMapVariances(x,y,ap_idx-1);
                EstimatedGlobalRadioMapVariancesInMeters_For_Search(Cell_Idx,2+ap_idx-1) = EstimatedGlobalRadioMapVariancesInMeters(x,y,ap_idx-1);
                
                Cell_Idx = Cell_Idx + 1;
            end
        end
    end
    
    %--> Localize the targeted mobile node using knn algorithm
    % record the index of anchors used in localization
    % and record the rss that the mobile node observes from anchors
    anchors_to_use_in_localization = [];
    average_rss_vector = [];
    for j = 1:number_of_wireless_nodes
        if j ~= mobile_node_index
            anchors_to_use_in_localization = [anchors_to_use_in_localization j];
            average_rss_vector = [average_rss_vector online_measurements_table{mobile_node_index+1,1+j}.RSS];
        end
    end
    
    %--> select the "healthy" radiomap records. These are the radiomap 
    %    records  that have minimum confidence level (variance)
    variance_vector = EstimatedGlobalRadioMapVariances_For_Search(:,2+anchors_to_use_in_localization)';
    MeanVarOfEstimatedGlobalRadioMapIndBm = mean(variance_vector);
    HealthyPointsIndeciesInTheEstimatedGlobalRadioMap = MeanVarOfEstimatedGlobalRadioMapIndBm <= rss_var_threshold;
    HealthyPointsIndecies = 1:length(HealthyPointsIndeciesInTheEstimatedGlobalRadioMap);
    HealthyPointsIndecies = HealthyPointsIndecies.*HealthyPointsIndeciesInTheEstimatedGlobalRadioMap;
    HealthyPointsIndecies(HealthyPointsIndecies == 0) = [];

    %--> calculate euclidean distances in rss domain between rss vector 
    %   of the mobile node and the rss vectors in the radiomap
    DistancesInRSS = zeros(length(HealthyPointsIndecies),1);
    LoopCounter = 1;
    for LoopCounter = 1:length(HealthyPointsIndecies)
        Cell_Idx = HealthyPointsIndecies(LoopCounter);
        RM_RSSI_Vect = EstimatedGlobalRadioMap_For_Search(Cell_Idx,2+anchors_to_use_in_localization);
        DistancesInRSS(LoopCounter) = sqrt(sum(( RM_RSSI_Vect-average_rss_vector).^2 ));
    end
    [ sorted_dist , IDX ]= sort(DistancesInRSS);
    
    %--> knn-algorithm
    EstimatedEP = 0;
    EstimatedNP = 0;
    AvgLen = length(sorted_dist)*(length(sorted_dist) < knn_num_of_records)+knn_num_of_records*(length(sorted_dist) > knn_num_of_records);
    Totalweights = 0;
    Estimatedvariance = 0;
    for x = 1:AvgLen
        weight = exp(-1*(sorted_dist(x)/sum(sorted_dist(1:AvgLen))));
        Totalweights = Totalweights + weight;
        EstimatedEP = EstimatedEP + EstimatedGlobalRadioMap_For_Search(HealthyPointsIndecies(IDX(x)),1)*weight;
        EstimatedNP = EstimatedNP + EstimatedGlobalRadioMap_For_Search(HealthyPointsIndecies(IDX(x)),2)*weight;
        Estimatedvariance = Estimatedvariance + mean(EstimatedGlobalRadioMapVariancesInMeters_For_Search(HealthyPointsIndecies(IDX(x)),2+anchors_to_use_in_localization))*weight;
    end
    EstimatedEP = EstimatedEP/Totalweights;
    EstimatedNP = EstimatedNP/Totalweights;
    Estimatedvariance = Estimatedvariance/Totalweights;
    HDOP = Estimatedvariance^.5;
    
    %--> display the estimated location and calculate the positioning error
    ap_anchor_locations{mobile_node_index+1}.est_x = EstimatedEP;
    ap_anchor_locations{mobile_node_index+1}.est_y = EstimatedNP;
    mobile_node_est_x = EstimatedEP;
    mobile_node_est_y = EstimatedNP;
    x_error = mobile_node_est_x - mobile_node_ref_x;
    y_error = mobile_node_est_y - mobile_node_ref_y;
    x_error_vector(length(x_error_vector) + 1) = x_error;
    y_error_vector(length(y_error_vector) + 1) = y_error;
    x_error_histogram_vector = histcounts(x_error_vector,error_bins);
    y_error_histogram_vector = histcounts(y_error_vector,error_bins);
    refreshdata(h_x_histogram);
    refreshdata(h_y_histogram);
    refreshdata(h_mobile_node_plot);

    %--> update the processing iteration
    processing_cycle_iteration = processing_cycle_iteration + 1;

end

disp('[PROCESSING FINISHED]');