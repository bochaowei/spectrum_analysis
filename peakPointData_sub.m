
% function for extracting data points in pd_current_plot.m using peak val
% alex crawford
% 12/30/2021

% isolates peak of 87 2-2 transition from data
% extracts data of height and chosen array for clean side of curve

function [dataPt, chooseTime, chooseA, ptm, ptam, F2] = peakPointData_sub(ChA,ChB,testTime, promVal)
% insert function explaination

    % smooth incoming data to help isolate peaks
    order = 2;
    framelen = 51;
    sgf = sgolayfilt(ChB, order, framelen);

    
    % isolating 87 peak by prominence of local max
    abs_max = max(sgf);
    [TF1max,Pmax] = islocalmax(sgf, 'MinProminence',promVal);
    prom_chb_max = ChB(TF1max);
    prom_cha_max = ChA(TF1max);
    ptam = prom_cha_max;
    prom_t_max = testTime(TF1max);
    ptm = prom_t_max;

    if isempty(Pmax)
        disp('ZERO PEAKS ... data lost');
    end

    % finding third peak from global max of sat spec
    ind_max = find(sgf(TF1max) == abs_max,1);
    
    % find min of chA by pulling index of chB max
%     minVal = min(ChA);    
    findMin = find(ChB == max(ChB),1);

    minVal = ChA(findMin);
%     findMin = find(ChA==minVal,1);


    try
        % Left side of max
        F2_ind_pmax = ind_max - 3;
        % reverse search index of isolated peak
        % NOTE: indication that probe beam is not perfectly perpendicular
        % to atom beam by the peak and current max not aligning perfectly

        rev_index_L = find(testTime == prom_t_max(F2_ind_pmax),1);
        indBuff = rev_index_L - 600; 
        pkMax = find(ChA == max(ChA(indBuff:rev_index_L+600)));
        ind_start = pkMax;
        % starting at peak of 2-2 transition ... max voltage val needed
        chooseTime = testTime(ind_start:findMin);
        chooseA = ChA(ind_start:findMin);
        maxVal = max(chooseA);
        dataPt = maxVal-minVal;
    catch
        disp('Try right peak');
        try
            % Right Side of max
            F2_ind_pmax = ind_max + 3;
            % check that index selected is not another max val
            if prom_chb_max(F2_ind_pmax)-abs_max == 0
                disp('Sub nan ... filter max');
                % last resort fill with nan ... usually doesnt happen
                chooseTime = testTime;
                chooseA = ChA;
                prom_t_max(F2_ind_pmax) = NaN;
                prom_cha_max(F2_ind_pmax) =NaN;
                dataPt = NaN;
            else
                % reverse search index of isolated peak
                disp('Use right peak');
                rev_index_R = find(testTime == prom_t_max(F2_ind_pmax),1);
                ind_start = rev_index_R +600;
%                 indBuff = rev_index_R + 300;
%                 pkMax = find(ChA == max(ChA(rev_index_R-300:indBuff)));
%                 ind_start = pkMax;
                % starting at peak of 2-2 transition ... max voltage val needed
                chooseTime = testTime(findMin:ind_start);
                chooseA = ChA(findMin:ind_start);
                maxVal = max(chooseA);
                dataPt = maxVal-minVal;
            end

        catch
            disp('Substitute NaN');
            % last resort fill with nan ... full data set plotted
            chooseTime = testTime;
            chooseA = ChA;
            prom_t_max(F2_ind_pmax) = NaN;
            prom_cha_max(F2_ind_pmax)= NaN;
            dataPt = NaN;
        end

    end
    F2 = F2_ind_pmax;
%        hold on
% %       % optional plot to visualize chosen peaks
%  %        plot(chooseTime,chooseA,prom_t_max(F2_ind_pmax),prom_cha_max(F2_ind_pmax),'ro');
%         plot(testTime, ChA, prom_t_max(F2_ind_pmax),prom_cha_max(F2_ind_pmax),'ro');
% % 
%        hold off
end