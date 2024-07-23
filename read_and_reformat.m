for func_num=1:29
    func = string(func_num);

    RDE = load('BC-SOPs\RDE\RDE_F'+func+'_Min_EV.mat').data;
    BlockEA = load('BC-SOPs\BlockEA\BlockEA_F'+func+'_Min_EV.mat').combinedMatrix;
    IEACOP = load('BC-SOPs\IEACOP\#1570992402_F'+func+'_Min_EV.mat').Min_EV;
    if(func_num > 1)
        func2 = string(func_num+1);
        jSOa = dlmread('BC-SOPs\jSOa\jSOa_D30_S'+func2+'.txt');
    else
        jSOa = dlmread('BC-SOPs\jSOa\jSOa_D30_S'+func+'.txt');
    end
    L_SRTDE = dlmread('BC-SOPs\L-SRTDE\L-SRTDE_F'+func+'_D30.txt');
    mLSHADE_LR = dlmread('BC-SOPs\mLSHADE_LR\mLSAHDE_LR_F#'+func+'_D#30.mat');

    dlmwrite('RDE_'+func+'.txt',RDE,'delimiter','\t','precision',20)
    dlmwrite('BlockEA_'+func+'.txt',BlockEA,'delimiter','\t','precision',20)
    dlmwrite('IEACOP_'+func+'.txt',IEACOP,'delimiter','\t','precision',20)
    dlmwrite('jSOa_'+func+'.txt',jSOa,'delimiter','\t','precision',20)
    dlmwrite('L_SRTDE_'+func+'.txt',L_SRTDE,'delimiter','\t','precision',20)
    dlmwrite('mLSHADE_LR_'+func+'.txt',mLSHADE_LR,'delimiter','\t','precision',20)
end

% avg = zeros(6,1000);
% for i=1:999
%     avg(1,i) = mean(RDE(i,:));
%     avg(2,i) = mean(BlockEA(i,:));
%     avg(3,i) = mean(IEACOP(i,:));
%     avg(4,i) = mean(jSOa(i,:));
%     avg(5,i) = mean(L_SRTDE(i,:));
%     avg(6,i) = mean(mLSHADE_LR(i,:));
% end

% plot(avg')
% set(gca, 'YScale', 'log');