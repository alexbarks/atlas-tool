function wss_batch(folder_name, TimeFlag)

% % MARFAN patients
% PATHNAME{1} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\AB20120210';
% PATHNAME{2} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\AM20150325';
% PATHNAME{3} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\AR20150428';
% PATHNAME{4} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\AS20140723';
% PATHNAME{5} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\AS20140814';
% PATHNAME{6} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\BK20140617';
% PATHNAME{7} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\BT20150804';
% PATHNAME{8} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\CB20120210';
% PATHNAME{9} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\CR20130924';
% PATHNAME{10} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\DB20141212';
% PATHNAME{11} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\DC20150624';
% PATHNAME{12} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\EM20150714';
% PATHNAME{13} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\ER20141014';
% PATHNAME{14} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\JB20150701';
% PATHNAME{15} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\JC20150819';
% PATHNAME{16} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\JJ20141211';
% PATHNAME{17} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\MH20150128';
% PATHNAME{18} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\MK20150414';
% PATHNAME{19} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\MK20150415';
% PATHNAME{20} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\MM20130102';
% PATHNAME{21} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\MS20140205';
% PATHNAME{22} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\PC20130220';
% PATHNAME{23} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\PS20150701';
% PATHNAME{24} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\RP20141003';
% PATHNAME{25} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\RS20150701';
% PATHNAME{26} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\SH20140627';
% PATHNAME{27} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\SM20150616';
% PATHNAME{28} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\ST20120817';
% PATHNAME{29} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\SW20140613';
% PATHNAME{30} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\TC20150611';
% PATHNAME{31} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Marfans\TT20150409';

% % CONTROLS
% PATHNAME{1} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\AC20141024';
% PATHNAME{2} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\AD20150106';
% PATHNAME{3} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\AH20141224';
% PATHNAME{4} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\CB20150213';
% PATHNAME{5} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\CF20150304';
% PATHNAME{6} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\CH20140924';
% PATHNAME{7} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\CI20150624';
% PATHNAME{8} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\CP20140716';
% PATHNAME{9} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\DA20150714';
% PATHNAME{10} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\DC20150212';
% PATHNAME{11} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\DF20150128';
% PATHNAME{12} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\EP20141024';
% PATHNAME{13} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\GR20140917';
% PATHNAME{14} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\JM20141014';
% PATHNAME{15} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\JS20150313';
% PATHNAME{16} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\KK20140821';
% PATHNAME{17} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\KM20150730';
% PATHNAME{18} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\MO20141219';
% PATHNAME{19} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\MW20141031';
% PATHNAME{20} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\NL20150721';
% PATHNAME{21} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\PR20140715';
% PATHNAME{22} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\RM20140723';
% PATHNAME{23} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\TC20150506';
% PATHNAME{24} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\TT20141029';
% PATHNAME{25} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\Normal\ZM20150415';

% % LDS
% PATHNAME{1} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\LDS\CG20140930';
% PATHNAME{2} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\LDS\JL20150508';
% PATHNAME{3} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\LDS\LB20150206';
% PATHNAME{4} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\LDS\MM20150225';
% PATHNAME{5} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\LDS\SW20141111';

% PATHNAME{1} = 'C:\Users\Emily\Desktop\Marfan_MRstructs\JB';

% scan all folders of the present directory
cd(folder_name)
folders = ls;
for i=1:size(folders,1)-2
    PATHNAME{i}=fullfile(pwd,folders(i+2,:));
end

currDir = pwd;
% Calculate WSS
for n = 1:size(PATHNAME,2)
    
    cd(PATHNAME{n})
    folders = ls;
    if exist('mrstruct','dir')==7
        MrstructPath = strcat(PATHNAME{n},'\mrstruct\');
    else
        for i=3:size(folders,1)
            [a,b]=find(folders(i,1:8)=='mrstruct');
            if sum(a)==8
                MrstructPath=strcat(PATHNAME{n},'\',folders(i,:),'\');
                break
            end
        end
    end
    cd(currDir);
    
    MimicsSegPath = strcat(PATHNAME{n},'\');
    
    % TimeFlag: 1 (5 systolic phases);  0 (peak systole)
    %TimeFlag=1;
    
    tic
    mimics_to_Wss([MrstructPath],[MimicsSegPath],1,1,1,1,TimeFlag,0,0,0,1);
    toc
    
%     close all
    
end