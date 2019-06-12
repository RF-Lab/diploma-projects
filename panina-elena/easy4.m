%EASY4   The master reciever position is computed like in EASY3.
%        Next the observations taken by the rover receiver are
%        introduced and the function baseline returns the baseline
%        components epoch by epoch.
%	     Note that the sequence of satellites in the stored data is
%        not the same at master and rover receivers. Therefore we
%        must introduce a matching mechanism.

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');
%----
sr_kv = 0; %��� ���������� �������� ��.��.

mean_dif_me = [];%Receiver - ������ ��.��. ��� N ���-�� �������� �� ������� ����� 
mean_dif_me1 = [];%Baseline - ������ ��.��. ��� N ���-�� �������� �� ������� �����
zn_sigm = []; %��� ������ �������� ����� (��� ������ �� ������)
    
sigm = 5;
for sigm = 5:1:8 % ���� �� ��������� �����(�� 5 �� 8)

    mean_sigm = [];    
    dif_me = [];%Receiver - ���������� ��� XYZ �������� ��� ������ �������� ������� �������� [3*gend]
    dif_me1 = [];%Baseline - ���������� ��� XYZ �������� ��� ������ �������� ������� �������� [3*gend]
    
    gend = 100; %���� �� ���������� �������� ��������
    for g = 1:gend
        % We identify the master observation file and open it
        ofile1 = 'SITE24~1.01O';
        fid1 = fopen(ofile1,'rt');
        [Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
        NoObs_types1 = size(Obs_types1,2)/2;
        Pos = [];
        rand_Pos =[];
        Gdop = [];
        
        error = randn(7,1)*sigm; %��� ��� ����������!!!
                
        % There are 22 epochs of data in ofile1
        qend = 22;
        for q = 1:qend %�� �����
            [time1, dt1, sats1, eof1] = fepoch_0(fid1);
            NoSv1 = size(sats1,1);
            % We pick the observed C1 pseudoranges
            obs1 = grabdata(fid1, NoSv1, NoObs_types1);
            i = fobs_typ(Obs_types1,'C1');
            
            %without error ��� ����
            [pos, el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);%obs1(:,i)+randn(7,1)*10
            Gdop = [Gdop gdop];
            Pos = [Pos pos];
            %with error  ��������� ���
            rand_obs=obs1(:,i)+error;
            [rand_pos, el, gdop] = recpo_ls(rand_obs(:,i),sats1,time1,Eph);
            rand_Pos = [rand_Pos rand_pos];
        end
        
        %without error
        me = mean(Pos,2); % �������� ���������� xyz        
        %with error
        rand_me = mean(rand_Pos,2);% �������� ���������� xyz
        
        %���������� �������� ��� ������ �������� ������� �������� [3*gend]
        dif_me = [dif_me rand_me];        
        %------------        
        % we need to close all open files and then open to read from the beginning
        fclose all;
        
        ofile1 = 'SITE24~1.01O';
        fid1 = fopen(ofile1,'rt');
        [Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
        NoObs_types1 = size(Obs_types1,2)/2;
        % Next we include the rover and identify the rover
        % observation file and open it
        ofile2 = 'SITE247j.01O';
        fid2 = fopen(ofile2,'rt');
        [Obs_types2, ant_delta2, ifound_types2, eof12] = anheader(ofile2);
        NoObs_types2 = size(Obs_types2,2)/2;
        master_pos = me;  % best possible estimate of master position
        bases = [];
        Omc = [];
        
        rand_bases = [];
        rand_Omc = [];
        
        for q = 1:qend
            [time1, dt1, sats1, eof1] = fepoch_0(fid1);
            [time2, dt2, sats2, eof2] = fepoch_0(fid2);
            if time1 ~= time2
                disp('Epochs not corresponding')
                break
            end;
            NoSv1 = size(sats1,1);
            NoSv2 = size(sats2,1);
            % We pick the observations
            obsm = grabdata(fid1, NoSv1, NoObs_types1);
            obsr = grabdata(fid2, NoSv2, NoObs_types2);
            i = fobs_typ(Obs_types1,'C1');
            obs1 = obsm(:,i);
            for s = 1:NoSv1
                ind = find(sats1(s) == sats2(:));
                obs2(s,1) = obsr(ind,1);
            end
            %master observations: obs1, and rover observations: obs2
            
            %without error ��� ����
            [omc,base] = baseline(master_pos,obs1,obs2,sats1,time1,Eph);%obs1+randn(7,1)*10
            Omc = [Omc, omc];
            bases = [bases base];
            %with error ��������� ���
            rand_obs1=obs1+error;
            [rand_omc,rand_base] = baseline(master_pos,rand_obs1,obs2,sats1,time1,Eph);
            rand_Omc = [rand_Omc, rand_omc];
            rand_bases = [rand_bases rand_base];
        end
        %without error
        me1 = mean(bases,2);%�������� ���������� xyz        
        %with error
        rand_me1 = mean(rand_bases,2);%�������� ���������� xyz
        
        %���������� �������� ��� ������ �������� ������� �������� [3*gend]
        dif_me1 = [dif_me1 rand_me1];
                
    end
        %Receiver - ������� ������� �����. ��� N ���-�� �������� �� ������� �����
        sr_kv = std(dif_me,1,2);        
        %Receiver - ������ ��.��. ��� N ���-�� �������� �� ������� �����
        mean_dif_me = [mean_dif_me sr_kv];
    
        %Baseline - ������� ������� �����. ��� N ���-�� �������� �� ������� �����
        sr_kv1 = std(dif_me1,1,2);        
        %Baseline - ������ ��.��. ��� N ���-�� �������� �� ������� �����
        mean_dif_me1 = [mean_dif_me1 sr_kv1];    
        
        %����� �������� ����� (����� ������ �������)
        zn_sigm = [zn_sigm sigm];
end
    %------�������
     figure(1);
     % �� ��� � �������� �����, �� � ��.��.
     plot(zn_sigm, mean_dif_me(1:3,:),'linewidth',2)
     title('Variation of Receiver ','fontsize',16)
     legend('X','Y','Z')
     xlabel('Sigma','fontsize',16)
     ylabel('[m]','fontsize',16)
     print -depsc easy4_1
     
     figure(2);
     % �� ��� � �������� �����, �� � ��.��.
     plot(zn_sigm, mean_dif_me1(1:3,:),'linewidth',2)
     title('Variation of Baseline Components ','fontsize',16)
     legend('X','Y','Z')
     xlabel('Sigma','fontsize',16)
     ylabel('[m]','fontsize',16)
     print -depsc easy4_2
    

%%%%%%%%%%%%%%%%%%% end easy4.m %%%%%%%%%%%%%%%