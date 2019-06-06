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
Smean_Mean_pos=0;
Smean_Base_comp=0;


gend = 5; %цикл по повторяемости

for g = 1:gend
    
    % We identify the master observation file and open it
    ofile1 = 'SITE24~1.01O';
    fid1 = fopen(ofile1,'rt');
    [Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
    NoObs_types1 = size(Obs_types1,2)/2;
    Pos = [];
    rand_Pos =[];
    Gdop = [];   
    dif_xyz =[];
    xyz_Mean_pos=[];
    % There are 22 epochs of data in ofile1
    qend = 22;
    
    for q = 1:qend
        [time1, dt1, sats1, eof1] = fepoch_0(fid1);
        NoSv1 = size(sats1,1);
        % We pick the observed C1 pseudoranges
        obs1 = grabdata(fid1, NoSv1, NoObs_types1);
        i = fobs_typ(Obs_types1,'C1');
        
        %without error
        [pos, el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);%obs1(:,i)+randn(7,1)*10
        Gdop = [Gdop gdop];
        Pos = [Pos pos];
        %with error
        rand_obs=obs1(:,i)+randn(7,1)*10;
        [rand_pos, el, gdop] = recpo_ls(rand_obs(:,i),sats1,time1,Eph);
        rand_Pos = [rand_Pos rand_pos];           
    end
    %without error
    me = mean(Pos,2);
    spread = std(Pos,1,2)
    fprintf('\nMean position as computed from %2.0f epochs:',qend)
    fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n\n', me(1,1), me(2,1), me(3,1))
    %with error
    rand_me = mean(rand_Pos,2);
    rand_spread = std(rand_Pos,1,2)
    xyz_Mean_pos(:,g)=rand_me;%Mean pos from
    
    %difference x-x_with_error...
    dif_me(1,g)=abs(me(1,1)-rand_me(1,1));
    dif_me(2,g)=abs(me(2,1)-rand_me(2,1));
    dif_me(3,g)=abs(me(3,1)-rand_me(3,1)); 

    
    
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
        
        %without error
        [omc,base] = baseline(master_pos,obs1,obs2,sats1,time1,Eph);%obs1+randn(7,1)*10
        Omc = [Omc, omc];
        bases = [bases base];                
        %with error
        rand_obs1=obs1+randn(7,1)*10;        
        [rand_omc,rand_base] = baseline(master_pos,rand_obs1,obs2,sats1,time1,Eph);
        rand_Omc = [rand_Omc, rand_omc];
        rand_bases = [rand_bases rand_base];
        
    end
    %without error
    me1 = mean(bases,2);   
    spread1 = std(bases,1,2)
    fprintf('\nBaseline Components as Computed From %2.0f Epochs:',qend)
    fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f', me1(1,1),me1(2,1),me1(3,1))
    %with error
    rand_me1 = mean(rand_bases,2);   
    rand_spread1 = std(rand_bases,1,2)    
    xyz_Base_comp(:,g) = rand_me1;
    
    %difference x-x_with_error...
    dif_me1(1,g)=abs(me1(1,1)-rand_me1(1,1));
    dif_me1(2,g)=abs(me1(2,1)-rand_me1(2,1));
    dif_me1(3,g)=abs(me1(3,1)-rand_me1(3,1)); 
        
end
    %------------
    
    %mean difference cooordinat's ME
    dif_xyz = mean(dif_me')  %mean array diff coordinat x-x_error   
    %print difference x-x_with_error...
    figure(1);
    plot(dif_me(:,:),'linewidth',2)
    title(['Differens Coordinates Over ',int2str(qend),' Epochs'],'fontsize',16)
    legend('X','Y','Z')
    xlabel('Epochs [1 s interval]','fontsize',16)
    ylabel('[m]','fontsize',16)
    print -depsc easy4_1
    
        %mean difference cooordinat's ME
    dif1_xyz = mean(dif_me1')  %mean array diff coordinat x-x_error   
    %print difference x-x_with_error...
    figure(2);
    plot(dif_me1(:,:),'linewidth',2)
    title(['Differens Coordinates Over ',int2str(qend),' Epochs'],'fontsize',16)
    legend('X','Y','Z')
    xlabel('Epochs [1 s interval]','fontsize',16)
    ylabel('[m]','fontsize',16)
    print -depsc easy4_2





%%%%%%%%%%%%%%%%%%% end easy4.m %%%%%%%%%%%%%%%



