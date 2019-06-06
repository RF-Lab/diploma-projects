%EASY4   The master reciever position is computed like in EASY3.
%        Next the observations taken by the rover receiver are
%        introduced and the function baseline returns the baseline
%        components epoch by epoch.
%	     Note that the sequence of satellites in the stored data is
%        not the same at master and rover receivers. Therefore we
%        must introduce a matching mechanism.


% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

%----
Smean_Mean_pos=0;
Smean_Base_comp=0;


gend = 5; %ciclpo povtor

for g = 1:gend
    
    % We identify the master observation file and open it
    ofile1 = 'SITE24~1.01O';
    fid1 = fopen(ofile1,'rt');
    [Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
    NoObs_types1 = size(Obs_types1,2)/2;
    Pos = [];
    Gdop = [];
    % There are 22 epochs of data in ofile1
    qend = 22;
    
    for q = 1:qend
        [time1, dt1, sats1, eof1] = fepoch_0(fid1);
        NoSv1 = size(sats1,1);
        % We pick the observed C1 pseudoranges
        obs1 = grabdata(fid1, NoSv1, NoObs_types1);
        i = fobs_typ(Obs_types1,'C1');
        
        rand_obs=obs1(:,i)+randn(7,1)*10;
        
        [pos, el, gdop] = recpo_ls(rand_obs(:,i),sats1,time1,Eph);
        %[pos, el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);%obs1(:,i)+randn(7,1)*10
        Gdop = [Gdop gdop];
        Pos = [Pos pos];
    end
    me = mean(Pos,2);
    spread = std(Pos,1,2)
    fprintf('\nMean position as computed from %2.0f epochs:',qend)
    fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n\n', me(1,1), me(2,1), me(3,1))
    
    xyz_Mean_pos(:,gend)=me;%mean pos 22 epochs
    
    
    % figure(1);
    % plot((Pos(1:3,:)-Pos(1:3,1)*ones(1,q))','linewidth',2)
    % title(['Variation of Receiver Coordinates Over ',int2str(qend),' Epochs'],'fontsize',16)
    % legend('X','Y','Z')
    % xlabel('Epochs [1 s interval]','fontsize',16)
    % ylabel('[m]','fontsize',16)
    % print -depsc easy4_1
    
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
        % master observations: obs1, and rover observations: obs2
        
        rand_obs1=obs1+randn(7,1)*10;
        
        [omc,base] = baseline(master_pos,rand_obs1,obs2,sats1,time1,Eph);%obs1+randn(7,1)*10
        %[omc,base] = baseline(master_pos,obs1,obs2,sats1,time1,Eph);%obs1+randn(7,1)*10
        Omc = [Omc, omc];
        bases = [bases base];
    end
    me1 = mean(bases,2);   
    spread1 = std(bases,1,2)
    fprintf('\nBaseline Components as Computed From %2.0f Epochs:',qend)
    fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f', me1(1,1),me1(2,1),me1(3,1))
    
    xyz_Base_comp(:,gend) = me1;%Baseline Components  22 Epochs
    
    %figure(2);
    %plot((bases-bases(:,1)*ones(1,q))','linewidth',2)
    %title(['Variation of Baseline Components Over ',int2str(qend),' Epochs'],'fontsize',16)
    %legend('X','Y','Z')
    %xlabel('Epochs [1 s interval]','fontsize',16)
    %ylabel('[m]','fontsize',16)
    %set(gca,'fontsize',16)
    %legend 
    
    
    %break
    
    % figure(3);
    % plot(Gdop,'linewidth',2)
    % axis([1 length(Gdop) 0 5])
    % title('GDOP')
    
end
    %mean value from xyz_Mean_pos
    Smean_Mean_pos = mean(xyz_Mean_pos') % ==? spread
    %mean value from по xyz_Base_comp
    Smean_Base_comp = mean(xyz_Base_comp')   %==?  spread1

    figure(1);
    plot(rand_obs(:,1))  %  rand_obs=obs1+randn(7,1)*10
    
    figure(2);
    plot(rand_obs1(:,1))  %  rand_obs1=obs1+randn(7,1)*10
   
    
print -deps easy4_2




%%%%%%%%%%%%%%%%%%% end easy4.m %%%%%%%%%%%%%%%

