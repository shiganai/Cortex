function [sdata,ocf]=autosmooth(rdata,np,nFr,time_int,g)

% residual analysis (Wells and Winter's(1980)Method)
%
% Usage
%   [sdata,ocf]=autosmooth(rdata,np, nFr, time_int, cfmin,cfmax)
%
% Imput
%   rdata         : raw data
%   np            : number of point
%   nFr           : number of frame
%   time_int      : time interval
%   cfmin         : Minimum cut-off frequency for regression line 
%   cfmax         : Maximum  cut-off frequency for regression line
%   g             :graphic on='y'/ off='n'
%
% Output
%   sdata         : smoothed data
%   ocf           : optimal cut-off frequency
%
% -----------------------------------------------------------
% (�c�������Ȃ���cfmin��cfmax�����肷��)

fs = 1/time_int; %1�b������̃R�}��
cfmin = fs*7/100;
cfmax = fs*10/100;
sdata = zeros(size(rdata));       %rdata�̍s�񐔂�0�s��
sd = zeros(size(rdata,1),cfmax);  %�s : �R�}��, �� : cfmax
dim = size(rdata,2)/np;           %������

for iaxis = 1:dim
    
    for ip = 1:np
        
        for icf = 1:cfmax
            
                [B,A] = butter(2,icf/(fs/2),'low'); %�o�^���[�X���[�p�X�t�B���^�[�̐݌v
                sd(:,icf) = filtfilt(B,A,rdata(:,(iaxis+(ip-1)*dim))); %�t�B���^�[�ɂ��f�[�^�̏���

                dif = sd(:,icf)-rdata(:,(iaxis+(ip-1)*dim));  %���f�[�^�Ƃ̍���(�c������)
                
                R(icf,:) = [icf,sqrt((dif'*dif)/nFr)];   %�e���g���ɂ����鐶�f�[�^�ƕ�������f�[�^�Ƃ�RMS
             
        end
        
        x(:,1) = R(cfmin:cfmax,1);
        x(:,2) = ones;
        y = R(cfmin:cfmax,2);
        
        answer = x\y;  %�����̌X���c���Ƃ̌�_�Z�o
        a = answer(1);
        Ts = answer(2);% Threshold of Residual
        
        if a > 0
            errordlg('Impossible to compute the optimal cut-off frequency!')
         return
         
        end
        
        cf = 1;  %�ŏ���1�Ɛݒ肵�Ă���
        
         while R(cf,2) > Ts
             cf = cf+1;
         end
         
         ocf(ip,iaxis) = cf;
         sdata(:,(iaxis+(ip-1)*dim)) = sd(:,cf);
        
         ax = [0,cfmax,0,R(1,2)];
                
         % graphic
         if g == 'y'
             figure
                set(gcf,'name',[num2str(ip) ' point     : Cut-off Frequency =  ' num2str(ocf(ip,iaxis))],'DoubleBuffer','on');
                lab(1,:) = line([0 cfmax],[Ts Ts+cfmax*a],'linewidth',1,'linestyle','-','color',[0 0 0]);
                lab(2,:) = line([0 cfmax],[Ts Ts],'linewidth',1,'linestyle','--','color',[0 .5 1]);
                lab(3,:) = line(R(:,1),R(:,2),'linewidth',3,'linestyle','-','color',[0 0 0]);
                lab(4,:) = line([ocf(ip,iaxis) ocf(ip,iaxis)],[0 R(1,2)],'linewidth',1,'linestyle',':','color',[1 .5 0]);
                
                xlabel('Frequency [Hz]','fontsize',30,'FontName','Times');
                ylabel('RMS [m]','fontsize',30,'FontName','Times');
                
                axis('normal');axis(ax);
                set(gca,'xtick',0:fs/100:fs*20/100);
                
                drawnow
         end
    end
end

end