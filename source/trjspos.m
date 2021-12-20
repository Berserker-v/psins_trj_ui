function pos = trjspos(avp0, wat, ts, repeats)
% Trajectory simulator (Ref. See my master's dissertation p63).
%
% Prototype: trj = trjspos(avp0, wat, ts, repeats)
% Inputs: avp0 - initial attitude, velocity and position
%         wat - segment parameter, each column described as follows:
%            wat(:,1) - period lasting time
%            wat(:,2) - period initial velocity magnitude
%            wat(:,3:5) - angular rate in trajectory-frame
%            wat(:,6:8) - acceleration in trajectory-frame
%         ts - simulation time step
%         repeats - trajectory repeats using the same 'wat' parameter. 
% Output: pos - trajectory position

% global glv
    if nargin<4, repeats = 1; end
    wat1 = repmat(wat, repeats, 1);
    att = avp0(1:3); vn = avp0(4:6); pos = avp0(7:9);
    eth = earth(pos, vn);
    len = fix(sum(wat1(:,1))/ts);
    avp = zeros(len,10);
    ki = timebar(1, len, 'Trajectory Preview.');
    for k=1:size(wat1,1)
        lenk = round(wat1(k,1)/ts);  % the lenght at k phase
        wt = wat1(k,3:5)'; at = wat1(k,6:8)';
        ufflag = 0;
        if norm(wt)==0 && norm(at)==0  % uniform phase flag
            ufflag = 1; 
%             vnr = a2mat(att)*[0;wat1(k,2);0];  % velocity damping reference
        end
        for j=1:lenk
            sa = sin(att); ca = cos(att);
            si = sa(1); sk = sa(3); ci = ca(1); ck = ca(3); 
            Cnt = [ ck, -ci*sk,  si*sk; 
                    sk,  ci*ck, -si*ck; 
                    0,   si,     ci ];
            att = att + wt*ts;
            if ufflag==1  % damping
                vn1 = vn;
                vn01 = (vn+vn1)/2;  vn = vn1;
            else
                an = Cnt*at;
                vn1 = vn + an*ts;  vn01 = (vn+vn1)/2;  vn = vn1;  % velocity
            end
            dpos01 = [vn01(2)/eth.RMh;vn01(1)/eth.clRNh;vn01(3)]*ts;
%             eth = earth(pos+dpos01, vn01);
            pos = pos+dpos01;      % position
            
            kts = ki*ts;
            avp(ki,:) = [att; vn; pos; kts]';
            ki = timebar;
            if(-1==ki)
                pos=-1; return; 
            end
        end
    end
    avp(ki:end,:) = [];
    avp = iatt2c(avp);
    pos = avp(:,7:10);

